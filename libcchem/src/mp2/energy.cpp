#include "mp2/mp2.hpp" 

#include "core/wavefunction.hpp"
#include "core/integral.hpp"
#include "core/transform/transform.hpp"
//#include "core/transform/lambda.hpp"
#include "core/submatrix.hpp"
#include "runtime.hpp"

#include "foreach.hpp"
#include "utility/timer.hpp"

#include <algorithm>
#include <iostream>

#include <boost/array.hpp>
#include <boost/progress.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/home/phoenix/operator.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/scope/let.hpp>
#include <boost/spirit/home/phoenix/scope/lambda.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement/sequence.hpp>
#include <boost/spirit/include/phoenix_scope.hpp>

#include "phoenix/thread.hpp"
#include "phoenix/ublas.hpp"
#include "phoenix/array.hpp"
#include "phoenix/bind_name.hpp"
#include "phoenix/tensor_function.hpp"
#include "phoenix/iterator_function.hpp"

#include "file/array.hpp"
#include "utility/array.hpp"
#include "boost/utility/profiler.hpp"

namespace mp2 {
namespace detail {

    struct energy_impl {
	template<class A, typename T, class V1, class V2>
	struct result { typedef double type; };
	template<class A, typename T, class V1, class V2>
	double operator()(const A &I, const T &e0, const V1 &e1, const V2 &e2) const {
	    double value = 0;
	    for (size_t b = 0; b < I.size2(); ++b) {
    	    	for (size_t a = 0; a < I.size1(); ++a) {
    	    	    double ab = I(a,b);
    	    	    double ba = I(b,a);
    		    double e = e0 + e1[a] + e2[b];
    	    	    // E += ab*(2*ab - ba)/de;
    		    value += ab*(2*ab - ba)/e;
    	    	}
    	    }
	    // std::cout << value << std::endl;
	    return value;
	}
    };

    static const boost::phoenix::function<energy_impl> energy;

}
}

struct increment {
    template<typename T>
    struct result { typedef T& type; };
    template<typename T>
    T& operator()(T &t) const { ++t;  return t; }
};

//namespace ublas = boost::numeric::ublas;
double mp2::energy(const Wavefunction &wf, const Integral::Screening &screening) {
    const Wavefunction::Basis& basis = wf.basis();

    const Wavefunction::Orbitals& active = wf.active();
    const Wavefunction::Orbitals& virtuals =  wf.virtuals();

    const Wavefunction::matrix_type& C = wf.C();

    namespace ublas = boost::numeric::ublas;
    ::Transform::Matrix C1(ublas::trans(submatrix(C, basis, virtuals)));
    ::Transform::Matrix C2(ublas::trans(submatrix(C, basis, active)));
    ::Transform::Matrix &C3 = C1;
    ::Transform::Matrix &C4 = C2;

    double E = 0;

    typedef Wavefunction::vector_type::const_iterator vector_iterator;
    vector_iterator e = wf.energies().begin();

    vector_iterator ea(e + virtuals.start()), eb(e + virtuals.start());
    vector_iterator ei(e + active.start()), ej(e + active.start());
    size_t nij = (active.size()*active.size() + active.size())/2;



    using utility::make_array;

    // size_t chunk[3] = { 50, nij, C1.size1() };
    // //size_t chunk[3] = { C1.size1(), nij, 1 };
    // {
    // 	const size_t max = (1 << 25);
    // 	size_t size = chunk[1];
    // 	size_t n = (chunk[0]*chunk[1]*sizeof(double))/max;
    // 	if (!n) n = 1;
    // 	chunk[1] = chunk[1]/n + bool(chunk[1]%n); 
    // 	std::cout << "size[1]/chunk[1]: " << size << "/" << chunk[1] << std::endl;
    // }

    std::cout << "Secondary storage: "
	      << (C1.size1()*nij*basis.size())/1e9 << " Gbytes"
	      << std::endl;

    std::cout << "Memory: "
	      << (nij*basis.size()*basis.max().size())/1e9 << " Gbytes"
	      << std::endl;

    File file = File::create("data.h5");
    //File::Array<3,double> F(file, "oopr",
    File::Array<3,double> T(file, "aijr", 
			    make_array(C1.size1(), nij, basis.size()));

    Transform::Triangular T_(Transform::Order(1,3,0,2), runtime("MP2"));

    {
	BOOST_AUTO(const &size, make_array(nij, basis.size(), basis.max().size()));
	// File::Array<3,double> F(file, "oopr",
	tensor_array<3,double> F(size);

	size_t N = basis.shells().size();

	typedef BOOST_TYPEOF(F) FType;

	namespace phoenix = boost::phoenix;

	BOOST_AUTO(i, phoenix::local_names::_i);
	BOOST_AUTO(T2, phoenix::local_names::_c);
	BOOST_AUTO(T3, phoenix::local_names::_d);

					    
	boost::progress_display display(basis.size());


	boost::phoenix::function<Integral::Function>
	    integral(Integral::function(basis, screening));
	// boost::phoenix::function<transform::apply> transform;



	utility::timer t;


	for (size_t j = 0; j < N; ++j) {
	    BOOST_AUTO(R, basis.shell(j));

	    BOOST_AUTO(array, phoenix::array<size_t>());

	    using namespace phoenix;

	    {
		BOOST_AUTO(T2, arg_names::_1);
		BOOST_AUTO(P, arg_names::_2);
		T_(C2, C4, Integral::function(basis, screening)(j),
		   phoenix::put(phoenix::ref(F), T2,
				array(0, phoenix::start(P), 0),
				array(nij, phoenix::stop(P), R.size())));
	    }

	    for (size_t i = 0; i < R.size(); ++i) {
	    	size_t start[] = {0, 0, i };
	    	size_t stop[] = { nij, basis.size(), i+1 };
	    	Transform::Matrix T2(nij, basis.size());
	    	F.get(&T2.data()[0], start, stop);
	    	BOOST_AUTO(T3, arg_names::_1);
	    	T_(C1, ublas::trans(T2),
	    	   phoenix::put(ref(T), &phoenix::ublas::data(T3)[0],
	    			array(0, 0, R.start() + i),
	    			array(C1.size1(), nij, R.start() + i+1)));
	    }

	    display += R.size();
	}

	std::cout << "T(p,q,r,s) -> T(a,i,j,r): " << t << std::endl;
	//std::cout << boost::utility::global_profiler() << std::endl;
	
    }

    {
	namespace phoenix = boost::phoenix;
	// BOOST_AUTO(i, phoenix::local_names::_a);
	// BOOST_AUTO(j, phoenix::local_names::_b);
	// BOOST_AUTO(T3, phoenix::local_names::_c);
	// BOOST_AUTO(T4, phoenix::local_names::_d);
	// BOOST_AUTO(E_, phoenix::local_names::_e);

    	// typedef BOOST_TYPEOF(F) FType;
    	size_t N = basis.size();
    	// BOOST_AUTO(start, phoenix::array<size_t>()(i, j, 0, 0));
    	// BOOST_AUTO(stop, phoenix::array<size_t>()(i+1, j+1, N, N));
	// using namespace phoenix;
	boost::progress_display progress(nij);
	phoenix::function<increment> increment;

	using detail::energy;

	phoenix::iterator_function<vector_iterator> eo(e + active.start());
	vector_iterator ev(e + virtuals.start());


	// boost::phoenix::function<transform::apply> transform;
	BOOST_AUTO(array, phoenix::array<size_t>());

	utility::timer t;

	size_t  na = C1.size1();

	struct triangular_number {
	    triangular_number(int n) {
      	        index2_ = int((1 + sqrt(8*n +1))/2 - 1);
		index1_ = n - (index2_*index2_ + index2_)/2;
	    }
	    int index1() const { return index1_; }
	    int index2() const { return index2_; }
	private:
	    int index1_, index2_;
	};

	// size_t no = active.size();
	size_t m = (nij*N*basis.max().size());
	m = std::max<size_t>(m/(na*N), 16);
	Transform::Matrix T3(na, N*m);
	for (size_t ij = 0; ij < nij; ij += m) {
	    size_t n = std::min(nij - ij, m);

	    size_t size[] = { na, n, N, 1 };
	    T_.array1().resize(size);

	    size_t start[] = { 0, ij, 0 }, stop[] = { na, ij+n, N };
	    T.get(T_.array1().data(), start, stop);

	    BOOST_AUTO(i, phoenix::arg_names::_1);
	    BOOST_AUTO(j, phoenix::arg_names::_2);
	    BOOST_AUTO(T4, phoenix::arg_names::_3);
	    BOOST_AUTO(E_, phoenix::local_names::_a);
	    BOOST_AUTO(scale, (2-(i==j)));
	    T_(ij, ij+n, T_.array1(), ublas::trans(C3),
	       ( phoenix::let(E_ = double(0))
		 [ E_ = scale*energy(T4, -(eo(i)+eo(j)), ev, ev),
		   phoenix::thread::critical[phoenix::ref(E) += E_],
		   phoenix::thread::critical[increment(phoenix::ref(progress))]
		 ]));

	}

	std::cout << "T(a,i,j,r) -> T(a,b,i,j): " << t << std::endl;
    }

    std::cout << "write/read times: "
	      << T.counters().write() << "/"
	      << T.counters().read() << std::endl;
	
    //std::cout << "E(mp2): "  << E << std::endl;
    return E;

}



