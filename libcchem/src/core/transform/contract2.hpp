#ifndef TRANSFORM_CONTRACT2_HPP
#define TRANSFORM_CONTRACT2_HPP

#include "multi_array/multi_array.hpp"

#include "core/transform/transform.hpp"
#include "blas.hpp"
#include "core/transform/contract1.hpp"
#include "core/transform/packing.hpp"

#include "core/integral/generator.hpp"
#include "core/submatrix.hpp"
#include "foreach.hpp"

#include <algorithm>

#include <boost/ref.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/mpl/has_xxx.hpp>

#include "multi_array/algorithm.hpp"
#ifdef HAVE_GPU_TRANSFORM
#include <boost/numeric/bindings/cublas/cublas.hpp>
#include "multi_array/gpu/algorithm.hpp"
#endif

#include "boost/utility/profiler.hpp"

namespace transform {
namespace detail {


    template<class E>
    void contract2(Transform<true>::Array1 &T1,
		   const ublas::matrix_expression<E> &C2,
		   Transform<true>::Array2 &T2,
		   Transform<true>::Matrix &T2_) {

	BOOST_PROFILE_FUNCTION(); 

	size_t m = C2().size1();
	size_t k = C2().size2();
	assert(T2.size<0>() == (m*m + m)/2);
	assert(T1.size<0>() == m);
	assert(!(T1.size<1>() < k));
	T2_.resize(m, m, false);
	// std::cout << m << " " << C2().size1() << " " << k << std::endl;

	using ublas::make_matrix;
	using ublas::upper;
	using ublas::column_major;
	typedef ublas::symmetric_adaptor<Transform<true>::Matrix, upper> S;

	for (size_t d = 0; d < T2.size<2>(); ++d) {
	    for (size_t c = 0; c < T2.size<1>(); ++c) {
		BOOST_AUTO(T1_, make_matrix<column_major>(m, k, T1[d][c].origin()));
		blas::gemm(1, T1_, trans(C2()), 0, T2_);
		BOOST_AUTO(t2, (make_matrix<upper, column_major>(m, T2[d][c].origin() )));
		ublas::noalias(t2) += S(T2_);
	    }
	}
    }

#ifdef HAVE_GPU_TRANSFORM
    void contract2(Transform<true>::Gpu::Array1 &T1,
		   const Transform<true>::Gpu::Matrix &C2,
		   Transform<true>::Array2 &T2,
		   Transform<true>::Gpu::Matrix &G2_,
		   Transform<true>::Matrix &T2_) {

	BOOST_PROFILE_FUNCTION();

	size_t m = C2.size1();
	size_t k = C2.size2();
	assert(T2.size<0>() == (m*m + m)/2);
	assert(T1.size<0>() == m);
	    
	assert(!(T1.size<1>() < k));
	G2_.resize(m, m, false);
	T2_.resize(m, m, false);
	// std::cout << m << " " << C2().size1() << " " << k << std::endl;

	// typedef ublas::column_major layout;
	// typedef ublas::matrix_adapter<double, layout> matrix_adapter;
	// typedef ublas::symmetric_matrix_adapter<double, ublas::upper, layout>
	//     symmetric_matrix_adapter;
	// typedef ublas::symmetric_adaptor<Transform<true>::Matrix, ublas::upper>
	//     symmetric_adaptor;

	for (size_t d = 0; d < T2.size<2>(); ++d) {
	    for (size_t c = 0; c < T2.size<1>(); ++c) {
		namespace cublas = boost::numeric::bindings::cublas;

		// {
		//     BOOST_PROFILE_FUNCTION("*gemm");
		//     cublas::matrix_adaptor<double, layout>
		// 	t1(m, k, T1[d][c].origin());
		//     cublas::gemm(1, t1, cublas::trans(C2), 0, G2_);
		// }
		// {
		//     BOOST_PROFILE_FUNCTION("*host");
		//     cublas::host(T2_) = G2_;
		// }
		// symmetric_matrix_adapter t2(m, T2[d][c].origin());
		// ublas::noalias(t2) += symmetric_adaptor(T2_);
	    }	        
	}
    }
#endif

    BOOST_MPL_HAS_XXX_TRAIT_DEF(gpu_tag)

    template<int Order, class M, class A>
    void copy(size_t start, const size_t (&dims)[4], const M &input, A &output) {
	BOOST_PROFILE_FUNCTION();

	multi_array_view<const double,4> I(dims, &input.begin2().begin()[0]);

	using boost::multi_array_types::index_range;		
	BOOST_AUTO(range, boost::indices
		   [index_range(0, output.shape()[0])]
		   [index_range(0, output.shape()[1])]
		   [index_range(start, output.shape()[2])]
		   [index_range(0, output.shape()[3])]);
	multi_array_view<double,4> O(output[range]);

#ifdef HAVE_GPU_TRANSFORM
	if (has_gpu_tag<A>::value) {		
	    //multi_array::gpu::algorithm::copy<Order>(I, O)();
	    return;
	}
#endif
	multi_array::algorithm::copy<Order>(I, O)();
    }
    

} // namespace detail


    /// T2(k,a,b) = C2(j,b)*C1(i,a)*Q(k,i,j)
    template<bool Tri>
    void contract2(const Matrix &C1, const Matrix &C2,
		   integral::generator<2> generator,
		   Transform<Tri> &T) {
	BOOST_PROFILE_FUNCTION();

	BOOST_VERIFY((T.order().switch12() || T.order().switch13()));

	typedef integral::generator<2>::Shell Shell;
	const Shell& C = generator.get<0>();
	const Shell& D = generator.get<1>();
	size_t max = generator.max_shell().size();
	int ni = C1.size1();//, nj = C2.size1();

	T.matrix1().resize(ni, max*C.size()*D.size());
	size_t size[] = { ni, generator.basis().size(), C.size(), D.size() };

	T.array1().resize(size);
	BOOST_AUTO(G, T.gpu());
#ifdef HAVE_GPU_TRANSFORM
	if (G) G->array1().resize(size);
#endif
	// T.resize1(size);

	struct {
	    void contract(const Matrix &C2, Transform<Tri> &T) {
		if (!size()) return;
		BOOST_AUTO(&C_, T.matrix3());
		C_.resize(C2.size1(), C2.size2(), false);
		size_t K = 
		index_.pack(&std::copy<const double*,double*>,
				       C2.size1(),
				       C2.data().begin(),
				       C_.data().begin());
		
#ifdef HAVE_GPU_TRANSFORM
		BOOST_AUTO(G, T.gpu());
		if (G) {
		    G->matrix3() = subrange(C_, 0, C_.size1(), 0, K);
		    detail::contract2(G->array1(),
		    		      G->matrix3(),//subrange(C_, 0, C_.size1(), 0, K),
		    		      T.array2(),
		    		      G->matrix2(),
		    		      T.matrix2());
		}
		else
#else
#endif
		{
		    detail::contract2(T.array1(),
				      subrange(C_, 0, C_.size1(), 0, K),
				      T.array2(),
				      T.matrix2());
		}
		index_.data.clear();
	    }
	    size_t size() const { return index_.size(); }
	    void push(const Basis::Shell &shell) {
		index_.push(shell);
	    }
	    void push(const std::vector<size_t> &count, const Basis::Block &block) {
		BOOST_AUTO(shells, block.begin());
		foreach(size_t n, count) {
		    if (n > 0) index_.push(*shells);
		    ++shells;
		}
	    }
	private:
	    detail::packing::index index_;
	} block;


	typedef Basis::Block Block;
	BOOST_AUTO(const &basis, generator.basis());
	//BOOST_AUTO(const &shells, basis.shells());

	//GPU_ARGUMENT(detail::Gpu1 gpu1(T.gpu().integrals()), {});

	foreach (const Block& B, basis.blocks()) {

	    integral::block_generator<3>::array QRS =
	    	{{ B, Block(&C), Block(&D) }};

	    T.order().change(QRS);

	    const Block& Q = QRS[0];
	    const Block& R = QRS[1];
	    const Block& S = QRS[2];
	    size_t N = (Q.functions().size()*R.functions().size()*S.functions().size());

	    T.matrix1().resize(ni, N, false);
	    T.matrix2().resize(N, basis.size(), false);
	    T.matrix3().resize(C1.size1(), C1.size2(), false);
	    using detail::host_gpu_tuple;
	    host_gpu_tuple<const Matrix&, const double*> C_(C1);
	    host_gpu_tuple<Matrix&, double*> T1(T.matrix1());
	    host_gpu_tuple<double*> T0(T.matrix2().data().begin());
	    host_gpu_tuple<Matrix&, double*> C0(T.matrix3());

#ifdef HAVE_GPU_TRANSFORM
	    struct matrix {
		typedef typename Transform<Tri>::Gpu Gpu;
		static size_t size(const Matrix &m) { return m.size1()*m.size2(); }
		static void resize(typename Transform<Tri>::Gpu::Matrix &G,
				   const Matrix &m) {
		    G.resize(m.size1(), m.size2(), false);
		}
	    };

	    if (G) {
		matrix::resize(G->matrix1(), T.matrix1());
		matrix::resize(G->matrix2(), T.matrix2());
		matrix::resize(G->matrix3(), T.matrix3());

		C_.device = G->C1().data().begin();
		T1.device = G->matrix1().data().begin();
		T0.device = G->matrix2().data().begin();
		C0.device = G->matrix3().data().begin();
	    }
#endif

	    integral::block_generator<3> generator3(generator, QRS);
	    BOOST_AUTO(const &index,
		       detail::contract1(generator3, C_, T1, T0, C0, G));
	    if (index.empty()) continue;

	    for (size_t i = 0; i < index.size(); ++i) {

		BOOST_PROFILE_FUNCTION("*copy");
		size_t dims[] = { B.shell().size(),
				  C.size(), D.size(),
				  T.matrix1().size1() };
		T.order().change(dims[0], dims[1], dims[2]);

		size_t N = B.shell().size()*C.size()*D.size();	

#define TRANSFORM_COPY_HOST(MASK) {					\
		    using boost::numeric::ublas::subrange;		\
		    BOOST_AUTO(const &T1,				\
			       subrange(T.matrix1(),			\
					0, dims[3], i*N, (i+1)*N));	\
		    detail::copy<(MASK)>(block.size(), dims,		\
					 T1, T.array1());		\
		}

#ifdef HAVE_GPU_TRANSFORM
#define TRANSFORM_COPY(MASK) {						\
		    if (G) {						\
			using boost::numeric::bindings::cublas::subrange; \
			BOOST_AUTO(const &T1,				\
			    subrange(G->matrix1(),			\
				     0, dims[3], i*N, (i+1)*N));	\
			detail::copy<(MASK)>(block.size(), dims,	\
					     T1, G->array1());		\
		    }							\
		    else TRANSFORM_COPY_HOST(MASK);			\
		}
#else
#define TRANSFORM_COPY(MASK) TRANSFORM_COPY_HOST(MASK)
#endif
	    
		if (T.order().switch12()) { TRANSFORM_COPY(0231); }
		if (T.order().switch13()) { TRANSFORM_COPY(0123); }
		block.push(B.begin()[index.at(i)]);
	    }
	    //block.contract(C2, T);
	}
	block.contract(C2, T);


    }

} // namespace transform


#endif // TRANSFORM_CONTRACT2_HPP
