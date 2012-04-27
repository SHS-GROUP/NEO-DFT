#include "boost/cuda/runtime.hpp"
#include "boost/cuda/mapped_ptr.hpp"

#include <rysq/cuda.hpp>
#include "cuda/detail.hpp"
#include <boost/utility/profiler.hpp>

namespace host = boost::cuda;

rysq::cuda::Eri::Eri(const rysq::Quartet<rysq::Shell> &quartet,
		     const cuda::stream &stream)
    : stream_(stream)
{
    Transpose transpose;
    try {
	this->kernel_.reset(new cuda::detail::Eri(quartet, transpose));
    }
    catch (const std::exception&) { this->kernel_.reset(NULL); }
    ///std::cout << this->kernel_.get() << std::endl;
}

rysq::cuda::Eri::~Eri() {}

void rysq::cuda::Eri::operator()(const Centers &centers,
				 const Quartets &quartets,
				 const std::vector<double*> &integrals,
				 const Parameters &parameters) {

    BOOST_PROFILE_FUNCTION();

	// for (int i = 0; i < quartets.size(); ++i)
	//     std::cout << quartets[i] << "\n";
    // host::device_ptr<double> eri =
    // 	host::malloc<double>(quartets.size()* kernel_->quartet().size());
    // cuda::vector<Center> c(centers.size());
    // boost::cuda::vector<Int4> q(quartets.size());

    // std::cout << stream << std::endl;

    //cudaThreadSynchronize();
    // c.assign(centers, stream);
    // q.assign(quartets);//, stream_);
    // std::cout << quartets.size() << " " << quartets.begin()[0] << std::endl;
    stream_.synchronize();
    quartets_.assign(quartets);
    integrals_.assign(integrals);

    using boost::cuda::device_ptr;
    
    // struct {
    // 	device_ptr<const Int4> quartets;
    // 	device_ptr<double*> integrals;
    // } data = {
    // 	device_ptr<const Int4>::reinterpret(quartets_.begin()),
    // 	device_ptr<double*>::reinterpret(&integrals_[0]) };

    // boost::cuda::device_ptr<const Int4>  = 
    BOOST_VERIFY(quartets.size() <= integrals.size());
    // 	boost::cuda::mapped_ptr<const Int4>::reinterpret(&quartets[0]);



    (*kernel_)(detail::Centers(centers.begin()),
	       detail::Quartets(quartets_), detail::Eri::List(integrals_.begin()),
	       //detail::Quartets(quartets_.begin(), quartets_.size()),
	       parameters, stream_);
    
}


void rysq::cuda::Eri::operator()(const Centers &centers, size_t size,
				 boost::cuda::device_ptr<const Int4> quartets,
				 boost::cuda::device_ptr<double*> I,
				 const Parameters &parameters) {
    //BOOST_PROFILE_FUNCTION();
    //std::cout << size << std::endl;
    (*kernel_)(detail::Centers(centers.begin()),
	       detail::Quartets(quartets, size),
	       detail::Eri::List(I),
	       parameters, stream_);
}

// void rysq::cuda::Eri::operator()(const  Centers &centers,
// 				  const  Quartets &quartets,
// 				 rysq::cuda::vector <double> &eri,
// 				 const rysq::Parameters &parameters) {
// 	// for (int i = 0; i < quartets.size(); ++i)
// 	//     std::cout << quartets[i] << "\n";
//     // host::device_ptr<double> eri =
//     // 	host::malloc<double>(quartets.size()* kernel_->quartet().size());
//     (*kernel_)(centers.data(), detail::Quartets(quartets.data(), quartets.size()),
// 	       eri.data(), parameters);
// }
