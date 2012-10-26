#include "rysq-cuda.hpp"
#include "assert.hpp"
#include "cuda/detail.hpp"
#include <boost/utility/profiler.hpp>

namespace rysq {
namespace cuda {

    Eri::Eri(Shell &a, Shell &b, Shell &c, Shell &d,
	     const Context &context) {
	bool t0 = (!a.L() && b.L());
	bool t1 = (!c.L() && d.L());
	bool t01 = (!(a.L() + b.L()) && (c.L() + d.L()));
	Transpose transpose(t0, t1, t01);
	detail::Shell::Quartet quartet(a, b, c, d);
	try {
	    this->kernel_.reset(new cuda::detail::Eri(quartet, context, transpose));
	}
	catch (...) {
	    this->kernel_.reset();
	}
    }

    Eri::~Eri() {}

    void Eri::operator()(const Centers &centers,
			 const Quartets &quartets,
			 double *I,
			 const Parameters &parameters,
			 Stream s) {
	RYSQ_ASSERT(kernel_.get());
	(*kernel_)(detail::Centers(centers),
		   detail::Quartets(quartets),
		   detail::Eri::List(I, kernel_->quartet().size()),
		   parameters,
		   static_cast<cudaStream_t>(s.data()));
    }

} // namespace cuda
} // namespace rysq

// Eri::~Eri() {}

// void Eri::operator()(const Centers &centers,
// 				 const Quartets &quartets,
// 				 const std::vector<double*> &integrals,
// 				 const Parameters &parameters) {

//     BOOST_PROFILE_FUNCTION();

// 	// for (int i = 0; i < quartets.size(); ++i)
// 	//     std::cout << quartets[i] << "\n";
//     // host::device_ptr<double> eri =
//     // 	host::malloc<double>(quartets.size()* kernel_->quartet().size());
//     // cuda::vector<Center> c(centers.size());
//     // boost::cuda::vector<Int4> q(quartets.size());

//     // std::cout << stream << std::endl;

//     //cudaThreadSynchronize();
//     // c.assign(centers, stream);
//     // q.assign(quartets);//, stream_);
//     // std::cout << quartets.size() << " " << quartets.begin()[0] << std::endl;
//     stream_.synchronize();
//     quartets_.assign(quartets);
//     integrals_.assign(integrals);

//     using boost::cuda::device_ptr;
    
//     // struct {
//     // 	device_ptr<const Int4> quartets;
//     // 	device_ptr<double*> integrals;
//     // } data = {
//     // 	device_ptr<const Int4>::reinterpret(quartets_.begin()),
//     // 	device_ptr<double*>::reinterpret(&integrals_[0]) };

//     // boost::cuda::device_ptr<const Int4>  = 
//     BOOST_VERIFY(quartets.size() <= integrals.size());
//     // 	boost::cuda::mapped_ptr<const Int4>::reinterpret(&quartets[0]);

//     (*kernel_)(detail::Centers(centers.begin()),
// 	       detail::Quartets(quartets_), detail::Eri::List(integrals_.begin()),
// 	       //detail::Quartets(quartets_.begin(), quartets_.size()),
// 	       parameters, static_cast<cudaStream_t>(stream_.data()));
    
// }


// void Eri::operator()(const Centers &centers, size_t size,
// 				 boost::cuda::device_ptr<const Int4> quartets,
// 				 boost::cuda::device_ptr<double*> I,
// 				 const Parameters &parameters) {
//     //BOOST_PROFILE_FUNCTION();
//     //std::cout << size << std::endl;
//     (*kernel_)(detail::Centers(centers.begin()),
// 	       detail::Quartets(quartets, size),
// 	       detail::Eri::List(I),
// 	       parameters, static_cast<cudaStream_t>(stream_.data()));
// }

// // void Eri::operator()(const  Centers &centers,
// // 				  const  Quartets &quartets,
// // 				 vector <double> &eri,
// // 				 const rysq::Parameters &parameters) {
// // 	// for (int i = 0; i < quartets.size(); ++i)
// // 	//     std::cout << quartets[i] << "\n";
// //     // host::device_ptr<double> eri =
// //     // 	host::malloc<double>(quartets.size()* kernel_->quartet().size());
// //     (*kernel_)(centers.data(), detail::Quartets(quartets.data(), quartets.size()),
// // 	       eri.data(), parameters);
// // }
