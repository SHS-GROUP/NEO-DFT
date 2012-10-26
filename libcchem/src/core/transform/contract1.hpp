#ifndef TRANSFORM_CONTRACT1_HPP
#define TRANSFORM_CONTRACT1_HPP

#include "core/transform/transform.hpp"
#include "blas.hpp"
#include "core/transform/packing.hpp"
#include "core/submatrix.hpp"
#include "foreach.hpp"

#include "core/integral.hpp"
#include <algorithm>
#include <boost/ref.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/numeric/ublas/adaptors.hpp>

#ifdef HAVE_GPU_TRANSFORM
#include <boost/numeric/bindings/cublas/cublas.hpp>
#include "core/integral/gpu.hpp"
#include "externals/gpu/gpu.hpp"
// #include "externals/gpu/blas.hpp"
#define GPU_ARGUMENT(arg, ...) arg
#else
#define GPU_ARGUMENT(arg, ...) __VA_ARGS__
#endif

#define BOOST_PROFILE_ENABLE
#include "boost/utility/profiler.hpp"

namespace transform {
namespace detail {

    template<class H, class D = H>
    struct host_gpu_tuple {
	typedef H host_type;
	typedef D device_type;
	host_gpu_tuple()
	    : host(), device() {}
	explicit host_gpu_tuple(H host)
	    : host(host), device() {}
	host_gpu_tuple(H host, D device)
	    : host(host), device(device) {}
	host_type host;
	device_type device;
    };
    
    packing::index2 integrals1(integral::block_generator<3> generator,
			       host_gpu_tuple<double*> G,
			       Gpu *gpu = 0) {

	typedef integral::block_generator<3>::Shell Shell;
	typedef integral::block_generator<3>::Block Block;

	host_gpu_tuple<packing::index2> index2;

	BOOST_AUTO(const &basis, generator.basis());
	BOOST_AUTO(const &centers, basis.centers());

	BOOST_PROFILE_FUNCTION();

	const Basis::Block &Q = generator.get<0>();
	const Basis::Block &R = generator.get<1>();
	const Basis::Block &S = generator.get<2>();

	size_t size = S.size()*R.size()*Q.size();
	index2.host.data().resize(size);
	index2.device.data().resize(size);
	
	typedef boost::array<Integral::index,4> quartet;
	struct {
	    struct {
		std::vector<double*> data;
		std::vector<double*> block0, block;
		std::vector<size_t> count;
		std::vector<quartet> quartets;
		boost::ptr_vector<Integral> eri;
	    } host;
	    struct {
		std::vector<double*> data;
		std::vector<double*> block;
		std::vector<size_t> count;
		std::vector<quartet> quartets;
		typedef GPU_ARGUMENT(Integral::Gpu::Quartet, bool) eri_type;
		std::vector<eri_type> eri;
	    } device;
	} integral;

	{
	    size_t size = (S.shell().size()*R.shell().size()*Q.shell().size());
	    double *h = G.host;
	    double *d = G.device;
	    for (size_t s = 0; s < S.size(); ++s) {
		for (size_t r = 0; r < R.size(); ++r) {
		    for (size_t q = 0; q < Q.size(); ++q) {
			integral.host.block.push_back(h);
			h += size*basis.size();
			integral.device.block.push_back(d);
			d += size*basis.size();
		    }
		}
	    }
	    integral.host.block0 = integral.host.block;
	}

	// std::cout << bool(gpu) << std::endl;

	foreach (const Basis::Block &P, basis.blocks()) {
	    packing::index2 block_index;

	    size_t size = (P.shell().size()*Q.shell().size()*
			   R.shell().size()*S.shell().size());
	    
	    bool use_gpu = false;
#ifdef HAVE_GPU_TRANSFORM
	    if (gpu) {
		BOOST_AUTO(eri, gpu->eri()(S, R, Q, P));
		use_gpu = true;//eri;
		if (use_gpu) integral.device.eri.push_back(eri);
	    }
#endif
	    if (!use_gpu)
		integral.host.eri.push_back(new Integral(S, R, Q, P));

	    BOOST_AUTO(block, ((use_gpu) ?
			       integral.device.block.begin() :
			       integral.host.block.begin()));
	    BOOST_AUTO(&quartets, ((use_gpu) ?
				  integral.device.quartets :
				  integral.host.quartets));
	    BOOST_AUTO(&pointers, ((use_gpu) ?
				   integral.device.data :
				   integral.host.data));

	    size_t count = 0;
	    foreach (const Shell& s, S) {
		foreach (const Shell& r, R) {
		    foreach (const Shell& q, Q) {
			block_index.push2();
			int l = basis.index(s);
			int k = basis.index(r);
			int j = basis.index(q);

			foreach (const Shell& p, P) {
			    if (!generator.test(p,q,r,s)) continue;
			    block_index.back2().push(p);
			    int i = basis.index(p);
			    quartet v = {{ l, k, j, i }};
			    quartets.push_back(v);
			    pointers.push_back(*block);
			    *block += size;
			    ++count;
			}
			++block;
		    }
		}
	    }
	    if (!count) {
		if (use_gpu) integral.device.eri.pop_back();
		else integral.host.eri.pop_back();
		continue;
	    }
	    if (use_gpu) {
		integral.device.count.push_back(count);
		index2.device.extend(block_index);
	    }
	    else if (!gpu) {
		integral.host.count.push_back(count);
		index2.host.extend(block_index);
	    }
	}

#ifdef HAVE_GPU_TRANSFORM
	if (gpu) {
	    gpu->quartets.assign(integral.device.quartets);
	    gpu->data.assign(integral.device.data);
	    for (size_t i = 0, block = 0; i < integral.device.count.size(); ++i) {
		size_t size = integral.device.count.at(i);
		BOOST_AUTO(&eri, integral.device.eri.at(i));
		//eri(size, gpu->quartets.begin() + block, gpu->data.begin() + block);
		block += size;
	    }
	}
#endif

	for (size_t i = 0; i < integral.host.block.size(); ++i) {
	    std::fill(integral.host.block0.at(i), integral.host.block.at(i), 0);
	}

	if (!gpu)
	for (size_t i = 0, block = 0; i < integral.host.count.size(); ++i) {
	    size_t size = integral.host.count.at(i);
	    BOOST_AUTO(&eri, integral.host.eri.at(i));
	    eri(centers, size,
	    	&integral.host.quartets.at(block),
	    	&integral.host.data.at(block));
	    block += size;
	}

#ifdef HAVE_GPU_TRANSFORM
	if (gpu) {
	    for (size_t i = 0; i < integral.host.block.size(); ++i) {
		typedef boost::cuda::device_ptr<double> device_ptr;
		// boost::cuda::copy(integral.host.block0.at(i),
		// 		  integral.host.block.at(i),
		// 		  device_ptr::reinterpret(integral.device.block.at(i)));
	    }
	    index2.device.extend(index2.host);
	    return index2.device;
	}
#endif

	return index2.host;
    }


    std::vector<size_t>
    contract1(integral::block_generator<3> generator,
	      host_gpu_tuple<const Matrix&, const double*> C,
	      host_gpu_tuple<Matrix&, double*> T,
	      host_gpu_tuple<double*> T0,
	      host_gpu_tuple<Matrix&, double*> C0,
	      Gpu *gpu = 0) {

	BOOST_PROFILE_FUNCTION();

	typedef integral::block_generator<3>::Shell Shell;
	typedef integral::block_generator<3>::Block Block;

	BOOST_AUTO(const &basis, generator.basis());
	// BOOST_AUTO(const &shells, basis.shells());
	BOOST_AUTO(const &Q, generator.get<0>()); 
	BOOST_AUTO(const &R, generator.get<1>()); 
	BOOST_AUTO(const &S, generator.get<2>()); 

	size_t M = C.host.size1();
	size_t N = (Q.shell().size()*R.shell().size()*S.shell().size());

	{
	    // size_t size1 = (T0.host.shape()[0]*T0.host.shape()[1]*T0.host.shape()[2]);
	// assert(T0.host.size1() == N && T0.host.size2() >= generator.basis().size());
	}
	BOOST_VERIFY(C0.host.size1() == M && C0.host.size2() >= C.host.size2());
	BOOST_VERIFY(T.host.size1() >= M && T.host.size2() >= N);

	BOOST_AUTO(const &index2, integrals1(generator, T0, gpu));
	// {
	//     BOOST_PROFILE_FUNCTION("synchronize");
	//     gpu.synchronize();
	// }

	BOOST_PROFILE_FUNCTION(" - integrals");

	std::vector<size_t> index;
	for (size_t i = 0, j = 0; i < index2.data().size(); ++i) {
	    double *H = T0.host + i*N*basis.size();

	    BOOST_AUTO(const &index1, index2.data().at(i).data);
	    BOOST_VERIFY(index1.size()%2 == 0);
	    if (index1.empty()) continue;

	    if (gpu) {
#ifdef HAVE_GPU_TRANSFORM
		size_t K;
		{
		    BOOST_PROFILE_FUNCTION("*pack");
		    // BOOST_AUTO(index, packing::index::pairs(index1));
		    // foreach (const BOOST_TYPEOF(index.front()) &r, index) {
		    // 	K += r.second - r.first;
		    // }
		    // gpu.pack(index, M, C.device, M, C0.device, M);
		    BOOST_AUTO(copy_device_to_device,
			       (&gpu::copy<gpu::device_to_device,
				const double*, double*>));
		    K = pack(copy_device_to_device,
			     M, C.device, C0.device, index1);
		}

		namespace cublas =  boost::numeric::bindings::cublas;
		// double *D = T0.device + i*N*basis.size();
		typedef cublas::matrix_adaptor<double> A;

		A C_(M, K, C0.device);
		A T_(N, K, T0.device + i*N*basis.size());
		A Q_(M, N, T.device + j*N*M);
		cublas::gemm(1, C_, cublas::trans(T_), 0, Q_);
		// // cublas::gemm(blas::no_transpose, blas::transpose,
		// //  	   M, N, K, C0.device, M, D, N, T.device + j*N*M, M);
#endif
	    }
	    else {
		size_t K = pack(&std::copy<const double*,double*>,
				M, C.host.data().begin(), C0.host.data().begin(),
				index1);

		using ublas::subrange;
		using ublas::trans;

		BOOST_AUTO(T_, subrange(T.host, 0, M, j*N, (j+1)*N));
		BOOST_AUTO(I_, ublas::make_matrix<ublas::column_major>(N, K, H));

		blas::gemm(1, subrange(C0.host, 0, M, 0, K), trans(I_), 0, T_);
	    }

	    index.push_back(i);
	    ++j;

	} // foreach
	return index;
    }

} // namespace detail
} // namespace transform


#endif // TRANSFORM_CONTRACT1_HPP
