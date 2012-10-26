#ifndef ARRAY_HDF5_HPP
#define ARRAY_HDF5_HPP

#include <hdf5.h>

#include "exception.hpp"
#include "array/array.hpp"

#include <string>
#include <vector>
#include <fstream>

#include <boost/array.hpp>
#include <boost/multi_array/index_range.hpp>

struct Array<>::HDF5 : Array<>::Factory {
    struct File;
    struct Implementation;
    typedef Array<>::Type Type;

    Array<>::Implementation*
    allocate(Type t, size_t N, const size_t *dims, const char *name) const;

    std::string str() const { return "HDF5"; }

    HDF5(const char *name);
    HDF5(const HDF5 &f);
    ~HDF5();
    hid_t file();
    hid_t dcpl() const { return this->dcpl_; } 
    bool parallel() const { return this->parallel_; }

    void set_chunk(size_t N, const size_t *chunk) {
	hsize_t hchunk[N];
	std::copy(chunk, chunk+N, hchunk);
	std::reverse(hchunk, hchunk+N);
	//std::cout << chunk[0] << " " << chunk[1] << std::endl;
	CCHEM_ASSERT(H5Pset_chunk(dcpl_, N, hchunk) > -1);
    }

    void set_deflate(int n) {
	CCHEM_ASSERT(H5Pset_deflate(dcpl_, n) > -1);
    }

#ifdef H5_HAVE_PARALLEL
    void set_plist(MPI_Comm comm) {
	CCHEM_ASSERT(H5Pset_fapl_mpio(plist_, comm, MPI_INFO_NULL) > -1);;
	parallel_ = true;
    }
#endif

    static hid_t type(Type t) {
	if (t.id == Array<>::Type::DOUBLE) return H5T_NATIVE_DOUBLE;
	if (t.id == Array<>::Type::INT) return H5T_NATIVE_INT;
	throw Error(__FILE__ ": unsupported HDF5 type");
    }

private:
    std::string name_;
    bool parallel_;
    hid_t file_, plist_, dcpl_;
    void operator=(const HDF5 &f);
};


struct Array<>::HDF5::Implementation : Array<>::Implementation {

    typedef Array<>::Type Type;
    const size_t N;

    Implementation(Type t, size_t N, const size_t *dims, const char *name,
		   const Array<>::HDF5 &f)
	: N(N), type_(Array<>::HDF5::type(t)), f_(f)
    {
	std::vector<hsize_t> hdims;
	for (int i = N-1; i >= 0; --i) {
	    hdims.push_back(dims[i]);
	}
	this->name_ = name;
	this->dataspace_ = H5Screate_simple(N, &hdims[0], NULL);
	this->dataset_ =
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	    H5Dcreate(f_.file(), name, type_, dataspace_, f_.dcpl());
#else
	H5Dcreate1(f_.file(), name, type_, dataspace_, f_.dcpl());
#endif
    }

    ~Implementation() {
	H5Dclose(dataset_);
	H5Sclose(dataspace_);
    }

    hid_t data() { return dataset_; }
    hid_t data() const { return dataset_; }

    bool parallel() const { return this->f_.parallel(); }

    void put(const void *buffer, const size_t *start, const size_t *stop) {
	apply(&H5Dwrite, *this, H5P_DEFAULT, buffer, start, stop);
    }

    void get(void *buffer, const size_t *start, const size_t *stop) const {
	apply(&H5Dread, *this, H5P_DEFAULT, buffer, start, stop);
    }

#ifdef H5_HAVE_PARALLEL
    void get_all(void *buffer, const size_t *start, const size_t *stop) const {
	hid_t plist = H5Pcreate(H5P_DATASET_XFER);
	CCHEM_ASSERT(H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE) > -1);
	apply(&H5Dread, *this, plist, buffer, start, stop);
	H5Pclose(plist);
    }
#endif

    void wait() {
	H5Fflush(dataset_, H5F_SCOPE_GLOBAL);
    }

private:
    std::string name_;
    hid_t type_;
    HDF5 f_;
    hid_t dataspace_, dataset_;

private:

    template<class F, class T>
    static void apply(F f,
		      const Implementation &a, hid_t plist,
		      T *buffer,
		      const size_t *start, const size_t *stop) {
	int N = a.N;

	// // mpio layer uses in to count bytes, split into 2GB chunks
	// if (a.parallel()) {
	//     size_t t = H5Tget_size(a.type_);
	//     size_t max = ((1LU << 31)-1)/t;
	//     size_t size = count(N, start, stop);
	//     if (size > max) {
	// 	if (count(N-1, start+1, stop+1) > 1) {
	// 	    CCHEM_ERROR("not implemented");
	// 	}
	// 	int dim = 0;
	// 	size_t begin[N], end[N];
	// 	std::copy(start, start+N, begin);
	// 	std::copy(stop, stop+N, end);
	// 	for (size_t i = start[dim]; i < stop[dim]; i += max) {
	// 	    begin[dim] = i;
	// 	    end[dim] = std::min(i+max, stop[dim]);
	// 	    apply(f, a, plist,
	// 		  (char*)buffer+t*(i-start[dim]), begin, end);
	// 	}
	// 	return;
	//     }
	// }

#ifndef H5_HAVE_THREADSAFE
#pragma omp critical(cchem_array_hdf5)
	//#warning "HDF5 not thread safe: !defined(H5_HAVE_THREADSAFE)"
#endif
	{
	    hid_t fspace = H5Dget_space(a.data());
	    hsize_t dims[N];
	    select(fspace, N, start, stop, dims);
	    hid_t mspace = H5Screate_simple(N, dims, NULL);
	    CCHEM_ASSERT
		(H5Sselect_all(mspace) >= 0);
	    CCHEM_ASSERT
		(f(a.data(), a.type_, mspace, fspace, plist, buffer) >= 0);
	    // {
	    // 	H5D_mpio_actual_io_mode_t mode;
	    // 	H5Pget_mpio_actual_io_mode(plist, &mode);
	    // 	std::cout << mode << " "
	    // 	    //<< H5D_MPIO_NO_COLLECTIVE_IO << " "
	    // 		  << H5D_MPIO_CONTIGUOUS_COLLECTIVE << " "
	    // 		  << H5D_MPIO_CHUNK_INDEPENDENT << " "
	    // 		  << H5D_MPIO_CHUNK_COLLECTIVE << " "
	    // 		  << H5D_MPIO_CHUNK_MIXED << " "
	    // 		  << std::endl;
	    // }
	    H5Sclose(mspace);
	    H5Sclose(fspace);
	}
    }

    static size_t count(int N, const size_t *start, const size_t *stop) {
	size_t size = N ? 1 : 0;
	for (int i = 0; i < N; ++i) {
	    size *= (stop[i] - start[i]);
	}
	return size;
    }

    static
    size_t select(hid_t space, int N, const size_t *start, const size_t *stop,
		  hsize_t *dims) {
	hsize_t offset[N];
	size_t size = 1;
	for (int i = 0, j = N-1; i < N; ++i, --j) {
	    offset[i] = start[j];
	    dims[i] = stop[j] - start[j];
	    size *= dims[i];
	}
	CCHEM_ASSERT
	    (H5Sselect_hyperslab(space, H5S_SELECT_SET,
				 offset, NULL, dims, NULL)
	     >= 0);
	return size;
    }    

};

inline Array<>::Implementation*
Array<>::HDF5::allocate(Type t, size_t N, const size_t *dims, const char *name) const {
    return new Implementation(t, N, dims, name, *this);
}


#endif // ARRAY_HDF5_HPP
