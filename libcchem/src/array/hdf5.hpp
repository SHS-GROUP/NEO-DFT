#ifndef ARRAY_H5_HPP
#define ARRAY_H5_HPP

#include "array/array.hpp"

#include <H5Cpp.h>
#include <string>
#include <vector>
#include <fstream>

#include <boost/array.hpp>
#include <boost/multi_array/index_range.hpp>

#ifndef H5_HAVE_THREADSAFE
#warning "HDF5 not thread safe: !defined(H5_HAVE_THREADSAFE)"
#endif

struct Array<>::HDF5 : Array<>::Factory {
    struct File;
    struct Implementation;
    typedef Array<>::TYPE TYPE;

    Array<>::Implementation*
    allocate(TYPE t, size_t N, const size_t *dims, const char *name) const;

    typedef H5::H5File H5File;

    struct Properties : H5::DSetCreatPropList {
	template<size_t M>
	Properties& set_chunk(const size_t (&chunk)[M]) {
	    boost::array<hsize_t,M> chunk_;
	    std::copy(chunk, chunk+M, chunk_.begin());
 	    return set_chunk(chunk_);
	}
	template<size_t M>
	Properties& set_chunk(const boost::array<hsize_t,M> &chunk) {
	    this->setChunk(M, chunk.data());
	    return *this;
	}
    };

    HDF5(const char *file = "data.h5", const Properties &properties = Properties()) {
	this->file = get(file);
	this->properties = properties;
    }

    static H5::PredType type(TYPE t) {
	switch (t) {
	case Array<>::DOUBLE: return H5::PredType::NATIVE_DOUBLE;
	default: throw;
	}
	throw;
    }

protected:

    static H5File create(const char *name) {
	return H5File(name, H5F_ACC_TRUNC);
    }
    static H5File get(const char *name) {
	H5File file;
	if (!std::ifstream(name)) {
	    file = H5File(name, H5F_ACC_EXCL);
	}
	file.openFile(name, H5F_ACC_RDWR);
	return file;
    }
    static H5File open(const char *name) {
	H5File file;
	file.openFile(name, H5F_ACC_RDWR);
	return file;
    }

private:
    friend class Array<>::HDF5::Implementation;
    H5File file;
    Properties properties;
};


struct Array<>::HDF5::Implementation : Array<>::Implementation {

    typedef Array<>::TYPE TYPE;

    Implementation(TYPE t, size_t N, const size_t *dims, const char *name,
		   const HDF5 &f)
	: type_(Array<>::HDF5::type(t))
    {
	std::vector<hsize_t> hdims;
	for (int i = N-1; i >= 0; --i) {
	    hdims.push_back(dims[i]);
	}
	H5::FloatType datatype(this->type_);
	H5::DataSpace fspace = H5::DataSpace(hdims.size(), &hdims[0]);
	name_ = name;
	factory_ = f;
	dataset_ = factory_.file.createDataSet(name, datatype, fspace, f.properties);
    }

    const Array<>::Factory& factory() const { return factory_; }

    H5::DataSet& data() { return dataset_; }
    const H5::DataSet& data() const { return dataset_; }

    void put(const void *buffer, const size_t *start, const size_t *stop) {
	put(*this, buffer, start, stop);
    }
    void get(void *buffer, const size_t *start, const size_t *stop) const {
	get(*this, buffer, start, stop);
    }
    void synchronize() {
	dataset_.flush(H5F_SCOPE_GLOBAL);
    }

    void open() {
	dataset_ = factory_.file.openDataSet(name_);
    }

    void close() {
	dataset_.close();
    }

private:
    Array<>::HDF5 factory_;
    std::string name_;
    H5::DataSet dataset_;
    H5::PredType type_;

private:

    static
    void put(Implementation &a, const void *buffer,
	     const size_t *start, const size_t *stop) {
	H5::DataSpace fspace = a.data().getSpace();
	size_t size = select(fspace, start, stop);
	hsize_t mdims[] = { size };
	H5::DataSpace mspace(1, mdims);
	mspace.selectAll();
	a.data().write(buffer, a.type_, mspace, fspace);
    }

    static
    void get(const Implementation &a, void *buffer,
	     const size_t *start, const size_t *stop) {
	H5::DataSpace fspace = a.data().getSpace();
	size_t size = select(fspace, start, stop);
	hsize_t mdims[] = { size };
	H5::DataSpace mspace(1, mdims);
	mspace.selectAll();
	a.data().read(buffer, a.type_, mspace, fspace);
    }

    static
    size_t select(const H5::DataSpace &space,
		  const size_t *start, const size_t *stop) {
	size_t N = space.getSimpleExtentNdims();
	hsize_t fstart[N];
	hsize_t fstride[N]; // Stride of hyperslab
	hsize_t fcount[N];  // Block count
	hsize_t fblock[N];  // Block sizes
	size_t size = 1;
	for (size_t i = 0, j = N-1; i < N; ++i, --j) {
	    fstart[i] = start[j];
	    fcount[i] = stop[j] - start[j];
	    fstride[i] = 1;
	    fblock[i] = 1;
	    size *= fcount[i];
	}
	space.selectHyperslab(H5S_SELECT_SET, fcount, fstart, fstride, fblock);
	return size;
    }    

};

inline Array<>::Implementation*
Array<>::HDF5::allocate(TYPE t, size_t N, const size_t *dims, const char *name) const {
    return new Implementation(t, N, dims, name, *this);
}


#endif // ARRAY_H5_HPP
