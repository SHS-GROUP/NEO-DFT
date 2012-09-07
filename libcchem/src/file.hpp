#ifndef FILE_HPP
#define FILE_HPP

#include <string>
#include <fstream>

#include <hdf5.h>

#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/type_traits.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/at_key.hpp>

#include "exception.hpp"

struct File {

    struct Attribute;
    struct Dataset;
    struct Group;

    template<typename T>
    static hid_t h5_predtype() {
	using boost::fusion::make_map;
	BOOST_AUTO(pt,
		   (make_map<int,double>
		    (H5T_NATIVE_INT,
		     H5T_NATIVE_DOUBLE)));
	typedef typename boost::remove_const<T>::type U;
	return boost::fusion::at_key<U>(pt);
    }

    struct Attribute {
	hid_t& h5() { return handle_; }
	template<typename T>
	operator T() const { return read<T>(); }
	template<typename T>
	void operator=(const T &value) { write(value); }
	Attribute(const Attribute &a) {
	    handle_ = a.handle_;
	    H5Iinc_ref(handle_);
	}
	~Attribute() { H5Aclose(handle_); }
    private:
	template<typename T>
	T read() const {
	    T value;
	    H5Aread(handle_, h5_predtype<T>(), &value);
	    return value;
	}
	template<typename T>
	void write(const T &value) {
	    H5Awrite(handle_, h5_predtype<T>(), &value);
	}	
    private:
	friend class Dataset;
	explicit Attribute(const hid_t &data)
	    : handle_(data)
	{
	    assert(handle_ > 0);
	}
	void operator=(const Attribute&) {}
    private:
	hid_t handle_;
    };


    struct Dataset {
    	Dataset() { this->handle_ = 0; }
	Dataset(const Dataset &d) {
	    this->handle_ = d.handle_;
	    H5Iinc_ref(this->handle_);
	}
    	~Dataset() { this->close(); }

    	void close() const { H5Dclose(handle_); }
    	hid_t h5() { return handle_; }

    	template<typename T, size_t K>
    	static Dataset create(const std::string &name,
    			      const size_t *dims, hid_t id) {
    	    hsize_t hdims[K];
    	    std::reverse_copy(dims, dims+K, hdims);
	    hid_t fspace = H5Screate_simple(K, hdims, NULL);
	    return Dataset
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	    (H5Dcreate(id, name.c_str(), h5_predtype<T>(), fspace, H5P_DEFAULT));
#else
	    (H5Dcreate1(id, name.c_str(), h5_predtype<T>(), fspace, H5P_DEFAULT));
#endif
    	}
    	static Dataset open(const std::string &name, hid_t id) {
    	    return Dataset(H5Dopen1(id, name.c_str()));
    	}

    private:
    	hid_t handle_;
    	explicit Dataset(hid_t data) : handle_(data) {}
	void operator=(const Dataset&) {}

    public:

	template<typename T, typename S>
	void put(const T *buffer, const S &start, const S &stop) {
	    apply(&H5Dwrite, this->handle_, buffer, start, stop);
	}

	template<typename T, typename S>
	void get(T *buffer, const S &start, const S &stop) {
	    apply(&H5Dread, this->handle_, buffer, start, stop);
	}

    	template<typename T>
    	T get(int i0, int i1, int i2, int i3, int i4) const {
    	    T v;
    	    int start[5] = { i0, i1, i2, i3, i4 };
    	    int stop[5];
    	    for (int i = 0; i < 5; ++i) stop[i] = start[i] + 1;
    	    get<T>(&v, start, stop);
    	    return v;
    	}

    	template<typename T>
    	Attribute create_attribute(const std::string &name) {
	    hid_t space = H5Screate(H5S_SCALAR);
	    hid_t attr =
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	    H5Acreate(handle_, name.c_str(), h5_predtype<T>(), space, H5P_DEFAULT);
#else
	    H5Acreate1(handle_, name.c_str(), h5_predtype<T>(), space, H5P_DEFAULT);
#endif
	    H5Sclose(space);
	    return Attribute(attr);	    
    	}
    	Attribute attribute(const std::string &name) {
	    return Attribute(H5Aopen_name(handle_, name.c_str()));
    	}

    private:

	template<class F, typename T, typename S>
	static void apply(F f, hid_t dset, T *buffer,
			  const S &start, const S &stop) {
	    hid_t fspace = H5Dget_space(dset);
	    size_t size = select(fspace, start, stop);
	    hsize_t mdims[] = { size };
	    hid_t mspace = H5Screate_simple(1, mdims, NULL);
	    CCHEM_ASSERT
		(H5Sselect_all(mspace) >= 0);
	    CCHEM_ASSERT
		(f(dset, h5_predtype<T>(), mspace, fspace, H5P_DEFAULT, buffer)
	     >= 0);
	    H5Sclose(mspace);
	    H5Sclose(fspace);
	}

    	template<typename S>
    	static size_t select(hid_t space, const S &start, const S &stop) {
    	    size_t N = H5Sget_simple_extent_ndims(space);
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
	    CCHEM_ASSERT
		(H5Sselect_hyperslab(space, H5S_SELECT_SET,
				     fstart, fstride, fcount, fblock)
		 >= 0);
    	    return size;
    	}    

    };

    struct Group {
	Group(const Group& group) {
	    this->handle_ = group.handle_;
	    H5Iinc_ref(handle_);
	}
	~Group() { H5Gclose(handle_); }
	template<typename T, size_t N>
	Dataset create_dataset(const std::string &name, const size_t (&dims)[N]) {
	    return Dataset::create<T,N>(name, dims, handle_);
	}
	Dataset dataset(const std::string &name) {
	    return Dataset::open(name, handle_);
	}
    private:
	hid_t handle_;
	friend class File;
	explicit Group(hid_t handle) : handle_(handle) {}
	void operator=(const Group&) {}
    };

private:
    hid_t handle_;
    boost::ptr_map<std::string, Group*> groups_;
    File(hid_t file) : handle_(file) {}
    void operator=(const File&) {}

    static hid_t get(const char *name) {
	if (!std::ifstream(name)) {
	    return H5Fcreate(name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	}
	return H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT);;
    }

public:

    explicit File(const std::string &name) {
	handle_ = get(name.c_str());
    }

    File(const File &file) {
	this->handle_ = file.handle_;
	H5Iinc_ref(this->handle_);
    }

    ~File() {
	H5Fclose(handle_);
    }

    Group create_group(const std::string &name) {
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	return Group(H5Gcreate(handle_, name.c_str(), 0));
#else
	return Group(H5Gcreate1(handle_, name.c_str(), 0));
#endif
    }
    Group group(const std::string &name) {
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	return Group(H5Gopen1(handle_, name.c_str()));
#else
	return Group(H5Gopen2(handle_, name.c_str(), H5P_DEFAULT));
#endif
    }

    template<typename T, size_t N>
    Dataset dataset(const std::string &name, boost::array<size_t,N> dims);


};

#endif // FILE_HPP
