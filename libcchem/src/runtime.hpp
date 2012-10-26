#ifndef RUNTIME_HPP
#define RUNTIME_HPP

#include "config.h"

#include "parallel.hpp"
#include "array/hdf5.hpp"
#include "array/memory.hpp"
#if defined(HAVE_GA)
#include "array/ga.hpp"
#endif

#include "dictionary.hpp"
#include "file.hpp"


#include <string>
#include <vector>
#include <map>

#include <cstdlib>
#include <cstring>

#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/thread.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/progress.hpp>


struct runtime {

    static std::string default_value(const std::string &key) {
	std::map< std::string, std::string> default_value;
	//default_value["CC_GPU_DEVICES"] = "-1";
	default_value["MP2_GPU_DEVICES"] = "-1";
	return default_value[key];
    }

    runtime(const char *prefix = "") {
	std::string num_threads = variable(prefix, "NUM_THREADS");

	put<int>("team::count", 1);
	if (num_threads.empty()) {
	    put<int>("team::size", runtime::num_threads());
	}
	else {
	    size_t pos = num_threads.find("x");
	    put<int>("team::size", num_threads.substr(0, pos));
	    if (pos != std::string::npos) {
		put<int>("team::count", num_threads.substr(pos+1));
	    }
	}
	put<int>("threads", (get<int>("team::size")*get<int>("team::count")));
    }

    static runtime& rt() {
	if (!runtime::rt_.get()) rt_.reset(new runtime());
	return *rt_;
    }

    template<typename T>
    void put(const std::string &key, const T &val) {
	data_[key] = boost::any(val);
    }

    template<typename T>
    void put(const std::string &key, const std::string &val) {
	data_[key] = boost::any(boost::lexical_cast<T>(val));
    }

    struct no_such_key : std::runtime_error {
	no_such_key(const char *what) : runtime_error(what) {}
    };

    template<typename T>
    T get(const char *key) const {
	if (!data_.count(key)) throw no_such_key(key);
	return try_get<T>(key);
    }
    
    template<typename T>
    T get(const char *key, const T &value) const {
	if (!data_.count(key)) return value;
	return try_get<T>(key);	
    }

private:
    std::map<std::string, boost::any> data_;
    static std::auto_ptr<runtime> rt_;

    template<typename T>
    T try_get(const char *key) const {
	try {
	    return boost::any_cast<T>(data_.find(key)->second);
	}
	catch (boost::bad_any_cast) {
	    throw std::runtime_error(std::string(BOOST_CURRENT_FUNCTION) + ": " + key);
	}	
    }

public:

    static std::string variable(const char *prefix, const char *variable) {
	std::string var(prefix);
	if (!var.empty()) var += "_";
	var += variable;
	//std::cout << var << std::endl;

	char *value = std::getenv(var.c_str());

	if (value) return std::string(value);
	if (!default_value(var).empty()) return default_value(var);
	if (strlen(prefix)) return runtime::variable("", variable);
	return std::string("");
    }

    static int profile() {
	char *v = std::getenv("CCHEM_PROFILE");
	return (v ? atoi(v) : 0);
    }

    template<typename T>
    static
    T num_threads(const T &default_value, const char *variable = "NUM_THREADS") {
	char *v = std::getenv(variable);
	std::string value((v) ? v : "");
	if (value.empty()) return default_value;
	return boost::lexical_cast<T>(value);
    }

    static size_t num_threads() {
	// return 1;
	return num_threads<size_t>(boost::thread::hardware_concurrency());
    }

    static size_t num_threads(const char *method) {
	std::string name = std::string(method) + "_NUM_THREADS";
	char *v = 0;
	v = std::getenv(name.c_str());
	std::string value((v) ? v : "");
	if (value.empty()) return num_threads();
	return boost::lexical_cast<size_t>(value);
    }

    template<typename T>
    static
	std::vector<T> gpu_devices(std::string values,
				   const std::vector<T> &default_value) {
	if (values.empty()) return default_value;

	std::vector<T> devices;
	while (!values.empty()) {
	    size_t token = values.find(",");
	    T t = boost::lexical_cast<T>(values.substr(0,token));
	    if (t < T(0)) return std::vector<T>();
	    devices.push_back(t);
	    token -= (token == std::string::npos);
	    values.erase(0, token+1);
	}
	return devices;
    }

//     static std::vector<int> gpu_devices(std::string values) {
// 	std::vector<int> devices;
// 	return devices;
// #if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
// 	try {
// 	    devices = gpu_devices<int>(values, boost::cuda::devices(1.3));
// 	}
// 	catch (std::exception &e) {
// 	    std::cout << e.what() << std::endl;
// 	}
// #endif
// 	return devices;
//     }

};

struct Runtime : Dictionary {

    struct Exception {};

    File file(const std::string &name = "data.h5") {
	return File(name);
    }

    struct Array : boost::noncopyable {

	bool contains(const std::string &name) const {
	    return data_.count(name);
	}

	template<typename T, size_t N, typename U>
	::Array<T>* distributed(const char *name, const U &dims,
				const Parallel &pe,
				const size_t *chunk = NULL) {
	    typedef ::Array<T> A;
#ifdef HAVE_GA
	    try {
		typename A::GA ga;
		if (chunk) ga.set_chunk(N, chunk);
		data_[name] = new A(N, dims, name, typename A::GA());
		return find<A>(name);
	    }
	    catch (typename A::Error e) {
		size_t size = 1;
		for (size_t i = 0; i < N; ++i) size *= dims[i];
		pe.cout() << "Not enough distributed memory to allocate "
			  << name << ", " << size*sizeof(T)/1e9 << " GB"
			  << std::endl;
		//CCHEM_MESSAGE("%s\n", e.what());
	    }
#endif
#if defined(HAVE_MPI) && defined(H5_HAVE_PARALLEL)
	    try {
		typename A::HDF5 hdf5("arrays.h5");
		if (chunk) {
		    hdf5.set_chunk(N, chunk);
		    //hdf5.set_deflate(5);
		}
		hdf5.set_plist(parallel::mpi_comm());
		data_[name] = new A(N, dims, name, hdf5);
		return find<A>(name);
	    }
	    catch (typename A::Error e) {
		//CCHEM_MESSAGE("%s\n", e.what());
	    }
#endif
	    throw Exception();
	}

	template<typename T, size_t N>
	::Array<T>* file(const char *name, const size_t *dims) {
	    typedef ::Array<T> A;
	    typename A::HDF5 hdf5("arrays.h5");
	    data_[name] = new A(N, dims, name, hdf5);
	    return find<A>(name);
	}

	template<typename T, size_t N>
	void core(const char *name, const size_t *dims) {
	    data_[name] = new ::Array<T>(N, dims, name, ::Array<>::Memory());
	}

	template<typename T, size_t N>
	void allocate(const char *name, boost::array<size_t,N> dims,
		      const Parallel &pe, const size_t *chunk = NULL) {
	    allocate<T,N>(name, &dims[0], pe, chunk);
	}

	template<typename T, size_t N>
	void allocate(const char *name, const size_t *dims,
		      const Parallel &pe, const size_t *chunk = NULL) {
// 	    for (int i = 0; i < N; ++i) {
// 		std::cout << dims[i] << ", ";
// 	    }
// 	    std::cout << std::endl;
	    try {
		distributed<T,N>(name, dims, pe, chunk);
	    }
	    catch (Exception e) {
		if (pe.size() > 1) {
		    pe.cout() << "Failed to allocate distributed "
			      << name << " array" << std::endl;
		    throw cchem::exception();
		}
		file<T,N>(name, dims);
	    }
	}

	template<class A>
	void erase(const std::string &name) {
	    if (contains(name)) delete find<A>(name);
	    data_.erase(name);
	}
	template<class A>
	A* find(const std::string &name) {
	    BOOST_AUTO(it, data_.find(name));
	    if (it == data_.end())
		throw std::runtime_error(name + ": no such key");
	    return boost::any_cast<A*>(it->second);
	}
    private:
	std::map<std::string, boost::any> data_;
    };
    Array& arrays() { return arrays_; }

public:
    static Runtime& rt();

    struct Stream {
	explicit Stream(std::ostream *os = &std::cout)
	    : os_(os) {}
	operator bool() const { return (os_ != NULL); }
	Stream& operator<<(std::ostream& (*f)(std::ostream&)) {
	    if (os_) f(*os_);
	    return *this;
	}
	template<class O>
	Stream& operator<<(const O &o) {
	    if (os_) *os_ << o;
	    return *this;
	}
    private:
	std::ostream *os_;
    };

    void set_cout(std::ostream *os) { cout_ = os; }
    Stream cout() { return Stream(cout_); }

    struct Progress : boost::noncopyable {
	explicit Progress(const Stream &stream = Stream())
	    : stream_(stream) {}
	void reset(size_t size) {
	    if (stream_) display_.reset(new boost::progress_display(size));
	}
	void operator+=(int step) {
	    if (display_.get()) (*display_.get()) += step;
	}
	void operator++(int) { (*this) += 1; }
    private:
	Stream stream_;
	std::auto_ptr<boost::progress_display> display_;
    };

    Runtime();

    size_t num_threads() const { return num_threads_; }
    const std::vector<int>& devices() const { return devices_; }

    std::vector<int> devices(const parallel::Comm &comm) const {
	std::vector<int> devices, all = this->devices();
	for (size_t i = comm.rank(); i < all.size(); i += comm.size()) {
	    devices.push_back(all.at(i));
	}
	return devices;
    }

    struct Memory : boost::noncopyable {
	~Memory() { clear(); }
	size_t available() const { return available_; }
	size_t used() const { return used_; }
	template<typename T>
	T* realloc(T* ptr, size_t size) {
	    return (T*)this->realloc((void*)ptr, size*sizeof(T));
	}
	template<typename T>
	T* malloc(size_t size) {
	    return realloc<T>(NULL, size);
	}
	void free(void* ptr) { this->realloc((void*)ptr, 0); }
	void* malloc(size_t size) { return this->realloc(NULL, size); };
	void* realloc(void* ptr, size_t size);
	void clear();
    private:
	size_t available_, used_;
	std::map<void*,size_t> data_;
	friend class Runtime;
	Memory() {
	    available_ = 0;
	    used_ = 0;
	}
    };

    Memory& memory() {
	return memory_;
    }

private:
    Array arrays_;
    std::ostream *cout_;
    size_t available_;
    size_t num_threads_;
    std::vector<int> devices_;
    Memory memory_;

private:
    static std::auto_ptr<Runtime> rt_;
    void set_options(std::string argv);
    void set_devices(std::string value);
};

std::ostream& operator<<(std::ostream &os, const Runtime &rt);

#endif /* RUNTIME_HPP */
