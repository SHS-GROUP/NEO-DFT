#ifndef RUNTIME_HPP
#define RUNTIME_HPP

#include <string>
#include <vector>
#include <cstdlib>

#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/thread.hpp>

#include "config.h"
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
#include "boost/cuda/runtime.hpp"
#endif

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

	put("gpu::devices", std::vector<int>());
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
	put("gpu::devices",
	    runtime::gpu_devices(variable(prefix, "GPU_DEVICES")));
#endif
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

    static std::vector<int> gpu_devices(std::string values) {
	std::vector<int> devices;
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
	try {
	    devices = gpu_devices<int>(values, boost::cuda::devices(1.3));
	}
	catch (std::exception &e) {
	    std::cout << e.what() << std::endl;
	}
#endif
	return devices;
    }

};

#endif /* RUNTIME_HPP */
