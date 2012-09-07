#include "config.h"

#include "runtime.hpp"
#include "exception.hpp"
#include "foreach.hpp"

#include "omp.hpp"
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
#include "cuda.hpp"
#endif

#include <cstdlib>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

std::auto_ptr<runtime> runtime::rt_;

std::auto_ptr<Runtime> Runtime::rt_;

Runtime& Runtime::rt() {
    if (!Runtime::rt_.get()) {
	Runtime::rt_.reset(new Runtime());
    }
    return *rt_.get();
}

namespace {

    size_t page_size() {
	return sysconf(_SC_PAGE_SIZE);
    }

    size_t avphys_pages() {
	return sysconf(_SC_AVPHYS_PAGES);
    }

    size_t phys_pages() {
	return sysconf(_SC_PHYS_PAGES);
    }

    struct options {
    private:
	std::map<std::string,std::string> data_;

	std::vector<std::string>
	split(std::string string, const char *chars) {
	    std::vector<std::string> strings;
	    if (!string.empty())
		boost::split(strings, string, boost::is_any_of(chars));
	    return strings;
	}

    public:
	explicit options(const char *string) {
	    const std::vector<std::string> &argv =
		split((string ? string : ""), ";");

	    foreach (std::string arg, argv) {
		boost::trim(arg);
		if (arg.empty()) continue;
		size_t pos = arg.find("=");
		bool key = (pos != std::string::npos);
		std::string k = arg.substr(0,pos);
		std::string v = (key) ? arg.substr(pos+1) : "";
		boost::trim(k);
		boost::trim(v);
		if (k.empty()) throw CCHEM_EXCEPTION("invalid key/value: " + arg);
		data_[k] = v;
	    }

	}

	bool test(const std::string &k) const {
	    return (data_.count(k) > 0);
	}

	template <typename T>
	T get(const std::string &k, const T& v) const {
	    if (!test(k)) return v;
	    return boost::lexical_cast<T>(data_.find(k)->second);
	}

	size_t memory() const {
	    if (!test("memory")) {
		size_t memory = (size_t)(0.9*avphys_pages()*page_size());
		Parallel().broadcast(&memory, 1, 0);
		return memory;
	    }
	    std::string v = data_.find("memory")->second;
	    BOOST_AUTO(last, v.end()-1);
	    size_t scale = 1;
	    char suffix = tolower(*last);
	    if (suffix == 'k') scale *= (1<<10);		
	    if (suffix == 'm') scale *= (1<<20);		
	    if (suffix == 'g') scale *= (1<<30);		
	    if (scale != 1) v.erase(last);
	    return (size_t)(scale*boost::lexical_cast<float>(v));
	}

	size_t num_threads() const {
	    return get("num_threads", 0);
	}

	std::vector<int> devices() {
	    std::vector<int> devices;

	    if (!test("devices")) {
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
		try {
		    int N = cuda::get_device_count();
		    for (int i = 0; i < N; ++i) {
			BOOST_AUTO(const &p, cuda::get_device_properties(i));
			if (p.major > 1 || (p.major == 1 && p.minor == 3)) {
			    devices.push_back(i);
			}
		    }
		}
		catch (cuda::error e) {
		    //std::cerr << e.what() << std::endl;
		}
#endif
	    }
	    else {
		// std::cout << get<std::string>("devices", "") << std::endl;
		const std::vector<std::string> &argv =
		    split(get<std::string>("devices", ""), ",");
		foreach (const std::string &arg, argv) {
		    // std::cout << arg << std::endl;
		    devices.push_back(boost::lexical_cast<int>(arg));
		}
	    }

	    return devices;
	}

    };

}

Runtime::Runtime() {
    cout_ = &std::cout;
    options p(std::getenv("CCHEM"));

    num_threads_ = p.num_threads();
    devices_ = p.devices();
    memory_.available_ = p.memory();

    if (Parallel().rank() == 0) {
	std::cout << "libcchem: " << std::endl;
	std::cout << "memory: " << memory_.available()/double(1<<20)
		  << " MB" << std::endl;
	std::cout << "threads: " << omp_get_max_threads() << std::endl;
	std::cout << "gpu: { ";
	for (size_t i = 0; i < devices_.size(); ++i) {
	    if (i) std::cout << ", ";
	    std::cout << devices_.at(i);
	}
	std::cout << " }" << std::endl;
	std::cout << *this << std::endl;
    }
}

void Runtime::Memory::clear() {
    BOOST_AUTO(it, data_.begin());
    while (it != data_.end()) {
	this->free((it++)->first);
    }
}

void* Runtime::Memory::realloc(void *ptr, size_t size) {
    if (!ptr && !size) return NULL;
#pragma omp critical(cchem_memory)
    {
	BOOST_AUTO(it, data_.find(ptr));
	size_t allocated = 0;
	if (ptr) {
	    if (it == data_.end()) CCHEM_ERROR("pointer %p not mapped\n", ptr);
	    allocated = it->second;
	}
	if (!size || (ptr && size > allocated)) { // free || reallocate
	    available_ += allocated;
	    used_ -= allocated;
	    data_.erase(it);
	    // CCHEM_MESSAGE("free(%lu)\n", ptr);
	    ::free(ptr);
	    ptr = NULL;
	}
	// size <= allocated bypasses this
	if (size > allocated) { // need to allocate
	    if (available_ < size) {
		CCHEM_ERROR("not enough memory: %lu bytes requested, %lu available\n",
			    size, available_);
	    }
	    size_t align = getpagesize();
	    int err = ::posix_memalign(&ptr, align, size);
	    if (err) CCHEM_ERROR("posix_memalign(%p, %lu, %lu) returned %i\n",
				 &ptr, align, size, err);
	    // CCHEM_MESSAGE("posix_memalign(%p, %lu, %lu)\n", ptr, align, size);
	    data_[ptr] = size;
	    available_ -= size;
	    used_ += size;
	}
    }
    return ptr;
}

std::ostream& operator<<(std::ostream &os, const Runtime &rt) {
    os << "sysconf:" << std::endl;
    os << "    page size: " << page_size() << std::endl;
    os << "    phys pages: " << phys_pages() << std::endl;
    os << "    avphys pages: " << avphys_pages() << std::endl;
    return os;
}
