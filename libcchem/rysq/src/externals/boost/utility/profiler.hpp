#ifndef BOOST_UTILITY_PROFILER_HPP
#define BOOST_UTILITY_PROFILER_HPP

// #include "boost/utility/timer.hpp"

#include "boost/date_time/posix_time/posix_time_types.hpp"

#include <string>
#include <map>
#include <typeinfo>
#include <boost/thread/mutex.hpp>
#include <boost/preprocessor/cat.hpp>

#define BOOST_PROFILE__(tag) BOOST_PP_CAT(boost_profile__,  tag)

#ifdef BOOST_PROFILE_ENABLE
#define BOOST_PROFILE_FUNCTION(...)					\
    boost::utility::profiler::event BOOST_PROFILE__(__LINE__)		\
    (boost::utility::global_profiler()					\
     [boost::utility::detail::profiler_event				\
      (__PRETTY_FUNCTION__)(__VA_ARGS__)])

#define BOOST_PROFILE_TYPE(T, ...)					\
    boost::utility::profiler::event BOOST_PROFILE__(__LINE__)		\
    (boost::utility::global_profiler()					\
     [boost::utility::detail::profiler_event				\
      (typeid(T).name())(__VA_ARGS__)])
	
#define BOOST_PROFILE(...)					\
    boost::utility::profiler::event BOOST_PROFILE__(__LINE__)	\
    (boost::utility::global_profiler()				\
     [boost::utility::detail::profiler_event(__VA_ARGS__)])

#else
#define BOOST_PROFILE_FUNCTION(...)
#define BOOST_PROFILE_TYPE(T, ...)
#define BOOST_PROFILE(...)
#endif

namespace boost {
namespace utility {
namespace detail {

    struct profiler_event {
	profiler_event(const std::string &key) : data_(key) {}
	operator const std::string&() const { return data_; }
	profiler_event& operator()(const std::string &key) {
	    data_ += (":" + key);
	    return *this;
	}
	profiler_event& operator()() { return *this; }
    private:
	std::string data_;
    };

    template<class>
    struct profiler {

	typedef std::string event_key;

	struct event_data {
	    event_data(): size_(0), value_(0) {}
	    event_data(const event_data &e)
		: size_(e.size_), value_(e.value_) {} 
	    event_data& operator+=(double t) {
		boost::lock_guard<boost::mutex> lock(mutex_);
		++size_;
		value_ += t;
		return *this;
	    }
	    event_data& operator++() { return (*this += 1); }
	    std::ostream& to_stream(std::ostream &ostream) const {
		boost::lock_guard<boost::mutex> lock(mutex_);
		ostream << value_ << "/" << size_;
		return ostream;
	    }
	private:
	    size_t size_;
	    double value_;
	    mutable boost::mutex mutex_;
	};

	struct event {
	    event(event_data &data)
		: data_(data),
		  start_(microsec_clock::universal_time()) {}
	    ~event() {
		boost::posix_time::time_duration t = 
		    (microsec_clock::universal_time() - start_);
		// std::cout << timer_ << std::endl;
		data_ += t.total_microseconds()/1.0e6;
	    }
	private:
	    typedef boost::posix_time::microsec_clock microsec_clock;
	    event_data &data_;
	    boost::posix_time::ptime start_;
	    // utility::timer timer_;
	};

	event_data& operator[](const event_key &key) {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    return events_[key];
	}
	std::ostream& to_stream(std::ostream &ostream) const {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    typename std::map<event_key, event_data>::const_iterator
		it = events_.begin();
	    while (it != events_.end()) {
		ostream << it->first << ": ";
		it->second.to_stream(ostream);
		ostream << std::endl;
		++it;		
	    }
	    return ostream;
	}
	void clear() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    events_.clear();
	}
	// static profiler global;
    private:
	std::map<event_key, event_data> events_;
	mutable boost::mutex mutex_;
    };

    // template<typename T>
    // profiler<T> profiler<T>::global;

    template<typename T>
    inline std::ostream& operator<<(std::ostream &ostream, const profiler<T> &p) {
	return p.to_stream(ostream);
    }
    
} // namespace detail

    typedef detail::profiler<void> profiler;
    inline profiler& global_profiler() {
	static profiler p;
	return p;
    }
    

}
}


#endif // BOOST_UTILITY_PROFILER_HPP
