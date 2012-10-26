#ifndef BOOST_UTILITY_PROFILER_HPP
#define BOOST_UTILITY_PROFILER_HPP

// #include "boost/utility/timer.hpp"

#include "boost/date_time/posix_time/posix_time_types.hpp"

#include <string>
#include <sstream>
#include <map>
#include <typeinfo>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/typeof/typeof.hpp>

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

#define BOOST_PROFILE_LINE							\
    boost::utility::profiler::event BOOST_PROFILE__(__LINE__)			\
    (boost::utility::global_profiler()						\
     [boost::utility::detail::profiler_event(std::string(__FILE__) + ":" +	\
					     (BOOST_PP_STRINGIZE(__LINE__)))])

#define BOOST_PROFILE_REGISTER_THREAD			\
    boost::utility::global_profiler().register_thread()

#else
#define BOOST_PROFILE_FUNCTION(...)
#define BOOST_PROFILE_TYPE(T, ...)
#define BOOST_PROFILE_LINE
#define BOOST_PROFILE(...)
#define BOOST_PROFILE_REGISTER_THREAD
#endif

#define BOOST_PROFILE_DUMP(stream)			\
    boost::utility::global_profiler().dump((stream))

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
		++size_;
		value_ += t;
		return *this;
	    }
	    event_data& operator++() { return (*this += 1); }
	    template<class O>
	    O& to_stream(O &ostream) const {
		ostream << value_ << "/" << size_;
		return ostream;
	    }
	private:
	    size_t size_;
	    double value_;
	};

	struct event {
	    event(event_data *data)
		: data_(data),
		  start_(microsec_clock::universal_time()) {}
	    ~event() {
		boost::posix_time::time_duration t = 
		    (microsec_clock::universal_time() - start_);
		// std::cout << timer_ << std::endl;
		if (data_) {
		    *data_ += t.total_microseconds()/1.0e6;
		}
	    }
	private:
	    typedef boost::posix_time::microsec_clock microsec_clock;
	    event_data *data_;
	    boost::posix_time::ptime start_;
	    // utility::timer timer_;
	};

	void register_thread() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    //std::cout << "register thread " << thread();
	    events_[thread()];
	}

	event_data* operator[](const event_key &key) {
	    //std::cout << key << std::endl;
	    std::string thread = this->thread();
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    if (events_.find(thread) == events_.end()) return NULL;
	    //std::cout << thread << ": " << key << std::endl;
	    return &events_[thread][key];
	}

	void clear() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    BOOST_AUTO(it, events_.begin());
	    while (it != events_.end()) {
		it++->second.clear();
	    }
	}

	template<class S>
	S& to_stream(S &ostream) const {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    typedef std::map<event_key, event_data> events_t;
	    typename std::map<std::string, events_t>::const_iterator
		e = events_.begin();
	    //std::cout << events_.size() << std::endl;
	    while (e != events_.end()) {
		ostream << "Thread " << e->first << ":" << std::endl;
		typename events_t::const_iterator it = e->second.begin();
		while (it != e->second.end()) {
		    ostream << "    " << it->first << ": ";
		    it->second.to_stream(ostream);
		    ostream << std::endl;
		    ++it;
		}
		++e;
	    }
	    return ostream;
	}

	template<class S>
	void dump(S &stream) {
	    this->to_stream(stream);
	    clear();
	}

    private:
	mutable boost::mutex mutex_;
	std::map<
	std::string, std::map<event_key, event_data>
	> events_;

	static std::string thread() {
	    std::stringstream ss;
	    ss << boost::this_thread::get_id();
	    return ss.str();
	}
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
