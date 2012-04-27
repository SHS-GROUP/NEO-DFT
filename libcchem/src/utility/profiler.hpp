#ifndef BOOST_UTILITY_PROFILER_HPP
#define BOOST_UTILITY_PROFILER_HPP

#include "utility/timer.hpp"

#include <string>
#include <map>
#include <boost/thread.hpp>

#define BOOST_PROFILE_FUNCTION(...)					\
    utility::profiler::event					\
    event__(utility::profiler::global[				\
        utility::detail::profiler_event(__PRETTY_FUNCTION__)(__VA_ARGS__)])
	

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
	    typedef boost::tuple<profiler&, const event_key&> constructor;
	    size_t size_;
	    double value_;
	    mutable boost::mutex mutex_;
	};

	struct event {
	    event(event_data &data) : data_(data) {}
	    ~event() {
		// std::cout << timer_ << std::endl;
		data_ += double(timer_);
	    }
	    event_data &data_;
	    utility::timer timer_;
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
	static profiler global;
    private:
	std::map<event_key, event_data> events_;
	mutable boost::mutex mutex_;
    };

    template<typename T>
    profiler<T> profiler<T>::global;

    template<typename T>
    inline std::ostream& operator<<(std::ostream &ostream, const profiler<T> &p) {
	return p.to_stream(ostream);
    }
    
} // namespace detail

    typedef detail::profiler<void> profiler;

}
}


#endif // BOOST_UTILITY_PROFILER_HPP
