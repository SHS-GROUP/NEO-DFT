#ifndef UTILITY_TIMER_HPP
#define UTILITY_TIMER_HPP

#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace utility {

    struct timer {
    private:
	boost::posix_time::ptime start_;
	boost::posix_time::time_duration duration_;
	bool running_;
    public:
	typedef boost::posix_time::time_duration value_type;
	typedef boost::posix_time::microsec_clock microsec_clock;
	timer(bool running = true) {
	    running_ = running;
	    if (running_) start_ = microsec_clock::universal_time();
	}
	boost::posix_time::time_duration duration() const {
	    value_type d;
	    if (running_) d = (microsec_clock::universal_time() - start_);
	    return duration_ + d;
	}
	void reset() {
	    duration_ = value_type();
	    start_ = microsec_clock::universal_time();
	}
	void stop() {
	    duration_ = duration();
	    start_  = microsec_clock::universal_time();
	    running_ = false;
	}
	void start() {
	    start_ = microsec_clock::universal_time();
	    running_ = true;
	}
	operator value_type() const {
	    return duration();
	}
	operator double() const {
	    return duration().total_microseconds()/1e6;
	}    
	long total_seconds() const {
	    return duration().total_seconds();
	}
    };
    inline std::ostream& operator << (std::ostream &os, const timer &t) {
	return os << boost::posix_time::to_simple_string(t.duration());
    }

}

#endif // UTILITY_TIMER_HPP
