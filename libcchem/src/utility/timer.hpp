#ifndef UTILITY_TIMER_HPP
#define UTILITY_TIMER_HPP

#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace utility {

    struct timer {
	boost::posix_time::ptime start_;
	typedef boost::posix_time::microsec_clock microsec_clock;
	timer() : start_(microsec_clock::universal_time()) {}
	boost::posix_time::time_duration duration() const {
	    return microsec_clock::universal_time() - start_;
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
