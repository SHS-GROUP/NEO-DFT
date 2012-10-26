#ifndef CCHEM_UTILITY_PROGRESS_HPP
#define CCHEM_UTILITY_PROGRESS_HPP

#include <memory>
#include <boost/progress.hpp>

namespace utility {

    struct Progress {
	std::auto_ptr<boost::progress_display> display_;
	void reset(int n) {
	    display_.reset(new boost::progress_display(n));
	}
	void operator+=(int n) {
	    if (display_.get()) (*display_) += n;
	}
	void operator++() {
	    this->operator+=(1);
	}
	size_t count() const {
	    if (display_.get()) return display_->count();
	    return 0;
	}
	void jump(int pos) {
	    this->operator+=(pos - this->count());
	}
    };

} // utility

#endif /* CCHEM_UTILITY_PROGRESS_HPP */
