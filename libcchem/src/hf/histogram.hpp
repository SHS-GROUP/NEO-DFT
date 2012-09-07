#ifndef HF_HISTOGRAM_HPP
#define HF_HISTOGRAM_HPP

#include <map>
#include <iostream>
#include <boost/typeof/typeof.hpp>

namespace hf {

    struct Histogram {
	std::map<int, double> gpu, cpu;
	double total() const {
	    double total;
	    for (int i = 0; i < 2; ++i) {
		BOOST_AUTO(const &data, (i == 0) ? this->gpu : this->cpu);
		BOOST_AUTO(it, data.begin());
		while (it != data.end()) {
		    total += it->second;
		    ++it;
		}
	    }
	    return total;
	}
    };

    std::ostream& operator << (std::ostream &os, const Histogram &h) {
	for (int i = 0; i < 2; ++i) {
	    BOOST_AUTO(const &data, (i == 0) ? h.gpu : h.cpu);
	    BOOST_AUTO(it, data.begin());
	    while (it != data.end()) {
		os << it->first << ": "
		   << it->second << " "
		   << ((i == 0) ? "gpu" : "cpu")
		   << "\n";
		++it;
	    }
	}
	return os;
    }

} // namespace hf

#endif /* HF_HISTOGRAM_HPP */
