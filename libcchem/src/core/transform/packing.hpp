#ifndef TRANSFORM_PACKING_HPP
#define TRANSFORM_PACKING_HPP

#include <vector>
#include <algorithm>
#include <boost/typeof/typeof.hpp>

namespace transform {
namespace detail {

    template<class R, class Input, class Output, class Index>
    static size_t pack(R (*copy)(Input, Input, Output),
		       size_t N, Input input, Output output,
		       const Index &index) {
	size_t k = 0;
	for (size_t i = 0; i < index.size(); i += 2) {
	    size_t start = index.at(i), stop = index.at(i+1);
	    copy(input + start*N, input + stop*N, output + k*N);
	    k += (stop - start);
	}
	return k;
    }

    struct packing {

	struct index {
	    index(size_t reserve = 64) {
		data.reserve(reserve);
	    }
	    template<class R>
	    void push(const R& range) {
		size_t stop = range.start() + range.size();
		if ((data.empty()) || (size_t(range.start()) != data.back())) {
		    data.push_back(range.start());
		    data.push_back(stop);
		}
		else {
		    data.back() = stop;
		}
	    }
	    void extend(const index &index) {
		BOOST_AUTO(begin, index.data.begin());
		if (index.data.empty()) return;
		if ((!data.empty()) && (index.data.front() == data.back())) {
		    data.pop_back();
		    ++begin;
		}
		data.insert(data.end(), begin, index.data.end());
	    }
	    size_t size() const {
		size_t size = 0;
		for (size_t i = 0; i < data.size(); i += 2) {
		    size += data.at(i+1) - data.at(i);
		}
		return size;
	    }
	    static std::vector<std::pair<size_t,size_t> >
	    pairs(const std::vector<size_t> &data) {
		typedef std::pair<size_t,size_t> pair;
		std::vector<pair> pairs;
		for (size_t i = 0; i < data.size(); i += 2) {
		    pairs.push_back(pair(data.at(i), data.at(i+1)));
		}
		return pairs;
	    }

	    template<class R, class Input, class Output>
	    size_t pack(R (*copy)(Input, Input, Output),
			size_t N, Input input, Output output) const {
		return detail::pack(copy, N, input, output, data);
	    }

	    std::vector<size_t> data;
	};

	struct index2 {
	    typedef std::vector<packing::index> index_type;
	    void push2() {
		data_.push_back(packing::index());
	    }
	    packing::index& back2() { return data_.back(); }
	    index_type& data() { return data_; }
	    const index_type& data() const { return data_; }
	    void extend(const index2 &other) {
		data_.resize(other.data_.size());
		for (size_t i = 0; i < other.data().size(); ++i) {
		    data_[i].extend(other.data_[i]);
		}
	    }
	private:
	    index_type data_;
	};

    };

}
}

#endif // TRANSFORM_PACKING_HPP
