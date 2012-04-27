#ifndef _UTIL_SLICE_HPP_
#define _UTIL_SLICE_HPP_

#include "util/range.hpp"

namespace util {

    struct range_bounds {
	typedef int index_type;
	range_bounds() : mask_(0) {}
	range_bounds(index_type from) : mask_(1) {
	    bounds_[0] = from;
	}
	range_bounds(index_type from, index_type to) : mask_(3) {
	    bounds_[0] = from; bounds_[1] = to;
	}
	range_bounds from(index_type index) const {
	    range_bounds rb(*this);
	    rb.mask_ |= 1; rb.bounds_[0] = index;
	    return rb;
	}
	range_bounds to(index_type index) const {
	    range_bounds rb(*this);
	    rb.mask_ |= 2; rb.bounds_[1] = index;
	    return rb;
	}
	index_type apply_lower(index_type from) const {
	    return ((mask_ & 1) ? bounds_[0] : from);
	}
	index_type apply_upper(index_type to) const {
	    return ((mask_ & 2) ? bounds_[1] : to);
	}
    private:
	int bounds_[2], mask_;
    };

    using boost::array;

    template<size_t N>
    struct range_n {

	typedef int index_type;

	range_n() : bounds_(N, range_bounds()) {}

	range_n<N+1> operator()(const index_type (&index)[2]) const {
	    return range_n<N+1>(this->bounds_, range_bounds(index));
	}

	range_n<N> operator()(index_type index) const {
	    return range_n<N>(this->bounds_, range_bounds(index, index));
	}

	range_n<N+1> operator()(index_type from, index_type to) const {
	    return range_n<N+1>(this->bounds_, range_bounds(from, to));
	}

	range_n<N+1> operator()() const {
	    return range_n<N+1>(this->bounds_, range_bounds());
	}


	template<size_t M>
	array<index_type,M> lower(index_type from) const {
	    index_type a[M] = { from };
	    return lower(a);
	}

	template<size_t M>
	array<index_type,M> upper(index_type index) const {
	    index_type a[M] = { index };
	    return upper(a);
	}

	template<size_t M>
	array<index_type,M> lower(const array<index_type,M> &index) const {
	    assert(bounds_.size() == M);
	    array<index_type,M> a;
	    for (XRANGE(i, M)) { a[i] = bounds_[i].apply_lower(index[i]); }
	    return a;
	}

	template<size_t M>
	array<index_type,M> upper(const array<index_type,M> &index) const {
	    assert(bounds_.size() == M);
	    array<index_type,M> a;
	    for (XRANGE(i, M)) { a[i] = bounds_[i].apply_upper(index[i]); }
	    return a;
	}

    protected:
	friend class range_n<N-1>;
	typedef std::vector<range_bounds> bounds_vector;
	bounds_vector bounds_;
	range_n(const bounds_vector &B, const range_bounds &b) : bounds_(B) {
	    bounds_.push_back(b);
	    assert(bounds_.size() >= N);
	}
	    
    };


    struct range : range_n<1> {
	typedef range_n<1> base;
	typedef base::index_type index_type;
	typedef index_type index_array[2];
	static const range i, j, k, l;
	range() {}
	range(index_type index) { set(range_bounds(index, index)); }
	range(index_type from, index_type to) {  set(range_bounds(from, to)); }
	range(const index_array &r) { set(range_bounds(r[0], r[1])); }
	range from(index_type index) const { return range(bounds().from(index)); }
	range to(index_type index) const { return range(bounds().to(index)); }
    private:
	range(const range_bounds &b) { set(b); }
	void set(const range_bounds &b) { base::bounds_.back() = b; }
	range_bounds bounds() const { return base::bounds_.back(); }
    };


    range_n<1> operator<(int i, const range &r) {
	return range(r).from(i-1);
    }

    range_n<1> operator<=(int i, const range &r) {
	return range(r).from(i);
    }

    range_n<1> operator<(const range &r, int i) {
	return range(r).to(i-1);
    }

    range_n<1> operator<=(const range &r, int i) {
	return range(r).to(i);
    }

    template<size_t N>
    range_n<N+1> operator,(const range_n<N> &r, const range &r1);

}

#endif /* _UTIL_SLICE_HPP_ */
