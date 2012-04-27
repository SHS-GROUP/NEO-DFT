#ifndef CC_TRIPLES_TRIPLES_HPP
#define CC_TRIPLES_TRIPLES_HPP

#include "cc/cc.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace cc {
namespace triples {
namespace detail {

    typedef boost::numeric::ublas::column_major Layout;
    typedef boost::numeric::ublas::matrix<double, Layout> Matrix;
    typedef boost::numeric::ublas::vector<double> Vector;

    template<size_t N>
    struct index_base : boost::array<size_t,N> {
	operator const size_t*() const { return this->data(); }
    protected:
	void assign(const size_t (&data)[N]) {
	    std::copy(data, data + N, this->begin());
	}
    };

    template<size_t N>
    struct index;

    template<>
    struct index<2> : index_base<2> {
	index(const size_t (&data)[2]) {
	    this->assign(data);
	}
	index(size_t i, size_t j) {
	    size_t data[] = { i, j };
	    this->assign(data);
	}
    };

    template<>
    struct index<4> : index_base<4> {
	index(const size_t (&data)[4]) {
	    this->assign(data);
	}
	index(size_t i, size_t j, size_t k, size_t l) {
	    size_t data[] = { i, j, k, l };
	    this->assign(data);
	}
    };
	    
    typedef cc::Triples::Correction Correction;

    struct Data {
	const size_t no, nv;
	Vector eh, ep;
	Data(size_t no, size_t nv)
	    : no(no), nv(nv), eh(no), ep(nv) {}
	const Array<double>*& operator[](const std::string &key) {
	    return array_[key];
	}
	const Array<double>& operator[](const std::string &key) const {
	    if (!array_.count(key)) throw std::range_error(key);
	    return *(array_.find(key)->second);
	}
	
    private:
	std::map< std::string, const Array<double>* > array_;
    };

    struct Result {
	Correction C;
	Matrix u1;
    };


    struct Task {
	struct Tuple : boost::array<int,3> {
	    Tuple(bool state = false) {
		this->assign(0);
		state_ = state;
	    }
	    operator bool&() { return state_; }
	private:
	    bool state_;
	};
	Task(size_t N) :
	    size_((N*(N+4)*(N-1))/6), // A005581
	    index_(0) {}
	size_t size() const { return size_; }
	Tuple next() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    Tuple t;
	    try { t = operator[](index_); }
	    catch (std::range_error) { return Tuple(false); }
	    index_++;
	    return t;
	}
	Tuple operator[](int position) const {
	    if (!(position >= 0 && position < int(size_)))
		throw std::range_error(BOOST_CURRENT_FUNCTION);
	    Tuple t(true);
	    for (int i = 0; i <= position; ++i) {
		++t[2];
		if (t[2] > std::min(t[0]-1, t[1])) { t[2] = 0; ++t[1]; }
		//if (t[2] > t[1]) { t[2] = 0; ++t[1]; }
		if (t[1] > t[0]) { t[1] = 0; ++t[0]; }
	    }
	    return t;
	}
    private:
	size_t size_, index_;
	boost::mutex mutex_;
    };


    template<class T, size_t N>
    struct Permutations;


    template<class T>
    struct Permutations<T,3> {
      typedef typename boost::remove_reference<T>::type U;
	T ij, ik, jk;
	Permutations(const U &ij, const U &ik, const U &jk)
	    : ij(ij), ik(ik), jk(jk) {}
    };

    template<class T>
    struct Permutations<T,6> {
        typedef typename boost::remove_reference<T>::type U;
	T ij, ik, jk, ji, ki, kj;
	Permutations(const U &ij, const U &ik, const U &jk,
		     const U &ji, const U &ki, const U &kj)
	    : ij(ij), ik(ik), jk(jk), ji(ji), ki(ki), kj(kj) {}
    };


    struct Queue {
	struct Tuple {
	    int size, begin;
	};
	explicit Queue(int size = 0)
	    : size_(size), begin_(0) {}
	Tuple advance(int n = 1) {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    n = std::min(n, (size_ - begin_));
	    n *= (n > 0);
	    Tuple t = { n, begin_ };
	    begin_ += n;
	    return t;
	}
	void reset(int size) {
	    size_ = size;
	    begin_ = 0;
	}
    private:
	int size_, begin_;
	boost::mutex mutex_;
    };


} // namespace detail
} // namespace triples
} // namespace cc

#endif // CC_TRIPLES_TRIPLES_HPP
