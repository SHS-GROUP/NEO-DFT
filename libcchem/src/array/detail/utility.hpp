#ifndef ARRAY_DETAIL_UTILITY_HPP 
#define ARRAY_DETAIL_UTILITY_HPP

#include <algorithm>
#include <boost/config.hpp>

namespace array {
namespace detail {

    struct range {
	BOOST_GPU_ENABLED
	range(int start, int finish)
	    : start_(start), size_(finish - start) {}
	BOOST_GPU_ENABLED
	int start() const { return start_; }
	BOOST_GPU_ENABLED
	int size() const { return size_; }
    private:
	int start_, size_;
    };

    struct serial {
	BOOST_GPU_ENABLED static int x() { return 0; }
	BOOST_GPU_ENABLED static int y() { return 0; }
	BOOST_GPU_ENABLED static int nx() { return 1; }
	BOOST_GPU_ENABLED static int ny() { return 1; }	    
    };

    struct thread {
	thread(int id, int count) : id(id), count(count) {}
	const int id, count;
    };

    template<typename T>
    BOOST_GPU_ENABLED
    void swap(T &a, T &b) {
	T _ = a;
	a = b;
	b = _;
    }

    template<typename T>
    BOOST_GPU_ENABLED
    void swap(T begin, T end, T other) {
	while (begin != end) {
	    swap(*begin++, *other++);
	}
    }

    template<typename T, class B>
    BOOST_GPU_ENABLED
    void swap(T begin, T end, T other, const B &b) {
	begin += b.x();
	other += b.x();
	while (begin < end) {
	    swap(*begin, *other);
	    begin += b.nx();
	    other += b.nx();
	}
    }
    
    template<typename T>
    BOOST_GPU_ENABLED
    const T& min(const T &a, const T &b) {
	return (a < b) ? a : b;
    }

}
}

#endif // ARRAY_DETAIL_UTILITY_HPP
