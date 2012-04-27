#ifndef ARRAY_DETAIL_GRID_HPP 
#define ARRAY_DETAIL_GRID_HPP

#include "array/detail/utility.hpp"
#include <boost/config.hpp>
#include <boost/typeof/typeof.hpp>

namespace array {
namespace detail {

    struct grid1 {
	grid1(int id, int size) : id_(id), size_(size) {}
	int id() const { return id_; }
	int size() const { return size_; }
    private:
	int id_, size_;
    };


    struct grid {

	struct serial {
	    static int x() { return 0; }
	    static int y() { return 0; }
	    static int nx() { return 1; }
	    static int ny() { return 1; }
	};

	struct X : grid1 {
	    X(int id, int size) : grid1(id, size) {}
	    int x() const { return grid1::id(); }
	    int nx() const { return grid1::size(); }
	    static int y() { return 0; }
	    static int ny() { return 1; }
	};

	struct Y : grid1 {
	    Y(int id, int size) : grid1(id, size) {}
	    int y() const { return grid1::id(); }
	    int ny() const { return grid1::size(); }
	    static int x() { return 0; }
	    static int nx() { return 1; }
	};

	struct Z : grid1 {
	    Z(int id, int size) : grid1(id, size) {}
	    int z() const { return grid1::id(); }
	    int nz() const { return grid1::size(); }
	    static int x() { return 0; }
	    static int y() { return 0; }
	    static int nx() { return 1; }
	    static int ny() { return 1; }
	};

	grid(int x, int y, int nx, int ny)
	    : x_(x), y_(y), nx_(nx), ny_(ny) {}
	int x() const { return x_; }
	int y() const { return y_; }
	int nx() const { return nx_; }
	int ny() const { return ny_; }
	static int z() { return 0; }
	static int nz() { return 1; }
    private: 
	int x_, y_;
	int nx_, ny_;
    };

}
}

#endif // ARRAY_DETAIL_GRID_HPP
