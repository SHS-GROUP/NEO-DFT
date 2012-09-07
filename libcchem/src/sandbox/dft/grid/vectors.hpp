#ifndef DFT_GRID_VECTORS_HPP
#define DFT_GRID_VECTORS_HPP

#include <boost/typeof/typeof.hpp>

#include "blas.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/mpl/if.hpp>
namespace dft {
namespace grid {

    struct Vectors {

	typedef boost::mpl::true_ ao_major;

	typedef boost::numeric::ublas::column_major layout;
	
	typedef boost::numeric::ublas::matrix<double, layout> matrix_type;
	typedef boost::numeric::ublas::matrix_row<matrix_type> row_type;
	typedef boost::numeric::ublas::matrix_column<matrix_type> column_type;
	typedef BOOST_TYPEOF(boost::numeric::ublas::trans(matrix_type())) transpose_type;

	typedef boost::mpl::if_<
	    ao_major, const matrix_type&, transpose_type>::type type;
	typedef boost::mpl::if_<
	    ao_major, column_type, row_type>::type vector_type;

	void resize_ao(size_t f, size_t p) {
	    resize(ao_, f, p);
	}
	void resize_dq(size_t f, size_t p) {
	    resize(dx_, f, p);
	    resize(dy_, f, p);
	    resize(dz_, f, p);
	}

	matrix_type& mo() { return mo_; }
	column_type mo(int i) { return column(mo_,i); }

	type ao() { return transform(ao_); }
	type dx() { return transform(dx_); }
	type dy() { return transform(dy_); }
	type dz() { return transform(dz_); }

	vector_type ao(int i) { return vector(ao_,i); }
	vector_type dx(int i) { return vector(dx_,i); }
	vector_type dy(int i) { return vector(dy_,i); }
	vector_type dz(int i) { return vector(dz_,i); }

	matrix_type& tmp() { return tmp_; }

    private:
	matrix_type ao_, mo_;
	matrix_type dx_, dy_, dz_;
	matrix_type tmp_;

	static void resize(matrix_type& m, size_t f, size_t p) {
	    if (ao_major::value) m.resize(f,p);
	    else m.resize(p,f);
	}
	static column_type column(matrix_type& m, int i) {
	    return boost::numeric::ublas::column(m,i);
	}
	static row_type row(matrix_type& m, int i) {
	    return boost::numeric::ublas::row(m,i);
	}

	static type transform(const matrix_type &m) {
	    return transform(m, ao_major());
	}
	static const matrix_type& transform(const matrix_type &m, boost::mpl::true_) {
	    return m;
	}
	static transpose_type transform(const matrix_type &m, boost::mpl::false_) {
	    return boost::numeric::ublas::trans(m);
	}

	static vector_type vector(matrix_type& m, int i) {
	    return vector(m, i, ao_major());
	}
	static column_type vector(matrix_type &m, size_t i, boost::mpl::true_) {
	    return boost::numeric::ublas::column(m,i);
	}
	static row_type vector(matrix_type &m, size_t i, boost::mpl::false_) {
	    return boost::numeric::ublas::row(m,i);
	}



    };
    
}
}

#endif // DFT_GRID_VECTORS_HPP
