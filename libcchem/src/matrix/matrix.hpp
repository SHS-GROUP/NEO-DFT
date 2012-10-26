#ifndef QC_UTIL_MATRIX_H
#define QC_UTIL_MATRIX_H


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <boost/mpl/void.hpp>

namespace ublas = boost::numeric::ublas;

namespace  matrix {

    namespace mpl = boost::mpl;


    template<typename T>
    struct matrix_semantics {
	typedef int index_type;
	typedef T& reference;
	typedef const T& const_reference;
        virtual ~matrix_semantics() {}
	virtual T& operator()(index_type i, index_type j) = 0;
	virtual const T& operator()(index_type i, index_type j) const = 0;
    };

    template<typename T>
    struct abstract_matrix {
	virtual T* operator()(int i, int j) = 0;
	virtual const T* operator()(int i, int j) const = 0;
    };

    template<typename T> class block_matrix;
    template<class Matrix> class meta_matrix;
    template<class Matrix> class block_meta_matrix;

}

#include "matrix/rectangular.hpp"
#include "matrix/block.hpp"
#include "matrix/meta.hpp"
#include "matrix/symmetric.hpp"

#include "matrix/ops.hpp"
#include "matrix/permute.hpp"

//typedef ublas::matrix<double> Matrix;
typedef ublas::symmetric_matrix<double> SymmetricMatrix;
typedef matrix::block_matrix<double> BlockMatrix;

#endif // _MATRIX_H
