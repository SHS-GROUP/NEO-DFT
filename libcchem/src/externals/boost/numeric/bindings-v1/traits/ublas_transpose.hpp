//
// Copyright Kresimir Fresl, Toon Knapen, and Karl Meerbergen 2002, 2003
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_BINDINGS_UBLAS_TRANSPOSE_HPP
#define BOOST_BINDINGS_UBLAS_TRANSPOSE_HPP

#include <boost/numeric/bindings/traits/transpose.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/mpl/bool.hpp>
#include "typename.hpp"

namespace boost {
namespace numeric {
namespace bindings {
namespace traits {

    namespace transpose {

	template<class E, class M = E>
	struct transpose_detail : mpl::false_ {
	    static const char option = NO_TRANSPOSE;
	    typedef E& no_transpose_type;
	    static no_transpose_type no_transpose(E &e) { return e; }
	};

	template<class E, class M, class F>
	struct transpose_detail<E, ublas::matrix_unary2<M, F> > : mpl::true_ {
	    typedef E expression;
	    static const char option = TRANSPOSE;
	    typedef typename E::expression_closure_type no_transpose_type;
	    static no_transpose_type
	    no_transpose(E &e) {
		return e.expression();
	    }
	};

	struct transpose_option {
	    transpose_option(char option) : data_(option) {}
	    operator char() const { return data_; }
#ifdef BOOST_NUMERIC_BINDINGS_CBLAS_ENUM_HPP
	    operator CBLAS_TRANSPOSE() const {
		if (data_ == TRANSPOSE) return CblasTrans;
		if (data_ ==  CONJUGATE) return CblasConjTrans;
		return CblasNoTrans;
	    }
	    operator bool() const { return (data_ != NO_TRANSPOSE); }
#endif
	private: const char data_;
	};


	template<class E>
	transpose_option option(const E&) {
	    return transpose_detail<E>::option;
	}

	template<class E>
	typename transpose_detail<E>::no_transpose_type
	no_transpose(E &e) {
	    return transpose_detail<E>::no_transpose(e);
	}

	template<class E>
	typename transpose_detail<const E, E>::no_transpose_type
	no_transpose(const E &e) {
	    return transpose_detail<const E, E>::no_transpose(e);
	}


	template<class E>
	transpose_option option(const ublas::matrix_expression<E>&) {
	    return transpose_detail<E>::option;
	}

	template<class E>
	typename transpose_detail<E>::no_transpose_type
	no_transpose(ublas::matrix_expression<E> &e) {
	    return  transpose_detail<E>::no_transpose(e());
	}

	template<class E>
	typename transpose_detail<const E, E>::no_transpose_type
	no_transpose(const ublas::matrix_expression<E> &e) {
	    return  transpose_detail<const E, E>::no_transpose(e());
	}
	    
    }	
		

}}}}

#endif 
