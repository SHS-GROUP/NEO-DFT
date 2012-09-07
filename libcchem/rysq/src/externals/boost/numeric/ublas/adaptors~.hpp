#ifndef UBLAS_ADAPTERS_HPP
#define UBLAS_ADAPTERS_HPP

// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/symmetric.hpp>

// namespace boost {
// namespace numeric {
// namespace ublas {

//     template<typename T, typename L>
//     struct matrix_adapter
// 	: ublas::matrix<T, L, ublas::array_adaptor<T> >
//     {
// 	typedef ublas::array_adaptor<T> array_type;
// 	typedef ublas::matrix<T, L, array_type> base;

// 	matrix_adapter(size_t size1, size_t size2, T *data)
// 	    : base(size1, size2, array_type(0, pointer()))
// 	{
// 	    array_type data_(size1*size2, data);
// 	    this->data().swap(data_);	    
// 	}

// 	base& operator()() { return *this; }
// 	const base& operator()() const { return *this; }

// 	template<class E>
// 	matrix_adapter& operator+=(const ublas::matrix_expression<E> &e) {
// 	    base::plus_assign(e);
// 	    return *this;
// 	}

// 	template<class E>
// 	matrix_adapter& operator-=(const ublas::matrix_expression<E> &e) {
// 	    base::minus_assign(e);
// 	    return *this;
// 	}

//     private:
// 	typedef typename array_type::pointer pointer;
//     };

//     template<typename T, class O, class L>
//     struct symmetric_matrix_adapter
// 	: ublas::symmetric_matrix<T, O, L, ublas::array_adaptor<T> >
//     {
// 	typedef ublas::array_adaptor<T> array_type;
// 	typedef ublas::symmetric_matrix<T, O, L, array_type> base;

// 	symmetric_matrix_adapter(size_t size1, T *data)
// 	    : base(size1, array_type(0, pointer()))
// 	{
// 	    array_type data_((size1*size1 + size1)/2, data);
// 	    this->data().swap(data_);
// 	    assert(this->data().begin() == data);
// 	}

// 	template<class E>
// 	symmetric_matrix_adapter& operator+=(const ublas::matrix_expression<E> &e) {
// 	    base::plus_assign(e);
// 	    return *this;
// 	}

// 	template<class E>
// 	symmetric_matrix_adapter& operator-=(const ublas::matrix_expression<E> &e) {
// 	    base::minus_assign(e);
// 	    return *this;
// 	}

//     private:
// 	typedef typename array_type::pointer pointer;
//     };

// }
// }
// }


#endif // UBLAS_ADAPTERS_HPP
