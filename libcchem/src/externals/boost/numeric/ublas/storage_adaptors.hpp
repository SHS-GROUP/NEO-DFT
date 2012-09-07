//
//  Copyright (c) 2000-2006
//  Joerg Walter, Mathias Koch, Michael Stevens, Gunter Winkler
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_STORAGE_ADAPTORS_H
#define BOOST_UBLAS_STORAGE_ADAPTORS_H

#include <algorithm>

#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp>


namespace boost {
namespace numeric {
namespace ublas {


  /** \brief gives read only access to a chunk of memory.
   *
   * This class partially models the storage concept. Only
   * the immutable interface is provided.
   *
   * example:
   * <code>
   *    T data[dim];
   *    // ... fill data ...
   *    typedef readonly_array_adaptor<T> a_t;
   *    typedef vector<const T, a_t>      v_t;
   *    a_t pa(dim, &(data[0]));
   *    v_t v(dim, pa);
   *
   * </code>
   */

    template<class T, class P = const T*>
    class const_array_ref_adaptor:
        public storage_array<const_array_ref_adaptor<T> > {
        typedef const_array_ref_adaptor<T> self_type;
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        typedef const T &const_reference;
        typedef const T *const_pointer;
    public:

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_array_ref_adaptor(P data = 0) {
	    resize(0, data);
        }
        BOOST_UBLAS_INLINE
        const_array_ref_adaptor(size_type size, P data) {
	    resize(size, data);
        }
        BOOST_UBLAS_INLINE
        ~const_array_ref_adaptor () {
        }

        const_array_ref_adaptor (const const_array_ref_adaptor& rhs) {
	    resize(rhs.size_, rhs.data_);
	}

	typedef const_array_ref_adaptor constructor;

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            BOOST_UBLAS_CHECK (capacity_ < size, bad_size ());
            size_ = size;
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size, P data) {
	    capacity_ = size;
            size_ = size;
            data_ = data;
        }

        // Random Access Container
        BOOST_UBLAS_INLINE
        size_type max_size () const {
            return capacity_;
        }
        
        BOOST_UBLAS_INLINE
        bool empty () const {
            return size_ == 0;
        }
            
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Iterators simply are pointers.
        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        // this typedef is used by vector and matrix classes
        typedef const_pointer iterator;

        // Reverse iterators
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

        BOOST_UBLAS_INLINE
	void swap(const_array_ref_adaptor &other) {
	    std::swap(size_, other.size_);
	    std::swap(capacity_, other.capacity_);
	    std::swap(data_, other.data_);
	}

    protected:
        size_type size_;
	size_type capacity_;
	P data_;
    private:
	const_array_ref_adaptor& operator=(const const_array_ref_adaptor &a) {
	}
    };


    template<class T>
    class array_ref_adaptor:
	public const_array_ref_adaptor<T, T*>,
	public storage_array<array_ref_adaptor<T> >
    {
        typedef array_ref_adaptor<T> self_type;
	typedef const_array_ref_adaptor<T, T*> base_type;
    public:
	typedef typename base_type::size_type size_type;
        typedef T &reference;
        typedef T *pointer;
        typedef const T *const_pointer;
    public:

        // Construction and destruction
        BOOST_UBLAS_INLINE
        array_ref_adaptor (T* data) : base_type(data) {}

        BOOST_UBLAS_INLINE
        array_ref_adaptor (size_type size, pointer data)
	    : base_type(size, data) {}
        BOOST_UBLAS_INLINE
        ~array_ref_adaptor () {}

        array_ref_adaptor (const array_ref_adaptor& rhs)
	    : base_type(rhs.size_, rhs.data_) {}

        // Element access
	using base_type::operator[];

        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            BOOST_UBLAS_CHECK (i < this->size_, bad_index ());
            return this->data_ [i];
        }

        // Iterators simply are pointers.
        typedef pointer iterator;
        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return this->data_;
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return this->data_ + this->size_;
        }

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return this->data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return this->data_ + this->size_;
        }

        // Reverse iterators
        typedef std::reverse_iterator<iterator> reverse_iterator;

        BOOST_UBLAS_INLINE
        reverse_iterator rbegin () {
            return reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator rend () {
            return reverse_iterator (begin ());
        }
    private:
	array_ref_adaptor& operator=(const array_ref_adaptor &a) {
	}
    };


    template<class T>
    class reserve_array :
	public storage_array<reserve_array<T> >,
	public unbounded_array<T>
    {
        typedef unbounded_array<T> base_type;
    public:
        typedef std::size_t size_type;
        typedef T value_type;
        typedef const T *const_pointer;
    public:
        // Construction and destruction
        reserve_array() {}
	reserve_array(size_type size) : base_type(size) {}
        reserve_array(const reserve_array& rhs) {
	    *this = rhs;
	}
	reserve_array& operator=(const reserve_array &a) {
	    if (this != &a) {
		resize(a.size());
		std::copy(a.begin(), a.end(), this->begin());
	    }
	    return *this;
	}

        // Resizing
        BOOST_UBLAS_INLINE
        void resize(size_type size) {
	    if (this->size() < size) base_type::resize(size);
        }
    };


    // // template<typename T>
    // // vector<T, array_adaptor<T> >
    // // make_vector(size_t size, T *data) {
    // //     typedef array_adaptor<T> a_t;
    // // 	return vector<T, a_t>(size, a_t(size, data));
    // // }

    // /** \brief converts a chunk of memory into a (readonly) usable ublas vector.
    //  *
    //  * <code>
    //  *   double data[10]
    //  *   vector<double> v(5);
    //  *   matrix<double> m(5,10);
    //  *   v = prod(m, make_vector_from_pointer(10, &(data[0])));
    //  * </code>
    //  */
    // template <class T>
    // vector<const T, const_array_ref_adaptor<T> >
    // make_vector(const size_t size, const T * data)
    // {
    //     typedef const_array_ref_adaptor<T> a_t;
    //     typedef vector<const T, a_t>      v_t;
    //     return v_t(size, a_t(size, data));
    // }

    // /** \brief converts a chunk of memory into a (readonly) usable dense matrix.
    //  *
    //  * <code>
    //  *   double data[50]
    //  *   vector<double> v(5);
    //  *   vector<double> x(10);
    //  *   matrix<double> m(5,10);
    //  *   v = prod(make_matrix_from_pointer(5, 10, &(data[0])), x);
    //  * </code>
    //  */
    // template <class LAYOUT, class T>
    // matrix<const T, LAYOUT, const_array_ref_adaptor<T> >
    // make_matrix(const size_t size1, const size_t size2, const T * data)
    // {
    //     typedef const_array_ref_adaptor<T> a_t;
    //     typedef matrix<const T, LAYOUT, a_t>      m_t;
    //     return m_t(size1, size2, a_t(size1*size2, data));
    // }

    // template <class F1, class LAYOUT, class T>
    // symmetric_matrix<const T, F1, LAYOUT, const_array_ref_adaptor<T> >
    // make_matrix(const size_t size1, const T * data)
    // {
    //     typedef const_array_ref_adaptor<T> a_t;
    //     typedef symmetric_matrix<const T, F1, LAYOUT, a_t> m_t;
    //     return m_t(size1, a_t((size1*size1 + size1)/2, data));
    // }


    // // default layout: row_major
    // template <class T>
    // matrix<const T, row_major, const_array_ref_adaptor<T> >
    // make_matrix(const size_t size1, const size_t size2, const T * data)
    // {
    // 	  return make_matrix<row_major>(size1, size2, data);
    // }

    // /** \brief converts a C-style 2D array into a (readonly) usable dense matrix.
    //  *
    //  * <code>
    //  *   double data[5][10];
    //  *   vector<double> v(5);
    //  *   vector<double> x(10);
    //  *   matrix<double> m(5,10);
    //  *   v = prod(make_matrix_from_pointer(data), x);
    //  * </code>
    //  */
    // template <class T, size_t M, size_t N>
    // matrix<const T, row_major, const_array_ref_adaptor<T> >
    // make_matrix(const T (&array)[M][N])
    // {
    //     typedef const_array_ref_adaptor<T> a_t;
    //     typedef matrix<const T, row_major, a_t>      m_t;
    //     return m_t(M, N, a_t(M*N, array[0]));
    // }
    // template <class T, size_t M, size_t N>
    // matrix<const T, row_major, const_array_ref_adaptor<T> >
    // make_matrix(const T (*array)[M][N])
    // {
    //     typedef const_array_ref_adaptor<T> a_t;
    //     typedef matrix<const T, row_major, a_t>      m_t;
    //     return m_t(M, N, a_t(M*N, (*array)[0]));
    // }

}}}

#endif
