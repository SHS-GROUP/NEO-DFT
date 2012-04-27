#ifndef TENSOR_VIEW_HPP
#define TENSOR_VIEW_HPP

#include "tensor/forward.hpp"
#include "tensor/storage/storage.hpp"
#include "tensor/traits.hpp"
#include "tensor/functional.hpp"

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>


namespace tensor {

    template<class A, class R>
    struct const_tensor_view :
	detail::const_tensor_base<A>,
	expression<const_tensor_view<A,R> >
    {
    private:
	typedef const_tensor_view self_type;
	typedef detail::traits<A> traits;
	typedef detail::const_tensor_base<A> base;

    public:
	typedef R indices_type;
	typedef typename detail::result_of::keys<R>::type keys_type;
	TENSOR_ASSERT_UNIQUE_KEYS(keys_type);
	typedef A array_type;

	struct const_iterator: detail::traits<A>::const_iterator {
	    explicit
	    const_iterator(typename detail::traits<A>::const_iterator it) :
		detail::traits<A>::const_iterator(it) {}
	};

	typedef typename base::const_result_type const_result_type;

	template<class O>
	struct result : base::template result<O> {};

	using base::operator();

	const_tensor_view(const array_type &data)
	    : base(data), ranges_() {}

	const keys_type& keys() const { return keys_; }
	const indices_type& indices() const { return ranges_; }

	const_iterator begin() const {
	    return const_iterator(storage::begin(this->data()));
	}
	const_iterator end() const {
	    return const_iterator(storage::end(this->data()));
	}


    protected:
	friend class detail::functional;
	template<class M>
	const_result_type operator[](const M &indices) const {
	    return (traits::template element<const_result_type>
		    (base::data_, keys(), indices));
	}
    private:
	keys_type keys_;
	indices_type ranges_;
    };


    template<class A, class R>
    struct tensor_view :
	detail::tensor_base<A>,
	expression<tensor_view<A,R> >
    {
	typedef tensor_view self_type;
	typedef detail::tensor_base<A> base_type;

	typedef A array_type;
	typedef typename detail::traits<A>::value_type value_type;
	typedef typename detail::traits<A>::reference reference;
	typedef typename detail::traits<A>::const_reference const_reference;

	struct iterator: detail::traits<A>::iterator {
	    explicit iterator(typename detail::traits<A>::iterator it) :
		detail::traits<A>::iterator(it) {}
	};

	typedef R indices_type;
	typedef typename detail::result_of::keys<R>::type keys_type;
	TENSOR_ASSERT_UNIQUE_KEYS(keys_type);


	tensor_view(const array_type &data) : base_type(data) {
	    //ranges_();
	}

	void operator=(const value_type &s) {
	    assign(*this, s);
	}
	void operator*=(const value_type &s) {
	    times_assign(*this, s);
	}

	template<class E>
	void operator=(const expression<E> &e) {
	    assign(*this, e);
	}
	template<class E>
	void operator+=(const expression<E> &e) {
	    plus_assign(*this, e);
	}
	template<class E>
	void operator-=(const expression<E> &e) {
	    minus_assign(*this, e);
	}


	template<class O>
	struct result : base_type::template result<O> {};

	using base_type::operator();

	const keys_type& keys() const { return keys_; }
	const indices_type& indices() const { return ranges_; }


	iterator begin() {
	    return iterator(storage::begin(this->data()));
	}
	iterator end() {
	    return iterator(storage::end(this->data()));
	}

	using base_type::operator[];

    protected:
	friend class detail::functional;

	template<class M>
	const_reference operator[](const M &indices) const {
	    return storage::array<array_type>::template
		element<const_reference>(base_type::data_, keys(), indices);
	}

	template<class M>
	reference operator[](const M &indices) {
	    return storage::array<array_type>::template
		element<reference>(base_type::data_, keys(), indices);
	}

    private:
	keys_type keys_;
	indices_type ranges_;

    };

}


#endif // TENSOR_VIEW_HPP
