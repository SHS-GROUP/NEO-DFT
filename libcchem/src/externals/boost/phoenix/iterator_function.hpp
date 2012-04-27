#ifndef PHOENIX_ITERATOR_FUNCTION_HPP
#define PHOENIX_ITERATOR_FUNCTION_HPP


#include <boost/spirit/home/phoenix/function/function.hpp>

namespace boost { namespace phoenix {

template<class E, size_t N>
struct iterator_adapter;

template<class E>
struct iterator_adapter<E, 1> {
    typedef typename E::reference reference;
    template<typename Arg> struct result { typedef reference type; };
    iterator_adapter(E &e): e_(e) {}
    template<typename T1>
    reference operator()(const T1 &_1) const {
	return e_[_1];
    }
    E e_;
};

template<class E, size_t N = 1>
struct iterator_function : function<iterator_adapter<E,N> > {
    iterator_function(E e) : function<iterator_adapter<E,N> >(e) {}
};


} }


#endif // PHOENIX_ITERATOR_FUNCTION_HPP
