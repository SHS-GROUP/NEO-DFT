#ifndef PHOENIX_TENSOR_FUNCTION_HPP
#define PHOENIX_TENSOR_FUNCTION_HPP

#include <boost/spirit/home/phoenix/function/function.hpp>

namespace boost { namespace phoenix {

template<class E, size_t N>
struct tensor_adapter;

template<class E>
struct tensor_adapter<E, 1> {
    typedef typename E::reference reference;
    template<typename>
    struct result { typedef reference type; };
    tensor_adapter(E &e): e_(e) {}
    template<typename T1>
    reference operator()(const T1 &_1) const {
	return e_(_1);
    }
    E &e_;
};

template<class E>
struct tensor_adapter<E, 2> {
    typedef typename E::reference reference;
    template<typename, typename>
    struct result { typedef reference type; };
    tensor_adapter(E &e): e_(e) {}
    template<typename T1, typename T2>
    reference operator()(const T1 &_1, const T2 &_2) const {
	return e_(_1,_2);
    }
    E &e_;
};

template<class E>
struct tensor_adapter<E, 4> {
    typedef typename E::reference reference;
    template<typename, typename, typename, typename>
    struct result { typedef reference type; };
    tensor_adapter(E &e): e_(e) {}
    template<typename T1, typename T2, typename T3, typename T4>
    reference operator()(const T1 &_1, const T2 &_2, const T3 &_3, const T4 &_4) const {
	return e_(_1,_2,_3,_4);
    }
    E &e_;
};


template<class E, size_t N, class enable = void>
struct tensor_function : function<tensor_adapter<E,N> > {
    tensor_function(E &e) : function<tensor_adapter<E,N> >(e) {}
};


} }

#endif // PHOENIX_TENSOR_FUNCTION_HPP
