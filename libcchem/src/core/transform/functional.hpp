#ifndef CORE_TRANSFORM_FUNCTIONAL_HPP
#define CORE_TRANSFORM_FUNCTIONAL_HPP

#include "core/transform/transform.hpp"
#include "core/transform/contract1.hpp"
#include "core/integral/generator.hpp"

namespace transform {

    struct apply {

    	template<class T0, class T1, class T2, class T3 = void>
    	struct result;

    	template<class E1, class E2, class T> 
    	struct result<E1, E2, integral::generator<2>, T> {
	    typedef void type;
	};

    	template<class E1, class E2> 
	void operator()(const E1 &C1, const E2 &C2,
			integral::generator<2> integral,
			Transform<true> &T) const {
	    size_t size[] = { C1.size1(), C2.size1(),
			      integral.at(0).size(), integral.at(1).size() };
	    T.resize(size);
	    detail::apply(C1, C2, integral, T);
    	}

    	// template<class E1, class E2, class A1, class A2>
    	// struct result { typedef A2& type; };

    	// template<class E1, class A1, class A2>
    	// struct result<E1,A1,A2> { typedef A2& type; };

    	// template<class E, class A>
    	// A& operator()(const E &C1, const A &m1, A &m2) const {
    	//     using namespace boost::numeric::ublas;
    	//     m2 = prod(m1, trans(C1));
	//     return m2;
    	// }

    	// template<class E, class A>
    	// A& operator()(const E &C1, const E &C2, A &m1, A &m2) const {
    	//     //std::cout << "apply-2" << std::endl;
    	//     size_t size1 = m2.size1();
    	//     size_t size2 = m2.size2();
    	//     using namespace boost::numeric::ublas;
    	//     m2 = prod(C1, m1);
    	//     m1 = prod(m2, trans(C2));
    	//     m2.resize(size1, size2, false);
    	//     return m1;
    	// }

    };


}


#endif // CORE_TRANSFORM_FUNCTIONAL_HPP
