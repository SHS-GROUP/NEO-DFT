#ifndef RYSQ_ERI_NEW_HPP_
#define RYSQ_ERI_NEW_HPP_

#include <sstream>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>

#include "rysq-core.hpp"
#include "rysq-types.hpp"
#include "kernel/eri.hpp"


namespace rysq {
namespace kernel {

  struct invalid_quartet : std::runtime_error {
    static std::string str(const Quartet<Shell> &quartet) {
      std::ostringstream os;
      os << "invalid quartet " << quartet;
      return os.str();
    }
    invalid_quartet(const Quartet<Shell> &quartet)
      : std::runtime_error(str(quartet)) {}
  };

template< template<class, class> class T >
kernel::Eri<>* new_(const Quartet<Shell> &quartet) {
    type a = type(quartet[0]);
    type b = type(quartet[1]);
    type c = type(quartet[2]);
    type d = type(quartet[3]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types) &&		\
			  c == BOOST_PP_SEQ_ELEM(2, types) &&		\
			  d == BOOST_PP_SEQ_ELEM(3, types)) {		\
									\
	typedef typename						\
	    kernel::find<BOOST_PP_SEQ_ENUM(types)>::type kernel;	\
	return new kernel(quartet, new T<kernel::bra, kernel::ket>());	\
    } else								\

    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RYSQ_TYPES)(RYSQ_TYPES)(RYSQ_TYPES)(RYSQ_TYPES)) {
	throw invalid_quartet(quartet);
    }
#undef ERI
#undef TYPES

}


} // namespace kernel
} // namespace rysq

#endif /* RYSQ_ERI_NEW_HPP_ */
