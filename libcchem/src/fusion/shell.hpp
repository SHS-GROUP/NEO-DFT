#ifndef CCHEM_FUSION_SHELL_HPP
#define CCHEM_FUSION_SHELL_HPP

#include "basis/shell.hpp"

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include <boost/fusion/include/vector.hpp>

namespace cchem {
namespace fusion {

    template<int L, int M, int N>
    struct function {
	static const int l = L;
	static const int m = M;
	static const int n = N;
    };

    namespace detail {

	template<int X>
	struct hex_to_function {
	    typedef function<
		(X >> 8) & 0xF, (X >> 4) & 0xF, (X >> 0) & 0xF> type;
	};

	template<int I>
	struct index_to_shell {
	    static const basis::shell::type type =
		basis::shell::type(I-2 
#ifdef BASIS_SHELL_SP
				   -1
#endif
				   );
	};

    }
	
    template<basis::shell::type>
    struct shell;

#define CCHEM_FUSION_SHELL_FUNCTION(r, data, T)		\
    detail::hex_to_function<BOOST_PP_CAT(0x,T)>::type   

#define CCHEM_FUSION_SHELL(R, SP, FUNCTIONS)			\
    template<>							\
    struct shell<detail::index_to_shell<R>::type> {		\
	static const basis::shell::type type =			\
	    detail::index_to_shell<R>::type;			\
	static const int L = (type > 0) ? type : -type;		\
	typedef boost::fusion::vector				\
	<BOOST_PP_SEQ_ENUM					\
	 (BOOST_PP_SEQ_TRANSFORM(CCHEM_FUSION_SHELL_FUNCTION, 	\
				 (), FUNCTIONS))		\
	 > functions;						\
    };

    // BASIS_SHELL_FUNCTIONS
    BOOST_PP_SEQ_FOR_EACH(CCHEM_FUSION_SHELL, (), BASIS_SHELL_FUNCTIONS); 


} // namespace fusion
} // namespace cchem

#endif // CCHEM_FUSION_SHELL_HPP
