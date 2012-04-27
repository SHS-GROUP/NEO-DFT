// #include "tensor/tensor.hpp"
#include "typename.hpp"

#include <boost/spirit/home/phoenix/core/argument.hpp>
#include <boost/spirit/home/phoenix/statement/for.hpp>
#include <boost/spirit/home/phoenix/core/compose.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/home/phoenix/scope/let.hpp>
#include <boost/spirit/home/phoenix/scope/local_variable.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>

#include <boost/fusion/include/cons.hpp>
#include <boost/fusion/include/make_cons.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/as_vector.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_same.hpp>

namespace interval {


    namespace detail {

	template<int N> struct variable_key {};

#define RANGE_VARIABLE(N)						\
    boost::phoenix::actor<boost::phoenix::local_variable<detail::variable_key<N> > >

	template<int N>
	struct variable_step {
	    typedef BOOST_TYPEOF((++(RANGE_VARIABLE(N)()))) type;
	};

	template<int N, class F = typename variable_step<N>::type>
	struct variable {
	    variable(const F & f) : f(f) {}
	    F f;
	};

	template<int N>
	struct variable<N> {
	    typename variable_step<N>::type f;
	};
	    

    }

    namespace variables {
	RANGE_VARIABLE(0) i;
	RANGE_VARIABLE(1) j;
	RANGE_VARIABLE(2) k;
    }

    namespace detail {

	namespace mpl = boost::mpl;
	namespace phoenix = boost::phoenix;
	namespace fusion = boost::fusion;

	/// undefined or implicit object
	struct undefined {
	    // typedef fusion::result_of::vector_tie<>::type sequence_type;
	    struct sequence { typedef mpl::void_ type; };
	    sequence::type vector_;
	};

	template<typename T>
	struct is_implicit : mpl::false_ {};

	template<>
	struct is_implicit<undefined> : mpl::true_ {};

	template<class T, class U>
	struct bound;

	template<class T, class U>
	struct bound<phoenix::actor<T>, U> {
	    typedef phoenix::actor<T> type;
	    template<class  E>
	    struct result {
		typedef typename T::template result<E>::type type;
	    }; 
	    static type apply(const type &t) {
		return t;
	    }
	};


	template<class T, class U>
	struct bound {
	    typedef phoenix::actor<phoenix::value<T> > type;
	    template<class E>
	    struct result {
		typedef typename
		bound<typename bound::type,U>::template result<E>::type type;
	    }; 
	    static type apply(const T& t){
		return phoenix::val(t);
	    }
	};

	template<class  U>
	struct bound<undefined, U> {
	    typedef typename bound<U,U>::type type;
	    template<class E>
	    struct result {
		typedef typename
		bound<typename bound::type,U>::template result<E>::type type;
	    }; 
	    static type apply(const undefined&){
		return bound<U,U>::apply(U());
	    }
	};

	template<typename T>
	struct is_value : mpl::false_ {};

	template<typename T>
	struct is_value<boost::phoenix::value<T> > : mpl::true_ {};

	template<class I, class C, class E>
	struct deduce_variable_type {
	    typedef typename boost::remove_reference<
		typename fusion::result_of::back<C>::type>::type cond_type;
	    typedef typename mpl::if_<is_implicit<I>, cond_type, I>::type init_type;
	    typedef typename mpl::eval_if<is_value<init_type>,
					  mpl::identity<init_type>,
					  phoenix::as_actor_base<init_type>
					  >::type init_actor_base;
	    typedef typename init_actor_base::template result<E>::type type;
	};

	template<class E, class S, typename T>
	struct expression_environment {
	    typedef typename boost::remove_reference<
		typename fusion::result_of::front<S>::type>::type F;
	    typedef typename F::variable K;
	    typedef phoenix::scoped_environment<
		E, E, fusion::vector<T>, 
		phoenix::detail::map_local_index_to_tuple<K> > environment_;
	    typedef typename expression_environment<
		environment_, typename fusion::result_of::pop_front<S>::type, T
		>::type type;
	};

	struct expression_eval {

	    template<class E, class S>
	    struct result {
		// typedef typename boost::remove_reference<
		//     typename fusion::result_of::front<
		// 	typename S::template result<E>::type
		// 	>::type>::type G;
		// typedef typename G::template result<E>::type type;
		typedef void type;
	    };

	    template<class S, class F>
	    struct actor {
		typedef typename phoenix::as_composite<
		    expression_eval,
		    typename fusion::result_of::push_back<
			const typename S::cons_type, F>::type
		    >::type composite;
		typedef phoenix::actor<composite> type;
	    };

	    template<class S, class F>
	    static typename actor<S, F>::type
	    compose(const S &s, const F& f) {
		return phoenix::compose<expression_eval>(fusion::push_back(s.cons, f));
	    }

	    template<typename R, class E, class S>
	    static void eval(const E &e, const phoenix::value<S> &v){
		BOOST_AUTO(s, fusion::as_vector(v.val));
		eval_(e, fusion::pop_back(s), fusion::back(s));
	    } 

	    template<class E, class S, class G>
	    static void eval(const E &e, const phoenix::value<S> &v, const G& g) {
		BOOST_AUTO(s, fusion::as_vector(v.val));
		// std::cout << TYPENAME(g(fusion::back(s))) << std::endl;
		eval_(e, fusion::pop_back(s), g(fusion::back(s)));
	    } 

	private:
	
	    template<class E, class S, class F>
	    static void eval_(const E &e, const S &s, const F &f) {
		typedef BOOST_TYPEOF(fusion::front(s)) S0;
		typedef typename deduce_variable_type<
		typename S0::init_type, typename S0::cond_type, E>::type T;
		eval_<T>(e, s, f);
	    }	    

	    template<class T, class E, class S, class F>
	    static typename boost::disable_if<fusion::result_of::empty<S> >::type
	    eval_(const E &e, const S &s, const F &f) {
	    	BOOST_AUTO(t, fusion::back(s));
		using namespace phoenix;
	    	eval_<T>(e, fusion::pop_back(s),
	    	      let(t.var = t.template init<T>())
	    		 [ for_(t.var, t.template cond<T>(), t.step()) [ f ] ]);
	    }

	    template<class T, class E, class S, class F>
	    static typename boost::enable_if<fusion::result_of::empty<S> >::type
	    eval_(const E &e, const S&, const F &f) {
		// std::cout << TYPENAME(f) << std::endl;
		f.eval(e);
	    }
	};

	template<int N, class L, class R, class Step>
	struct expression_node {

	    typedef L init_type;
	    typedef R cond_type;
	    typedef Step step_type;

	    L init_;
	    R cond_;
	    Step step_;

	    typedef RANGE_VARIABLE(N) variable;
	    variable var;

	    expression_node(const L &lower, const R& upper, const Step &step)
		: init_(lower), cond_(upper), step_(step) {}

	    template<typename T>
	    typename bound<L,T>::type init() const {
		return bound<L,T>::apply(init_);
	    }

	    template<typename T>
	    typename bound<R,T>::type cond() const {
		return bound<R,T>::apply(cond_);
	    }

	    Step step() const { return step_; }

	};

	struct innermost {
	    typedef fusion::nil cons_type;
	    cons_type cons;
	};

	template<int N, class Init, class Cond, class Step, class Inner = innermost>
	struct expression {
	    typedef expression self_type;
	    typedef expression_node<N, Init, Cond, Step> node;
	    typedef fusion::cons<node, typename Inner::cons_type> cons_type;
	    cons_type cons;

	    template<class E, class F>
	    struct result {
		typedef typename expression_environment<E,self_type,int>::type E_;
		typedef typename F::template result<E_>::type type;
	    };

	    template<class E>
	    struct extend {
		// assert E is expression
		typedef expression<N, Init, Cond, Step, E> type;
	    };

	    expression(const Init& init, const Cond& cond,
		       const Step& step, const Inner &inner)
		: cons(node(init, cond, step), inner.cons) {}

	    expression(const expression<N, Init, Cond, Inner> &e)
		: cons(e.cons) {}

	    const Init& init() const { return cons.car.init_; }
	    const Cond& cond() const { return cons.car.cond_; }
	    const Step& step() const { return cons.car.step_; }

	    template<class E>
	    typename extend<E>::type
	    operator()(const E &e) const {
		return typename extend<E>::type(init(), cond(), step(), e);
	    }


	    template<class F>
	    typename expression_eval::actor<expression, F>::type
	    operator[](const F &f) const {
	    	return expression_eval::compose(*this, f);
 	    }

	};

	template<int N, class I, class C, class S>
	expression<N, I, C, S>
	make_expression(const I &init, const C& cond, const S &step) {
	    return expression<N, I, C, S>(init, cond, step, innermost());
	}

	template<int N, class  I, class C, class S>
	expression<N, I, C, S>
	make_expression(const expression<N, I, undefined, S> &e, const C& cond) {
	    return expression<N, I, C, S>(e.init(), cond, e.step(), innermost());
	}

	template<class E, int N, class L, class R, class S>
	typename expression<N, L, R, S>::template extend<E>::type
	operator,(const E &e, const expression<N, L, R, S> &d)  {
	    return d(e);
	}

// #define RANGE_VARIABLE(N) phoenix::actor<phoenix::argument<N> >
#define ACTOR_COMPOSITE(O,L,R)						\
	phoenix::actor<typename phoenix::as_composite<O, L, R>::type>      
	

#define LOWER_BOUND_OPERATOR(op, eval, param)				\
	template <typename T, int N, class F>				\
	expression<N, param, undefined, F>				\
	operator op (const param& left, const variable<N,F> &i) {	\
	    return make_expression<N>(left, undefined(), i.f);		\
	}


#define UPPER_BOUND_OPERATOR_VARIABLE(op, eval, param)			\
	template <typename T, int N, class F>				\
	expression<N, undefined,					\
		   ACTOR_COMPOSITE(eval, RANGE_VARIABLE(N), param), F>	\
	operator op (const variable<N,F> &i, const param& e) {		\
 	    return make_expression<N>(undefined(),			\
 				      (RANGE_VARIABLE(N)() op e), i.f);	\
 	}

#define UPPER_BOUND_OPERATOR_EXPRESSION(op, eval)			\
	template <typename T, int N, class E, class S>			\
	expression<N, E, ACTOR_COMPOSITE(eval, RANGE_VARIABLE(N), T), S> \
	operator op (const expression<N, E, undefined, S> &left,	\
		     const T& right) {					\
	    return make_expression(left,				\
				   (RANGE_VARIABLE(N)() op right));	\
	}

#define UPPER_BOUND_OPERATOR(op, eval)					\
	UPPER_BOUND_OPERATOR_VARIABLE(op, eval, T)			\
	UPPER_BOUND_OPERATOR_VARIABLE(op, eval, phoenix::actor<T>)	\
	UPPER_BOUND_OPERATOR_EXPRESSION(op, eval)

	LOWER_BOUND_OPERATOR( <= , phoenix::less_equal_eval, T)
	LOWER_BOUND_OPERATOR( <= , phoenix::less_equal_eval, phoenix::actor<T>)

	UPPER_BOUND_OPERATOR( < , phoenix::less_eval)
	UPPER_BOUND_OPERATOR( <= , phoenix::less_equal_eval)

    }

    namespace detail {


	struct reduce_eval {
	    template<class B, class A, typename T>
	    struct tuple {
		tuple(const B& b, const A& a, const T& t) : b(b), a(a), t(t) {}
		B b; A a; T t;
	    };    

	    template<class B, class E, class T>
	    struct actor {
		typedef phoenix::actor<
		    typename phoenix::as_composite<
			reduce_eval, reduce_eval::tuple<B,E,T>
			>::type> type;
	    };

	    template<class B, class E, class T>
	    static typename actor<B, E, T>::type
	    compose(const B& b, const E& e, const T& t) {
		return phoenix::compose<reduce_eval>(tuple<B, E, T>(b,e,t));
	    }

	    // typedef typename fusion::result_of::front<A>::type expression; 

	    template<class B, class R, class F>
	    struct curry_eval {
		typedef mpl::true_ no_nullary;
		template<class E>
		struct result { typedef void type; };
	    	curry_eval(const B& b, const R &r, const F &f)
		    : b(b), r(r), f(f) {}
	    	template<class E>
	    	void eval(const E &e) const {
		    BOOST_AUTO(const &v, f.eval(e));
		    b(r(), v);
		}
	    	B b; R r; F f;
	    };

	    template<class B, typename R>
	    struct curry_gen {
	    	curry_gen(const B& b, const R &r) : b(b), r(r) {}
	    	template<class F>
	    	phoenix::actor<curry_eval<B, R, F> >
		operator()(const F &f) const {
		    return curry_eval<B, R, F>(b, r, f);
		}
	    	B b; R r;
	    };

	    template<class B, typename R>
	    curry_gen<B, R> static curry(const B& b, const R &r) {
		return curry_gen<B, R>(b,r);
	    }

	    template<class E, class T>
	    struct result { typedef double type; };

	    template<typename R, class E, class T>
	    static R eval(const E& e, const phoenix::value<T>& v) {
		const T& tuple = v.val;
		R r(v.val.t);
		expression_eval::eval(e, fusion::front(tuple.a),
				      curry(tuple.b, phoenix::ref(r)));
		return r;
	    }
	};

	template<class B, class T, class E = void>
	struct reduce_gen;

	template<class B, class T>
	struct reduce_gen<B,T> {
	    
	    reduce_gen(const B &b, const T &t) : b(b), t(t){}

	    template<class E>
	    typename reduce_eval::actor<B,E,T>::type
	    operator[](const E& e) const {
		return reduce_eval::compose(b, e, t);
	    }

	    template<int N, class L, class R, class S>
	    reduce_gen<B, T, expression<N, L, R, S> >
	    operator()(const expression<N, L, R, S>& e) const {
		return reduce_gen<B, T, expression<N, L, R, S> >(b, t, e);
	    }
	    B b; T t;
	};


	template<class B, class T, int N, class L, class R, class C> 
	struct reduce_gen<B, T, expression<N, L, R, C> > {
	    typedef expression<N, L, R, C> E;
	    reduce_gen(const B &b, const T &t, const E &e)
		: b(b), t(t), e(e) {}
	    B b; T t; E e;
	    template<class E2>
	    reduce_gen<B, T, typename E::template extend<E2>::type>
	    operator()(const E2 &e2) const {
		typedef typename E::template extend<E2>::type E12;
		return reduce_gen<B, T, E12>(e(e2));
	    }

	    template<class F>
	    typename reduce_eval::actor<
		B, typename expression_eval::actor<E,F>::type, T>::type
	    operator[](const F &f) const {
		return reduce_gen<B,T>(b, t)[e[f]];
	    }
	};


	typedef BOOST_TYPEOF(phoenix::arg_names::_1 += phoenix::arg_names::_2)
	    sum_functor;
	typedef reduce_gen<sum_functor, int> sum_gen;

    }

    template<typename T, class G, class V>
    detail::reduce_gen<G,V> reduce(const G &g, const V &v = V()) {
	return detail::reduce_gen<G,V>(g,v);
    }

    template<class G, class V>
    detail::reduce_gen<G,V> reduce(const G &g, const V &v = V()) {
	return detail::reduce_gen<G,V>(g,v);
    }
    
    detail::sum_gen sum(detail::sum_functor(), 0);


    template<int N>
    detail::variable<N> range( const RANGE_VARIABLE(N)& ) {
	return detail::variable<N>();
    }

    template<int N, class F>
    detail::variable<N,F> range(const RANGE_VARIABLE(N)&, const F &f) {
	return detail::variable<N,F>(f);
    }


}    


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/spirit/include/phoenix_function.hpp>


template<typename T>
struct matrix {
    typedef boost::numeric::ublas::matrix<T> ublas_matrix;
    matrix(ublas_matrix &m) : m(m){}
    ublas_matrix & m;
    template<typename, typename>
    struct result { typedef T& type; };
    T& operator()(int i, int j) const { return m(i,j); }
};

int main(){
    

    // formula domain language rough prototype
    // stuff in brackets can be any valid phoenix lambda expression

    // physicist notation, lower bound may be implicit
    // (range(i) <= j, 3 <= range(j) < 4)[std::cout << i << " " << j << std::endl];

    // implicit physicist notation, not done yet
    //(range(i) < range(j) < 4)[...]

    using interval::range;
    using namespace interval::variables;


    // programmer notation, same as above , but order is reversed
    // (0 < range(j) < 4)[std::cout << boost::phoenix::val(1) << std::endl];
    (range(i) < 2, 1 <= range(j) < 4)[std::cout << i << " " << j << std::endl]();
    

    using namespace boost::phoenix;
    using namespace boost::phoenix::arg_names;
    using interval::reduce;
    using interval::sum;

    // std::cout << interval::sum<double>(1 <= range(i) < 4000)[1.0/i]() << std::endl;
    // (range(i) < 2)[std::cout << i<< std::endl]();
    std::cout << reduce(_1 /= _2, 1)[(1 <= range(i) < 40)[i] ]() << std::endl;

    std::cout << sum(1.0 <= range(i,i += 0.2) < 40.0)[0.2/i]() << std::endl;

    
    actor<local_variable<int> > W;

    size_t N = 4;
    namespace ublas = boost:: numeric::ublas;
    ublas::matrix<double> a(N,N), b(N,N), c(N,N);
    function<matrix<double> > A(a), B(b), C(c);

    (range(i) < N)(range(j) < N)
	[ A(i,j) = sum(range(k) < N)[B(k,i)*C(k,j)/A(k,j)] ]();

    

    // ignore, some early prototype for lambda tensor
    // size_t N = 4;
    // tensor::tensor<4> T(N,N,N,N);

     // tensor::function<tensor::tensor<4> > T_(T);
    // (range(j) < 4)(range(i) <= j)[std::cout << T_(i,j,0,0)];

    // (range(i) < j, range(j) < N)[T_(i,j,0,0) = T_(j,i,0,0)];
    // sum(j < val(N))[T_(i,j,0,0)];

}
