#ifndef PHOENIX_FOR_RANGE_HPP
#define PHOENIX_FOR_RANGE_HPP


#include <boost/spirit/home/phoenix/core.hpp>
#include <boost/spirit/home/phoenix/scope/local_variable.hpp>
#include <boost/spirit/home/phoenix/scope/let.hpp>
#include <boost/spirit/home/phoenix/operator.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement.hpp>

#include <boost/fusion/include/vector_tie.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/sequence/intrinsic.hpp>
#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/array.hpp>
#include <boost/fusion/include/io.hpp>

#include <stdexcept>
#include <boost/array.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/ref.hpp>
#include "typename.hpp"

namespace boost {
namespace phoenix {


    template<class Env, class Vars, class Map>
    struct let_gen_scope {
	
	typedef typename 
	boost::fusion::result_of::as_vector<
	    typename boost::fusion::result_of::transform<
		Vars, boost::phoenix::detail::initialize_local<Env>
		>::type
	    >::type locals_type;

	typedef scoped_environment<Env, Env, locals_type, Map> environment_type;

	static locals_type locals(const Env &env,
				  const let_actor_gen<Vars,Map> &gen) {
	    typedef phoenix::detail::initialize_local<Env> initialize;
	    return fusion::as_vector(fusion::transform(gen.vars, initialize(env)));
	}

	let_gen_scope(const Env &env, const let_actor_gen<Vars,Map> &gen)
	    :  locals_(locals(env, gen)),
	       environment_(env, env, locals_) {}

	environment_type environment() const {
	    return environment_;
	}
    private:
	locals_type locals_;
	environment_type environment_;
    };

    template<class Env, class Vars, class Map>
    let_gen_scope<Env,Vars,Map>
    let_scope(const Env &env, const let_actor_gen<Vars,Map> &gen) {
	return let_gen_scope<Env,Vars,Map>(env, gen);
    }


    struct for_range_eval {

	template<class G, class F>
	struct actor {
	    typedef boost::fusion::vector<G&,F> vector_type;
	    typedef boost::phoenix::actor<
		typename boost::phoenix::as_composite<
		    for_range_eval, vector_type>::type> type;
	};

	template<class G, class F>
	static typename actor<G,F>::type
	compose(G &g, const F &f) {
	    typename actor<G,F>::vector_type vector(g, f);
	    return boost::phoenix::compose<for_range_eval>(vector);
	}

	template<typename Env, class V>
	struct result {typedef void type; };

	template <typename RT, typename Env, class V>
	static void
	eval(const Env &env, const V &v) {
	    using boost::phoenix::as_actor_base;
	    using boost::fusion::at_c;
	    BOOST_AUTO(&generator, at_c<0>(v.val));
	    BOOST_AUTO(const &f, at_c<1>(v.val));

	    BOOST_AUTO(it, generator.begin(env));

	    while (true) {
	    	try {
	    	    eval_(env, as_actor_base<BOOST_TYPEOF(f)>::convert(f),
	    		  it.local(), it.next());
	    	}
	    	catch (std::out_of_range) { break; }
	    }
	}

    private:

	template<class Env, class L, class V, class F>
	static void
	eval_(const Env &env, const F &f, const L &locals, const V &values,
	      boost::mpl::true_ = boost::mpl::true_()) {
	    using namespace boost::fusion;
	    using boost::phoenix::let;
	    typename boost::phoenix::actor<BOOST_TYPEOF(back(locals))> k;
	    BOOST_AUTO(const &v, back(values));
	    static const size_t size = boost::fusion::result_of::size<L>::value;
	    eval_(
		  let_scope(env, let(k = val(v))).environment(), f,
		  pop_back(locals), pop_back(values),
		  boost::mpl::bool_<size-1>());
	}

	template<class Env, class L, class V, class F>
	static void
	eval_(const Env &env, const F &f, const L &locals, const V &values,
	      boost::mpl::false_) {
	    f.eval(env);
	}

	template<class V, class M, class B>
	static boost::phoenix::let_actor<B, V, M>
	make_let_actor(const boost::phoenix::let_actor_gen<V,M> &let, const B &base) {
	    return boost::phoenix::let_actor<B, V, M>(base, let.vars);
	}

    };


    template<class C>
    struct for_range_expression {
	typedef typename C::base_type base_type;
	typedef typename C::eval_policy_type eval_policy_type;
	typedef typename boost::mpl::if_<
	    boost::is_same<eval_policy_type, boost::phoenix::sequence_eval>,
	    base_type, boost::fusion::vector<base_type> >::type vector_type;
	static const size_t size = 
	    boost::fusion::result_of::size<vector_type>::value;

	template<size_t I>
	struct values {

	    struct at {
		template<class> struct result;
		template<class F, class A>
		struct result<F(A)> {
		    typedef typename boost::remove_reference<A>::type S;
		    BOOST_MPL_ASSERT((boost::fusion::traits::is_sequence<S>));
		    typedef typename
		    boost::fusion::result_of::value_at_c<S,I>::type type;
		};

		template<class V>
		typename result<at(V)>::type
		operator()(const V &v) const {
		    return boost::fusion::at_c<I>(v);
		}

	    };

	    typedef typename boost::fusion::result_of::as_vector<
		typename boost::fusion::result_of::transform<
		    vector_type, at>::type >::type type;

	    static type convert(const C &expression) {
		namespace fusion = boost::fusion;
		BOOST_MPL_ASSERT((fusion::traits::is_sequence<vector_type>));
		return fusion::transform(vector_type(expression),at());
	    }
	};

    };


    template<class G, class Env>
    struct for_range_iterator {

	typedef typename boost::phoenix::as_actor_base<
	    typename G::expression_type>::type expression_type;
	typedef for_range_expression<expression_type> expression_traits;
	typedef typename expression_traits::template values<0>::type local_types;
	typedef typename expression_traits::template values<1> ranges;
	typedef typename ranges::type range_types;

	typedef typename
	boost::fusion::result_of::value_at_c<range_types,0>::type range0_type;
	typedef typename range0_type::template value<Env>::type value0_type;
	typedef boost::array<value0_type, expression_traits::size> value_type;

	for_range_iterator(G &generator, const Env &env)
	    : generator_(generator), env_(env),
	      range_(ranges::convert(generator_.expression()))
	{
	    state_ = true;
	    position_ = 0;
	    initialize(env_, zero);
	}

	const local_types& local() const { return local_; }

	value_type next() {
	    size_t forward =  generator_.advance() - position_;
	    position_ += forward;
	    while (state_ && forward--) {
		bool state = advance(env_, zero);
		state_ = state && state_; 
	    }
	    if (!state_) throw std::out_of_range("");
	    return value_;
	}

    private:

	template<class Index>
	typename boost::fusion::result_of::at<value_type, Index>::type
	value(Index) {	
	    return boost::fusion::at<Index>(value_);
	}

	template<class Index>
	typename boost::fusion::result_of::at<range_types, Index>::type
	range(Index) {	
	    return boost::fusion::at<Index>(range_);
	}

	template<class Index>
	boost::phoenix::actor<
	    typename boost::fusion::result_of::value_at<local_types, Index>::type>
	local(Index) const {	
	    typedef BOOST_TYPEOF(boost::fusion::at<Index>(local_)) local_type;
	    return boost::phoenix::actor<local_type>();
	}


	boost::mpl::int_<0> zero;
	typedef boost::mpl::int_<expression_traits::size> stop_index; 

	template<class OuterEnv, int I>
	void initialize(const OuterEnv &outer, boost::mpl::int_<I> index) {
	    boost::mpl::int_<I+1> next;
	    BOOST_AUTO(&v, value(index));
	    v = range(index).start.eval(outer);
	    using namespace boost::phoenix;
	    BOOST_AUTO(k, local(index));
	    initialize(let_scope(outer, let(k = v)).environment(), next);
	}

	template<class OuterEnv>
	void initialize(const OuterEnv &outer, stop_index) {}

	template<class OuterEnv, int I>
	bool advance(const OuterEnv &outer, boost::mpl::int_<I> index) {
	    boost::mpl::int_<I+1> next;

	    namespace fusion = boost::fusion;
	    BOOST_AUTO(&v, value(index));
	    BOOST_AUTO(k, local(index));

	    using namespace boost::phoenix;
	    // std::cout << v << std::endl;
	    BOOST_AUTO(const &env, let_scope(outer, let(k = val(v))).environment());
	    BOOST_AUTO(const &stop, range(index).stop.eval(env));

	    bool state = advance(env, next);
	    if (!state) {
		//std::cout << "state "  << v << " " << stop << std::endl;
		if (v < stop) v += range(index).step.eval(env);
		if (!(v < stop)) return false;
		initialize(env, next);
	    }
	    return true;
	}

	template<class OuterEnv>
	bool advance(const OuterEnv &outer, stop_index) {
	    return false;
	}

    private:
	G &generator_;
	Env env_;
	range_types range_;
	local_types local_;

	value_type value_;
	size_t position_;
	bool state_;
    };


    template<class E, class S= void>
    struct for_range_generator {

	typedef E expression_type;
	typedef for_range_generator self;

	for_range_generator(S &subclass, const E &expression)
	    : subclass_(subclass), expression_(expression) {}

	size_t advance() { return subclass_.advance(); }

	template<class Env>
	for_range_iterator<self,Env> begin(const Env &env) {
	    return for_range_iterator<self,Env>(*this, env);
	}
    
	template<class F>
	typename for_range_eval::actor<for_range_generator, F>::type
	operator[](const F &f) {
	    return for_range_eval::compose(*this, f);
	}
	const E& expression() const { return expression_; }
    private:
	S &subclass_;
	E expression_;
    };


    template<class E>
    struct for_range_generator<E>
	: for_range_generator<E,for_range_generator<E> > {
	typedef for_range_generator<E> base_type;
	for_range_generator(const E &expression)
	    : base_type(*this, expression), value_(0) {}
	size_t advance() { return ++value_; }
    private:
	size_t value_;
    };


    template<class Start, class Stop, class Step = boost::phoenix::value<size_t> >
    struct range_actor {
	typedef boost:: mpl::false_ no_nullary;
	range_actor(const Start &start, const Stop &stop, const Step &step)
	    : start(start), stop(stop), step(step) {}

	template<class Env>
	struct result {
	    typedef struct do_not_instantiate type;
	};

	template<class Env>
	struct value {
	    typedef typename Stop::template result<Env>::type type;
	};

	const Start start;
	const Stop stop;
	const Step step;
    };
	    

    template<class Stop>
    boost::phoenix::actor<
	range_actor<
	    boost::phoenix::value<size_t>,
	    typename boost::phoenix::as_actor_base<Stop>::type> >
    range(Stop stop) {
	typedef boost::phoenix::value<size_t> value;
	typedef boost::phoenix::as_actor_base<Stop> base;
	return range_actor<value, typename base::type>(0, base::convert(stop), 1);
    }



} // namespace phoenix
} // namespace boost


#endif // PHOENIX_FOR_RANGE_HPP
