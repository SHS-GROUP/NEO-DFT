#ifndef PHOENIX_THREAD_HPP
#define PHOENIX_THREAD_HPP

#include "phoenix/for_range.hpp"

#include <boost/typeof/typeof.hpp>
#include <boost/thread.hpp>
#include <boost/ref.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/spirit/home/phoenix/core.hpp>
#include <boost/spirit/home/phoenix/scope/local_variable.hpp>
#include <boost/spirit/home/phoenix/scope/let.hpp>
#include <boost/spirit/home/phoenix/operator.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement.hpp>

#include <boost/fusion/include/vector_tie.hpp>
#include <boost/fusion/sequence/intrinsic.hpp>

#include <boost/mpl/bool.hpp>


namespace boost {
namespace phoenix {
namespace thread {


    template<class F>
    struct critical_actor {

	template<class Env>
	struct result { typedef void type; };

	typedef boost::mpl::false_ no_nullary;

	template<class Env>
	void eval(const Env &env) const {
	    boost::mutex::scoped_lock lock(*mutex_);
	    function.eval(env);
	}

	critical_actor(const F &f, boost::mutex *mutex)
	    : function(f), mutex_(mutex) {}
    private:
	F function;
	boost::shared_ptr<boost::mutex> mutex_;
    };


    struct critical_gen {
	template<class F>
	boost::phoenix::actor<critical_actor<F> >
	operator[](const F &f) const {
	    return critical_actor<F>(f, new boost::mutex());
	}
    };
    
    const critical_gen critical = critical_gen();


    struct throttle {

	template <typename Env, class V>
	struct result { typedef void type; };

	template <typename RT, typename Env, class V>
	static void eval(const Env &env, const V &v) {
	    using boost::fusion::at_c;
	    throttle& throttle = at_c<0>(v.val);
	    BOOST_AUTO(const &f, at_c<1>(v.val));
	    throttle.eval_(env, f);
	}

	explicit throttle(size_t throttle)
	    : throttle_(throttle_), size_(0) {}

	template <typename F>
	boost::phoenix::actor<
	    typename boost::phoenix::as_composite<
		throttle,
		boost::fusion::vector<throttle&, const F&> >::type>
	operator[](const F &f) {
	    using boost::fusion::vector_tie;
	    return boost::phoenix::compose<throttle>(vector_tie(*this, f));
	}

    private:
	boost::mutex mutex_;
	boost::condition_variable condition_;
	const size_t throttle_;
	size_t size_;
	bool wait_;
	template <typename Env, class F>
	void eval_(const Env &env, const F &f) {
	    {	
		boost::unique_lock<boost::mutex> lock(mutex_);
		size_ = std::min(size_+1, throttle_);
		while (throttle_ <= size_) condition_.wait(lock);
	    }
	    f.eval(env);
	    {
		boost::lock_guard<boost::mutex> lock(mutex_);
		--size_; 
	    }
	    condition_.notify_one();
	}
    };


    template<class E>
    struct parallel_for_generator
	: for_range_generator<E, parallel_for_generator<E> > {

	typedef for_range_generator<E, parallel_for_generator> base_type;

	parallel_for_generator(const E &expression, boost::mutex *mutex)
	    :  base_type(*this, expression),
	       mutex_(mutex), value_(0) {}

	size_t advance() {
	    boost::lock_guard<boost::mutex> lock(*mutex_);
	    return value_++;
	}
    private:
	boost::shared_ptr<boost::mutex> mutex_;
	size_t value_;
    };


    template<class E>
    static parallel_for_generator<E>
    parallel_for(const E &expression) {
	return parallel_for_generator<E>(expression, new boost::mutex());
    }


    struct group : boost::thread_group {
	explicit group(size_t size) : size_(size) {}
	template<class F>
	boost::thread_group& operator()(F f) {
	    for (size_t i = 1; i < size_; ++i) {
		std::cout << "new thread" << std::endl;
		this->add_thread(new boost::thread(f));
	    }
	    f();
	    this->join_all();
	    return *this;
	}
	~group() {
	    this->join_all();
	}
    private:
	size_t size_;
    };


} // namespace thread
} // namespace phoenix
} // namespace boost


#endif // PHOENIX_THREAD_HPP
