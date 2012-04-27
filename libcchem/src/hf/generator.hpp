#ifndef _GENERATOR_HPP_
#define _GENERATOR_HPP_

#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/noncopyable.hpp>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/erase.hpp>

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/include/as_vector.hpp>

#include <memory>

/**
   for loop generator which can use lambda expressions.
  
   For example:
   @code
   using namespace generator;
   using namespace boost::lambda;
   make_for(N, N, range(bind(std::max<int>, _1, _2), N), range(_2, _3+1));
   // equivalent to  pseudocode
   // for l=0,N: for k=0,N: for j=max(l,k),N: for i=k,j
   @endcode

   If range is given as upper bound only,
   lower bound is assumed to be default constructed
   Lambda placeholders may only reference first three indices.
*/

namespace generator {
    namespace detail {

	using boost::lambda::constant_type;
	using boost::lambda::constant;

	/// lambda expression identity
	template<class E, class enable = void>
	struct lambda {
	    typedef E type;
	};

	/// transform/construct constant lambda expression from non-lambda
	template<class E>
	struct lambda<E, typename boost::disable_if<
			     boost::lambda::is_lambda_functor<E> >::type>
	{
	    struct constant : boost::lambda::constant_type<E>::type {
		typedef typename boost::lambda::constant_type<E>::type base_type;
		constant() : base_type(boost::lambda::constant(E())) {}
		constant(const E &e) : base_type(boost::lambda::constant(e)) {}
	    };
	    typedef constant type;
	};

	/// range functor
	template<class L, class U>
	struct range_ {
	    typedef boost::array<int,4> index_type;
	    range_(U upper) : bounds_(typename lambda<L>::type(), upper) {}
	    range_(L lower, U upper) : bounds_(lower, upper) {}

	    template< typename T, size_t N>
	    T lower(const boost::array<T,N> &index) {
		return bound<0>(index);
	    }

	    template< typename T, size_t N>
	    T upper(const boost::array<T,N> &index) {
		return bound<1>(index);
	    }

	private:
	    template<bool b, typename T>
	    T bound(const boost::array<T,1> &index) {
		return (boost::fusion::at_c<b>(bounds_))(index[0]);
	    }

	    template<bool b, typename T>
	    T bound(const boost::array<T,2> &index) {
		return (boost::fusion::at_c<b>(bounds_))(index[0], index[1]);
	    }

	    template<bool b, typename T, size_t N>
	    T bound(const boost::array<T,N> &index) {
		using boost::fusion::at_c;
		return (at_c<b>(bounds_))(index[0], index[1], index[2]);
	    }

	    boost::fusion::vector<typename lambda<L>::type,
				  typename lambda<U>::type> bounds_;
	};

	template<typename T, size_t N>
	struct for_base {
	    typedef boost::array<T,N> value_type;
	    virtual ~for_base() {}
	    virtual value_type next() = 0;
	    virtual operator bool() const = 0;
	    virtual for_base* new_() const = 0;
	};
	    
	/// N-index generator
	template<typename T, size_t N, class R, class I>
	struct for_ : for_base<T,N> {
	    typedef typename for_base<T,N>::value_type value_type;
	    typedef R range_tuple;
	    for_(const range_tuple &r) : r_(r), state_(true) {
		boost::fusion::for_each(r_, initialize(index));
	    }
	    /// @return new generator
	    for_* new_() const { return new for_(r_); }
	    /// @return next index value and increment
	    value_type next() {
		value_type next;
		namespace lambda = boost::lambda;
		using lambda::var;
		typename value_type::iterator n = next.begin();
		typename value_type::iterator i = index.begin();
		boost::mpl::for_each<I>(*(var(n))++ = var(i)[lambda::_1]);

		state_ = advance<N>(r_, index);
		return next;
	    }
	    /// @return false if out of bounds, true otherwise
	    operator bool() const { return state_; }

	private:
	    /// initialize indices
	    struct initialize {
		value_type &index_;
		mutable size_t i_;
		initialize(value_type &index) : index_(index), i_(0) {}
		template<class R_> void operator()(R_& r) const {
		    index_[i_++] = r.lower(index_);
		}
	    };

	    /// advance index[0:M)
	    template<size_t M>
	    struct advance {
		/// stop recursion
		struct stop {
		    stop(R r, value_type &index) {}
		};
		/// advance index
		/// @param r range tuple
		/// @param index  index array
		advance(R &r, value_type &index) : index_(index), i_(0) {
		    namespace fusion = boost::fusion;
		    index[M-1] += 1; // increment index
		    fusion::for_each(r, *this); // update indices
		    state_ = index[M-1] >= fusion::at_c<M-1>(r).upper(index);
		    if (state_) { // out of bounds
			typename boost::mpl::if_c<(M > 1),
			    advance<M-1>, stop>::type(r, index);
		    }
		}
		/// apply lower bound of range to index
		template<typename R_> void operator()(R_& r) const {
		    if (i_ >= M) index_[i_] = r.lower(index_);
		    ++i_;
		}
		/// @return false if out of bounds, true otherwise
		operator bool() const { return state_; }
	    private:
		value_type &index_; ///< index array reference
		mutable size_t i_; ///< running index
		bool state_;	///< out of bounds state
	    };	    

	    value_type index;
	    range_tuple r_;
	    bool state_;
	};


	/// polymorphic generator template base
	template<typename T,size_t N>
	struct For {
	    For(const For &f) : for_(f.for_->new_()) {}
	    For(const for_base<T,N> &f) : for_(f.new_()) {}
	    typedef boost::array<T,N> value_type;
	    /// @return next index value and increment
	    value_type next() { return for_->next(); }
	    /// @return false if out of bounds, true otherwise
	    operator bool() const { return for_; }
	protected:
	    For& operator=(const For &f) {}
	    /// reset smart pointer
	    void reset(for_base<T,N> *f) { for_.reset(f); }
	    std::auto_ptr<for_base<T,N> > for_;
	};

	/// range [T,R) type
	template<typename T, typename R>
	struct range_type {
	    typedef range_<T,R> type;
	};

	/// range identity specialization
	template<typename T, class L, class U>
	struct range_type<T, range_<L,U> > {
	    typedef range_<L,U> type;
	};

	namespace fusion = boost::fusion;
	namespace mpl = boost::mpl;

	template<typename T, size_t N, class R1, class R2, class R3, class R4>
	struct range_tuple {
	    // full range vector
	    typedef typename mpl::vector<R1,R2,R3,R4> v;
	    typedef typename mpl::end<v>::type end;
	    typedef typename mpl::advance_c<typename mpl::begin<v>::type, N>::type pos;
	    // [0:N) range vector
	    typedef typename mpl::erase<v, pos, end>::type t;
	    // transform into proper range fusion::vector
	    typedef typename fusion::result_of::as_vector<
		typename mpl::transform<t,range_type<T, mpl::_1> >::type
		>::type type;
	};


	template<typename T, size_t N,
		 class R1, class R2, class R3, class R4,
		 class O>
	struct for_type {
	    typedef typename range_tuple<T,N,R1,R2,R3,R4>::type range_tuple;
	    typedef for_<T, N, range_tuple, O> type;
	};

    } // namespace detail


    /// default index order, [0:N)
    template<size_t  N>
    struct order {
	typedef boost::mpl::range_c<size_t,0, N> type;
    };

    /// N-loop generator, 0 < N <= 5
    /// @tparam T index type
    /// @tparam N number of indices/loops
    /// @tparam R1,... range types
    /// @tparam O index order
    template<typename T, size_t N,
	     class R1, class R2 = void, class R3 = void, class R4 = void,
	     class O = typename order<N>::type>
    struct for_ : detail::for_type<T, N, R1, R2, R3, R4, O>::type {
    	typedef typename detail::for_type<T, N, R1, R2, R3, R4, O>::type base_type;
    	typedef typename base_type::range_tuple range_tuple;
    	for_(const range_tuple &range) : base_type(range) {}
    };

    /// loop range [L:U)
    /// @tparam L lower bound type
    /// @tparam U upper bound type
    /// @return range
    template<class L, class U>
    detail::range_<L,U> range(L lower, U upper) {
	return detail::range_<L,U>(lower, upper);
    }

    /// make 4-loop generator with specified index ordering
    template<typename T, class R1, class R2, class R3, class R4, class O>
    for_<T, 4, R1, R2, R3, R4, O>
    make_for(R1 r1, R2 r2, R3 r3, R4 r4, const O&) {
	typedef for_<T, 4, R1, R2, R3, R4, O> F;
    	return F(typename F::range_tuple(r1, r2, r3, r4));
    }

    /// polymorphic generator template forward declaration
    template<typename T,size_t N>
    struct For;

    /// polymorphic 4-loop generator
    template<typename T>
    struct For<T,4> : detail::For<T,4> {
	typedef detail::For<T,4> base_type;
	For(const For &F) : base_type(F) {}
	/// generator with default index ordering
	template<class R1, class R2, class R3, class R4>
	For(R1 r1, R2 r2, R3 r3, R4 r4)
	    : base_type(make_for<T>(r1, r2, r3, r4)) {}
	/// generator with specified index ordering
	template<class R1, class R2, class R3, class R4, class O>
	For(R1 r1, R2 r2, R3 r3, R4 r4, O o)
	    : base_type(make_for<T>(r1, r2, r3, r4, o)) {}
    };

}


#endif /* _GENERATOR_HPP_ */
