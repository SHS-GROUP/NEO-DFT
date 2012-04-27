#ifndef INTEGRAL_INTEGRAL_GPU_HPP
#define INTEGRAL_INTEGRAL_GPU_HPP

#include "core/integral.hpp"
#include "rysq/cuda.hpp"
#include "boost/utility/profiler.hpp"

#include <boost/array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

namespace integral {

    namespace detail {
	
	struct cache {
	    typedef Basis::Shell Shell;
	    typedef boost::array<Shell::Data::key_type, 4> key;
	    typedef boost::ptr_map<key, rysq::cuda::Eri> type;
	    rysq::cuda::Eri& find(const Shell::Data &P, const Shell::Data &Q,
				  const Shell::Data &R, const Shell::Data &S,
				  const rysq::cuda::stream &stream) 
	    {
		BOOST_PROFILE_FUNCTION();
	    	
		key k = {{ P, Q, R, S }};
		
		BOOST_AUTO(iterator, data_.find(k));
		if (iterator == data_.end()) {
		    BOOST_PROFILE_FUNCTION("create");
		    adapter::rysq::Quartet quartet(P,Q,R,S);
		    rysq::cuda::Eri *eri =
			new rysq::cuda::Eri(quartet, stream);
		    data_.insert(k, eri);
		    return *eri;
		}
		return *(iterator->second);
	    }
	    void clear() { data_.clear(); }
	private:
	    type data_;
	};

	template<class>
	struct static_cache_type {
	    static cache instance;
	};

	template<class T>
	cache static_cache_type<T>::instance;

	typedef static_cache_type<void> static_cache;
	
    }
    
    struct Integral::Gpu {
	typedef Basis::Shell Shell;
	typedef rysq::cuda::Centers Centers;
	typedef boost::array<index,4> index4;
	typedef std::vector<index4> Quartets;
	
	class Quartet {
	    friend class Integral::Gpu;
	    rysq::cuda::Eri *eri_;
	    const Centers *centers_;
	    Quartets *quartets_;
	    Quartet(rysq::cuda::Eri *eri, const Centers *centers, Quartets *quartets)
		: eri_(eri), centers_(centers), quartets_(quartets) {}
	public:
	    Quartet() : eri_(), centers_(), quartets_() {}
	    operator bool() const { return (eri_ != NULL && *eri_); }
	    void operator()(size_t size,
			    boost::cuda::device_ptr<const index4> quartets,
			    boost::cuda::device_ptr<double*> I) {
		// BOOST_PROFILE_FUNCTION();
		(*eri_)(*centers_, size, quartets, I);
	    }
	};


	void set(const std::vector<Basis::Center> &centers) {
	    centers_.assign(centers);
	}

	void synchronize() {
	    boost::cuda::thread::synchronize();
	}

	Quartet operator()(const Basis::Block &P, const Basis::Block &Q,
			   const Basis::Block &R, const Basis::Block &S) {
	    return (*this)(P.shell(), Q.shell(), R.shell(), S.shell());
	}

	Quartet operator()(const Shell::Data &P, const Shell::Data &Q,
			   const Shell::Data &R, const Shell::Data &S) {
	    ++data_;
	    return Quartet(&cache_.find(P,Q,R,S, data_.stream()),
			   &centers_, &data_.quartets());
	}

	detail::cache& cache() { return detail::static_cache::instance; }
    private:

	struct data {
	    data() {
		for (int i = 0; i < 1; ++i)
		    tuple_.push_back(new tuple());
		index_ = 0;
	    }
	    Quartets& quartets() { return tuple_[index_].quartets; }
	    boost::cuda::stream& stream() { return tuple_[index_].stream; }
	    data& operator++() {
		++index_;
		if (index_ == tuple_.size()) index_ = 0;
		return *this;
	    }
	private:
	    size_t index_;
	    struct tuple {
		tuple() : stream(0) {}
		Quartets quartets;
		boost::cuda::stream stream;
	    };
	    boost::ptr_vector<tuple> tuple_;
	};

	data data_;
	detail::cache cache_;
	rysq::cuda::Centers centers_;

	// struct thread {
	//     ~thread() { rysq::cuda::thread::exit(); }
	// } thread_;

    };



}

#endif // INTEGRAL_INTEGRAL_GPU_HPP
