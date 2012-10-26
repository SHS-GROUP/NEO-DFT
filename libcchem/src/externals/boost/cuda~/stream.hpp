#ifndef BOOST_CUDA_STREAM_HPP
#define BOOST_CUDA_STREAM_HPP

#include <boost/type_traits/remove_const.hpp>

namespace boost {
namespace cuda {

    struct synchronous_tag {};
    const synchronous_tag synchronous = synchronous_tag();
	
    struct stream {
	typedef size_t type;
    public:
	static stream create() {
	    return stream(stream::create_());
	}
	stream();
	stream(synchronous_tag)
	    : data_(0), synchronous_(1) {}
	explicit stream(type data);
	stream(const stream &s) { data_ = 0; copy(s); }
	stream& operator=(const stream &s) { copy(s); return *this; }
	~stream();
	type data() const;
	template<typename T>
	T data() const;
	bool synchronous() const { return synchronous_; }
	void synchronize();
    private:
	struct impl;
	impl *data_;
	bool synchronous_;
	void copy(const stream &s);
	static type create_();
	static void destroy_(type);
    };

}
}

#endif /* BOOST_CUDA_STREAM_HPP */
