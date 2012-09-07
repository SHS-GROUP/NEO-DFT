#ifndef THREADS_WORK_GENERATOR_HPP
#define THREADS_WORK_GENERATOR_HPP

#include <boost/array.hpp>

template<class B1, class B2 = void>
struct threads_work_generator;

template<>
struct threads_work_generator<size_t,size_t> {
    typedef boost::array<size_t,2> value_type;
    threads_work_generator(size_t m, size_t n) : m_(m), n_(n)
    {
	current_.assign(0);
    }
    bool empty() const { return current_[1] >= n_; }
    value_type next() {
	if (empty()) throw std::exception();
	value_type next = current_;
	
	bool advance;
	current_[0] += 1;
	advance = (current_[0] == m_);
	if (advance) {
	    current_[1] += 1;
	    current_[0] = 0;
	}

	//std::cout << next << std::endl;
	return next;
    }
private:
    size_t m_, n_;
    value_type current_;
};


template<>
struct threads_work_generator<size_t> {
    typedef boost::array<size_t,1> value_type;
    threads_work_generator(size_t m) : m_(m)
    {
	current_.assign(0);
    }
    bool empty() const { return current_[0] >= m_; }
    value_type next() {
	if (empty()) throw std::exception();
	value_type next = current_;
	
	current_[0] += 1;

	//std::cout << next << std::endl;
	return next;
    }
private:
    size_t m_;
    value_type current_;
};

#endif // THREADS_WORK_GENERATOR_HPP
