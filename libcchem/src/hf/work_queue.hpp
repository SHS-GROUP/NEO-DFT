#ifndef _HF_WORK_QUEUE_HPP_
#define _HF_WORK_QUEUE_HPP_


#include <queue>
#include <memory>
#include <boost/array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>

#include <boost/array.hpp>


#include "generator.hpp"
#include "parallel/counter.hpp"
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <algorithm>
#include "externals/cxx/array.hpp"

namespace hf {


    struct work_queue_lock {
	boost::mutex backlog;
	boost::condition_variable throttle_condition;
	boost::mutex queue;
    };    

    template<class G>
    struct threads_work_queue {

	typedef typename G::value_type value_type;

	threads_work_queue(work_queue_lock &lock, const G &generator,
			   size_t throttle = 8196)
	    : lock_(lock),
	      generator_(generator),
	      backlog_(lock_, throttle) { }
	bool empty() const { return generator_.empty(); }

	// value_type front() {
	//     //boost::mutex::scoped_lock(mutex_);
	//     boost::lock_guard<boost::mutex> lock(mutex_);
	//     return current_;
	// }

	value_type pop() {
	    boost::lock_guard<boost::mutex> lock(lock_.queue);
	    //std::cout << this << std::endl;
	    if (empty()) throw std::exception();
	    return generator_.next();
	}

	struct backlog_queue {
	    explicit backlog_queue(work_queue_lock &lock,
				   size_t throttle = (1 << 16))
		: lock_(lock)
	    {
		throttle_ = std::max<size_t>(throttle, 2);
		throttle_set_ = false;
	    }
	    bool empty () const { return queue_.empty(); }
	    void push(const value_type &work) {
		// boost::mutex::scoped_lock throttle_lock(throttle_mutex_);
		boost::mutex::scoped_lock lock(lock_.backlog);
		while (throttle_set_) lock_.throttle_condition.wait(lock);
		queue_.push(work);
		if (queue_.size() > throttle_) {
		    //std::cout << "throttle set" << std::endl;
		    throttle_set_ = true;
		}
	    }
	    value_type pop() {
		//boost::unique_lock<boost::mutex> lock(mutex_);
		boost::mutex::scoped_lock lock(lock_.backlog);
		if (queue_.empty()) {
		    throw std::exception();
		}
		value_type front = queue_.front();
		queue_.pop();
		if (queue_.size() < throttle_/2) {
		    throttle_set_ = false;
		    lock_.throttle_condition.notify_all();
		}
		//std::cout <<  queue_.size()<< front << std::endl;
		return front;
	    }
	private:
	    std:: queue< value_type> queue_;
	    size_t throttle_;
	    work_queue_lock &lock_;
	    bool throttle_set_;
	};
	backlog_queue& backlog() { return backlog_; }
    private:
	work_queue_lock &lock_;
	G generator_;
	backlog_queue backlog_;
    };


    struct work_queue_generator {
	typedef boost::array<int,4> value_type;
	int N;
	value_type current;
	parallel::counter<size_t> &counter_;
	size_t last_;
	work_queue_generator(size_t N,
			     parallel::counter<size_t> &counter)
	    : N(N), counter_(counter)
	{
	    current.assign(0);
	    last_ = 0;
	}
	bool empty() const { return current[3] >= N; }
	value_type front() {
	    return current;
	}
	value_type next() {
	    value_type next = current;

	    size_t end = (counter_++)+1;
	    for (size_t i = last_; i < end; ++i) {

		// std::cout << i << std::endl;

		if (empty()) throw std::exception();
		next = current;

		// for (RANGE(lb, basis.blocks())) {
		// 	for (RANGE(jb, basis.blocks())) {
		// 	    for (RANGE(kb, std::max(lb, jb), basis. blocks().end())) {
		// 		for (RANGE(ib, jb, kb+1)) {

		bool advance;
		current[0] += 1;

		advance = (current[0] >= (current[1]+1));
		if (advance) {
		    current[1] += 1;
		    current[0] = current[2];
		}

		advance = advance && (current[1] == N);
		if (advance) {
		    current[2] += 1;
		    current[1] = std::max(current[2], current[3]);
		    current[0] = current[2];
		}

		advance = advance && (current[2] == N);
		if (advance) {
		    current[3] += 1;
		    current[2] = 0;
		    current[1] = std::max(current[2], current[3]);
		    current[0] = current[2];
		}

	    }
	    last_ = end;

	    //std::cout << next << std::endl;
	    return next;
	}
    };

    struct work_queue :
	threads_work_queue<work_queue_generator> {
	work_queue(hf::work_queue_lock& lock, size_t size,
		   parallel::counter<size_t> & counter)
	    : threads_work_queue<work_queue_generator>
	      (lock, work_queue_generator(size, counter)) {}
    };

}


// for (RANGE(lb, N)) {
// 	for (RANGE(jb, N)) {
// 	    for (RANGE(kb, std::max(lb, jb), N)) {
// 		for (RANGE(ib, jb, kb+1)) {



#endif /* _HF_WORK_QUEUE_HPP_ */
