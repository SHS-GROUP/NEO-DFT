#ifndef HF_WORK_QUEUE_HPP
#define HF_WORK_QUEUE_HPP


#include "parallel.hpp"
#include "omp.hpp"
#include "generator.hpp"

#include <queue>
#include <algorithm>
#include <boost/array.hpp>


namespace hf {


    struct work_queue_lock {
	omp::lock backlog, queue;
    	// boost::mutex backlog;
    	// boost::condition_variable throttle_condition;
    	// boost::mutex queue;
	struct guard : omp::scoped_lock {
	    guard(omp::lock &lock) : omp::scoped_lock(lock) {}
	};
    };    

    template<class G>
    struct threads_work_queue {

	typedef typename G::value_type value_type;

	threads_work_queue(work_queue_lock &lock, const G &generator,
			   size_t throttle = 8196)
	    : lock_(lock),
	      generator_(generator) {}

	bool empty() const { return generator_.empty(); }

	value_type pop() {
	    work_queue_lock::guard lock(lock_.queue);
	    //std::cout << this << std::endl;
	    return generator_.next();
	}

    private:
	work_queue_lock &lock_;
	G generator_;
    };


    struct work_queue_generator {
	typedef boost::array<int,4> value_type;
	int N;
	value_type current;
	Parallel::Task &counter_;
	size_t last_;
	work_queue_generator(size_t N, Parallel::Task &counter)
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
		   Parallel::Task &counter)
	    : threads_work_queue<work_queue_generator>
	      (lock, work_queue_generator(size, counter)) {}
    };

}


// for (RANGE(lb, N)) {
// 	for (RANGE(jb, N)) {
// 	    for (RANGE(kb, std::max(lb, jb), N)) {
// 		for (RANGE(ib, jb, kb+1)) {



#endif /* HF_WORK_QUEUE_HPP */
