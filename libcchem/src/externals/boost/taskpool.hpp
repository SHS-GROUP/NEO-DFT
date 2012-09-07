#ifndef BOOST_TASKPOOL_HPP
#define BOOST_TASKPOOL_HPP

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace boost {

    template<class T>
    struct taskpool {

	struct task {
	    virtual void operator()(T &t) const = 0;
	};

	struct noop : task {
	    void operator()(T &t) const { }
	};

	struct thread : T {
	    thread(const T &t, const taskpool &pool) : T(t), pool_(pool) {}
	    // thread(const thread &thread) ; 
	    void operator()() {
		//std::cout << "while" << std::endl;
		while (true) {

		    task *task = 0;
		    {
			std::cout  << "wait "
				   << this_thread::get_id() << std::endl;
			boost::unique_lock<boost::mutex> lock(pool_);
			if (pool_.tasks().empty()) {
			    boost::system_time timeout =
				boost::get_system_time() +
				boost::posix_time::milliseconds(1000);
			    pool_.cond_.timed_wait(lock, timeout);
			}
			//std::cout << "xxx " << pool_.tasks().size() << std::endl;
			if (pool_.tasks().empty()) continue;
			task = pool_.tasks().front();
			pool_.tasks().erase(pool_.tasks().begin());
		    }
		    std::cout  << "task "
			       << this_thread::get_id() << " " 
			       << pool_.tasks().size() << std::endl;
		    if (!task) break;
		    (*task)(*this);
		}
		std::cout << "break " << pool_.tasks().size() << std::endl;
	    }
	private:
	    const taskpool &pool_;
	};

	taskpool() : task_(NULL) {}
	taskpool(const T &t, size_t size = 1) {
	    for (size_t i = 0; i < size; ++i) {
	        add_thread(t);
	    }
	    barrier();
	}
	~taskpool() {
	    // std::cout << "destroy" << std::endl;
	    //barrier();
	    all(NULL);
	    threads_.join_all();
	}
	void add_thread(const T &t) {
	    using boost::cref;
	    threads_.add_thread(new boost::thread(thread(cref(t), cref(*this))));
	}

	void barrier() {
	    noop noop;
	    all(noop);
	}

	operator bool() const { return task_; }
	operator boost::mutex&() const { return mutex_; }
	template<class L>
	void wait(L &lock) const {
	    cond_.wait(lock);
	}
	std::vector<task*>& tasks() const { return tasks_; }
	void all(task &task) {
	    all(&task);
	}
	void all(task &task, T &t) {
	    all(&task, &t);
	}
    private:

	struct barrier_task : task {
	    barrier_task(size_t size, task *task)
		: barrier_(size), task_(task) {}
	    void  operator()(T &t) const {
		run(t);
		wait();
	    }
	    void run(T &t) const { (*task_)(t); }
	    void wait() const {
		std::cout << "block" << std::endl;
 barrier_.wait(); }
	private:
	    mutable boost::barrier barrier_;
	    task *task_;
	};

	void all(task *task, T *t = 0) {
	    size_t size = threads_.size();
	    barrier_task task_(size+1, task);
	    if (task) task = &task_;

	    {
		boost::lock_guard<boost::mutex> lock(mutex_);
		if (tasks_.size()) throw std::runtime_error("task queue not empty");
		tasks_.resize(tasks_.size() + size, task);
	    }
	    cond_.notify_all();
	    
	    // while(true) {
	    // 	cond_.notify_one(); 
	    // 	boost::lock_guard<boost::mutex> lock(mutex_);
	    // 	if (tasks_.empty()) break;
		
	    // }

	    if (task) {
		if (t) task_.run(*t);		
		task_.wait();
	    }
	}
    private:
	mutable boost::condition_variable cond_;
	mutable boost::mutex mutex_;
	task *task_;
	mutable std::vector<task*> tasks_;
	boost::thread_group threads_;
	// boost::barrier barrier_;
    };

}

#endif // BOOST_TASKPOOL_HPP
