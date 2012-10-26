// #ifndef PARALLEL_HPP
// #define PARALLEL_HPP

// #include <iostream>
// #include <memory>
// #include <boost/noncopyable.hpp>

// namespace parallel {

//     struct Comm : boost::noncopyable {
// 	virtual int rank() const = 0;
// 	virtual int size() const = 0;
// 	virtual void barrier() = 0;
//     };


//     template<typename T = int>
//     struct Task : boost::noncopyable {
// 	Task() : initial_(0) { reset(); }
// 	void reset();
// 	T next();
// 	T operator++(int) { return next(); }
//     private:
// 	T initial_;
//    };

//     // inline
//     // std::ostream& operator<<(std::ostream &ostream, const Environment &environment) {
//     // 	ostream << "parallel::Environment(";
//     // 	ostream << environment.rank() << ",";
//     // 	ostream << environment.size();
//     // 	ostream << ")";
//     // 	return ostream;
//     // }

// }

// // struct Parallel : boost::noncopyable {
// //     typedef parallel::Comm Comm;
// //     typedef parallel::Task<int> Task;
// //     Parallel();
// //     Comm& node() { return *node_; }
// // private:
// //     std::auto_ptr<Comm> node_;
// // };

// #endif // PARALLEL_HPP
