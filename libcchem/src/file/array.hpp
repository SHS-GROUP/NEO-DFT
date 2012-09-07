#ifndef FILE_ARRAY_HPP
#define FILE_ARRAY_HPP

#include <H5Cpp.h>
#include <string>

#include <boost/array.hpp>
#include <boost/multi_array/index_range.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include "boost/utility/timer.hpp"

// struct File {

//     template<size_t N, typename T = double>
//     struct Array;

//     struct Counters {
// 	typedef boost::utility::timer timer;
// 	typedef timer::value_type value_type;
// 	void clear() {
// 	    *this = Counters();
// 	}
// 	const value_type& read() const { return read_; }
// 	const value_type& write() const { return write_; }
//     private:
// 	template<size_t, typename T> friend class Array;
// 	value_type read_, write_;
//     };

//     File() {}

//     H5::H5File& data() { return file_; }

//     template<typename T>
//     Array<4,T> array(size_t N1, size_t N2, size_t N3, size_t N4);

//     static File create(const std::string &name) {
// 	return File(H5::H5File(name, H5F_ACC_TRUNC));
//     }

//     static File open(const std::string &name) {
// 	H5::H5File file;
// 	file.openFile(name, H5F_ACC_RDWR);
// 	return File(file);
//     }

// private:
//     H5::H5File file_;
//     File(const H5::H5File &file) : file_(file) {}
// };


// #ifndef H5_HAVE_THREADSAFE
// #warning "HDF5 not thread safe: !defined(H5_HAVE_THREADSAFE)"
// #endif


// template<size_t N, typename T>
// struct File::Array {

//     typedef hsize_t size_type;
//     typedef boost::detail::multi_array::index_range<size_t,size_t> index_range;

//     struct Properties : H5::DSetCreatPropList {
// 	Properties& set_chunk(const size_t (&chunk)[N]) {
// 	    boost::array<size_type,N> chunk_;
// 	    std::copy(chunk, chunk+N, chunk_.begin());
//  	    return set_chunk(chunk_);
// 	}
// 	Properties& set_chunk(const boost::array<size_type,N> &chunk) {
// 	    this->setChunk(N, chunk.data());
// 	    return *this;
// 	}
//     };

//     Array(File file, const std::string &name, const size_t (&dims)[N],
// 	  const Properties &properties = Properties()) {
// 	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
// 	H5::DataSpace fspace = H5::DataSpace(N, reverse(dims).data());
// 	name_ = name;
// 	file_ = file;
// 	dataset_ = file_.data().createDataSet(name, datatype, fspace, properties);
// 	mutex_.reset(new boost::mutex());
// 	counters_.reset(new File::Counters());
//     }

//     void put(const T *buffer, const size_t *start, const size_t *stop) {
// 	H5::DataSpace fspace = dataset_.getSpace();
// 	size_t size = select(fspace, start, stop);
// 	size_type mdims[] = { size };
// 	H5::DataSpace mspace(1, mdims);
// 	mspace.selectAll();
// 	boost::lock_guard<boost::mutex> lock(*mutex_);
// 	typename File::Counters::timer timer;
// 	dataset_.write(buffer, H5::PredType::NATIVE_DOUBLE, mspace, fspace);
// 	counters_->write_ += timer;
//     }

//     void get(T *buffer, const size_t *start, const size_t *stop) const {
// 	H5::DataSpace fspace = dataset_.getSpace();
// 	size_t size = select(fspace, start, stop);
// 	size_type mdims[] = { size };
// 	H5::DataSpace mspace(1, mdims);
// 	mspace.selectAll();
// 	boost::lock_guard<boost::mutex> lock(*mutex_);
// 	typename File::Counters::timer timer;
// 	dataset_.read(buffer, H5::PredType::NATIVE_DOUBLE, mspace, fspace);
// 	counters_->read_ += timer;
//     }

//     void open() {
// 	dataset_ = file_.data().openDataSet(name_);
//     }

//     void close() {
// 	dataset_.close();
//     }

//     void synchronize() {
// 	file_.data().flush(H5F_SCOPE_GLOBAL);
//     }

//     // template<class C>
//     // C& get(C &m, const size_t *start, const size_t *stop) const {
//     // 	H5::DataSpace fspace = dataset_.getSpace();
//     // 	size_t size = select(fspace, start, stop);
//     // 	size_type mdims[] = { size };
//     // 	h5::dataspace mspace(1, mdims);
//     // 	mspace.selectall();
//     // 	t* buffer = m.data().begin();
//     // 	dataset_.read(buffer, h5::predtype::native_double, mspace, fspace);
//     // 	return m;
//     // }

//     const File::Counters& counters() const { return *counters_; }

// private:
//     std::string name_;
//     H5::DataSet dataset_;
//     File file_;
//     mutable boost::shared_ptr<boost::mutex> mutex_;
//     mutable boost::shared_ptr<File::Counters> counters_;

//     template<typename U, size_t M>
//     static boost::array<size_type,M> reverse(const U (&a)[M]) {
// 	boost::array<size_type,M> array;
// 	for (size_t i = 0, j = N-1; i < M; ++i, --j) {
// 	    array[i] = a[j];
// 	}
// 	return array;
//     }

//     static size_t select(const H5::DataSpace &space, size_t start, size_t stop) {
// 	size_t sstart[] = { start };
// 	size_t sstop[] = { stop };
// 	return select(space, sstart, sstop);
//     }

//     static size_t select(const H5::DataSpace &space,
// 			 const size_t *start, const size_t *stop) {
// 	size_type fstart[N];
// 	size_type fstride[N]; // Stride of hyperslab
// 	size_type fcount[N];  // Block count
// 	size_type fblock[N];  // Block sizes
// 	size_t size = 1;
// 	for (size_t i = 0, j = N-1; i < N; ++i, --j) {
// 	    fstart[i] = start[j];
// 	    fcount[i] = stop[j] - start[j];
// 	    fstride[i] = 1;
// 	    fblock[i] = 1;
// 	    size *= fcount[i];
// 	}
// 	space.selectHyperslab(H5S_SELECT_SET, fcount, fstart, fstride, fblock);
// 	return size;
//     }
// };

#endif // FILE_ARRAY_HPP
