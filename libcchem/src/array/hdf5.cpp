#include "runtime.hpp"
#include "array/hdf5.hpp"
#include <string>
#include <map>

namespace cchem {
namespace detail {

    std::map<std::string, hid_t> hdf5_;

}
}

Array<>::HDF5::HDF5(const char *name) {
    this->name_ = name;
    this->file_ = 0;
    this->plist_ = H5Pcreate(H5P_FILE_ACCESS);
    this->dcpl_ = H5Pcreate(H5P_DATASET_CREATE);
    this->parallel_ = false;
}

Array<>::HDF5::HDF5(const HDF5 &f) {
    this->name_ = f.name_;
    this->file_ = f.file_;
    this->plist_ = f.plist_;
    this->dcpl_ = f.dcpl_;
    this->parallel_ = f.parallel_;
    H5Iinc_ref(this->plist_);
    H5Iinc_ref(this->dcpl_);
    if (this->file_) H5Iinc_ref(this->file_);
}

Array<>::HDF5::~HDF5() {
    H5Pclose(this->plist_);
    H5Pclose(this->dcpl_);
    if (this->file_) {
	int ref = H5Iget_ref(this->file_);
	H5Fclose(this->file_);
	--ref;
	std::cout << "closed " << this->file_
		  << " count = " << ref
		  << std::endl;
	if (!ref) {
	    std::cout << "erase " << this->name_ << std::endl;
	    cchem::detail::hdf5_.erase(this->name_);
	}
    }
}

hid_t Array<>::HDF5::file() {
    if (!this->file_) {
	using cchem::detail::hdf5_;
	const char *name = this->name_.c_str();
	if (!hdf5_.count(name)) {
	    if (!std::ifstream(name)) {
		this->file_ = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, this->plist_);
	    }
	    else {
		this->file_ = H5Fopen(name, H5F_ACC_RDWR, this->plist_);
	    }
	    hdf5_[name] = this->file_;
	}
	else {
	    this->file_ = hdf5_[name];
	    H5Iinc_ref(this->file_);
	}
    }
    return this->file_;
}

