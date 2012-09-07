#ifndef UTILITY_ANY_HPP
#define UTILITY_ANY_HPP

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include <boost/any.hpp>

namespace utility {

    template<typename K>
    struct any {
	typedef K key_type;

	struct no_such_key : std::runtime_error {
	    no_such_key(const K& k) : std::runtime_error(str(k)) {}
	};

	struct bad_cast : std::runtime_error {
	    bad_cast(const K& k) : std::runtime_error(str(k)) {}
	};

	any(const K& k = K()) : key_(k) {}

	template<class O>
	key_type add(const O &o) {
	    key_type key = keys_.empty() ? key_++ : keys_.back();
	    if (!keys_.empty()) keys_.pop_back();
	    data_[key] = boost::any(o);
	    return key;
	}

	bool contains(const K &k) const {
	    return bool(data_.count(k));
	}

	void erase(const K& k) {
	    if (contains(k)) data_.erase(data_.find(k));
	}

	template<class O>
	O find(const K& key) const {
	    if (!data_.count(key)) throw no_such_key(key);
	    try {
		return boost::any_cast<O>(data_.find(key)->second);
	    }
	    catch (boost::bad_any_cast) {
		throw bad_cast(key);
	    }	
	}

    private:
	key_type key_;
	std::vector<key_type> keys_;
	std::map<key_type, boost::any> data_;

	static std::string str(const K& key) {
	    std::stringstream ss;
	    ss << key;
	    return ss.str();
	}

    };

}

#endif /* UTILITY_ANY_HPP */
