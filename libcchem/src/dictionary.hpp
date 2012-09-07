#ifndef DICTIONARY_HPP
#define DICTIONARY_HPP

#include <string>
#include <stdexcept>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/typeof/typeof.hpp>

struct Dictionary {
public:

    struct exception : std::runtime_error {
	exception(const std::string &what)
	    : std::runtime_error(what.c_str()) {}
    };

    bool has(const std::string &key) const {
	return data_.count(key);
    }

    template<typename T>
    void set(const std::string &key, const T &val) {
	data_[key] = boost::any(val);
    }

    template<typename T>
    T get(const char *key) const {
	return try_get<T>(key);
    }
    
    template<typename T>
    T get(const char *key, const T &value) const {
	if (!data_.count(key)) return value;
	return try_get<T>(key);	
    }

private:
    std::map<std::string, boost::any> data_;

    template<typename T>
    T try_get(const std::string &key) const {
	typename std::map<std::string, boost::any>::const_iterator it =
	    data_.find(key);
	if (it == data_.end())
	    throw exception(std::string(BOOST_CURRENT_FUNCTION) +
			    ": no such key: " + key);
	try {
	    return boost::any_cast<T>(it->second);
	}
	catch (boost::bad_any_cast) {
	    throw exception(std::string(BOOST_CURRENT_FUNCTION) +
			    ": bad cast: " + key);
	}	
    }

};

#endif /* DICTIONARY_HPP */
