#ifndef _UTIL_MAP_CACHE_HPP_
#define _UTIL_MAP_CACHE_HPP_

#include <list>
#include <vector>
#include <map>
#include <algorithm>

namespace util {

    template<class Key, class Object>
    class map_cache {
	std::list<Object> objects_;
	std::list<Object*> free_;
	std::list<Key> keys_; 
	std::map<Key, Object*> allocated_;

    public:

	map_cache() {}

	void add(Object object) {
	    objects_.push_back(object);
	    free_.push_back(&objects_.back());
	}

	const std::list<Object>& objects() { return objects_; }

	bool contains(const Key &key) const {
	    return (keys_.end() != std::find(keys_.begin(), keys_.end(), key));
	}

	bool full() const { return free_.empty(); }

	bool update(const Key &key) {
	    if (! contains(key)) return false;
	    keys_.remove(key);
	    keys_.push_front(key);
	    return true;
	}

	Key pop() {
	    Key k = pop(keys_);
	    Object *o = allocated_[k];
	    allocated_.erase(k);
	    free_.push_back(o);
	    return k;
	}

	Object& push(const Key &k) {
	    assert(!full() && !contains(k));
	    keys_.insert(keys_.begin(), k);
	    allocated_[k] = pop(free_);
	    return *allocated_[k];
	}

    private:

	template<class T>
	static T pop(std::list<T> &v) {
	    T t = v.back();
	    v.pop_back();
	    return t;
	}

    };

}

#endif /* _UTIL_MAP_CACHE_HPP_ */
