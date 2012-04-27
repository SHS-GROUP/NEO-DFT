#include <stdint.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/array.hpp>

#include "cxx/foreach.hpp"
#include "externals/cxx/array.hpp"
#include "externals/cxx/utility.hpp"
#include "cxx/iterator/product.hpp"
#include "cxx/sugar.hpp"

#include <rysq.hpp>


namespace rysq {

    template<class C>
    struct Index {
	typedef std::vector<C> vector_type;
	typedef std::map<uint32_t,size_t> map_type;

	static uint32_t key(const Quartet<rysq::type> &quartet) {
	    return cxx::utility::pack(QUARTET(quartet));
	}

	static vector_type generate(const Quartet<rysq::type> &quartet) {
	    rysq::shell shells[] = { QUARTET(quartet) };
	    vector_type index;
	    typedef boost::array<int,4> product;
	    boost::fortran_storage_order order;
	    // std::cout << quartet << std::endl;
	    foreach (const product p, (cxx::iterator::product(shells, order))) {
		// std::cout << p << "\n";
	        int x = 0, y = 0, z = 0;
		int m = 1;
	        for (int i = 0; i < 4; ++i) {
		    x += rysq::LX[p[i]]*m;
		    y += rysq::LY[p[i]]*m;
		    z += rysq::LZ[p[i]]*m;
		    m *= (shells[i].L+1);
	        }
		typename vector_type::value_type value = { x,y,z };
		// std::cout <<  p << value << std::endl;
	        index.push_back(value);
	    }
	    return index;
	}

	Index() {
	    std::vector<rysq::type> types;
	    foreach (rysq::type type, rysq::types) types.push_back(type);
	    initialize(types);
	}

	Index(const std::vector<rysq::type> &types) {
	    initialize(types);
	}

	void initialize(const std::vector<rysq::type> &types) {
	    typename vector_type::iterator it = index_.begin();

	    // std::ostream_iterator<ushort3, char> ushort3_output(std::cout, "\n");
	    // std::ostream_iterator< int > key_output(std::cout, " " );

	    typedef boost:: array< rysq:: type,4> product;
	    foreach (const product q,  (cxx::iterator::product<4>(types ))) {
		size_t size =  rysq::Quartet<rysq::Shell>::size(q.elems);
		if (size > 2400) continue;
		// std::cout << q << " " << size << "\n";
		assert(map_.size() < map_.max_size());

		if (!this->find(q)) {
		    vector_type index = generate(q);
		    map_[key(q)] = index_.size();
		    std::copy(index.begin(), index.end(), std::back_inserter(index_));
		}
	    }
	}
	const C* find(const boost::array<rysq::type,4> &quartet) const {
	    typename map_type::const_iterator it = map_.find(key(quartet));
	    if (it == map_.end()) return NULL;
	    return &index_[(*it).second];
	}
	const vector_type& data() const { return index_; }
	const map_type& map() const { return map_; }

    private:
	map_type map_;
	vector_type index_;
    };
}
