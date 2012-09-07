#ifndef _BASIS_BLOCK_HPP_
#define _BASIS_BLOCK_HPP_

#include "basis/forward.hpp"
#include "basis/shell.hpp"

#include <boost/ref.hpp>
#include <boost/assert.hpp>
#include <algorithm>
#include <iostream>

namespace basis {

    class Block {
    public:
	typedef std::vector<Shell> vector_type;
	typedef const Shell* const_iterator;

	explicit Block(const Shell *shell) :
	    shells_(shell), start_(shell->index()), size_(1) {}

	explicit Block(const Shell *shells, size_t start) :
	    shells_(shells+start), start_(start), size_(0) {}

	const Shell* shells() const { return shells_; }

	size_t size() const { return size_; }
	int start() const { return start_; }
	int stop() const { return start_ + size_; }

	const_iterator begin() const { return shells_; }
	const_iterator end() const { return begin() + size(); }

	const Shell::Data& shell() const { return *(begin()->data()); }

	Shell::Function::Range range() const {
	    size_t size = size_*shell().size();
	    size_t start = begin()->start();
	    return Shell::Function::Range(start, start+size);
	}

	// static bool compare_size(const Block &a, const Block &b) {
	//     return a.size() > b.size();
	// }

	bool operator==(const Block &b) const { return (this == &b); }

    private:
	friend class basis::Basis;
	const Shell *shells_;
	int start_, size_;
	bool operator!=(const Shell::Data &shell) const {
	    return shells_->data() != &shell;
	}
    };

    inline bool operator<(const Block&a, const Block &b) {
	return a.range().size() > b.range().size();
    }
    
    inline std::ostream& operator<<(std::ostream& os, const Block &block) {
	return (os << block.shell());
    }


}

#endif /* _BASIS_BLOCK_HPP_ */
