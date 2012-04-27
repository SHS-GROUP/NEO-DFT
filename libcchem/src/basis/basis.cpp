#include <functional>
#include <iostream>

#include <boost/typeof/typeof.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "basis/basis.hpp"
#include "util/util.hpp"
#include "iterator/iterator.hpp"


struct Basis::Mapping {
    typedef const Basis::Shell::Data* key_type;
    Mapping() {}
    Mapping(size_t shell, size_t function, key_type key = key_type())
	: shell(shell), function(function), key_(key) {}
    size_t shell, function;
    key_type key_;
    operator const key_type&() const { return key_; }
};

size_t Basis::add(Shell shell, Mapping mapping, bool block) {
    size_t index = shells_.size();
    shell.index_ = index;
    shell.start_ = num_functions_;
    shells_.push_back(shell);

    if ((!block) || (blocks_.empty()) || (blocks_.back() != shell)) {
	const Shell *shells = &shells_[0];
	blocks_.push_back(Block(shells, index));
	// vector may reallocate
	foreach (Block &b, blocks_) {
	    b.shells_ = shells + b.start();
	}
    }
    ++blocks_.back().size_;

    size_t size = shell.size();
    shell_permutation_.push_back(mapping.shell);
    std::generate_n(std::back_inserter(this->function_permutation_),
		    size, util::generator<int>(mapping.function));
    std::fill_n(std::back_inserter(this->function_to_shell_),
		size, index);
    num_functions_ += size;

    return index;

}

void Basis::add(const Shell::Data &shell, int atom) {
    
    const Shell::Data *shell_ = NULL;

    // search for equivalent shell
    foreach (const Shell::Data *s, this->shell_data_) {
        if (*s == shell) { shell_ = s; break; }
    }

    // append new shell
    if (!shell_) {
	{
	    typedef Shell::Data::key_type Key;
	    Key key;
	    if (this->shell_data_.empty()) key = Key('a');
	    else key = Key(*(this->shell_data_.back()));
	    Shell::Data *s = new Shell::Data(shell);
	    s->key_ = ++key;
	    shell_data_.push_back(s);
	}
	shell_ = shell_data_.back();
	if (!max_ || (shell_->size() > max_->size())) {
	    max_ = shell_;
	}
    }

    add(Shell(shell_, (*this->molecule_)(atom), atom),
	Mapping(shells_.size(), num_functions_),
	false);

}

template<class C>
typename C::difference_type distance(const C &range,
				     typename C::const_reference r0,
				     typename C::const_reference r1) {
    return (std::find(range.begin(), range.end(), r0) -
	    std::find(range.begin(), range.end(), r1));
}

void Basis::sort() {

    using namespace boost::lambda;

    std::vector<Mapping> index(shells_.size());
    for (size_t i = 0; i < shells_.size(); ++i) {
	BOOST_AUTO(const &shell, shells_[i]);
	index[i] = Mapping(i, shell.start(), shell.data());
    }
    
    // sorts shells by size
    std::sort(shell_data_.begin(), shell_data_.end(),
	      bind(&Shell::Data::size, _1) > bind(&Shell::Data::size, _2));
    
    // sort mapped shells by shell order
    std::sort(index.begin(), index.end(),
	      bind(&distance<BOOST_TYPEOF(shell_data_)>,
		   boost::ref(shell_data_), _1, _2) > 0);

    BOOST_TYPEOF(shells_) shells(shells_);
    
    this->clear();

    foreach (Mapping i, index) {
	add(shells.at(i.shell), i, true);
    }

    // foreach (const Shell * shell, shells_)
    // 	std::cout << * shell << std::endl;

}


std::ostream& operator<<(std::ostream &output, const basis::Basis &b) {
    for (unsigned int i = 0; i < b.shells().size(); ++i) {
        //output << *(b.shells_.at(i)) << std::endl;
    }

    size_t i = 0;
    foreach (const Basis::Block &block, b.blocks()) {
	foreach (const Basis::Shell &shell, block) {
	    output << i << " -> " << shell.data();

	    // output.width(13);
	    // output.precision(8);
	    output << " " << shell.start() << ':' << shell.stop()-1 << " ";
	    output << b.shell_permutations().at(i) << "/"
		   << b.function_permutations().at(shell.start());

	    // output.width(13);
	    // output.precision(8);
	    // output << shell.center(0) << "\t"
	    //        << shell.center(1) << "\t"
	    //        << shell.center(2) << std::endl;
	
	    output << std::endl;
	    ++i;
	}
    }

    return output;
}
