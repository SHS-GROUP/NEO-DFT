#ifndef BASIS_BASIS_HPP
#define BASIS_BASIS_HPP

/**
   @file
   @brief
*/

#include "basis/forward.hpp"
#include "molecule.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include "basis/shell.hpp"
#include "basis/block.hpp"

#include <boost/range/iterator_range.hpp>

/**
   @brief
*/
class basis::Basis {

    /**
       @brief
    */

public:

    typedef basis::Shell Shell;
    typedef basis::Block Block;
    typedef basis::Center Center;

    typedef std::vector<Shell>::const_iterator const_iterator;
    typedef std::vector<Shell>::iterator iterator;

    /**
       @brief Basis constructor
       @param molecule
    */
    Basis(const Molecule &molecule) :
	molecule_(&molecule), max_(0), 
	num_functions_(0), sorted_(false) {}

    std::vector<Shell>& shells() { return shells_; }
    const std::vector<Shell>& shells() const { return shells_; }
    const Shell& operator[](size_t i) const { return shells_.at(i); }

    const_iterator begin() const { return shells().begin(); }
    const_iterator end() const { return shells().end(); }

    boost::iterator_range<std::vector<Shell>::const_iterator>
    range(const basis::Shell* shell) const {
	size_t distance = (shell - &shells().front());
	if (distance > shells().size()) 
	    throw std::range_error("invalid shell range");
	return boost::make_iterator_range(shells().begin(),
					  shells().begin() + distance);
    }

    int index(const Shell &shell) const {
	return shell.index(shells());
    }

    size_t size() const { return num_functions_; }

    /** @brief Get number of shells */
    int num_shells() const { return shells_.size(); }

    /** @brief Get total number of functions */
    int num_functions() const { return num_functions_; }

    /** @brief Get the shells */
    const Shell& shell(int i) const { return shells_[i]; }

    /**
       @brief
       @param K
       @param exp
       @param C
       @param L
       @param center
    */
    void add(int K, double *exp, double *C, int L, int center) {
	add(Shell::Data(K, exp, C, L), center);
    }

    /**
       @brief
       @param K
       @param exp
       @param Cs
       @param Cp
       @param center
    */
    void add(int K, double *exp, double *Cs, double *Cp, int center) {
	add(Shell::Data(K, exp, Cs, Cp), center);
    }

    /**
       @brief
    */
    const Shell::Data& max() const { return *max_; }

    size_t max_block() const {
	size_t max_block = 0;
	for (size_t i = 0; i < blocks().size(); ++i) {
	    max_block = std::max(max_block, blocks()[i].range().size());
	}
	return max_block;
    }

    const std::vector<Center> centers() const {
	std::vector<Center> centers(shells().size());
	for (size_t i = 0; i < centers.size(); ++i) {
	    centers[i] = shells().at(i).center();
	}
	return centers;
    }

    const std::vector<const Shell::Data*>& data() const {
	return shell_data_;
    }

    void sort();
    void reverse_sort();
    void reverse();

    const std::vector<int>& function_permutations() const {
	return function_permutation_;
    }

    int function_permutation(int i) const { return function_permutation_.at(i); }

    const std::vector<int>& index() const {
	return function_permutation_;
    }

    template<class A>
    A index() const {
	A index(size());
	for (size_t i = 0; i < size(); ++i) {
	    index[i] = function_permutation_[i];
	}
	// std::copy(function_permutation_.begin(),
	// 	  function_permutation_.end(),
	// 	  &index[0]);
	return index;
    }


    const std::vector<int>& shell_index() const {
    	return shell_permutation_;
    }

    int shell_index(int i) const { return shell_permutation_.at(i); }

    template<class A>
    A shell_index() const {
	A index(num_shells());
	std::copy(shell_permutation_.begin(),
		  shell_permutation_.end(),
		  &index[0]);
	return index;
    }

    typedef std::vector<Block>::iterator block_iterator;
    typedef std::vector<Block>::const_iterator const_block_iterator;

    const Block& block(int i) const { return blocks_.at(i); }
    const std::vector<Block>& blocks() const { return blocks_; }
    const Block* firstBlock() const { return &(blocks_.front()); }
    const Block* lastBlock() const { return &(blocks_.back()); }
    std::vector<Block> blocks(size_t max) const;

    std::vector<double> N() const;

private:


    /** creates blocks based contiguous _shells of same type in same block
	@brief
	@param s
	@param center
    */
    //void blockShell::Datas();


    struct Mapping;

    void add(const Shell::Data &shell, int atom);
    size_t add(Shell::Data shell);
    size_t add(Shell shell, Mapping index, bool block);

    void clear(){
	shells_.clear();
	function_to_shell_.clear();
	shell_permutation_.clear();
	function_permutation_.clear();
	blocks_.clear();
	num_functions_ = 0;
    }

    /** Molecule Object
	@brief
	@param molecule
    */
    const Molecule *molecule_;
    /** @brief List of shells*/
    std::vector<const Shell::Data*> shell_data_;

    /** @brief List of mapped shells*/
    std::vector<Shell> shells_;

    std::vector<int> function_to_shell_;
    std::vector<int> function_permutation_;
    std::vector<int> shell_permutation_;

    std::vector<Block> blocks_;

    const Shell::Data* max_;

    size_t num_functions_;

    bool sorted_;
};

typedef basis::Basis Basis;

namespace basis {

    inline void recenter(Basis &basis,
			 const double* const (&r)[3], const int (&inc)[3]) {
	for (size_t i = 0; i < basis.shells().size(); ++i) {
	    int atom = basis.shells().at(i).atom();
	    Center c;
	    for (int j = 0; j < 3; ++j) {
		c[j] = *(r[j] + atom*inc[j]);
	    }
	    basis.shells().at(i).recenter(c);
	}
    }

    std::ostream& operator<<(std::ostream &os, const Basis &basis);

}

#endif // BASIS_BASIS_HPP
