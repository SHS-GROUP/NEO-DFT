#ifndef _MOLECULE_H
#define _MOLECULE_H

#include <string>
#include <vector>
#include <iostream>
#include <boost/array.hpp>

/**
   @class Molecule
*/
class Molecule {
    friend std::ostream& operator<<(std::ostream &output, const Molecule &m);

public:

    typedef boost::array<double,3> Center;

    /**
       @brief Molecule constructor
    */
    Molecule(int numAtoms, int *Z, double *r3, double *Q = NULL);

    /**
       @brief Number of atoms in Molecule
    */
    uint size() const { return _centers.size(); }

    const Center& operator()(int center) const {
	return _centers[center];
    }

    double operator()(int center, int coord) const {
	return (_centers[center])[coord];
 }

    const std::vector<Center>& centers() const { return _centers; }

    //double Vnuc();

    class Atom {
	friend std::ostream& operator<<(std::ostream& output, const Atom& a);

    public:
	/**
	   @brief Atom constructor
	*/
	Atom(int Z, double Q) : _Z(Z), _Q(Q) {}

	/**
	   @brief Nuclear charge of atom
	*/
	int Z() const { return _Z; }

	/**
	   @brief Charge of atom
	*/
	double Q() const { return _Q; }

    private:
	/** @brief Nuclear charge */
	int _Z;
	/** @brief Charge */
	double _Q;

    };

private:
    /** @brief List of Atoms in Molecule */
    std::vector<Atom> _atoms;
    std::vector<Center> _centers;

};


#endif // _MOLECULE_H
