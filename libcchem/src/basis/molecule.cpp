#include "molecule.hpp"

#include <math.h>



/** @see molecule.h */
Molecule::Molecule(int numAtoms, int *Z, double *r3, double *Q) {
    for (int i = 0; i < numAtoms; ++i) {
	double q = (Q) ? Q[i] : double(Z[i]);
	Center center = {{ r3[0+i*3], r3[1+i*3], r3[2+i*3] }};
	_atoms.push_back(Atom(Z[i], q));
	_centers.push_back(center);
    }
}

// /** @see molecule.h */
// double Molecule::Vnuc() {
//     double V = 0;
//     for (unsigned int j = 0; j < _atoms.size(); ++j) {
// 	double Qj = _atoms[j].Q();
// 	double xj = _atoms[j](0);
// 	double yj = _atoms[j](1);
// 	double zj = _atoms[j](2);
// 	for (unsigned int i = 0; i < j; ++i) {
// 	    double Qi = _atoms[i].Q();
// 	    double dx2 = (_atoms[i](0) - xj) * (_atoms[i](0) - xj);
// 	    double dy2 = (_atoms[i](1) - yj) * (_atoms[i](1) - yj);
// 	    double dz2 = (_atoms[i](2) - zj) * (_atoms[i](2) - zj);
// 	    double rij = sqrt(dx2 + dy2 + dz2);
// 	    V += (Qi*Qj) / rij;
// 	}
//     }
//     return V;
// }


std::ostream& operator<<(std::ostream& output, const Molecule& molecule) {
    for (uint i = 0; i < molecule.size(); ++i) {
	output << molecule._atoms[i] << "\t";

	output.width(20);
	output.precision(5);
	output << std::left << molecule(i, 0) << '\t';
	
	output.width(20);
	output.precision(5);
	output << std::left << molecule(i, 1) << '\t';
	
	output.width(20);
	output.precision(5);
	output << std::left << molecule(i, 2) << std::endl;
    }

    return output;
}

std::ostream& operator<<(std::ostream& output, const Molecule::Atom& a) {
    output.width(20);
    output.precision(5);
    output << std::left << a.Z();
    return output;
}

#undef OFFSET
