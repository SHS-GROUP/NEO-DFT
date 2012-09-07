#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "rysq-core.hpp"

#include "roots/roots.hpp"
// #include "quadrature.h"

using namespace rysq;

static bool initialized = false;

int rysq::initialize() {
    if (!initialized) {
	roots_initialize();
	initialized = true;
    }
    return 0;
}

int rysq::finalize() {
    roots_finalize();
    initialized = false;
    return 0;
}

const rysq::index_list::type rysq::index_list::list_[] = 
    {{0,1}, {2,3}, {0,2}, {0,3}, {1,2}, {1,3}};

// void  rysq::eri(const Quartet<const Shell> &shells,
// 		const Quartet<const Center> &centers,
// 		double *I,
// 		const Parameters &parameters) {
//     Quadrature q = Quadrature(shells, parameters);
//     q.evaluate(centers, I, parameters);
// }



// void  rysq::FockJK(const Quartet<const Shell> &shells,
// 		   const Quartet<const Center> &centers,
// 		   const Density &D, const Fock &F,
// 		   const Parameters &parameters) {
//     Quadrature q = Quadrature(shells, parameters);
//     q.evaluate(centers, D, F, parameters);
// }


// void rysq::fock2(const Quartet<const Shell> &shells,
// 		 const BlockMatrix &D, BlockMatrix &F,
// 		 const Array<int,4> &quartets, double *Eri,
// 		 const Parameters &parameters) {

//     for (uint q = 0; q < quartets.length; ++q) {

// 	int i = quartets[q][0];
// 	int j = quartets[q][1];
// 	int k = quartets[q][2];
// 	int l = quartets[q][3];

// 	Density D6 = {{ D(i,j), D(k,l), D(i,k), D(i,l), D(j,k), D(j,l) }};
// 	Fock F6 = {{ F(i,j), F(k,l), F(i,k), F(i,l), F(j,k), F(j,l) }};
		
// 	Parameters p = parameters;
// 	p.scale /= symmetry(i, j, k, l);

// 	Quadrature::evaluate(shells, D6, F6, &Eri[q*shells.size()], p);

//     }

// }


// void  rysq::FockJK(const Quartet<const Shell> &shells,
// 		   const std::vector<Center> &centers,
// 		   const BlockMatrix &D, BlockMatrix &F,
// 		   const std::vector<Int4> &quartets,
// 		   const Parameters & parameters) {

//     Quadrature q = Quadrature(shells, parameters);

//     // evaluate quartet vector
//     foreach (Int4 quartet, quartets) {

// 	int i = quartet[0];
// 	int j = quartet[1];
// 	int k = quartet[2];
// 	int l = quartet[3];

// 	Quartet<const Center> c(centers.at(i), centers.at(j),
// 				 centers.at(k), centers.at(l));

// 	Density D6 = {{ D(i,j), D(k,l), D(i,k), D(i,l), D(j,k), D(j,l) }};
// 	Fock F6 = {{ F(i,j), F(k,l), F(i,k), F(i,l), F(j,k), F(j,l) }};
		
// 	Parameters p = parameters;
// 	p.scale /= symmetry(i, j, k, l);

// 	q.evaluate(c, D6, F6, p);
    
//     }

// }
