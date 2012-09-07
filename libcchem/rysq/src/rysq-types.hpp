/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#ifndef RYSQ_TYPES_HPP
#define RYSQ_TYPES_HPP


namespace rysq {

    enum type {
	SP = -1, S, P, D, F, G
    };
    
    static const type types[] = {
	SP, S, P, D, F, G
    };

#ifdef __CUDACC__
    __constant__
#endif
    const int LX[] = {
	0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1
    };

#ifdef __CUDACC__
    __constant__
#endif
    const int LY[] = {
	0, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1
    };

#ifdef __CUDACC__
    __constant__
#endif
    const int LZ[] = {
	0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2
    };

#ifdef __CUDACC__

    template<int N>
    struct Orbitals {
	long int x_, y_;
	__host__ __device__
	int x(int i) const { return ((x_ >> N*i) & ((1 << N)-1)); }
	__host__ __device__
	int y(int i) const { return ((y_ >> N*i) & ((1 << N)-1)); }
    };

    __constant__
    const Orbitals<3> orbitals[] = {
	{ 8,
	  64 }
	,
	{ 0,
	  0 }
	,
	{ 1,
	  8 }
	,
	{ 4610,
	  33296 }
	,
	{ 136356867,
	  151585304 }
	,
	{ 5087659341316,
	  5583743582752 }
    };

#endif


}


#define RYSQ_LMAX 4
#define RYSQ_TYPES (rysq::SP)(rysq::S)(rysq::P)(rysq::D)(rysq::F)(rysq::G)

#endif // RYSQ_TYPES_HPP
