/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#ifndef _RYSQ_NORMALIZE_H
#define _RYSQ_NORMALIZE_H

#define RYSQ_NORMAL_000 1.0 /**< @brief (0 0 0) normalization constant */ 
#define RYSQ_NORMAL_100 1.0 /**< @brief (1 0 0) normalization constant */ 
#define RYSQ_NORMAL_200 1.0 /**< @brief (2 0 0) normalization constant */ 
#define RYSQ_NORMAL_110 1.73205080757 /**< @brief (1 1 0) normalization constant */ 
#define RYSQ_NORMAL_300 1.0 /**< @brief (3 0 0) normalization constant */ 
#define RYSQ_NORMAL_210 2.2360679775 /**< @brief (2 1 0) normalization constant */ 
#define RYSQ_NORMAL_111 3.87298334621 /**< @brief (1 1 1) normalization constant */ 
#define RYSQ_NORMAL_400 1.0 /**< @brief (4 0 0) normalization constant */ 
#define RYSQ_NORMAL_310 2.64575131106 /**< @brief (3 1 0) normalization constant */ 
#define RYSQ_NORMAL_220 3.31662479036 /**< @brief (2 2 0) normalization constant */ 
#define RYSQ_NORMAL_211 5.9160797831 /**< @brief (2 1 1) normalization constant */ 

namespace rysq {

    // static inline void normalize(const Shell::Impl &a, const Shell::Impl &b,
    // 				 const Shell::Impl &c, const Shell::Impl &d,
    // 				 double scale, double *I) {
    // 	for(int l = d.first, ijkl = 0; l <= d.last; ++l) {
    // 	    for(int k = c.first; k <= c.last; ++k) {
    // 		double qkl = scale*NORMALIZE[l]*NORMALIZE[k];
    // 		for(int j = b.first; j <= b.last; ++j) {
    // 		    for(int i = a.first; i <= a.last; ++i, ++ijkl) {
    // 			I[ijkl] *= qkl*NORMALIZE[j]*NORMALIZE[i];
    // 		    }
    // 		}
    // 	    }
    // 	}    
    // }

}

#endif // _RYSQ_NORMALIZE_H
