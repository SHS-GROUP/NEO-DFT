/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#include "normalize.h"


/**@defgroup const Constant
	@internal 
	@{
*/

/** @brief Cartesian normalization constants
	@see rysq.h
*/
const double rysq::NORMALIZE[] = {
	RYSQ_NORMAL_000,
	RYSQ_NORMAL_100, RYSQ_NORMAL_100, RYSQ_NORMAL_100,
	RYSQ_NORMAL_200, RYSQ_NORMAL_200, RYSQ_NORMAL_200,
	RYSQ_NORMAL_110, RYSQ_NORMAL_110, RYSQ_NORMAL_110,
	RYSQ_NORMAL_300, RYSQ_NORMAL_300, RYSQ_NORMAL_300,
	RYSQ_NORMAL_210, RYSQ_NORMAL_210, RYSQ_NORMAL_210, RYSQ_NORMAL_210, RYSQ_NORMAL_210, RYSQ_NORMAL_210,
	RYSQ_NORMAL_111,
	RYSQ_NORMAL_400, RYSQ_NORMAL_400, RYSQ_NORMAL_400,
	RYSQ_NORMAL_310, RYSQ_NORMAL_310, RYSQ_NORMAL_310, RYSQ_NORMAL_310, RYSQ_NORMAL_310, RYSQ_NORMAL_310,
	RYSQ_NORMAL_220, RYSQ_NORMAL_220, RYSQ_NORMAL_220,
	RYSQ_NORMAL_211, RYSQ_NORMAL_211, RYSQ_NORMAL_211};

/**  @brief Cartesian function l values, L = l+m+n 
	@see rysq.h
*/
const int rysq::LX[] = {
	0,
	1, 0, 0,
	2, 0, 0,
	1, 1, 0,
	3, 0, 0,
	2, 2, 1, 0, 1, 0,
	1,
	4, 0, 0,
	3, 3, 1, 0, 1, 0,
	2, 2, 0,
	2, 1, 1};

/**  @brief Cartesian function m values, L = l+m+n 
 	@see rysq.h
 */
const int rysq::LY[] = {
	0,
	0, 1, 0,
	0, 2, 0,
	1, 0, 1,
	0, 3, 0,
	1, 0, 2, 2, 0, 1,
	1,
	0, 4, 0,
	1, 0, 3, 3, 0, 1,
	2, 0, 2,
	1, 2, 1};

/**  @brief Cartesian function n values, L = l+m+n 
	 @see rysq.h
*/
const int rysq::LZ[] = {
	0,
	0, 0, 1,
	0, 0, 2,
	0, 1, 1,
	0, 0, 3,
	0, 1, 0, 1, 2, 2,
	1,
	0, 0, 4,
	0, 1, 0, 1, 3, 3,
	0, 2, 2,
	1, 1, 2};

/** @}*/
