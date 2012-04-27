#ifndef RYSQ_TYPES_HPP
#define RYSQ_TYPES_HPP

#include "rysq/config.hpp"

#define RYSQ_TYPES_0 (rysq::S)
#define RYSQ_TYPES_1 RYSQ_TYPES_0(rysq::P)(rysq::SP)
#define RYSQ_TYPES_2 RYSQ_TYPES_1(rysq::D)
#define RYSQ_TYPES_3 RYSQ_TYPES_2(rysq::F)
#define RYSQ_TYPES_4 RYSQ_TYPES_3(rysq::G)
#define RYSQ_TYPES_5 RYSQ_TYPES_4(rysq::H)
#define RYSQ_TYPES_6 RYSQ_TYPES_5(rysq::I)

#if RYSQ_LMAX == 0
#define RYSQ_TYPES RYSQ_TYPES_0
#elif RYSQ_LMAX == 1
#define RYSQ_TYPES RYSQ_TYPES_1
#elif RYSQ_LMAX == 2
#define RYSQ_TYPES RYSQ_TYPES_2
#elif RYSQ_LMAX == 3
#define RYSQ_TYPES RYSQ_TYPES_3
#elif RYSQ_LMAX == 4
#define RYSQ_TYPES RYSQ_TYPES_4
#elif RYSQ_LMAX == 5
#define RYSQ_TYPES RYSQ_TYPES_5
#elif RYSQ_LMAX == 6
#define RYSQ_TYPES RYSQ_TYPES_6
#else
#error "must enable L to " RYSQ_LMAX
#endif


#endif // RYSQ_TYPES_HPP
