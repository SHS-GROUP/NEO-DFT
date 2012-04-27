/**
@file
@brief librysq fortran bindings
@details The fortran binding provides lower case subroutines with and without trailing underscore for compatibility with different fortran compilers.
The default size of fortran integers is 4 bytes. If RYSQ_WITH_INTEGER8 is defined, the size is 8 bytes. 
*/

#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include <rysq.h>
#include "fortran.h"

extern "C" {

    /**@{*/

 	/**
 	@brief Initializes librysq
 	@see Rysq_init ()
 	*/
    void rysq_init() {
	Rysq_init();
#ifdef RYSQ_CUDA_ENABLE
	cuRysq_init(0);
#endif
	//    Rysq_roots_init();
    }
	
    /*@{*/
    /** @see rysq_init */   
    void rysq_init_();


#pragma weak rysq_init_ = rysq_init

/*@}*/

}
