/* ----- this sets up running GA over MPI: used only with LIBCCHEM ----- */

#if defined MPI
#include <mpi.h>
void gamess_mpi_init_() {
    int provided;
    int argc = 0;
    char **argv;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    // printf("MPI_Init_thread returned %i\n", provided);
}
#endif
