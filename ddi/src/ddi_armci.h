#ifndef DDI_ARMCI_H
#define DDI_ARMCI_H

/**
   @file
   @brief DDI ARMCI one-sided communication implementation.
   @note The ARMCI one-sided communication implementation is meant to be used alongside
   the MPI-1 point-to-point implementation.
   @note The distributed array communicator is the communicator in which the distributed
   array was created.
   @note The distributed array process id is the id of the process where the distributed
   array segment is local.  The distributed array process id is the distributed array
   communicator process id.
   @warning The implementation macros, fields, and functions should not be accessed by the
   users of the DDI library directly.
*/

#if !defined DDI_MPI
#error "DDI ARMCI implementation requires MPI"
#endif

#include <mpi.h>
#include <armci.h>

#ifdef DOXYGEN_DDI_ARMCI_LOCK /* define DDI_ARMCI_LOCK only for Doxygen */
/**
   @brief Implicitly lock distributed data arrays.
   @details If DDI_ARMCI_LOCK is defined,
   DDI_ARMCI_Acquire(int handle, int pid, int ltype, void **array, int *armciPid) and
   DDI_ARMCI_Release(int handle, int pid, int ltype) will call DDI_ARMCI_Lock(int mutex, int pid) and
   DDI_ARMCI_Unlock(int mutex, int pid) respectively.
   @note DDI_ARMCI_LOCK is not defined by default.
   @see DDI_ARMCI_Acquire(int handle, int pid, int ltype, void **array, int *armciPid)
   @see DDI_ARMCI_Release(int handle, int pid, int ltype)
   @see DDI_ARMCI_Lock(int mutex, int pid)
   @see DDI_ARMCI_Unlock(int mutex, int pid)
*/
#define DDI_ARMCI_LOCK
#endif /* DOXYGEN_DDI_ARMCI_LOCK */

#ifdef DOXYGEN_DDI_ARMCI_IMPLICIT_NBPUT /* define DDI_ARMCI_IMPLICIT_NBPUT only for Doxygen */
/**
   @brief Implicit non-blocking put operations.
   @details If DDI_ARMCI_IMPLICIT_NBPUT is defined and DDI_ARMCI_IMPLICIT_WAIT is undefined,
   DDI ARMCI put operations are implicitly non-blocking.
   If DDI_ARMCI_IMPLICIT_NBPUT and DDI_ARMCI_IMPLICIT_WAIT are defined,
   DDI_ARMCI_Put(DDI_Patch *patch, void *buf) will wait for all outstanding non-blocking
   operations to complete.
   @note DDI_ARMCI_IMPLICIT_NBPUT is not defined by default.
   @see DDI_ARMCI_Put(DDI_Patch *patch, void *buf)
   @see DDI_ARMCI_Put_domain_SMP(DDI_Patch *patch, void *buf, int domain)
   @see DDI_ARMCI_Put_proc(DDI_Patch *patch, void *buf, int pid)
   @see DDI_ARMCI_Barrier(MPI_Comm comm)
   @see DDI_ARMCI_IMPLICIT_WAIT
*/
#define DDI_ARMCI_IMPLICIT_NBPUT
#endif /* DOXYGEN_DDI_ARMCI_IMPLICIT_NBPUT */

#ifdef DOXYGEN_DDI_ARMCI_IMPLICIT_NBGET /* define DDI_ARMCI_IMPLICIT_NBGET only for Doxygen */
/**
   @brief Implicit non-blocking get operations.
   @details If DDI_ARMCI_IMPLICIT_NBGET is defined and DDI_ARMCI_IMPLICIT_WAIT is undefined,
   DDI ARMCI get operations are implicitly non-blocking.
   If DDI_ARMCI_IMPLICIT_NBGET and DDI_ARMCI_IMPLICIT_WAIT are defined,
   DDI_ARMCI_Get(DDI_Patch *patch, void *buf) will wait for all outstanding non-blocking
   operations to complete.
   @note DDI_ARMCI_IMPLICIT_NBGET is not defined by default.
   @see DDI_ARMCI_Get(DDI_Patch *patch, void *buf)
   @see DDI_ARMCI_Get_domain_SMP(DDI_Patch *patch, void *buf, int domain)
   @see DDI_ARMCI_Get_proc(DDI_Patch *patch, void *buf, int pid)
   @see DDI_ARMCI_Barrier(MPI_Comm comm)
   @see DDI_ARMCI_IMPLICIT_WAIT
*/
#define DDI_ARMCI_IMPLICIT_NBGET
#endif /* DOXYGEN_DDI_ARMCI_IMPLICIT_NBGET */

#ifdef DOXYGEN_DDI_ARMCI_IMPLICIT_NBGET /* define DDI_ARMCI_IMPLICIT_NBGET only for Doxygen */
/**
   @brief Implicit non-blocking acc operations.
   @details If DDI_ARMCI_IMPLICIT_NBACC is defined and DDI_ARMCI_IMPLICIT_WAIT is undefined,
   DDI ARMCI acc operations are implicitly non-blocking.
   If DDI_ARMCI_IMPLICIT_NBACC and DDI_ARMCI_IMPLICIT_WAIT are defined,
   DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf) will wait for all outstanding
   non-blocking operations to complete.
   @note DDI_ARMCI_IMPLICIT_NBACC is not defined by default.
   @see DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf)
   @see DDI_ARMCI_Acc_domain_SMP(DDI_Patch *patch, double alpha, void *buf, int domain)
   @see DDI_ARMCI_Acc_proc(DDI_Patch *patch, double alpha, void *buf, int pid)
   @see DDI_ARMCI_Barrier(MPI_Comm comm)
   @see DDI_ARMCI_IMPLICIT_WAIT
*/
#define DDI_ARMCI_IMPLICIT_NBACC
#endif /* DOXYGEN_DDI_ARMCI_IMPLICIT_NBACC */

#ifdef DOXYGEN_DDI_ARMCI_IMPLICIT_WAIT /* define DDI_ARMCI_IMPLICIT_WAIT only for Doxygen */
/**
   @brief Implicitly wait for non-blocking operation completion.
   @details DDI_ARMCI_Put(DDI_Patch *patch, void *buf), 
   DDI_ARMCI_Get(DDI_Patch *patch, void *buf), and
   DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf) may call multiple
   DDI_ARMCI_Put_proc(DDI_Patch *patch, void *buf, int pid),
   DDI_ARMCI_Get_proc(DDI_Patch *patch, void *buf, int pid), and
   DDI_ARMCI_Acc_proc(DDI_Patch *patch, double alpha, void *buf, int pid) operations respectively.
   DDI_ARMCI_IMPLICIT_WAIT ensures that all those multiple operations have completed.
   @note DDI_ARMCI_IMPLICIT_WAIT is not defined by default.
   @see DDI_ARMCI_Put(DDI_Patch *patch, void *buf)
   @see DDI_ARMCI_Get(DDI_Patch *patch, void *buf)
   @see DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf)
 */
#define DDI_ARMCI_IMPLICIT_WAIT
#endif /* DOXYGEN_DDI_ARMCI_IMPLICIT_WAIT */

/** @brief DDI_ARMCI_Init(int argc, char **argv) error. */
#define DDI_ARMCI_INIT_ERROR            2000
/** @brief DDI_ARMCI_MPI_Init(int argc, char **argv) error. */
#define DDI_ARMCI_MPI_INIT_ERROR        2001
/** @brief DDI_ARMCI_Memory_init(size_t size) error. */
#define DDI_ARMCI_MEMORY_INIT_ERROR     2002
/** @brief DDI_ARMCI_Memory_finalize() error. */
#define DDI_ARMCI_MEMORY_FINALIZE_ERROR 2003

/** @brief DDI_ARMCI_Init(int argc, char **argv) error message. */
#define DDI_ARMCI_INIT_ERROR_MESSAGE            NULL
/** @brief DDI_ARMCI_MPI_Init(int argc, char **argv) error message. */
#define DDI_ARMCI_MPI_INIT_ERROR_MESSAGE        NULL
/** @brief DDI_ARMCI_Memory_init(size_t size) error message. */
#define DDI_ARMCI_MEMORY_INIT_ERROR_MESSAGE     NULL
/** @brief DDI_ARMCI_Memory_finalize() error message. */
#define DDI_ARMCI_MEMORY_FINALIZE_ERROR_MESSAGE NULL

void* gv(armci_mem_addr)[MAX_PROCESSORS];     /**< @brief Remote ARMCI base addresses */
void* gv(armci_dlb_counter)[MAX_PROCESSORS];  /**< @brief Dynamic load balancer remote ARMCI addresses. */
void* gv(armci_gdlb_counter)[MAX_PROCESSORS]; /**< @brief Global dynamic load balancer remote ARMCI addresses. */

/** @brief Distributed array remote memory offset and mutex. */
typedef struct {
  size_t offset;  /**< @brief Distributed array offset. */
  int mutex;      /**< @brief Distributed array mutex. */
} DDA_Remote_Index;

/** @brief DDA_Remote_Index structures indexed by distributed array and distributed array process. */
DDA_Remote_Index gv(dda_remote_index)[MAX_DD_ARRAYS][MAX_PROCESSORS];

/**
   @brief Initializes DDI ARMCI implementation.
   @details Initializes ARMCI runtime and DDI structures such topology information, etc.
   @note This operation is collective in MPI_COMM_WORLD.
   @param[in] argc The number of arguments.
   @param[in] argv Pointer to the argument vector.
 */
void DDI_ARMCI_Init(int argc, char **argv);

/**
   @brief Initializes ARMCI MPI.
   @note This operation is collective in MPI_COMM_WORLD.
   @param[in] argc The number of arguments.
   @param[in] argv Pointer to the argument vector.
 */
void DDI_ARMCI_MPI_Init(int argc, char **argv);

/**
   @brief Determines MPI_COMM_WORLD rank.
   @param[out] rank MPI_COMM_WORLD rank of the calling process.
 */
void DDI_ARMCI_MPI_global_rank(int *rank);

/**
   @brief Initializes DDI communicator.
   @note This operation is collective in MPI_COMM_WORLD.
   @param[in] np The number of ARMCI processes.
   @param[in] pids The array of ARMCI process ids.
   @param[in] pid Calling process ARMCI process id.
   @param[in] domainNp The number of processes in the ARMCI_DOMAIN_SMP domain of the calling process.
   @param[in] domainPids The array of ARMCI process ids in the ARMCI_DOMAIN_SMP domain of the calling process.
   @param[in] domainPid Calling process id in the ARMCI_DOMAIN_SMP domain.
   @param[out] comm DDI communicator.
 */
void DDI_ARMCI_MPI_Comm_init(int np, int *pids, int pid, int domainNp,
			     int *domainPids, int domainPid, DDI_Comm *comm);
/**
   @brief Finalizes ARMCI runtime.
   @note This operation is collective in MPI_COMM_WORLD.
 */
void DDI_ARMCI_Finalize();

/**
   @brief Allocates ARMCI distributed memory, creates mutexes, and initializes DDI memory variables.
   @note This operation is collective in MPI_COMM_WORLD.
   @param[in] size Memory size in bytes
*/
void DDI_ARMCI_Memory_init(size_t size);

/**
   @brief Frees ARMCI distributed memory and destroys mutexes.
   @note This operation is collective in MPI_COMM_WORLD.
 */
void DDI_ARMCI_Memory_finalize();

/**
   @brief Initializes remote ARMCI counter addresses and resets local counters to 0.
*/
void DDI_ARMCI_Counters_init();

/**
   @brief Resets the dynamic load balancer to 0.
   @note This operation has effect only if called by the process 0 of DDI_WORKING_COMM.
 */
void DDI_ARMCI_DLBReset();

/**
   @brief Stores and then increments the dynamic load balancer.
   @param[out] counter Dynamic load balancer value.
 */
void DDI_ARMCI_DLBNext(size_t *counter);

/**
   @brief Resets the global dynamic load balancer to 0.
   @note This operation has effect only if called by the rank 0 of DDI_COMM_WORLD.
*/
void DDI_ARMCI_GDLBReset();

/**
   @brief Stores on all processes in DDI_WORKING_COMM and then increments
   the global dynamic load balancer.
   @note This operation is collective in DDI_WORKING_COMM.
   @param[out] counter Global dynamic load balancer value.
*/
void DDI_ARMCI_GDLBNext(size_t *counter);

/**
   @brief Acquires the distributed array on the specified process.
   @note If armciPid != NULL, ARMCI process id is returned.
   @note If pid < 0, id of the calling process is assumed.
   @note Lock type is ignored.
   @note DDI_ARMCI_Lock(int pid, int mutex) is called if DDI_ARMCI_LOCK is defined.
   @see DDI_ARMCI_LOCK
   @param[in] handle Distributed array handle.
   @param[in] pid Distributed array process id.
   @param[in] ltype Lock type.
   @param[out] array Pointer to the distributed array remote ARMCI address.
   @param[out] armciPid ARMCI process id.
 */
void DDI_ARMCI_Acquire(int handle, int pid, int ltype, void **array, int *armciPid);

/**
   @brief Releases the distributed array lock on the specified process.
   @note If pid < 0, id of the calling process is assumed.
   @note Lock type is ignored.
   @note DDI_ARMCI_Unlock(int pid, int mutex) is called if DDI_ARMCI_LOCK is defined.
   @see DDI_ARMCI_LOCK
   @param[in] handle Distributed array handle.
   @param[in] pid Distributed array process id.
   @param[in] ltype Lock type.
 */
void DDI_ARMCI_Release(int handle, int pid, int ltype);

/**
   @brief Acquires mutex on the ARMCI process.
   @param[in] mutex Mutex.
   @param[in] pid ARMCI process id.
 */
inline void DDI_ARMCI_Lock(int mutex, int pid);

/**
   @brief Releases mutex on the ARMCI process.
   @param[in] mutex Mutex.
   @param[in] pid ARMCI process id.
 */
inline void DDI_ARMCI_Unlock(int mutex, int pid);

/**
   @brief Synchronizes ARMCI processes and memory in the specified MPI communicator.
   @details Combines ARMCI_WaitAll() and MPI_Barrier(MPI_Comm comm);
   @note This operation is collective in the specified MPI communicator.
   @param[in] comm MPI communicator.
 */
void DDI_ARMCI_Barrier(MPI_Comm comm);

/**
   @brief One-sided put operation.
   @details Copies the local buffer to the distributed array segment.
   @warning This operation may or may not be blocking depending on the values of
   DDI_ARMCI_IMPLICIT_NBPUT and DDI_ARMCI_IMPLICIT_WAIT.
   @see DDI_ARMCI_IMPLICIT_NBPUT
   @see DDI_ARMCI_IMPLICIT_WAIT
   @param[in] patch Distributed array segment.
   @param[in] buf Local buffer.
 */
void DDI_ARMCI_Put(DDI_Patch *patch, void *buf);

/**
   @brief One-sided get operation.
   @details Copies the distributed array segment to the local buffer.
   @warning This operation may or may not be blocking depending on the values of
   DDI_ARMCI_IMPLICIT_NBGET and DDI_ARMCI_IMPLICIT_WAIT.
   @see DDI_ARMCI_IMPLICIT_NBGET
   @see DDI_ARMCI_IMPLICIT_WAIT
   @param[in] patch Distributed array segment.
   @param[out] buf Local buffer.
 */
void DDI_ARMCI_Get(DDI_Patch *patch, void *buf);

/**
   @brief One-sided accumulate operation.
   @details Adds the local buffer scaled by alpha to the distributed array segment.
   @note Local buffer is not changed by this operation.
   @warning This operation may or may not be blocking depending on the values of
   DDI_ARMCI_IMPLICIT_NBACC and DDI_ARMCI_IMPLICIT_WAIT.
   @see DDI_ARMCI_IMPLICIT_NBACC
   @see DDI_ARMCI_IMPLICIT_WAIT
   @param[in] patch Distributed array segment.
   @param[in] alpha Scale factor.
   @param[in] buf Local buffer.
 */
void DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf);

/**
   @brief One-sided atomic get/accumulate operation.
   @details Adds the local buffer to the distributed array segment and stores the original
   distributed array segment in the local buffer.
   @warning This operation should only be called in the context of
   DDI_ProcDLB_next(int handle, int pid, int *counter).
   @see DDI_ProcDLB_next(int handle, int pid, int *counter)
   @param[in] patch Distributed array segment.
   @param[in,out] buf Local buffer.
 */
void DDI_ARMCI_GetAcc(DDI_Patch *patch, void *buf);

/**
   @brief One-sided get operation from the specified ARMCI_DOMAIN_SMP domain.
   @details Copies the distributed array segment local to the specified ARMCI_DOMAIN_SMP
   domain to the local buffer.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBGET.
   @see DDI_ARMCI_IMPLICIT_NBGET
   @param[in] patch Distributed array segment.
   @param[out] buf Local buffer.
   @param[in] domain ARMCI_DOMAIN_SMP domain.
   @return Number of ARMCI get operations.
 */
inline int DDI_ARMCI_Get_domain_SMP(DDI_Patch *patch, void *buf, int domain);

/**
   @brief One-sided put operation to the specified ARMCI_DOMAIN_SMP domain.
   @details Copies the local buffer to the distributed array segment local to the specified
   ARMCI_DOMAIN_SMP domain.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBPUT.
   @see DDI_ARMCI_IMPLICIT_NBPUT
   @param[in] patch Distributed array segment.
   @param[in] buf Local buffer.
   @param[in] domain ARMCI_DOMAIN_SMP domain.
   @return Number of ARMCI put operations.
 */
inline int DDI_ARMCI_Put_domain_SMP(DDI_Patch *patch, void *buf, int domain);

/**
   @brief One-sided accumulate operation to the specified ARMCI_DOMAIN_SMP domain.
   @details Adds the local buffer scaled by alpha to the distributed array segment local
   to the specified ARMCI_DOMAIN_SMP domain.
   @note Local buffer is not changed by this operation.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBACC.
   @see DDI_ARMCI_IMPLICIT_NBACC
   @param[in] patch Distributed array segment.
   @param[in] alpha Scale factor.
   @param[in] buf Local buffer.
   @param[in] domain ARMCI_DOMAIN_SMP domain.
   @return Number of ARMCI accumulate operations.
 */
inline int DDI_ARMCI_Acc_domain_SMP(DDI_Patch *patch, double alpha, void *buf, int domain);

/**
   @brief One-sided get operation from the specified ARMCI process.
   @details Copies the distributed array segment local to the specified ARMCI process
   to the local buffer.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBGET.
   @see DDI_ARMCI_IMPLICIT_NBGET
   @param[in] patch Distributed array segment.
   @param[out] buf Local buffer.
   @param[in] pid ARMCI process id.
   @return Number of ARMCI get operations.
*/
inline int DDI_ARMCI_Get_proc(DDI_Patch *patch, void *buf, int pid);

/**
   @brief One-sided put operation to the specified ARMCI process.
   @details Copies the local buffer to the distributed array segment local to the specified
   ARMCI process.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBPUT.
   @see DDI_ARMCI_IMPLICIT_NBPUT
   @param[in] patch Distributed array segment.
   @param[in] buf Local buffer.
   @param[in] pid ARMCI process id.
   @return Number of ARMCI put operations.
 */
inline int DDI_ARMCI_Put_proc(DDI_Patch *patch, void *buf, int pid);

/**
   @brief One-sided accumulate operation to the specified ARMCI process.
   @details Adds the local buffer scaled by alpha to the distributed array segment local
   to the specified ARMCI process.
   @warning This operation may or may not be blocking depending on the value of
   DDI_ARMCI_IMPLICIT_NBACC.
   @see DDI_ARMCI_IMPLICIT_NBACC
   @param[in] patch Distributed array segment.
   @param[in] alpha Scale factor.
   @param[in] buf Local buffer.
   @param[in] pid ARMCI process id.
   @return Number of ARMCI accumulate operations.
 */
inline int DDI_ARMCI_Acc_proc(DDI_Patch *patch, double alpha, void *buf, int pid);

#endif /* DDI_ARMCI_H */
