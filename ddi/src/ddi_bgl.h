#ifndef DDI_BGL_H
#define DDI_BGL_H

#include <stdio.h>

/**
   @file
   @brief IBM Blue Gene/L DDI port specifics.
   @warning These macros, fields, and functions should not be accessed by the
   users of the DDI library directly.
 */

/**
   @brief Hooks DDI_BGL_File_IO_rank(DDI_Comm *comm) into DDI_File_IO_rank(DDI_Comm *comm).
   @see DDI_BGL_File_IO_rank(DDI_Comm *comm)
   @see DDI_File_IO_rank(DDI_Comm *comm)
 */
#define DDI_FILE_IO_RANK_(c)     DDI_BGL_File_IO_rank((c))

/**
   @brief Hooks DDI_BGL_File_IONode_rank(DDI_Comm *comm) into DDI_File_IONode_rank(DDI_Comm *comm).
   @see DDI_BGL_File_IONode_rank(DDI_Comm *comm)
   @see DDI_File_IONode_rank(DDI_Comm *comm)
 */
#define DDI_FILE_IONODE_RANK_(c) DDI_BGL_File_IONode_rank((c))

/**
   @brief Hooks DDI_BGL_Runtime_print(FILE *stream) into DDI_Runtime_print(FILE *stream).
   @see DDI_BGL_Runtime_print(FILE *stream)
   @see DDI_Runtime_print(FILE *stream)
 */
#define DDI_RUNTIME_PRINT_(s)    DDI_BGL_Runtime_print(s)

/** @brief BG/L MPI rank and pset */
typedef struct {
  int rank;  /**< @brief MPI rank */
  int pset;  /**< @brief BG/L pset */
} BGLRank;

/**
   @brief Compares psets of two BGLRanks.
   @param[in] a BGLRank a. 
   @param[in] b BGLRank b.
   @return 1 if pset of a is greater, 0 if two psets are equal, -1 if pset of b is greater.
*/
int BGLRank_pset_comparator(const void *a, const void *b);

/**
   @brief BG/L pset-aware DDI communicator division
   @details If groups can be aligned to BG/L psets, the function makes a call to
   Comm_create(int np, int *ids, int ngroups, int mygroup, int commId, int *newCommId) and
   returns 0.  If groups cannot be aligned to psets, either because comm is pset-incomplete
   or psets cannot be divided among groups evenly, function returns 1 and does not call
   Comm_create.
   @param[in] nGroups The number of groups.
   @param[in] nNodesByGroup Array of the number of nodes by group.
   @param[in] commId Communicator id.
   @param[out] newCommId New communicator id.
   @return 0 if groups are aligned to BG/L psets and the new communicator has been created,
   -1 otherwise.
 */
int DDI_BGL_Comm_divide_custom(int nGroups, int *nNodesByGroup, int commId, int *newCommId);

/**
   @brief Returns rank in pset which can be used to determine I/O topology.
   @param[in] comm DDI communicator pointer.  Ignored.
   @return Rank in pset.
*/
int DDI_BGL_File_IO_rank(DDI_Comm *comm);

/**
   @brief Returns pset rank which can be used to determine I/O topology. 
   @param[in] comm DDI communicator pointer.  Ignored.
   @return Pset rank.
 */
int DDI_BGL_File_IONode_rank(DDI_Comm *comm);

/**
   @brief Prints DDI Blue Gene/L runtime.
   @details Prints POSIX runtime information and Blue Gene/L specific details, such as
   topology, I/O configuration, and environmental variables.
   @param[in] stream File stream.
*/
void DDI_BGL_Runtime_print(FILE *stream);

#endif /* DDI_BGL_H */
