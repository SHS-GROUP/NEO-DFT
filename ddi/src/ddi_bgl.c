#include "ddi_base.h"

#include <rts.h>
#include <bglpersonality.h>

/** @see ddi_bgl.h */
int BGLRank_pset_comparator(const void *a, const void *b) {
    int pset_a = ((BGLRank*)a)->pset;
    int pset_b = ((BGLRank*)b)->pset;

    if (pset_a > pset_b) return 1;
    if (pset_a == pset_b) return 0;
    if (pset_a < pset_b) return -1;
}

/** @see ddi_bgl.h */
int DDI_BGL_Comm_divide_custom(int nGroups, int *nNodesByGroup, int commId, int *newCommId) {
    const DDI_Comm *comm = (const DDI_Comm*)Comm_find(commId);
    
    int i, j, rankIndex;
    BGLPersonality personality;
    BGLRank bglRank;
    BGLRank bglRanks[comm->np];
    int numNodesInPset, nPsets;
    int aligned;

    int *ids = NULL;
    int group, nGroupRanks;
    int nodeIndex,groupIndex, nextGroupIndex;

    bglRank.rank = comm->me;
    rts_get_personality(&personality, sizeof(personality));
    bglRank.pset = BGLPersonality_psetNum(&personality);
    numNodesInPset = BGLPersonality_numNodesInPset(&personality);

    /* Gather and sort ranks by pset */
    MPI_Allgather(&bglRank, sizeof(BGLRank), MPI_CHAR, bglRanks, sizeof(BGLRank), MPI_CHAR, comm->compute_comm);
    qsort(bglRanks, comm->np, sizeof(BGLRank), BGLRank_pset_comparator);
    
    /* determine if groups can be aligned to psets. */
    aligned = 1;
    nPsets = 0;
    for (i = 0; i < comm->np;) {
	for (j = i + 1; j < comm->np; ++j)
	    if (bglRanks[i].pset != bglRanks[j].pset) break;
	if (j - i != numNodesInPset) {
	    aligned = 0;
	    break;
	}
	++nPsets;
	i = j;
    }
    if (nPsets % nGroups != 0) aligned = 0;
    if (aligned == 0) return 1;

    /* locate rank */
    for (i = 0; i < comm->np; ++i) {
	if (bglRanks[i].rank == comm->me) {
	    nodeIndex = i;
	    break;
	}
    }

    /* determine the group and group boundaries */
    groupIndex = 0;
    nextGroupIndex = 0;
    for (i = 0; i < nGroups; ++i) {
	nextGroupIndex += nNodesByGroup[i];
	if (groupIndex  <= nodeIndex && nodeIndex < nextGroupIndex) {
	    group = i;
	    break;
	}
	groupIndex = nextGroupIndex;
    }
    /* No SMP in BG/L */
    nGroupRanks = nNodesByGroup[group];

    /* build list of group ids */
    ids = (int*)malloc(nGroupRanks * sizeof(int));
    for (i = 0, j = groupIndex; j < nextGroupIndex; ++i, ++j)
	ids[i] = bglRanks[j].rank;

    /*
    for (i = 0; i < comm->np; ++i) {
	if (comm->me == i) {
	    printf("%s: ids = [", DDI_Id());
	    for (j = 0; j < nGroupRanks; ++j) {
		if (j > 0) printf(", ");
		printf("%i", ids[j]);
	    }
	    printf("]\n");
	    fflush(stdout);
	}
	Comm_sync(126,comm);
    }
    */

    Comm_create(nGroupRanks, ids, nGroups, group, commId, newCommId);

    free(ids);

    return 0;
}

/** @see ddi_bgl.h */
int DDI_BGL_File_IO_rank(DDI_Comm *comm) {
    BGLPersonality personality;
    int numComputePerIO;

    rts_get_personality(&personality, sizeof(personality));
    return BGLPersonality_rankInPset(&personality);
}

/** @see ddi_bgl.h */
int DDI_BGL_File_IONode_rank(DDI_Comm *comm) {
    BGLPersonality personality;
    
    rts_get_personality(&personality, sizeof(personality));
    return BGLPersonality_psetNum(&personality);
}

/** @see ddi_bgl.h */
void DDI_BGL_Runtime_print(FILE *stream) {
    BGLPersonality personality;
    int rank,nprocs;
    int dim = 0,torus_dim = 0;
    char topology[] = "torus";
    char *topology_axis, mesh_axis[] = "XYZ", torus_axis[] = "XYZ";
    char *bglmpi_eager = NULL;
    char *bglmpi_mapping = NULL;
    char *bglmpi_pacing = NULL;
    
    rts_get_personality(&personality, sizeof(personality));
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /* determine mesh */
    dim = 0;
    strcpy(topology,"mesh");
    strcpy(mesh_axis,"");
    if (BGLPersonality_xSize(&personality) > 1) { ++dim; strcat(mesh_axis,"X"); }
    if (BGLPersonality_ySize(&personality) > 1) { ++dim; strcat(mesh_axis,"Y"); }
    if (BGLPersonality_zSize(&personality) > 1) { ++dim; strcat(mesh_axis,"Z"); }
    if (dim == 0) { dim = 1; strcpy(mesh_axis,"X"); }
    topology_axis = mesh_axis;
    
    /* determine torus */
    torus_dim = 0;
    strcpy(torus_axis,"");
    if (BGLPersonality_isTorusX(&personality)) { ++torus_dim; strcat(torus_axis,"X"); }
    if (BGLPersonality_isTorusY(&personality)) { ++torus_dim; strcat(torus_axis,"Y"); }
    if (BGLPersonality_isTorusZ(&personality)) { ++torus_dim; strcat(torus_axis,"Z"); }
    if (torus_dim > 0) { dim = torus_dim; strcpy(topology,"torus"); topology_axis = torus_axis; }

    /* determine BGLMPI_MAPPING */
    bglmpi_eager = getenv("BGLMPI_EAGER");
    bglmpi_mapping = getenv("BGLMPI_MAPPING");
    bglmpi_pacing = getenv("BGLMPI_PACING");
    
    /* print DDI Posix runtime */
    DDI_POSIX_Runtime_print(stream);

    /* print DDI Blue Gene/L runtime */
    fprintf(stream,"%i compute nodes, %s mode, %i I/O nodes\n",
	    BGLPersonality_numComputeNodes(&personality),
	    BGLPersonality_virtualNodeMode(&personality) ? "VN" : "CO",
	    BGLPersonality_numIONodes(&personality));
    fprintf(stream,"%i-D %s(%s) <%i,%i,%i>\n",
	    dim,topology,topology_axis,
	    BGLPersonality_xSize(&personality),
	    BGLPersonality_ySize(&personality),
	    BGLPersonality_zSize(&personality));
    if (bglmpi_eager) fprintf(stream,"BGLMPI_EAGER=%s\n", bglmpi_eager);
    if (bglmpi_mapping) fprintf(stream,"BGLMPI_MAPPING=%s\n", bglmpi_mapping);
    if (bglmpi_pacing) fprintf(stream,"BGLMPI_PACING=%s\n", bglmpi_pacing);
    fprintf(stream,"MPI %i/%i <%i,%i,%i> %iMHz %iMB\n",
	    rank,nprocs,
	    BGLPersonality_xCoord(&personality),
	    BGLPersonality_yCoord(&personality),
	    BGLPersonality_zCoord(&personality),
	    BGLPersonality_clockHz(&personality)/1000000,
	    BGLPersonality_DDRSize(&personality)/(1024*1024));

}
