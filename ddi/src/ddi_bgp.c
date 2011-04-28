/*   GDF 10.08.09  new source code for BG/P, we are indebted to Vitali Morozov (LCF)   */

#include "ddi_bgp.h"
#include "ddi_runtime.h"

#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>

/*   GDF 10.08.09  DDI_BGL_File_IO_rank and DDI_BGL_File_IONode_rank were not used??  */

/** print DDI Blue Gene/P runtime */
void DDI_BGP_Runtime_print(FILE *stream) {
    _BGP_Personality_t personality;
    int rank,nprocs;
    int dim = 0,torus_dim = 0;
    char topology[] = "torus";
    char *topology_axis, mesh_axis[] = "XYZ", torus_axis[] = "XYZ";
    char *bgpmpi_eager = NULL;
    char *bgpmpi_mapping = NULL;
    
    Kernel_GetPersonality(&personality, sizeof(personality));
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /* determine mesh */
    dim = 0;
    strcpy(topology,"mesh");
    strcpy(mesh_axis,"");
    if (BGP_Personality_xSize(&personality) > 1) { ++dim; strcat(mesh_axis,"X"); }
    if (BGP_Personality_ySize(&personality) > 1) { ++dim; strcat(mesh_axis,"Y"); }
    if (BGP_Personality_zSize(&personality) > 1) { ++dim; strcat(mesh_axis,"Z"); }
    if (dim == 0) { dim = 1; strcpy(mesh_axis,"X"); }
    topology_axis = mesh_axis;
    
    /* determine torus */
    torus_dim = 0;
    strcpy(torus_axis,"");
    if (BGP_Personality_isTorusX(&personality)) { ++torus_dim; strcat(torus_axis,"X"); }
    if (BGP_Personality_isTorusY(&personality)) { ++torus_dim; strcat(torus_axis,"Y"); }
    if (BGP_Personality_isTorusZ(&personality)) { ++torus_dim; strcat(torus_axis,"Z"); }
    if (torus_dim > 0) { dim = torus_dim; strcpy(topology,"torus"); topology_axis = torus_axis; }

    /* determine BGPMPI_MAPPING */
    bgpmpi_eager = getenv("DCMF_EAGER");   /* GDF 12.09.09 what if DCMF not used ? */
    bgpmpi_mapping = getenv("BG_MAPPING");
    
    /* print DDI Posix runtime */
    DDI_POSIX_Runtime_print(stream);

    /* print DDI Blue Gene/P runtime */
    fprintf(stream,"%i compute nodes, %s mode, %i I/O nodes\n",
            BGP_Personality_numComputeNodes(&personality),

/* GDF 12.08.09  cannot find this, replaced with 'processConfig' ?
            BGP_Personality_virtualNodeMode(&personality) ? "VN" : "CO",
 */
            BGP_Personality_processConfig(&personality),

            BGP_Personality_numIONodes(&personality));
    fprintf(stream,"%i-D %s(%s) <%i,%i,%i>\n",
    dim,topology,topology_axis,
            BGP_Personality_xSize(&personality),
            BGP_Personality_ySize(&personality),
            BGP_Personality_zSize(&personality));
    if (bgpmpi_eager) fprintf(stream,"BGPMPI_EAGER=%s\n", bgpmpi_eager);
    if (bgpmpi_mapping) fprintf(stream,"BGPMPI_MAPPING=%s\n", bgpmpi_mapping);
    fprintf(stream,"MPI %i/%i <%i,%i,%i> %iMHz %iMB\n",
    rank,nprocs,
            BGP_Personality_xCoord(&personality),
            BGP_Personality_yCoord(&personality),
            BGP_Personality_zCoord(&personality),

/* GDF 12.08.09  these two have been re-named and converted
            BGP_Personality_clockHz(&personality)/1000000,
            BGP_Personality_DDRSize(&personality)/(1024*1024));
 */
            BGP_Personality_clockMHz(&personality),
            BGP_Personality_DDRSizeMB(&personality));

}
