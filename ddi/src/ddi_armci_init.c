#include "ddi_base.h"

/** @see ddi_armci.h */
void DDI_ARMCI_Init(int argc,char **argv) {
  int code;
  int pid;
  int i,j,k;
  int nprocs,*pids;
  int domain_my_id;
  int domainNp, *domainPids, domainPid, domainCount;
  int myds;
  DDI_Comm *comm = (DDI_Comm*)&gv(ddi_base_comm);

  /* initialize ARMCI message passing */
  DDI_ARMCI_MPI_Init(argc, argv);
  DDI_ARMCI_MPI_global_rank(&pid);

  code = ARMCI_Init();
  if (code != 0) {
    fprintf(DDI_STDERR, "MPI process %i: ARMCI_Init() returned %i\n", pid, code);
    DDI_Error(DDI_ARMCI_INIT_ERROR, DDI_ARMCI_INIT_ERROR_MESSAGE);
  }
  if (pid == 0) fprintf(DDI_STDERR, "DDI ARMCI initialized\n");
  fflush(DDI_STDERR);

  Init_scratch(argc, argv);
  
  /* ARMCI_DOMAIN_SMP domain */
  domainNp = armci_domain_nprocs(ARMCI_DOMAIN_SMP, -1);
  domain_my_id = armci_domain_my_id(ARMCI_DOMAIN_SMP);
  domainPids = (int*)Malloc(sizeof(int)*domainNp);
  for (i = 0; i < domainNp; ++i) {
    domainPids[i] = armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP,domain_my_id,i);
    if (domainPids[i] == pid) domainPid = i;
  }

  /* number of ARMCI_DOMAIN_SMP domains and total number of processes */
  domainCount = armci_domain_count(ARMCI_DOMAIN_SMP);
  nprocs = 0;
  for (i = 0; i < domainCount; ++i) nprocs += armci_domain_nprocs(ARMCI_DOMAIN_SMP, i);
  
  /* compile time limits */
#if FULL_SMP
  if (domainNp > MAX_SMP_PROCS) DDI_Error(DDI_MAX_SMP_PROCS_ERROR, NULL);
  if (domainCount > MAX_NODES) DDI_Error(DDI_MAX_NODES_ERROR, NULL);
#else
  if (nprocs > MAX_NODES) DDI_Error(DDI_MAX_NODES_ERROR, NULL);
#endif
  if (nprocs > MAX_PROCESSORS) DDI_Error(DDI_MAX_PROCESSORS_ERROR, NULL);

  /* process ids ordered by SMP domains */
  pids = (int*)Malloc(sizeof(int)*nprocs);
  for (i = 0, k = 0; i < domainCount; ++i)
    for (j = 0; j < armci_domain_nprocs(ARMCI_DOMAIN_SMP, i); ++j, ++k) {
      pids[k] = armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP, i, j);
      gv(ddiprocs)[k].node = i;
    }

  /* node information */
  for (i = 0; i < domainCount; ++i) {
    myds = (pid % armci_domain_nprocs(ARMCI_DOMAIN_SMP,i));    
    gv(ddinodes)[i].cpus       = armci_domain_nprocs(ARMCI_DOMAIN_SMP,i);
    gv(ddinodes)[i].myds       = armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP, i, myds);
    gv(ddinodes)[i].nodemaster = armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP, i, 0);
  }

  /* initialize communicator */
  DDI_ARMCI_MPI_Comm_init(nprocs, pids, pid, domainNp, domainPids, domainPid, comm);
}

/** @see ddi_armci.h */
void DDI_ARMCI_MPI_Init(int argc,char *argv[]) {
  int code;

  code = MPI_Init(&argc, &argv);
  if (code != MPI_SUCCESS) {
    fprintf(DDI_STDERR, "MPI_Init(%p, %p) returned %i\n", &argc, &argv);
    DDI_Error(DDI_ARMCI_MPI_INIT_ERROR, DDI_ARMCI_MPI_INIT_ERROR_MESSAGE);
  }
}

/** @see ddi_armci.h */
void DDI_ARMCI_MPI_global_rank(int *rank) {
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

/** @see ddi_armci.h */
void DDI_ARMCI_MPI_Comm_init(int nprocs, int *pids, int pid, int domainNp,
			     int *domainPids, int domainPid, DDI_Comm *comm) {
  int domain_master = 0;
  int domain_my_id = armci_domain_my_id(ARMCI_DOMAIN_SMP);
  int i;

  MPI_Group Comm_World_grp;
  MPI_Group SMP_World_grp;
  MPI_Group SMP_Compute_grp;
  MPI_Group DDI_World_grp;
  MPI_Group DDI_Compute_grp;
  MPI_Comm SMP_World_comm;
  MPI_Comm SMP_Compute_comm;
  MPI_Comm SMP_Masters_comm;
  MPI_Comm DDI_World_comm;
  MPI_Comm DDI_Compute_comm;

  /* SMP domain communicator */
  MPI_Comm_group(MPI_COMM_WORLD, &Comm_World_grp);
  MPI_Group_incl(Comm_World_grp, domainNp, domainPids, &SMP_World_grp);
  MPI_Comm_create(MPI_COMM_WORLD, SMP_World_grp, &SMP_World_comm);
  MPI_Barrier(MPI_COMM_WORLD);

  /* SMP domain masters communicator */
  if (pid == armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP, domain_my_id, 0)) domain_master = 1;
  MPI_Comm_split(MPI_COMM_WORLD, domain_master, 0, &SMP_Masters_comm);
  MPI_Barrier(MPI_COMM_WORLD);

  /* DDI_Compute_comm communicator */
  MPI_Group_incl(Comm_World_grp, nprocs, pids, &DDI_Compute_grp);
  MPI_Comm_create(MPI_COMM_WORLD, DDI_Compute_grp, &DDI_Compute_comm);
  MPI_Barrier(MPI_COMM_WORLD);

  /* DDI_World_comm communicator */
  MPI_Group_incl(Comm_World_grp, nprocs, pids, &DDI_World_grp);
  MPI_Comm_create(MPI_COMM_WORLD, DDI_World_grp, &DDI_World_comm);
  MPI_Barrier(MPI_COMM_WORLD);

  /* SMP_Compute_comm communicator */
  MPI_Group_intersection(DDI_Compute_grp, SMP_World_grp, &SMP_Compute_grp);
  MPI_Comm_create(MPI_COMM_WORLD, SMP_Compute_grp, &SMP_Compute_comm);
  MPI_Barrier(MPI_COMM_WORLD);

  /* populate comm */
  comm->np           = nprocs;
  comm->me           = pid;
  comm->id           = DDI_COMM_WORLD;
  MPI_Comm_dup(DDI_World_comm, &(comm->world_comm));
  MPI_Comm_dup(DDI_Compute_comm, &(comm->compute_comm));

#if FULL_SMP
  comm->nn           = armci_domain_count(ARMCI_DOMAIN_SMP);
  comm->my           = domain_my_id;
  comm->np_local     = domainNp;
  comm->me_local     = domainPid;
  MPI_Comm_dup(SMP_Compute_comm, &(comm->smp_comm));
  MPI_Comm_dup(SMP_Masters_comm, &(comm->node_comm));
#else
  comm->nn       = comm->np;
  comm->my       = comm->me;
  comm->np_local = 1;
  comm->me_local = 0;
  MPI_Comm_dup(MPI_COMM_SELF, &(comm->smp_comm));
  MPI_Comm_dup(DDI_Compute_comm, &(comm->node_comm));
#endif

  /* Communicators have been duplicated and can be freed */
  MPI_Comm_free(&SMP_World_comm);
  MPI_Comm_free(&SMP_Compute_comm);
  MPI_Comm_free(&SMP_Masters_comm);
  MPI_Comm_free(&DDI_World_comm);
  MPI_Comm_free(&DDI_Compute_comm);
}
