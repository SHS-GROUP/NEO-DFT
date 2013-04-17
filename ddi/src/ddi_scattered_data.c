/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with determining the distribution by node of
 * a given scattered data segment.
 *
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"


/* -------------------------------------------------------------------- *\
   DDI_Scattered_data(handle,scattered,ibuff,nsubs,ranks,subs)
   ===========================================
   [IN]  handle     - Handle of the DD array.
   [IN]  scattered  - Dimensions of full scatter.
   [IN]  ibuff      - coordinates of scattered data.
   [OUT] nsubs      - Number of subscatters.
   [OUT] ranks      - ranks of nodes for each subscatter.
   [OUT] subs       - Array of subscatters disjoint by node.
   
   A segment of scattered data may be spread across multiple nodes.
   This subroutine determines the set of subscatters disjoint by node
   needed to for the requested scattered.
\* -------------------------------------------------------------------- */
void DDI_Scattered_data(int handle,const DDI_Scattered* scattered,long *ibuff,int *nsubs,int *ranks,DDI_Scattered *subs) {
    DDI_Scattered_comm(handle,scattered,nsubs,ranks,subs,DDI_WORKING_COMM,ibuff);

   }

void DDI_Scattered_comm(int handle,const DDI_Scattered* scattered,int *nsubs,int *ranks,DDI_Scattered *subs,int commid,long * ibuff) {
    int i,j,k,np,me,nsub;
    DDI_Scattered s,*t;
    const DDI_Comm *comm = (const DDI_Comm *) Comm_find(commid);

      np = comm->nn;
      me = comm->my;

      for(i=0,nsub=0; i<np; i++) {

	  DDI_Scattered_NDistribP(handle,i,&s);
	  if(s.ihi > 0){printf("Error: DDI_Scatter_Acc is not implemented for arrays with more than one row.\n\n"); fflush(stdout); Fatal_error(911);}

	  for (j=0; j<scattered->nelem; j++) {

	      if(s.jlo > (ibuff[j]-1) || s.jhi < (ibuff[j]-1)) { } else {
		  ranks[nsub] = i;
		  t = &subs[nsub];
		  t->oper   = scattered->oper;
		  t->handle = handle;
		  t->nelem = scattered->nelem;
		  t->nlocal = 0;
		  
		  //		  		  t->start = j;
		  //		  		  t->end = 0;
		  for(k=j; k<scattered->nelem; k++){
		      if( s.jhi < (ibuff[k]-1) ) {break;} else {t->nlocal += 1;}
		  }
		  //		  		  t->end = t->start + t->nlocal -1;
		  t->size = t->nlocal;
		  t->size *= sizeof(double);
		  ++nsub;
		  break;
	      } //(s.jlo > (ibuff[j]-1) || s.jhi < (ibuff[j]-1))

	  } //(j=0; j<scattered->nelem; j++)

      }

      *nsubs = nsub;
}

/* -------------------------------------------------------------------- *\
   DDI_Quicksort(ibuffer,buffer,m,n)
   ===========================================
   [IN]  ibuff  - coordinates of scattered data.
   [IN]  buff   - scattered data to be operated on.
   [IN]  m      - sorting index
   [IN]  n      - sorting index

   Vanillia flavor of quicksort.
\* -------------------------------------------------------------------- */
void DDI_Quicksort(long *ibuff,double *buff,long m,long n){
    long key,i,j,k;
    if( m < n){
	k = DDI_Choose_Pivot(m,n);
	DDI_I_Swap(&ibuff[m],&ibuff[k]);
	DDI_D_Swap(&buff[m],&buff[k]);
	key = ibuff[m];
	i = m+1;
	j = n;
	while(i <= j){
	    while((i <= n) && (ibuff[i] <= key))i++;
	    while((j >= m) && (ibuff[j] > key) )j--;
	    if( i < j){
		DDI_I_Swap(&ibuff[i],&ibuff[j]);
		DDI_D_Swap(&buff[i],&buff[j]);
	    }
	}
	/* swap two elements */
	DDI_I_Swap(&ibuff[m],&ibuff[j]);
	DDI_D_Swap(&buff[m],&buff[j]);
	
	/* recursively sort the lesser list */
	DDI_Quicksort(ibuff,buff,m,j-1);
	DDI_Quicksort(ibuff,buff,j+1,n);
    }
}
void DDI_I_Swap(long *x,long *y){
    long temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
void DDI_D_Swap(double *x,double *y){
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
long DDI_Choose_Pivot(long i,long j ){
    return((i+j) /2);
}
