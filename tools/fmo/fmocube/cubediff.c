/*
 * Utility to subtract two Gaussian-like cube files
 * written by D.G. Fedorov, AIST, Japan
 * version 1.0 12.11.04
 *
 * Usage:
 * compile as "cc -o cubediff cubediff.c"
 *
 * cubediff input1.cube input2.cube >input1-2.cube
 * */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char*argv[])
{
	int i,j,m,n;
	double diff,maxd=0,v1,v2;
	FILE *f1,*f2;
   char title[80+1];
   double z,c[3];
   double origin[3],u1[3],u2[3],u3[3];
	int nat,n1,n2,n3,r;

	if(argc<3 || (f1=fopen(argv[1],"r"))==NULL || (f2=fopen(argv[2],"r"))==NULL) return; 

   fgets(title,80,f1);
   r=fscanf(f1,"%d%lf%lf%lf",&nat,origin,origin+1,origin+2)
    +fscanf(f1,"%d%lf%lf%lf",&n1,u1,u1+1,u1+2)
    +fscanf(f1,"%d%lf%lf%lf",&n2,u2,u2+1,u2+2)
    +fscanf(f1,"%d%lf%lf%lf",&n3,u3,u3+1,u3+2);
   if(r!=16) return 1;

   fgets(title,80,f2);
	r=fscanf(f2,"%d%lf%lf%lf",&nat,origin,origin+1,origin+2)
	 +fscanf(f2,"%d%lf%lf%lf",&n1,u1,u1+1,u1+2)
	 +fscanf(f2,"%d%lf%lf%lf",&n2,u2,u2+1,u2+2)
	 +fscanf(f2,"%d%lf%lf%lf",&n3,u3,u3+1,u3+2);
	if(r!=16) return 1;

   puts(title);
   printf("%5d%12.6f%12.6f%12.6f\n",nat,origin[0],origin[1],origin[2]);
   printf("%5d%12.6f%12.6f%12.6f\n",n1,u1[0],u1[1],u1[2]);
   printf("%5d%12.6f%12.6f%12.6f\n",n2,u2[0],u2[1],u2[2]);
   printf("%5d%12.6f%12.6f%12.6f\n",n3,u3[0],u3[1],u3[2]);

   for(i=0;i<nat;i++)
     if(fscanf(f1,"%d%lf%lf%lf%lf",&j,&z,c,c+1,c+2)+
		  fscanf(f2,"%d%lf%lf%lf%lf",&j,&z,c,c+1,c+2)!=10) return 2;
     else printf("%5d%12.6f%12.6f%12.6f%12.6f\n",j,z,c[0],c[1],c[2]);

	m=n2*n1; 
	n=n3;

	for(i=0;i<m;i++)	
	{
		for(j=0;j<n;j++)
		{
			if(fscanf(f1,"%lf",&v1)+fscanf(f2,"%lf",&v2)!=2) return 3;
/*			fprintf(stderr,"%f %f\n",v1,v2);	*/
			diff=v1-v2;
			if(fabs(diff)>maxd) maxd=fabs(diff);
			printf("%13.5e",diff); 
			if(j%6==5 || j==n-1) putchar('\n');
		}
	}
	fprintf(stderr,"Max diff=%.5e\n",maxd); 
	fclose(f2);
	fclose(f1);
/*	free(v2);
	free(v1); */
}
