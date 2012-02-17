/*
Sparse cube fie -> regular Gaussian-like cube file convertor
written by D.G. Fedorov, AIST, Japan
version 1.1  1.12.04 minor changes in scube format
version 1.0 12.11.04

Usage:
compile as "cc -o scubecon scubecon.c"

scubecon [n] <input.scube >output.cube
n is optional number of passes (Default: 1).
*/

#include <stdio.h>
#include <stdlib.h>

#define min(a,b) ((a)>(b)?(b):(a))

int main(int argc, char*argv[])
{
	int npass=argc<2?1:atoi(argv[1]);
	char title[80+1];
	double *v,*vi,val,z,c[3];
	double origin[3],u1[3],u2[3],u3[3]; 
	int nat,nfg,n1,n2,n3,n1p,ip,ifg,i1,i2,i3,mini1,maxi1,
	    min1,max1,min2,max2,min3,max3,r=0,iifg,
		 maxgx,maxgy,maxgz,ilay;
	long pos;

	fgets(title,80,stdin); puts(title);
	r+=scanf("%d%d%lf%lf%lf",&nfg,&nat,origin,origin+1,origin+2);
	r+=scanf("%d%d%lf%lf%lf",&maxgx,&n1,u1,u1+1,u1+2);
	r+=scanf("%d%d%lf%lf%lf",&maxgy,&n2,u2,u2+1,u2+2);
	r+=scanf("%d%d%lf%lf%lf",&maxgz,&n3,u3,u3+1,u3+2);
	if(r!=20) return 1; 

	printf("%5d%12.6f%12.6f%12.6f\n",nat,origin[0],origin[1],origin[2]);
	printf("%5d%12.6f%12.6f%12.6f\n",n1,u1[0],u1[1],u1[2]);
	printf("%5d%12.6f%12.6f%12.6f\n",n2,u2[0],u2[1],u2[2]);
	printf("%5d%12.6f%12.6f%12.6f\n",n3,u3[0],u3[1],u3[2]);

	for(i1=0;i1<nat;i1++)
	  if(scanf("%d%d%lf%lf%lf%lf",&i2,&i3,&z,c,c+1,c+2)!=6) return 2;
	  else printf("%5d%12.6f%12.6f%12.6f%12.6f\n",i3,z,c[0],c[1],c[2]); 
	
/*	if(feof(stdin) || ferror(stdin)) return 2; */

	pos=ftell(stdin);
	n1p=(n1-1)/npass+1;

	if((v=malloc(sizeof(*v)*n1p*n2*n3))==NULL) return 3; 

	for(ip=1;ip<=npass;ip++)
	{
		if(npass!=1) fprintf(stderr,"Doing pass %d\r\n",ip);
	   mini1=(ip-1)*n1p+1;
	   maxi1=min(ip*n1p,n1);
		for(vi=v,i1=0;i1<n1p*n2*n3;i1++) *vi++=0;
		for(ifg=1;ifg<=nfg;ifg++)
		{
			if(scanf("%d%d%d%d%d%d%d%d",
						&iifg,&ilay,&min1,&max1,&min2,&max2,&min3,&max3)!=8)
				return 4;
		   fprintf(stderr,"Doing frg %d\r",iifg);
			for(i1=min1;i1<=max1;i1++)
				for(i2=min2;i2<=max2;i2++)
					for(i3=min3;i3<=max3;i3++)
					{
			   		if(scanf("%lf",&val)!=1) return 5;
						if(i1>=mini1 && i1<=maxi1)
							v[((i1-mini1)*n2+i2-1)*n3+i3-1]+=val;
					}
		}
		if(npass!=1) fprintf(stderr,"Writing data...\r\n");
		for(i1=mini1;i1<=maxi1;i1++)
			for(i2=1;i2<=n2;i2++)
				for(i3=1;i3<=n3;i3++)
				{
					printf("%13.5e",v[((i1-mini1)*n2+i2-1)*n3+i3-1]);
			      if((i3-1)%6==5 || i3==n3) putchar('\n');
				}
		if(ip<npass) fseek(stdin,pos,SEEK_SET);
	}
	fprintf(stderr,"Converted %d fragments in %d pass(es).\r\n",nfg,npass);
	free(v);
}
