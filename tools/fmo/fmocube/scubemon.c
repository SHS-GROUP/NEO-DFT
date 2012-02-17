/*
Extract monomer densities in Gaussian-like cube file format from 
a sparse cube file
written by D.G. Fedorov, AIST, Japan
version 1.0 1.12.04

Usage:
compile as "cc -o scubemon scubemon.c"

scubemon file [n]
n is optional number specifying fragments
n=0 all fragments
n=ifg only fragment ifg (Default: 0).

file names: 
input: file.scube 
output: file.ifg.cube , where ifg is either a single number or 
                        running index for all 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define min(a,b) ((a)>(b)?(b):(a))
#define MAXNAME 1024

int main(int argc, char*argv[])
{
	char title[80+1];
	double val;
	double origin[3],u1[3],u2[3],u3[3],origini[3];
	int jfg=argc<3?0:atoi(argv[2]);
	int nat,nati,nfg,n1,n2,n3,ifg,i1,i2,i3,
	    min1,max1,min2,max2,min3,max3,iifg,
		 maxgx,maxgy,maxgz,ilay;
	FILE *f1,*f2;
	int *indat,*ian;
	double *fmozan,*fmoc;
	char fin[MAXNAME],fout[MAXNAME];
	char *scube=".scube",*cube=".cube";

	if(jfg<0) return 1;

	if(argc<2 || 
		(f1=fopen(strncat(strncpy(fin,argv[1],MAXNAME),scube,MAXNAME),"r"))==NULL)
	  	return 2;

	fgets(title,80,f1); title[strlen(title)-1]=0; 

	if(fscanf(f1,"%d%d%lf%lf%lf",&nfg,&nat,origin,origin+1,origin+2)+
	   fscanf(f1,"%d%d%lf%lf%lf",&maxgx,&n1,u1,u1+1,u1+2)+
	   fscanf(f1,"%d%d%lf%lf%lf",&maxgy,&n2,u2,u2+1,u2+2)+
	   fscanf(f1,"%d%d%lf%lf%lf",&maxgz,&n3,u3,u3+1,u3+2)!=20) return 3; 

	if(jfg>nfg) return 1;

	if((indat=malloc(sizeof(*indat)*nat))==NULL ||
		(ian=malloc(sizeof(*ian)*nat))==NULL ||
		(fmozan=malloc(sizeof(*fmozan)*nat))==NULL ||
		(fmoc=malloc(sizeof(*fmoc)*nat*3))==NULL) return 4;

	for(i1=0;i1<nat;i1++)
	  if(fscanf(f1,"%d%d%lf%lf%lf%lf",indat+i1,ian+i1,fmozan+i1,
				  fmoc+i1*3,fmoc+i1*3+1,fmoc+i1*3+2)!=6) return 5;

	for(ifg=1;ifg<=nfg;ifg++)
	{
		fprintf(stderr,"Doing frg %d\r",ifg);

		if(fscanf(f1,"%d%d%d%d%d%d%d%d",
					&iifg,&ilay,&min1,&max1,&min2,&max2,&min3,&max3)!=8) return 6;

		if(iifg==jfg || jfg==0)
		{
		   sprintf(fout,"%s.%d%s",argv[1],ifg,cube);
		   if((f2=fopen(fout,"w"))==NULL) return 7;

			for(i1=0;i1<3;i1++) origini[i1]=origin[i1]+
										(min1-1)*u1[i1]+(min2-1)*u2[i1]+(min3-1)*u3[i1];

			for(nati=0,i1=0;i1<nat;i1++)
				if(indat[i1]==iifg) nati++;

	   	fprintf(f2,"%s monomer %d\n\n",title,iifg);
			fprintf(f2,"%5d%12.6f%12.6f%12.6f\n",nati,
							origini[0],origini[1],origini[2]);
			fprintf(f2,"%5d%12.6f%12.6f%12.6f\n",(max1-min1+1),u1[0],u1[1],u1[2]);
			fprintf(f2,"%5d%12.6f%12.6f%12.6f\n",(max2-min2+1),u2[0],u2[1],u2[2]);
			fprintf(f2,"%5d%12.6f%12.6f%12.6f\n",(max3-min3+1),u3[0],u3[1],u3[2]);

			for(i1=0;i1<nat;i1++)
				if(indat[i1]==iifg) fprintf(f2,"%5d%12.6f%12.6f%12.6f%12.6f\n",
					ian[i1],fmozan[i1],fmoc[i1*3],fmoc[i1*3+1],fmoc[i1*3+2]); 

		}

		for(i1=min1;i1<=max1;i1++)
			for(i2=min2;i2<=max2;i2++)
				for(i3=min3;i3<=max3;i3++)
		   		if(fscanf(f1,"%lf",&val)!=1) return 8;
					else  if(iifg==jfg || jfg==0) 
							{
								fprintf(f2,"%13.5e",val); 
								if((i3-min3)%6==5 || i3==max3) fputc('\n',f2);
							}
	   if(iifg==jfg || jfg==0) fclose(f2);
	}
	free(fmoc); free(fmozan); free(ian); free(indat);
	fclose(f1);
	fprintf(stderr,"Scanned %d fragments.\r\n",nfg);
}
