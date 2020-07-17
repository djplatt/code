#include "stdio.h"
#include "stdlib.h"

#define ifname ("L_half_data.dat")
#define ofname ("L_half_data1.dat")


void printz(long unsigned int i, double *z)
{
  printf("%11lu - %f %f\n",i,z[0],z[1]);
} 

int main()
     {
       double z[2];
       long unsigned int P=(((long unsigned int) 1<<34) + 25); // a large prime
       long unsigned int PHI=P-1;
       
       FILE *outfile,*infiler,*infilei;
       if(!(outfile=fopen(ofname,"wb")))
	 {
	   printf("Error opening file %s for binary output. Exiting.\n",ofname);
	   exit(0);
	 }
       if(!(infiler=fopen(ifname,"rb")))
	 {
	   printf("Error opening file %s for binary input. Exiting.\n",ifname);
	   exit(0);
	 }
       if(!(infilei=fopen(ifname,"rb")))
	 {
	   printf("Error opening file %s for binary input. Exiting.\n",ifname);
	   exit(0);
	 }

       fread(z,sizeof(double),1,infiler);
       z[1]=0.0;
       fwrite(z,sizeof(double),2,outfile);

       printz(0,z);

       if(fseek(infilei,-sizeof(double)*3,SEEK_END)!=0)
	 printf("fseek failed\n");

       for(unsigned long int i=1;i<PHI/2;i++)
	 {

	   fread(z,sizeof(double),1,infiler);
	   fread(&z[1],sizeof(double),1,infilei);
	   fseek(infilei,-sizeof(double)*2,SEEK_CUR);
	   if((i&0x3fffff)==0)
	     printz(i,z);
	   fwrite(z,sizeof(double),2,outfile);
	 }

       fread(z,sizeof(double),1,infiler);
       z[1]=0.0;
       fwrite(z,sizeof(double),2,outfile);
       printz(PHI/2,z);
       return(0);
     };
