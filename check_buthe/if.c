/*2:*/
#line 92 "./if.w"

#include <time.h> 
#include <unistd.h> 
#include <sys/times.h> 
#include <sys/stat.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <stdarg.h> 
#include <string.h> 
#include <limits.h> 

#include "asm/config.h"
#include "if.h"

int verbose= 0;
static ulong used_cols,ncol= 80;

#ifdef BIGENDIAN
u32_t bswap32(u32_t);
u64_t bswap64(u64_t);
#endif

/*:2*//*3:*/
#line 115 "./if.w"

void*xmalloc(size_t size)
{
char*x;
if(size==0)return NULL;
if((x= malloc(size))==NULL)complain("xmalloc %llu: %m",(u64_t)size);
return x;
}

/*:3*//*4:*/
#line 125 "./if.w"

void*xvalloc(size_t size)
{
char*x;
if(size==0)return NULL;
if((x= valloc(size))==NULL)complain("xvalloc: %m");
return x;
}

/*:4*//*5:*/
#line 135 "./if.w"

void*
xcalloc(size_t n,size_t s)
{
void*p;
if(n==0||s==0)return NULL;
if((p= calloc(n,s))==NULL)complain("calloc: %m");
return p;
}

/*:5*//*6:*/
#line 146 "./if.w"

void*xrealloc(void*x,size_t size)
{
char*y;
if(size==0){
if(x!=NULL)free(x);
return NULL;
}
if((y= realloc(x,size))==NULL&&size!=0)complain("xrealloc: %m");
return y;
}

/*:6*//*7:*/
#line 159 "./if.w"

FILE*
xfopen(const char*fn,const char*format)
{
FILE*res;
res= fopen(fn,format);
if(res==NULL)complain("Cannot open %s in format %s: %m\n",fn,format);
return res;
}

/*:7*//*8:*/
#line 170 "./if.w"

FILE*
xyzopen4r(char*fn,int*is_pipe)
{
FILE*r;
struct stat sb;
char*fn1;

r= fopen(fn,"r");
if(r!=NULL){
*is_pipe= 0;
return r;
}
asprintf(&fn1,"%s.gz",fn);
if(stat(fn1,&sb)==0){
free(fn1);
asprintf(&fn1,"gunzip --stdout %s.gz",fn);
r= popen(fn1,"r");
free(fn1);
}else{
free(fn1);
asprintf(&fn1,"%s.bz2",fn);
if(stat(fn1,&sb)==0){
free(fn1);
asprintf(&fn1,"bunzip2 --stdout %s.bz2",fn);
r= popen(fn1,"r");
free(fn1);
}
}
if(r==NULL)
complain("Cannot open or decompress %s : %m\n",fn);
*is_pipe= 1;
return r;
}

/*:8*//*9:*/
#line 207 "./if.w"

FILE*logfile= NULL;

void complain(char*fmt,...)
{
va_list arglist;
va_start(arglist,fmt);
gmp_vfprintf(stderr,fmt,arglist);
if(logfile!=NULL)gmp_vfprintf(logfile,fmt,arglist);
exit(1);
}

/*:9*//*10:*/
#line 220 "./if.w"

void Schlendrian(char*fmt,...)
{
va_list arglist;
va_start(arglist,fmt);
gmp_vfprintf(stderr,fmt,arglist);
if(logfile!=NULL)gmp_vfprintf(logfile,fmt,arglist);
abort();
}

/*:10*//*11:*/
#line 231 "./if.w"

#define NEW_NUMBER -0x10000
void logbook(int l,char*fmt,...)
{
if(l==NEW_NUMBER){
used_cols= 0;
return;
}
if(l<verbose){
va_list arglist;
char*output_str;
ulong sl;
va_start(arglist,fmt);
vasprintf(&output_str,fmt,arglist);
sl= strlen(output_str);
if(used_cols+sl> ncol){
fprintf(stderr,"\n");
if(logfile!=NULL)fprintf(logfile,"\n");
used_cols= 0;
}
fputs(output_str,stderr);
if(logfile!=NULL)fputs(output_str,logfile);
if(output_str[sl-1]=='\n')used_cols= 0;
else used_cols+= sl;
free(output_str);
}
}

/*:11*//*12:*/
#line 261 "./if.w"

int
errprintf(char*fmt,...)
{
va_list arglist;
int res;

va_start(arglist,fmt);
if(logfile!=NULL)gmp_vfprintf(logfile,fmt,arglist);
res= gmp_vfprintf(stderr,fmt,arglist);
return res;
}

/*:12*//*13:*/
#line 275 "./if.w"

void adjust_bufsize(void**buf,size_t*alloc,size_t req,
size_t incr,size_t item_size)
{
if(req> *alloc){
size_t new_alloc;
new_alloc= *alloc+incr*((req+incr-1-*alloc)/incr);
if(*alloc> 0)*buf= xrealloc(*buf,new_alloc*item_size);
else*buf= xmalloc(new_alloc*item_size);
*alloc= new_alloc;
}
}

/*:13*//*14:*/
#line 291 "./if.w"

int yn_query(char*fmt,...)
{
va_list arglist;
char answer[10];

va_start(arglist,fmt);
if(logfile!=NULL)gmp_vfprintf(logfile,fmt,arglist);
gmp_vfprintf(stderr,fmt,arglist);
if(!isatty(STDIN_FILENO)||!isatty(STDERR_FILENO))return 0;

fflush(stderr);
while(scanf("%9s",answer)!=1||
(strcasecmp(answer,"yes")!=0&&strcasecmp(answer,"no")!=0)){
fprintf(stderr,"Please answer yes or no!\n");
gmp_vfprintf(stderr,fmt,arglist);
}
if(strcasecmp(answer,"yes")==0)return 1;
return 0;
}


/*:14*//*15:*/
#line 320 "./if.w"

#if 0
ssize_t
skip_blanks_comments(char**iline,size_t*iline_alloc,FILE*ifi)
{
while(getline(iline,iline_alloc,ifi)> 0){
if(**iline!='#'&&strspn(*iline,"\n\t ")<strlen(*iline))
return 1;
}
return 0;
}
#endif

/*:15*//*16:*/
#line 334 "./if.w"

#ifdef BIGENDIAN
/*17:*/
#line 341 "./if.w"

static u32_t
bswap_32(u32_t x)
{
return((x&0x000000ffUL)<<24)|((x&0x0000ff00UL)<<8)|
((x&0x00ff0000UL)>>8)|((x&0xff000000UL)>>24);
}

/*:17*//*18:*/
#line 350 "./if.w"

static u64_t
bswap_64(u64_t x)
{
return((x&0xffULL)<<56)|((x&0xff00ULL)<<40)|((x&0xff0000ULL)<<24)|
((x&0xff000000ULL)<<8)|((x&0xff00000000ULL)>>8)|
((x&0xff0000000000ULL)>>24)|((x&0xff000000000000ULL)>>40)|
((x&0xff00000000000000ULL)>>56);
}

/*:18*//*19:*/
#line 361 "./if.w"

size_t
write_i64(FILE*ofile,i64_t*buffer,size_t count)
{
size_t i,res;

for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:19*//*20:*/
#line 376 "./if.w"

size_t
write_u64(FILE*ofile,u64_t*buffer,size_t count)
{
size_t i,res;

for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:20*//*21:*/
#line 391 "./if.w"

size_t
write_i32(FILE*ofile,i32_t*buffer,size_t count)
{
size_t i,res;

for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:21*//*22:*/
#line 406 "./if.w"

size_t
write_u32(FILE*ofile,u32_t*buffer,size_t count)
{
size_t i,res;

for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:22*//*23:*/
#line 421 "./if.w"

size_t
write_u16(FILE*ofile,u16_t*buffer,size_t count)
{
size_t res;

swab(buffer,buffer,count*sizeof(*buffer));
res= fwrite(buffer,sizeof(*buffer),count,ofile);
swab(buffer,buffer,count*sizeof(*buffer));
return res;
}

/*:23*//*24:*/
#line 434 "./if.w"

size_t
read_i64(FILE*ifile,i64_t*buffer,size_t count)
{
size_t i,res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:24*//*25:*/
#line 447 "./if.w"

size_t
read_u64(FILE*ifile,u64_t*buffer,size_t count)
{
size_t i,res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:25*//*26:*/
#line 460 "./if.w"

size_t
read_i32(FILE*ifile,i32_t*buffer,size_t count)
{
size_t i,res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:26*//*27:*/
#line 473 "./if.w"

size_t
read_u32(FILE*ifile,u32_t*buffer,size_t count)
{
size_t i,res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}


/*:27*//*28:*/
#line 487 "./if.w"

size_t
read_u16(FILE*ifile,u16_t*buffer,size_t count)
{
size_t res;

res= fread(buffer,sizeof(*buffer),count,ifile);
swab(buffer,buffer,count*sizeof(*buffer));
return res;
}

/*:28*/
#line 336 "./if.w"

#endif 


/*:16*//*29:*/
#line 500 "./if.w"

#ifdef NEED_GETLINE
#define GETL_INCR 128
ssize_t
getline(char**lineptr,size_t*n,FILE*stream)
{
int rv;

if(*n==0){
*n= GETL_INCR;
*lineptr= xmalloc(*n);
}

rv= 0;
for(;;){
int m;
m= *n-rv;
if(fgets(*lineptr+rv,m-1,stream)==NULL)break;
rv= strlen(*lineptr);
if(rv==0||(*lineptr)[rv-1]=='\n')break;
*n+= GETL_INCR;
*lineptr= xrealloc(*lineptr,*n);
}
return rv;
}
#endif

/*:29*//*30:*/
#line 528 "./if.w"

#ifdef NEED_ASPRINTF
int
vasprintf(char**ptr,const char*template,va_list ap)
{

int n,size= 32;

while(1){
*ptr= xmalloc(size);

n= vsnprintf(*ptr,size,template,ap);
if(1+strlen(*ptr)<size)return n;
size*= 2;
free(*ptr);
}
}
#endif

/*:30*//*31:*/
#line 548 "./if.w"

#ifdef NEED_ASPRINTF
int asprintf(char**ptr,const char*template,...)
{
int rv;
va_list ap;

va_start(ap,template);
rv= vasprintf(ptr,template,ap);
va_end(ap);
return rv;
}
#endif

/*:31*//*32:*/
#line 563 "./if.w"

#ifdef NEED_FNMATCH
int fnmatch(char*s,char*fname,int dummy)
{
while(*fname)fname++;
if(*(fname--)!='z')return 0;
if(*(fname--)!='g')return 0;
if(*(fname--)!='.')return 0;
return 1;
}
#endif

/*:32*//*33:*/
#line 577 "./if.w"

int
u32_cmp012(const void*x,const void*y)
{
const u32_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return-1;
if(*xx> *yy)return 1;
return 0;
}

/*:33*//*34:*/
#line 590 "./if.w"

int
u32_cmp210(const void*x,const void*y)
{
const u32_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return 1;
if(*xx> *yy)return-1;
return 0;
}

/*:34*//*35:*/
#line 603 "./if.w"

int
u64_cmp012(const void*x,const void*y)
{
const u64_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return-1;
if(*xx> *yy)return 1;
return 0;
}

/*:35*//*36:*/
#line 616 "./if.w"

int
u64_cmp210(const void*x,const void*y)
{
const u64_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return 1;
if(*xx> *yy)return-1;
return 0;
}/*:36*/
