/*1:*/
#line 18 "./if.w"

#include <stdarg.h> 
#include <stdio.h> 
#include <gmp.h> 

void*xmalloc(size_t size);
void*xvalloc(size_t size);
void*xcalloc(size_t n,size_t s);
void*xrealloc(void*x,size_t size);
FILE*xfopen(const char*,const char*);
FILE*xyzopen4r(char*,int*);
void complain(char*fmt,...)__attribute__((noreturn));
void Schlendrian(char*fmt,...)__attribute__((noreturn));
void logbook(int l,char*fmt,...);
int errprintf(char*fmt,...);
void adjust_bufsize(void**,size_t*,size_t,size_t,size_t);
extern int verbose;
extern FILE*logfile;
#ifdef BIGENDIAN
size_t write_i64(FILE*,i64_t*,size_t);
size_t write_u64(FILE*,u64_t*,size_t);
size_t write_i32(FILE*,i32_t*,size_t);
size_t write_u32(FILE*,u32_t*,size_t);
size_t write_u16(FILE*,u16_t*,size_t);
size_t read_i64(FILE*,i64_t*,size_t);
size_t read_u64(FILE*,u64_t*,size_t);
size_t read_i32(FILE*,i32_t*,size_t);
size_t read_u32(FILE*,u32_t*,size_t);
size_t read_u16(FILE*,u16_t*,size_t);
#else
#define write_i64(ofile,buffer,count) fwrite((void*)(buffer),sizeof(i64_t),count,ofile)
#define write_u64(ofile,buffer,count) fwrite((void*)(buffer),sizeof(u64_t),count,ofile)
#define write_u32(ofile,buffer,count) fwrite((void*)(buffer),sizeof(u32_t),count,ofile)
#define write_u16(ofile,buffer,count) fwrite((void*)(buffer),sizeof(u16_t),count,ofile)
#define write_i32(ofile,buffer,count) fwrite((void*)(buffer),sizeof(i32_t),count,ofile)
#define read_i64(ofile,buffer,count) fread((void*)(buffer),sizeof(i64_t),count,ofile)
#define read_u64(ofile,buffer,count) fread((void*)(buffer),sizeof(u64_t),count,ofile)
#define read_u32(ofile,buffer,count) fread((void*)(buffer),sizeof(u32_t),count,ofile)
#define read_i32(ofile,buffer,count) fread((void*)(buffer),sizeof(i32_t),count,ofile)
#define read_u16(ofile,buffer,count) fread((void*)(buffer),sizeof(u16_t),count,ofile)
#endif 
int yn_query(char*fmt,...);
ssize_t skip_blanks_comments(char**,size_t*,FILE*);
int u32_cmp012(const void*,const void*);
int u32_cmp210(const void*,const void*);
int u64_cmp012(const void*,const void*);
int u64_cmp210(const void*,const void*);

#ifdef NEED_ASPRINTF
int asprintf(char**,const char*,...);
#ifndef NEED_VASPRINTF_DECL
#define NEED_VASPRINTF_DECL
#endif
#endif

#ifdef NEED_VASPRINTF_DECL
extern int vasprintf(char**strp,const char*fmt,va_list ap);
#endif

#ifdef NEED_GETLINE_DECL
ssize_t getline(char**lineptr,size_t*n,FILE*stream);
#endif

#ifdef NEED_GETLINE
ssize_t getline(char**,size_t*,FILE*);
#endif

#ifdef NEED_FNMATCH
int fnmatch(char*,char*,int)
#endif

typedef unsigned long long ullong;

/*:1*/
