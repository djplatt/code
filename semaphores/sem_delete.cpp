#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>

int main(int argc, char **argv)
{
if(argc!=2)
{
printf("Usage:- sem_delete <sem_name>\n");
exit(0);
}
sem_t *this_sem=sem_open(argv[1],0);
if(this_sem==SEM_FAILED)
printf("Failed to open semaphore %s.\n",argv[1]);
if(sem_close(this_sem)!=0)
printf("Failed to close semaphore %s.\n",argv[1]);
exit(0);
}
