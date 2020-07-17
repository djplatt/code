#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- sem_make <sem_name> <no_resources>\n");
      exit(1);
    }

  sem_t *this_sem=sem_open(argv[1],O_CREAT,S_IRUSR|S_IWUSR,atoi(argv[2]));

  if(this_sem==SEM_FAILED)
    {
      printf("Failed to create semaphore %s with %d resources.\n",argv[1],atoi(argv[2]));
      perror("Error reported was: ");
      exit(1);
    }

  exit(0);
}
