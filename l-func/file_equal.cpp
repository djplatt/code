#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{
  ifstream f1,f2;
  char ch1,ch2;

  if(argc!=3)
    {
      printf("Usage file_equal <f1> <f2>.\n");
      return(0);
    };

  f1.open(argv[1]);
  if(!f1.is_open())
    {
      printf("Failed to open %s.\n",argv[1]);
      return(0);
    };
  f2.open(argv[2]);
  if(!f2.is_open())
    {
      printf("Failed to open %s.\n",argv[2]);
      return(0);
    };

  while(true)
    {
      f1.read(&ch1,1);
      f2.read(&ch2,1);
      if(ch1!=ch2)
	{
	  printf("difference found.\n");
	  return(0);
	};
      if(f1.eof())
	{
	  if(f2.eof())
	    {
	      printf("they are the same.\n");
	      return(0);
	    }
	  else
	    {
	      printf("difference found.\n");
	      return(0);
	    }
	}
    };
  return(0);
}
