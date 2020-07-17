#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#define SUCCESS 0
#define ZEROS_FILE "zeros6"

using namespace std;

/* zeros stats for Andy */
/* 24/June/2008 */

int main(int argc, char **argv)
{
    ifstream zeros_file;
    int zero_ptr,last_one_zero;
    double this_zero,last_zero,diff,max_diff;

   if(argc==1)
      max_diff=1.0;
   if(argc==2)
   {
      stringstream(argv[1]) >> max_diff;
      if(max_diff<1.0)
      {
         cout << "Difference must be >= 1.0" << endl;
         max_diff=1.0;
      };
   };
   if(argc>2)
   {
      cout << "Too many arguments, setting difference to 1.0" << endl;
      max_diff=1.0;
   };

    zeros_file.open(ZEROS_FILE);
    if (!zeros_file.is_open())
    {
        cout << "Error opening " << ZEROS_FILE << endl;
        exit(0);
    };

   last_zero=0.0;
   last_one_zero=0;
   for(zero_ptr=1;!zeros_file.eof();zero_ptr++)
   {
      zeros_file >> this_zero;
      diff=this_zero-last_zero;
      if(diff>max_diff)
         last_one_zero=zero_ptr;
      last_zero=this_zero;
      if (!(zero_ptr&8191))    /* print a line evry so often */
          cout << "Zero Number " << zero_ptr << " " << diff << endl;
   };
   cout << "Last zero with gap > " << max_diff << " was at zero number " << last_one_zero << endl;
   cout << "Total zeros read: " << zero_ptr << endl;
   return SUCCESS;
};
