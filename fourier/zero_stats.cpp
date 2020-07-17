#include <iostream>
#include <ostream>
#include <fstream>

#define ZEROS_FILE "zeros6"

using namespace std;

/* zeros stats for Andy */
/* 24/June/2008 */

int main()
{
    ifstream zeros_file;

    zeros_file.open(ZEROS_FILE);
    if (!zeros_file.is_open())
    {
        cout << "Error opening " << ZEROS_FILE << endl;
        exit(0);
    };
   return 0;
};
