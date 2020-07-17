#include <iostream>
#include <cmath>

using namespace std;

double psi(double x) {
    if(x < 0) x = -x;
    if(x < .125) return 1;
    if(x > .25) return 0;
    x -= .125;
    x *= 8.0;
    return 1.0 - x;
}

int main(int argc, char ** argv) {
  if(argc!=6)
    {
      printf("Usage: %s <Q> <start> <end> <delta> <d>.\n",argv[0]);
      return 0;
    }
    int Q = atoi(argv[1]);
    double st = atof(argv[2]);
    double en = atof(argv[3]);
    double delta = atof(argv[4]);
    double d = atof(argv[5]);

    for(double xi = st; xi <= en; xi += delta) {
        double S = 0.0;
        for(int q = 1; q <= Q; q++) {
            int m = (int)floor(q*xi + .25);
            //if(m < 0 || m > q-1) continue;
            //for(int m = 0; m < q-1; m++) {
            //    S1 += psi(q*xi - m);
            //}
            S += psi(q*xi - m) * pow(q, -d/2.0);
        }
        cout << xi << " " << S << endl;
    }
    return 0;
}
