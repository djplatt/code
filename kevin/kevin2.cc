#include <iostream>
#include <cmath>

using namespace std;

double psi(double x) {
  printf("psi(%f) returning ",x);
    if(x < 0) x = -x;
    if(x <= .5) {printf("1\n");return 1;}
    printf("0\n");return 0;
    if(x > .25) {printf("0\n");return 0;}
    x -= .125;
    x *= 8.0;
    printf("%f\n",1.0-x);
    return 1.0 - x;
}

int main(int argc, char ** argv) {
    int J = atoi(argv[1]);
    double delta = atof(argv[2]);
    double d = atof(argv[3]);
    //printf("%d\n",sizeof(int));
    for(double xi = 0.0; xi < 1.0; xi += delta) {
        double S = 0.0;
        for(int j = 0; j <= J; j++) {
            for(int q = (1u << j); q < (1u << (j+1)); q++) {
                int m = (int)floor((1u<<(3*j))*xi + .25);
                //if(m < 0 || m > q-1) continue;
                //for(int m = 0; m < q-1; m++) {
                //    S1 += psi(q*xi - m);
                //}
		//printf("Adding j=%d q=%d\n",j,q);
                S += psi((1u<<(3*j))*xi - m) * pow(q, -d/2.0);
            }}
	    //printf("%e %e\n",xi,S);
	    cout << xi << " " << S << endl;
	    //}
    }
    return 0;
}
