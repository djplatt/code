#include "IntervalsSSE2.h"
#include "IntervalSuggested.h"
#include <iostream>
#include <iomanip>

using namespace RealLib;
using namespace std;

typedef MachineEstimate (*BinaryFuncME)(const MachineEstimate &a, const MachineEstimate &b);
typedef __m128d (*BinaryFuncM128D)(__m128d a, __m128d b);


bool AreSame(const MachineEstimate &a, const MachineEstimate &b)
{
	return _mm_movemask_pd(_mm_cmpeq_pd(a.GetInterval(), b.GetInterval())) == 3;
}

bool TestFunction(BinaryFuncME def, BinaryFuncM128D impl)
{
	const int cnt = 4;
	MachineEstimate a[cnt] =
	{		
		MachineEstimate(1, 2),
		MachineEstimate(3, 4),
		MachineEstimate(-6, -5),
		MachineEstimate(-7, 8)
	};

	for (int i=0;i<cnt;++i)
		for (int j=0;j<cnt;++j) {
			try {
				MachineEstimate v = def(a[i], a[j]);
				MachineEstimate w = impl(a[i].GetInterval(), a[j].GetInterval());
				if (!AreSame(v, w)) 
					return false;
			} catch (RealLibException &e) {
				try {
					MachineEstimate v = def(a[i], a[j]);
					return false;
				} catch (RealLibException &e1) {
					try {
						MachineEstimate w = impl(a[i].GetInterval(), a[j].GetInterval());
						return false;
					} catch (RealLibException &e2) {
						if (typeid(e1) != typeid(e2)) 
							return false;
					}
				}
			}
		}

	return true;
	
}

int main()
{
	MachineEstimate::Computation comp;

	MachineEstimate a(1, 2);
	MachineEstimate b(3, 4);
	MachineEstimate c(-6, -5);
	MachineEstimate d(-7, 8);

#define PRINT(z) try {\
		cout << #z ": "; \
		(z).PrintInterval(cout); \
		cout << endl;\
	} catch (RealLibException &e) { \
		cout << e << endl; \
	}

	PRINT(a);
	PRINT(b);
	PRINT(c);
	PRINT(d);

	PRINT(a+b);
	PRINT(c+d);
	PRINT(a-c);
	PRINT(d-c);

	PRINT(abs(a));
	PRINT(abs(b));
	PRINT(abs(c));
	PRINT(abs(d));

	PRINT(a*a);
	PRINT(sq(a));
	PRINT(a*b);
	PRINT(a*c);
	PRINT(a*d);

	PRINT(b*a);
	PRINT(b*b);
	PRINT(sq(b));
	PRINT(b*c);
	PRINT(b*d);

	PRINT(c*a);
	PRINT(c*b);
	PRINT(c*c);
	PRINT(sq(c));
	PRINT(c*d);

	PRINT(d*a);
	PRINT(d*b);
	PRINT(d*c);
	PRINT(d*d);
	PRINT(sq(d));

	PRINT(a/a);
	PRINT(a/c);
	PRINT(a/d);

	PRINT(c/a);
	PRINT(c/c);
	PRINT(c/d);

	PRINT(d/a);
	PRINT(d/c);
	PRINT(d/d);

	PRINT(recip(a));
	PRINT(recip(b));
	PRINT(recip(c));
	PRINT(recip(d));

	cout << boolalpha << endl;
	cout << "ivsub matches operator -: " << TestFunction(operator -, ivsub) << endl;
	cout << "IntervalMul matches operator *: " << TestFunction(operator *, IntervalMul) << endl;
	cout << "IntervalMul2 matches operator *: " << TestFunction(operator *, IntervalMul2) << endl;
	cout << "IntervalDiv matches operator /: " << TestFunction(operator /, IntervalDiv) << endl;
	cout << "IntervalDiv2 matches operator /: " << TestFunction(operator /, IntervalDiv2) << endl;

	return 0;
}
