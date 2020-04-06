#include <iostream>
#include <QuEST.h>
#include "randqalg.h"
#include <cmath>
#include <string>

#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#define LINES 3
#define COLS 4
#define DEPTH 56
#define STATISTICS 100
#define NOISEAVG 100

using namespace std;

int main()
{

	freopen("NewProgLog.txt", "w", stdout);
	srand( (unsigned)time(NULL) );
	QuESTEnv env = createQuESTEnv();
	RandQAlg alg(env, COLS, LINES, DEPTH);
	alg.setSingleErrRate(0);
	alg.setMultiErrRate(0);

	static constexpr int ampNum = (1 << (COLS * LINES));
	double amps[ampNum];

	cout << fixed << "{";
	for(int j = 0; j < STATISTICS; j++)
	{
		cerr << "Step " << j + 1 << " of "<< STATISTICS << "\n";

		alg.generate();
		alg.evaluate();
		for(int i = 0; i < ampNum; i++)
		{
			cout << alg.getSquaredAmp(i);
			if(i + 1 < ampNum)
				cout << ", ";
		}
	}
	cout << "}\n";

	return 0;
}
