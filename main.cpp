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
#define DEPTH 40
#define STATISTICS 50
#define AVG 100

using namespace std;

int main()
{

	freopen("NewProgLog.txt", "w", stdout);
	srand( (unsigned)time(NULL) );
	QuESTEnv env = createQuESTEnv();
	RandQAlg alg(env, COLS, LINES, DEPTH);

	alg.setSingleErrRate(0.0005);
	alg.setMultiErrRate(0.005);
	alg.setEnvCoupling(0);

	static constexpr int ampNum = (1 << (COLS * LINES));
	double amps[ampNum];

	cout << fixed << "{";
	for(int j = 0; j < STATISTICS; j++)
	{
		cerr << "Step " << j + 1 << " of "<< STATISTICS << "\n";

		alg.generate();

		for(int i = 0; i < ampNum; i++)
			amps[i] = 0;
		for(int k = 0; k < AVG; k++)
		{
			alg.evaluate();
			for(int i = 0; i < ampNum; i++)
				amps[i] += alg.getSquaredAmp(i);
		}
		for(int i = 0; i < ampNum; i++)
		{
			cout << amps[i] / AVG;
			if(i + 1 < ampNum)
				cout << ", ";
		}
	}
	cout << "}\n";

	return 0;
}
