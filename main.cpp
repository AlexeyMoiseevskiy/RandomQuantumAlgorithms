#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>

#define LINES 3
#define COLS 4
#define DEPTH 40
#define AVG 1

#define ERRDEPTH 2
#define ERRCORDS {0, 0}

using namespace std;

int main()
{
	freopen("NewProgLog.txt", "w", stdout);
	srand( (unsigned)time(NULL) );

	QubitArray regErr(COLS, LINES);
	regErr.setSingleErrRate(0.000);
	regErr.setMultiErrRate(0.00);
	regErr.setEnvCoupling(0);
	QubitArray regNoErr(COLS, LINES);
	regNoErr.setSingleErrRate(0.000);
	regNoErr.setMultiErrRate(0.00);
	regNoErr.setEnvCoupling(0);
	RandQAlg alg(COLS, LINES, DEPTH);

	const int ampsNum = (1 << (LINES * COLS));
	double ampsErr[DEPTH][ampsNum] = {{0.0}};
	double ampsNoErr[DEPTH][ampsNum] = {{0.0}};

	alg.generate();
	for(int j = 0; j < AVG; j++)
	{
		cerr << "Step " << j + 1 << " of "<< AVG << "\n";

		for(int d = 0; d < DEPTH; d++)
		{
			alg.evaluateLayer(regErr, d);
			alg.evaluateLayer(regNoErr, d);
			if(d == ERRDEPTH)
				regErr.XGate(ERRCORDS);
			for(int amp = 0; amp < ampsNum;	amp++)
			{
				ampsErr[d][amp] += regErr.getSquaredAmp(amp);
				ampsNoErr[d][amp] += regNoErr.getSquaredAmp(amp);
			}
		}
	}

	cout << fixed << "{";
	for(int d = 0; d < DEPTH; d++)
	{
		cout << "{";
		for(int amp = 0; amp < ampsNum;	amp++)
		{
			cout << "{" << ampsNoErr[d][amp] * ampsNum / AVG << ", " << ampsErr[d][amp] * ampsNum / AVG << "}";
			if(amp + 1 < ampsNum)
				cout << ", ";
		}
		cout << "}";
		if(d + 1 < DEPTH)
			cout << ", ";
		else
			cout << "\n\n";
	}
	cout << "}\n";

	return 0;
}
