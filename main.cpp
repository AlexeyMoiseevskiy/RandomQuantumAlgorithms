#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cmath>

#define LINES 4
#define COLS 4
#define DEPTH 100
#define AVG 1

#define ERRDEPTH 24
#define ERRCORDS {3, 3}

using namespace std;

int main()
{
	freopen("ErrAt(3,3).txt", "w", stdout);
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

	double probsErr[LINES][COLS][DEPTH] = {{{0.0}}};
	//double probsNoErr[LINES][COLS][DEPTH] = {{{0.0}}};

	alg.generate();

	for(int d = 0; d < DEPTH; d++)
	{
		//cerr << "Depth " << d << "\n";
		alg.evaluateLayer(regErr, d);
		alg.evaluateLayer(regNoErr, d);
		if(d == ERRDEPTH)
			regErr.XGate(ERRCORDS);
		for(int l = 0; l < LINES; l++)
			for(int c = 0; c < COLS; c++)
			{
				probsErr[l][c][d] = abs(regErr.calcProb({c, l}) - regNoErr.calcProb({c, l}));
			//	probsNoErr[l][c][d] = regNoErr.calcProb({c, l});
			}
	}

	cout << fixed << "{ ";

	for(int d = 0; d < DEPTH; d++)
	{
		cout << "{";
		for(int l = 0; l < LINES; l++)
			for(int c = 0; c < COLS; c++)
			{
				cout << "{" << l << ", " << c << ", " << probsErr[l][c][d] << "}";
				//cout << probsErr[l][c][d];
				if(l + 1 < LINES || c + 1 < COLS)
					cout << ", ";
			}
		cout << "}";
		if(d + 1 < DEPTH)
			cout << ", ";
	}
	cout << "}\n ";
	/*for(int d = 0; d < DEPTH; d++)
	{
		cout << "{";
		for(int l = 0; l < LINES; l++)
			for(int c = 0; c < COLS; c++)
			{
				cout << "{" << l << ", " << c << ", "<< probsErr[l][c][d] << "}";
				if(l + 1 < LINES || c + 1 < COLS)
					cout << ", ";
			}
		cout << "}";
		if(d + 1 < DEPTH)
			cout << ", ";
	}
	cout << "}\n} ";*/

	return 0;
}
