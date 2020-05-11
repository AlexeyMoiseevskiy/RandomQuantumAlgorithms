#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cmath>

#define LINES 4
#define COLS 4
#define DEPTH 40
#define AVR 50
#define ALGS 100

using namespace std;

int main()
{
	freopen("DistrAlgsErr16q.txt", "w", stdout);
	srand( (unsigned)time(NULL) );

	QubitArray noErrReg(COLS, LINES);
	noErrReg.setSingleErrRate(0.00);
	noErrReg.setMultiErrRate(0.0);
	noErrReg.setEnvCoupling(0.0);
	noErrReg.setSingleGateTime(0.01);
	noErrReg.setMultiGateTime(0.1);
	QubitArray errReg(COLS, LINES);
	errReg.setSingleErrRate(0.001);
	errReg.setMultiErrRate(0.01);
	errReg.setEnvCoupling(0.01);
	errReg.setSingleGateTime(0.01);
	errReg.setMultiGateTime(0.1);
	errReg.setLoseTime(1000);
	errReg.setDynamicNoise(0.04);
	errReg.setSpamErr0to1(0.01);
	errReg.setSpamErr1to0(0.03);

	RandQAlg alg(COLS, LINES, DEPTH);
	static constexpr int ampNum = (1 << LINES * COLS);
	double ampsErr[ampNum] = {0.0};
	double ampsNoErr[ampNum] = {0.0};

	cout << fixed << "{\n";
	for(int k = 0; k < ALGS; k++)
	{
		cerr << "Step " << k + 1 << " of " << ALGS << "\n";

		alg.generate();
		for(int i = 0; i < ampNum; i++)
			ampsErr[i] = ampsNoErr[i] = 0.0;
		for(int i = 0; i < AVR; i++)
		{
			alg.evaluate(noErrReg);
			alg.evaluate(errReg);
			for(int j = 0; j < ampNum; j++)
			{
				ampsErr[j] += errReg.getSquaredAmp(j);
				ampsNoErr[j] += noErrReg.getSquaredAmp(j);
			}
		}

		double xeb = 0;
		for(int i = 0; i < ampNum; i++)
		{
			xeb += (ampsNoErr[i] * ampsErr[i]) / (AVR * AVR);
		}
		cout << ampNum * xeb - 1;
		if(k + 1 < ALGS)
			cout << ", ";
	}
	cout << "\n}";

	return 0;
}
