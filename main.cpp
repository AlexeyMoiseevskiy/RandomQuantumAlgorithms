#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cmath>

#define LINES 3
#define COLS 4
#define DEPTH 40
#define AVG 100
#define ALGS 50

#define ERRDEPTH 24
#define ERRCORDS {3, 3}

using namespace std;

int main()
{
	freopen("LoseTime3.txt", "w", stdout);
	srand( (unsigned)time(NULL) );

	QubitArray qubits(COLS, LINES);
	qubits.setSingleErrRate(0.000);
	qubits.setMultiErrRate(0.00);
	qubits.setEnvCoupling(0);
	qubits.setLostProb(3);
	RandQAlg alg(COLS, LINES, DEPTH);
	alg.generate();

	static constexpr int ampNum = (1 << LINES * COLS);
	double amps[ampNum] = {0.0};

	cout << fixed << "{ ";
	for(int k = 0; k < ALGS; k++)
	{
		cerr << "Step " << k + 1 << " of " << ALGS << "\n";
		for(int i = 0; i < ampNum; i++)
			amps[i] = 0.0;
		alg.generate();
		for(int i = 0; i < AVG; i++)
		{
			alg.evaluate(qubits);
			for(int j = 0; j < ampNum; j++)
				amps[j] += qubits.getSquaredAmp(j);
		}

		for(int j = 0; j < ampNum; j++)
		{
			cout << amps[j] * ampNum / AVG;
			if(j + 1 < ampNum || k + 1 < ALGS)
				cout << ",\n";
		}
	}
	cout << "}";

	return 0;
}
