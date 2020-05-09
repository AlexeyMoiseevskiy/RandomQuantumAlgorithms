#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cmath>

#define LINES 1
#define COLS 2
#define DEPTH 40
#define AVG 100000
#define ALGS 50

using namespace std;

void czTomo(QubitArray &qubits, bool first, bool sec)
{
	double case00 = 0, case01 = 0, case10 = 0, case11 = 0;
	cords ctrl = {0, 0}, targ = {1, 0};
	for(int i = 0; i < AVG; i++)
	{
		qubits.init();
		if(first)
			qubits.pureX(ctrl);
		if(sec)
			qubits.pureX(targ);

		qubits.hadamardGate(targ);
		qubits.cz(ctrl, targ);
		qubits.hadamardGate(targ);

		/*int qub1 = qubits.meas(ctrl);
		int qub2 = qubits.meas(targ);

		if(qub1 == 0 && qub2 == 0)
			case00++;
		if(qub1 == 0 && qub2 == 1)
			case01++;
		if(qub1 == 1 && qub2 == 0)
			case10++;
		if(qub1 == 1 && qub2 == 1)
			case11++;
		*/
		case00 += qubits.getSquaredAmp(0);
		case10 += qubits.getSquaredAmp(1);
		case01 += qubits.getSquaredAmp(2);
		case11 += qubits.getSquaredAmp(3);
	}
	cout << "{"<< case00 * 100 / AVG << ", " << case01 * 100 / AVG << ", "
		 << case10 * 100 / AVG << ", " << case11 * 100 / AVG << "},\n";
}

int main()
{
	//freopen("CZLogs.txt", "w", stdout);
	srand( (unsigned)time(NULL) );

	QubitArray qubits(COLS, LINES);
	qubits.setSingleErrRate(0.001);
	qubits.setMultiErrRate(0.01);
	qubits.setEnvCoupling(0.01);
	qubits.setSingleGateTime(0.01);
	qubits.setMultiGateTime(0.1);
	qubits.setLoseTime(1000);
	qubits.setDynamicNoise(0.04);
	qubits.setSpamErr0to1(0.01);
	qubits.setSpamErr1to0(0.03);
	cout << "{\n";
	czTomo(qubits, 0, 0);
	czTomo(qubits, 0, 1);
	czTomo(qubits, 1, 0);
	czTomo(qubits, 1, 1);

	cout << "}\n";

	return 0;
}
