#include <iostream>
#include <QuEST.h>
#include "qubitarray.h"
#include <cmath>
#include <string>

#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#define LINES 3
#define COLS 4
#define GATESNUM 3
#define DEPTH 7
#define STATISTICS 10000

using namespace std;

enum gates
{
	TGate,
	sqrtX,
	sqrtY,
};

gates lastGate[COLS][LINES];
bool CZLast[COLS][LINES];
bool usedBefore[COLS][LINES];

void applyRandomGate(QubitArray &qubits, int x, int y)
{
	if(usedBefore[x][y])
	{
		qubits.TGate({x, y});
		usedBefore[x][y] = true;
		return;
	}
	if(!CZLast[x][y])
		return;

	int r;
	for (r = rand() % GATESNUM; r == lastGate[x][y]; r = rand() % GATESNUM);

	switch(r)
	{
	case gates::TGate :
		qubits.TGate({x, y});
		break;
	case gates::sqrtX :
		qubits.sqrtX({x, y});
		break;
	case gates::sqrtY :
		qubits.sqrtY({x, y});
		break;
	default :
		cout << "applyRandomGate : Unexpected error with random generator\n";
	}
	lastGate[x][y] = static_cast<gates>(r);
	CZLast[x][y] = false;
}

void fillLine(QubitArray &qubits, int line, int begin, int end, int gap)
{
	for(int i = 0; i < begin; i++)
		applyRandomGate(qubits, i, line);
	for(int i = begin; i + 1< end; i += gap + 1)
	{
		qubits.cz({i, line}, {i + 1, line});
		CZLast[i][line] = true;
		CZLast[i + 1][line] = true;
		for(int j = 0; j < gap && i + 2 + j < end; j++)
			applyRandomGate(qubits, i + 2 + j, line);
	}
	for(int i = end; i < qubits.getXSize(); i++)
		applyRandomGate(qubits, i, line);
}


void fillCol(QubitArray &qubits, int col, int begin, int end, int gap)
{
	for(int i = 0; i < begin; i++)
		applyRandomGate(qubits, col, i);
	for(int i = begin; i + 1< end; i += gap + 1)
	{
		qubits.cz({col, i}, {col, i + 1});
		CZLast[col][i] = true;
		CZLast[col][i + 1] = true;
		for(int j = 0; j < gap && i + 2 + j < end; j++)
			applyRandomGate(qubits, col, i + 2 + j);
	}
	for(int i = end; i < qubits.getYSize(); i++)
		applyRandomGate(qubits, col, i);
}

void applyLayersBlock(QubitArray &qubits)
{
	const int gap = 2;
	const int lines = qubits.getYSize(), cols = qubits.getXSize();
	//#1
	for(int i = 0; i < lines; i += 2)
		fillLine(qubits, i, 2, cols, gap);
	for(int i = 1; i < lines; i += 2)
		fillLine(qubits, i, 0, cols, gap);
	//#2
	for(int i = 0; i < lines; i += 2)
		fillLine(qubits, i, 0, cols, gap);
	for(int i = 1; i < lines; i += 2)
		fillLine(qubits, i, 2, cols, gap);
	//#3
	for(int i = 0; i < cols; i += 2)
		fillCol(qubits, i, 3, lines - 1, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 1, lines - 1, gap);
	//#4
	for(int i = 0; i < cols; i += 2)
		fillCol(qubits, i, 1, lines - 1, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 3, lines - 1, gap);
	//#5
	for(int i = 0; i < lines; i += 2)
		fillCol(qubits, i, 3, cols - 1, gap);
	for(int i = 1; i < lines; i += 2)
		fillCol(qubits, i, 1, cols - 1, gap);
	//#6
	for(int i = 0; i < lines; i += 2)
		fillCol(qubits, i, 1, cols - 1, gap);
	for(int i = 1; i < lines; i += 2)
		fillCol(qubits, i, 3, cols - 1, gap);
	//#7
	for(int i = 0; i < cols; i += 2)
		fillCol(qubits, i, 0, lines, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 2, lines, gap);
	//#8
	for(int i = 0; i < cols; i += 2)
		fillCol(qubits, i, 2, lines, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 0, lines, gap);
}

int main()
{

	freopen("3by4depth7statistics10k.txt", "w", stdout);
	srand( (unsigned)time(NULL) );
	QuESTEnv env = createQuESTEnv();
	QubitArray qubits(env, COLS, LINES);

	//H layer
	for(int i = 0; i < COLS; i++)
		for(int j = 0; j < LINES; j++)
			qubits.hadamardGate({i, j});

	qubits.setSingleErrRate(0);
	qubits.setMultiErrRate(0);

	for(int j = 0; j < STATISTICS; j++)
	{
		cerr << "Step " << j << " of "<< STATISTICS << "\n";
		for(int i = 0; i < DEPTH; i++)
			applyLayersBlock(qubits);
		//int outcome[COLS][LINES];

		cout << "{ ";
		for(int i = 0; i < COLS; i++)
		{
			for(int j = 0; j < LINES; j++)
				cout << qubits.meas({i, j}) << ", ";
		}
		cout << "},\n";
	}
	return 0;
}
