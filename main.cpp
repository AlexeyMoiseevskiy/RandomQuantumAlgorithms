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
#define STATISTICS 1000

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
	if(!usedBefore[x][y])
	{
		qubits.TGate({x, y});
		lastGate[x][y] = gates::TGate;
		CZLast[x][y] = false;
		usedBefore[x][y] = true;
		return;
	}
	if(!CZLast[x][y])
		return;

	int r = rand() % GATESNUM;
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
	int i;
	for(i = 0; i < begin; i++)
		applyRandomGate(qubits, i, line);
	for(i = begin; i + 1 < end; i += gap + 1)
	{
		qubits.cz({i, line}, {i + 1, line});
		//cout << "CZ for " << i << " " << line << "and" << i + 1 << " " << line << '\n';
		CZLast[i][line] = true;
		CZLast[i + 1][line] = true;
		for(int j = 0; j < gap && i + 2 + j < end; j++)
			applyRandomGate(qubits, i + 2 + j, line);
	}
	for( ; i < qubits.getXSize(); i++)
		applyRandomGate(qubits, i, line);
}


void fillCol(QubitArray &qubits, int col, int begin, int end, int gap)
{
	int i;
	for(i = 0; i < begin; i++)
		applyRandomGate(qubits, col, i);
	for(i = begin; i + 1< end; i += gap + 1)
	{
		qubits.cz({col, i}, {col, i + 1});
		//cout << "CZ for " << col << " " << i << "and" << col << " " << i + 1 << '\n';
		CZLast[col][i] = true;
		CZLast[col][i + 1] = true;
		for(int j = 0; j < gap && i + 2 + j < end; j++)
			applyRandomGate(qubits, col, i + 2 + j);
	}
	for(; i < qubits.getYSize(); i++)
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
		fillCol(qubits, i, 3, lines, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 1, lines, gap);
	//#4
	for(int i = 0; i < cols; i += 2)
		fillCol(qubits, i, 1, lines, gap);
	for(int i = 1; i < cols; i += 2)
		fillCol(qubits, i, 3, lines, gap);
	//#5
	for(int i = 0; i < lines; i += 2)
		fillLine(qubits, i, 3, cols, gap);
	for(int i = 1; i < lines; i += 2)
		fillLine(qubits, i, 1, cols, gap);
	//#6
	for(int i = 0; i < lines; i += 2)
		fillLine(qubits, i, 1, cols, gap);
	for(int i = 1; i < lines; i += 2)
		fillLine(qubits, i, 3, cols, gap);
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

	freopen("Hist3by4depth7stat1k.txt", "w", stdout);
	srand( (unsigned)time(NULL) );
	QuESTEnv env = createQuESTEnv();
	QubitArray qubits(env, COLS, LINES);

	qubits.setSingleErrRate(0);
	qubits.setMultiErrRate(0);

	//H layer
	for(int i = 0; i < COLS; i++)
		for(int j = 0; j < LINES; j++)
			qubits.hadamardGate({i, j});

	cout << fixed << "{";
	for(int j = 0; j < STATISTICS; j++)
	{
		qubits.init();
		cerr << "Step " << j << " of "<< STATISTICS << "\n";
		for(int i = 0; i < DEPTH; i++)
			applyLayersBlock(qubits);

		for(int i = 0; i < (1 << (COLS * LINES)); i++)
			cout << qubits.getSquaredAmp(i) << ", ";
	}
	cout << "}\n";
	return 0;
}
