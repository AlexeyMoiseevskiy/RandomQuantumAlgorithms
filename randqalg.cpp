#include "randqalg.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iostream>

RandQAlg::RandQAlg(QuESTEnv externEnv, int a, int b, int depth)
	: QubitArray(externEnv, a, b)
{
	givenDepth = depth + 1; // + 1 init layer
}

RandQAlg::~RandQAlg()
{
}

void RandQAlg::applyGate(Gate gate, cords qubit)
{
	switch(gate.name)
	{
	case gateName::TGate :
		TGate(qubit);
		//std::cerr << "Tgate for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::sqrtX :
		sqrtX(qubit);
		//std::cerr << "sqrtX for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::sqrtY :
		sqrtY(qubit);
		//std::cerr << "sqrtY for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::CZ :
		if(gate.isControl)
		{
			cz(qubit, gate.companion);
			//std::cerr << "CZ for (" << qubit.x << ", " << qubit.y << ") and (" << gate.companion.x << ", " << gate.companion.y << ")\n";
		}
		break;
	case gateName::wait :
		break;
	case gateName::init :
		break;
	}
}

void RandQAlg::applySingleGate(gateName name, cords qubit)
{
	Gate gate;
	gate.name = name;
	applyGate(gate, qubit);
}


void RandQAlg::addGate(Gate gate, cords qubit)
{
	if(layers.back()[getIndex(qubit)].name != gateName::wait)
		throw std::runtime_error("incorrect program generated");
	layers.back()[getIndex(qubit)] = gate;
	/*switch(gate.name)
	{
	case gateName::TGate :
		std::cerr << "Tgate for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::sqrtX :
		std::cerr << "sqrtX for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::sqrtY :
		std::cerr << "sqrtY for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::CZ :
		if(gate.isControl)
		{
			std::cerr << "CZ for (" << qubit.x << ", " << qubit.y << ") and (" << gate.companion.x << ", " << gate.companion.y << ")\n";
		}
		break;
	case gateName::wait :
		std::cerr << "wait for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	case gateName::init :
		std::cerr << "init lolwut for (" << qubit.x << ", " << qubit.y << ")\n";
		break;
	}*/
}

void RandQAlg::addSingleGate(gateName name, cords qubit)
{
	Gate gate;
	gate.name = name;
	addGate(gate, qubit);
}

gateName RandQAlg::randomSingleGate()
{
	return singleGates[static_cast<int>(rand() % singleGates.size())];
}

gateName RandQAlg::lastSingleGateApplyed(cords qubit)
{
	//std::cerr << "look for last single gate, layers.size = " << layers.size() << "\n";
	int index = getIndex(qubit);
	for(int i = layers.size() - 1; i >= 0; i--)
	{
		//std::cerr << "At depth " << i << " the gate is " << static_cast<int>(layers[i][index].name) << "\n";
		auto iterator = std::find(singleGates.begin(), singleGates.end(), layers[i][index].name);
		if(iterator != singleGates.end())
			return *iterator;
	}
	return gateName::init;
}

gateName RandQAlg::lastGateApplyed(cords qubit)
{
	int index = getIndex(qubit);
	for(int i = layers.size() - 1; i >= 0; i--)
		if(layers[i][index].name != gateName::wait)
			return layers[i][index].name;

	return gateName::init;
}

gateName RandQAlg::genSingleGate(cords target)
{
	gateName lastSingleGate = lastSingleGateApplyed(target);
	if(lastSingleGate == gateName::init)
	{
		return gateName::TGate;
	}

	gateName lastGate = lastGateApplyed(target);
	if(lastGate != gateName::CZ)
	{
		return gateName::wait;
	}

	gateName randGate = randomSingleGate();
	while (randGate == lastSingleGate)
		randGate = randomSingleGate();
	return randGate;
}

void RandQAlg::fillLine(int lineNum, int begin)
{
//	std::cerr << "Start filline\n";
	int count, end = getXSize();
	for(count = 0; count < begin; count++)
		addSingleGate(genSingleGate({count, lineNum}), {count, lineNum});

	//std::cerr << "filline main cycle, count + 1 = " << count + 1 << " end  = " << end << "\n";
	while(count + 1 < end)
	{
		cords control = {count, lineNum}, target = {count + 1, lineNum};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBtwCZ && count + i + 2 < end; i++)
			addSingleGate(genSingleGate({count + i + 2, lineNum}), {count + i + 2, lineNum});
		count += (gapBtwCZ + 2);
	}
	if(count < end)		//we couldn't apply last CZ because we reached the edge of array
	{
		//std::cerr << "speccial case\n";
		addSingleGate(genSingleGate({count, lineNum}), {count, lineNum});
	}
}

void RandQAlg::fillCol(int colNum, int begin)
{
	int count, end = getYSize();
	for(count = 0; count < begin; count++)
		addSingleGate(genSingleGate({colNum, count}), {colNum, count});

	while(count + 1 < end)
	{
		cords control = {colNum, count}, target = {colNum, count + 1};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBtwCZ && count + i + 2 < end; i++)
				addSingleGate(genSingleGate({colNum, count + i + 2}), {colNum, count + i + 2});
		count += (gapBtwCZ + 2);
	}
	if(count < end)		//we couldn't apply last CZ because we reached the edge of array
	{
		//std::cerr << "speccial case\n";
		addSingleGate(genSingleGate({colNum, count}), {colNum, count});
	}
}

bool RandQAlg::newLayer()
{
	if(static_cast<int>(layers.size()) == givenDepth)
		return false;
	Gate wait;
	wait.name = gateName::wait;
	layers.push_back(std::vector<Gate>(getXSize() * getYSize(), wait));
	return true;
}

void RandQAlg::addLayersBlock()
{
	int lines = getYSize(), cols = getXSize();
	if(!newLayer()) return;
	//std::cerr << "Mask 1\n";
	for(int i = 0; i < lines; i += lineSaplesInMask)	//#1
		fillLine(i, 2);
	for(int i = 1; i < lines; i += lineSaplesInMask)
		fillLine(i, 0);
	if(!newLayer()) return;
	//std::cerr << "Mask 2\n";
	for(int i = 0; i < lines; i += lineSaplesInMask)	//#2
		fillLine(i, 0);
	for(int i = 1; i < lines; i += lineSaplesInMask)
		fillLine(i, 2);
	if(!newLayer()) return;
	//std::cerr << "Mask 3\n";
	for(int i = 0; i < cols; i += lineSaplesInMask)		//#3
		fillCol(i, 3);
	for(int i = 1; i < cols; i += lineSaplesInMask)
		fillCol(i, 1);
	if(!newLayer()) return;
	//std::cerr << "Mask 4\n";
	for(int i = 0; i < cols; i += lineSaplesInMask)		//#4
		fillCol(i, 1);
	for(int i = 1; i < cols; i += lineSaplesInMask)
		fillCol(i, 3);
	if(!newLayer()) return;
	//std::cerr << "Mask 5\n";
	for(int i = 0; i < lines; i += lineSaplesInMask)	//#5
		fillLine(i, 3);
	for(int i = 1; i < lines; i += lineSaplesInMask)
		fillLine(i, 1);
	if(!newLayer()) return;
	//std::cerr << "Mask 6\n";
	for(int i = 0; i < lines; i += lineSaplesInMask)	//#6
		fillLine(i, 1);
	for(int i = 1; i < lines; i += lineSaplesInMask)
		fillLine(i, 3);
	if(!newLayer()) return;
	//std::cerr << "Mask 7\n";
	for(int i = 0; i < cols; i += lineSaplesInMask)		//#7
		fillCol(i, 0);
	for(int i = 1; i < cols; i += lineSaplesInMask)
		fillCol(i, 2);
	if(!newLayer()) return;
	//std::cerr << "Mask 8\n";
	for(int i = 0; i < cols; i += lineSaplesInMask)		//#8
		fillCol(i, 2);
	for(int i = 1; i < cols; i += lineSaplesInMask)
		fillCol(i, 0);
}

void RandQAlg::init()
{
	for(auto &layer: layers)
		layer.clear();
	layers.clear();

	Gate initGate;
	initGate.name = gateName::init;
	layers.push_back(std::vector<Gate>(getXSize() * getYSize(), initGate));
}

void RandQAlg::generate()
{
	init();
	while(static_cast<int>(layers.size()) < givenDepth)
		addLayersBlock();
}

void RandQAlg::evaluate()
{
	QubitArray::init();
	//H layer
	for(int i = 0; i < getXSize(); i++)
		for(int j = 0; j < getYSize(); j++)
			hadamardGate({i, j});

	for(auto &layer: layers)
		for(int i = 0; i < static_cast<int>(layer.size()); i++)
			applyGate(layer[i], {i % getXSize(), i / getXSize()});
}




