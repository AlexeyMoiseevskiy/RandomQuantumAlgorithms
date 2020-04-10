#include "randqalg.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm>

RandQAlg::RandQAlg(QuESTEnv externEnv, int a, int b, int depth)
	: QubitArray(externEnv, a, b)
{
	givenDepth = depth;
	layers.resize(givenDepth);
	for(auto &layer: layers)
		layer.resize(a * b, gateWait);
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
		break;
	case gateName::sqrtX :
		sqrtX(qubit);
		break;
	case gateName::sqrtY :
		sqrtY(qubit);
		break;
	case gateName::CZ :
		if(gate.isControl)
			cz(qubit, gate.companion);
		break;
	case gateName::wait :
		break;
	case gateName::init :
		break;
	}
}

void RandQAlg::addGate(Gate gate, cords qubit)
{
	layers[currentDepth][getIndex(qubit)] = gate;
}

Gate RandQAlg::lastSingleGateApplyed(cords qubit)
{
	int index = getIndex(qubit);
	for(int i = currentDepth - 1; i >= 0; i--)
	{
		auto wantedName = layers[i][index].name;
		auto iterator = std::find_if(singleGates.begin(), singleGates.end(),
			[wantedName](Gate a)
			{
				return(a.name ==  wantedName);
			}
		);
		if(iterator != singleGates.end())
			return *iterator;
	}
	return gateInit;
}

Gate RandQAlg::lastGateApplyed(cords qubit)
{
	int index = getIndex(qubit);
	for(int i = currentDepth - 1; i >= 0; i--)
		if(layers[i][index].name != gateName::wait)
			return layers[i][index].name;
	return gateInit;
}

Gate RandQAlg::genSingleGate(cords target)
{
	const Gate lastSingleGate = lastSingleGateApplyed(target);
	if(lastSingleGate.name == gateName::init)
		return Gate(gateName::TGate);

	const Gate lastGate = lastGateApplyed(target);
	if(lastGate.name != gateName::CZ)
		return Gate(gateName::wait);

	Gate randGate = singleGates[static_cast<int>(rand() % singleGates.size())];
	while (randGate.name == lastSingleGate.name)
		randGate = singleGates[static_cast<int>(rand() % singleGates.size())];
	return randGate;
}

void RandQAlg::fillLine(int lineNum, int begin)
{
	int count, end = getXSize();
	for(count = 0; count < begin; count++)
		addGate(genSingleGate({count, lineNum}), {count, lineNum});

	while(count + 1 < end)
	{
		cords control = {count, lineNum}, target = {count + 1, lineNum};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBeetwCZ && count + i + 2 < end; i++)
			addGate(genSingleGate({count + i + 2, lineNum}), {count + i + 2, lineNum});
		count += (gapBeetwCZ + 2);
	}
	if(count < end)		//we couldn't apply last CZ because we reached the edge of array
		addGate(genSingleGate({count, lineNum}), {count, lineNum});
}

void RandQAlg::fillCol(int colNum, int begin)
{
	int count, end = getYSize();
	for(count = 0; count < begin; count++)
		addGate(genSingleGate({colNum, count}), {colNum, count});

	while(count + 1 < end)
	{
		cords control = {colNum, count}, target = {colNum, count + 1};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBeetwCZ && count + i + 2 < end; i++)
			addGate(genSingleGate({colNum, count + i + 2}), {colNum, count + i + 2});
		count += (gapBeetwCZ + 2);
	}
	if(count < end)		//we couldn't apply last CZ because we reached the edge of array
		addGate(genSingleGate({colNum, count}), {colNum, count});
}

void RandQAlg::addLayer(const std::array<int, lineSamplesInMask> &beginPoints, int numOfLines, void (RandQAlg::*fillerFunc)(int, int))
{
	if(currentDepth >= givenDepth)
		return;

	for(int lineSample = 0; lineSample < lineSamplesInMask; lineSample++)
		for(int line = lineSample; line < numOfLines; line += lineSamplesInMask)
			(this->*fillerFunc)(line, beginPoints[lineSample]);
	currentDepth++;
}

void RandQAlg::generate()
{
	currentDepth = 0;
	int lines = getYSize(), cols = getXSize();
	while(currentDepth < givenDepth)
	{
		addLayer({2, 0}, lines, &RandQAlg::fillLine);
		addLayer({0, 2}, lines, &RandQAlg::fillLine);
		addLayer({3, 1}, cols, &RandQAlg::fillCol);
		addLayer({1, 3}, cols, &RandQAlg::fillCol);
		addLayer({3, 1}, lines, &RandQAlg::fillLine);
		addLayer({1, 3}, lines, &RandQAlg::fillLine);
		addLayer({0, 2}, cols, &RandQAlg::fillCol);
		addLayer({2, 0}, cols, &RandQAlg::fillCol);
	}
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




