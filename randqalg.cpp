#include "randqalg.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm>

#define _USE_MATH_DEFINES

RandQAlg::RandQAlg(int xSize, int ySize, int depth)
{
	cols = xSize;
	lines = ySize;
	givenDepth = depth;
	currentDepth = 0;
	layers.resize(givenDepth);
	for(auto &layer: layers)
		layer.resize(xSize * ySize, gateWait);
}

RandQAlg::~RandQAlg()
{
}

void RandQAlg::applyGate(QubitArray &qubits, Gate gate, cords qubit)
{
	switch(gate.name)
	{
	case gateName::TGate :
		qubits.applyRotation(qubit, {0, 0, 1}, M_PI_4);	//qubits.TGate(qubit); or PI/4 around Z-axis
		break;
	case gateName::sqrtX :
		qubits.applyRotation(qubit, {1, 0, 0}, M_PI_2); //qubits.sqrtX(qubit); or PI/2 around X-axis
		break;
	case gateName::sqrtY :
		qubits.applyRotation(qubit, {0, 1, 0}, M_PI_2); //qubits.sqrtY(qubit); or PI/2 around Y-axis
		break;
	case gateName::CZ :
		if(gate.isControl)
			qubits.cz(qubit, gate.companion);
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
		auto lastGateName = layers[i][index].name;
		auto iterator = std::find_if(singleGates.begin(), singleGates.end(),
			[lastGateName](Gate a)
			{
				return(a.name ==  lastGateName);
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
	if(begin > cols)
		return;

	int count;
	for(count = 0; count < begin; count++)
		addGate(genSingleGate({count, lineNum}), {count, lineNum});

	while(count + 1 < cols)
	{
		cords control = {count, lineNum}, target = {count + 1, lineNum};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBeetwCZ && count + i + 2 < cols; i++)
			addGate(genSingleGate({count + i + 2, lineNum}), {count + i + 2, lineNum});
		count += (gapBeetwCZ + 2);
	}
	if(count < cols)		//we couldn't apply last CZ because we reached the edge of array
		addGate(genSingleGate({count, lineNum}), {count, lineNum});
}

void RandQAlg::fillCol(int colNum, int begin)
{
	if(begin > lines)
		return;

	int count;
	for(count = 0; count < begin; count++)
		addGate(genSingleGate({colNum, count}), {colNum, count});

	while(count + 1 < lines)
	{
		cords control = {colNum, count}, target = {colNum, count + 1};
		addGate({gateName::CZ, true, target}, control);
		addGate({gateName::CZ, false, control}, target);
		for(int i = 0; i < gapBeetwCZ && count + i + 2 < lines; i++)
			addGate(genSingleGate({colNum, count + i + 2}), {colNum, count + i + 2});
		count += (gapBeetwCZ + 2);
	}
	if(count < lines)		//we couldn't apply last CZ because we reached the edge of array
		addGate(genSingleGate({colNum, count}), {colNum, count});
}

void RandQAlg::addLayer(const std::array<int, lineSamplesInMask> &beginPoints, int numOfLines, std::function<void (int, int)> fillerFunc)
{
	if(currentDepth >= givenDepth)
		return;

	for(int lineSample = 0; lineSample < lineSamplesInMask; lineSample++)
		for(int line = lineSample; line < numOfLines; line += lineSamplesInMask)
			fillerFunc(line, beginPoints[lineSample]);
	currentDepth++;
}

void RandQAlg::generate()
{
	currentDepth = 0;
	while(currentDepth < givenDepth)
	{
		addLayer({2, 0}, lines, lineFiller);
		addLayer({0, 2}, lines, lineFiller);
		addLayer({3, 1}, cols, colFiller);
		addLayer({1, 3}, cols, colFiller);
		addLayer({3, 1}, lines, lineFiller);
		addLayer({1, 3}, lines, lineFiller);
		addLayer({0, 2}, cols, colFiller);
		addLayer({2, 0}, cols, colFiller);
	}
}

void RandQAlg::init(QubitArray &qubits)
{
	qubits.reset();
	//H layer
	for(int i = 0; i < cols; i++)
		for(int j = 0; j < lines; j++)
			qubits.hadamardGate({i, j});
}

void RandQAlg::evaluate(QubitArray &qubits)
{
	if(qubits.getXSize() != cols || qubits.getYSize() != lines)
		throw std::length_error("Algorithm's layer size doesn't match qubit array size");

	init(qubits);
	for(int j = 0; j < givenDepth; j++)
		for(int i = 0; i <static_cast<int>(layers[j].size()); i++)
			applyGate(qubits, layers[j][i], {i % cols, i / cols});
}

void RandQAlg::evaluateLayer(QubitArray &qubits, int layerIndex)
{
	if(qubits.getXSize() != cols || qubits.getYSize() != lines)
		throw std::length_error("Algorithm's layer size doesn't match qubit array size");
	if(layerIndex >= givenDepth)
		throw std::out_of_range("layerNum must be less then algorithm depth");

	if(layerIndex == 0)
		init(qubits);
	for(int i = 0; i <static_cast<int>(layers[layerIndex].size()); i++)
		applyGate(qubits, layers[layerIndex][i], {i % cols, i / cols});
}


