#include "qubitarray.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

QubitArray::QubitArray(QuESTEnv externEnv, int a, int b)
{
	env = externEnv;
	setSize(a, b);

	envCoupling = 0.2;
	singleGateCoupling = 0.2;
	multiGateCoupling = 0.2;

	setSingleErrRate(0.01);
	setMultiErrRate(0.1);

	singleGateTime = 1;
	multiGateTime = 2;
}

void QubitArray::setSize(int a, int b)
{
	xSize = a;
	ySize = b;
	init();
}

int QubitArray::getIndex(cords c)
{
	if (c.x > xSize || c.y > ySize || c.y < 0 || c.x < 0)
		throw std::out_of_range("Given index is out of qubit register ranges");
	return c.y * xSize + c.x;
}

void QubitArray::updateTotalTime()
{
	if(singleGateInCurLayer && !multiGateInCurLayer)
		totalTime += singleGateTime;
	if(multiGateInCurLayer)
		totalTime += multiGateTime;
}

void QubitArray::init()
{
	qubits = createQureg(xSize * ySize, env);
	lastNoiseTime.resize(xSize * ySize);
	usedInCurLayer.resize(xSize * ySize);
	totalTime = 0;
	std::fill(lastNoiseTime.begin(), lastNoiseTime.end(), 0);
	singleGateInCurLayer = false;
	multiGateInCurLayer = false;
	startNewLayer();
}

void QubitArray::startNewLayer()
{
	updateTotalTime();
	//std::cerr << "New Layer, time is " << totalTime << "\n";
	singleGateInCurLayer = false;
	multiGateInCurLayer = false;
	std::fill(usedInCurLayer.begin(), usedInCurLayer.end(), false);
}

void QubitArray::generateBell(cords first, cords sec)
{
	hadamard(qubits, getIndex(first));
	controlledNot(qubits, getIndex(first), getIndex(sec));
	multiGateInCurLayer = true;
	usedInCurLayer[getIndex(first)] = usedInCurLayer[getIndex(sec)] = true;
	lastNoiseTime[getIndex(first)] = lastNoiseTime[getIndex(sec)] = totalTime + multiGateTime;
	//std::cerr << "Generated bell, time is " << totalTime << "\n";
}

void QubitArray::cz(cords control, cords target)
{
	if(usedInCurLayer[getIndex(control)] || usedInCurLayer[getIndex(target)])
		startNewLayer();
	usedInCurLayer[getIndex(control)] = usedInCurLayer[getIndex(target)] = true;

	applyNoise(control);
	applyNoise(target);
	applyMultiGateErr(control, target);

	hadamard(qubits, getIndex(target));
	controlledNot(qubits, getIndex(control), getIndex(target));
	hadamard(qubits, getIndex(target));

	multiGateInCurLayer = true;

	applyMultiGateErr(control, target);
}

void QubitArray::swap(cords first, cords sec)
{
	//std::cerr << "Moving " << getIndex(first) << " to " << getIndex(sec)
	//			 << " time is " << totalTime << "\n";
	if(usedInCurLayer[getIndex(first)] || usedInCurLayer[getIndex(sec)])
		startNewLayer();
	usedInCurLayer[getIndex(first)] = usedInCurLayer[getIndex(sec)] = true;

	applyNoise(first);
	applyNoise(sec);
	applyMultiGateErr(first, sec);

	swapGate(qubits, getIndex(first), getIndex(sec));
	multiGateInCurLayer = true;

	applyMultiGateErr(first, sec);
}

int QubitArray::move(cords init, cords dest)
{
	if(dest.x < 0 || dest.y < 0 || dest.x > xSize || dest.y > ySize || init.x < 0 || init.y < 0 || init.x > xSize || init.y > ySize)
		throw std::out_of_range("Given index is out of qubit register ranges");
	cords cur = init;
	int xStep = 1, yStep = 1;

	if (dest.x < init.x)
		xStep = -1;
	if(dest.y < init.y)
		yStep = -1;

	while (cur.y != dest.y)
	{
		swap(cur, {cur.x, cur.y + yStep});
		cur.y += yStep;
	}
	while (cur.x != dest.x)
	{
		swap(cur, {cur.x + xStep, cur.y});
		cur.x += xStep;
	}
	return abs(init.x - dest.x) + abs(init.y - dest.y);
}

void QubitArray::hadamardGate(cords target)
{
	if(usedInCurLayer[getIndex(target)])
		startNewLayer();
	usedInCurLayer[getIndex(target)] = true;

	applyNoise(target);
	applySingleGateErr(target);
	hadamard(qubits, getIndex(target));
	singleGateInCurLayer = true;

	applySingleGateErr(target);
}

double QubitArray::calcBellFidelity(cords first, cords sec)
{
	updateTotalTime();
	//std::cerr << "Calcing fidelity, time is " << totalTime << "\n";
	applyNoise();
	controlledNot(qubits, getIndex(first), getIndex(sec));
	hadamard(qubits, getIndex(first));
	Complex amp = getAmp(qubits, 0);
	//hadamard(qubits, getIndex(first));
	//controlledNot(qubits, getIndex(first), getIndex(sec));
	return amp.real * amp.real + amp.imag * amp.imag;
}

double QubitArray::calcBellFidelityDirect(cords first, cords sec)
{
	//std::cerr << "Calcing fidelity\n";
	updateTotalTime();
	applyNoise();

	//preparing ideal register
	Qureg noNoiseReg = createQureg(xSize * ySize, env);
	hadamard(noNoiseReg, 0);
	controlledNot(noNoiseReg, 0, 1);
	if(getIndex(first) != 0)
		swapGate(noNoiseReg, 0, getIndex(first));
	if(getIndex(sec) != 1)
		swapGate(noNoiseReg, 1, getIndex(sec));

	return calcFidelity(qubits, noNoiseReg);
}

void QubitArray::applyNoiseGate(int index, double coupling, double time)
{
	std::normal_distribution<double> rand(0, 1 * sqrt(time));	// \sigma = 1
	double w = rand(gen);
	rotateZ(qubits, index, 2 * sqrt(coupling) * w);
}

void QubitArray::applyNoise(int index)
{
	//std::cerr << "Noise to qubit " << index << " time is " << totalTime - lastNoiseTime[index] << '\n';
	applyNoiseGate(index, envCoupling, totalTime - lastNoiseTime[index]);
	lastNoiseTime[index] = totalTime;
}

void QubitArray::applyNoise()
{
	for(int i = 0; i < xSize * ySize; i++)
		applyNoise(i);
}

void QubitArray::applySingleGateErr(cords target)
{
	//std::cerr << "Single Gate error to qubit " << getIndex(target) << '\n';
	double effectiveTime = -(log(1 - 2 * singleErrRate) / (2 * singleGateCoupling));
	applyNoiseGate(getIndex(target), singleGateCoupling, effectiveTime);
	lastNoiseTime[getIndex(target)] = totalTime + singleGateTime;
}


void QubitArray::applyMultiGateErr(cords target)
{
	//std::cerr << "Multi Gate error to qubit " << getIndex(target) << '\n';
	double effectiveTime = -(log(1 - 2 * multiErrRate) / (2 * multiGateCoupling));
	applyNoiseGate(getIndex(target), multiGateCoupling, effectiveTime);
	lastNoiseTime[getIndex(target)] = totalTime + multiGateTime;
}
