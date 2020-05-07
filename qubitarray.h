#ifndef QUBITARRAY_H
#define QUBITARRAY_H

#include <QuEST.h>
#include <random>
#include <vector>
#include <functional>

#define _USE_MATH_DEFINES

struct cords
{
	int x;
	int y;
};

class QubitArray
{
public:
	QubitArray(int a, int b);
	~QubitArray(){}

	void setSize(int a, int b);
	void setEnvCoupling(int val){ envCoupling = val; }
	void setSingleGateCoupling(double val){ singleGateCoupling = val; }
	void setMultiGateCoupling(double val){ multiGateCoupling = val; }
	void setSingleErrRate(double val){ singleErrRate = val / 2; }
	void setMultiErrRate(double val){ multiErrRate = val / 2; }
	void setSingleGateTime(double val){ singleGateTime = val; }
	void setMultiGateTime(double val){ multiGateTime = val; }
	void setLoseTime(double val){ loseTime = val; }
	void setDynamicNoise(double val){ dynamicNoise = val; }

	int getXSize(){ return xSize; }
	int getYSize(){ return ySize; }

	void init();
	void generateBell(cords first, cords sec);
	void cz(cords control, cords target);
	void swap(cords first, cords sec);
	int move(cords init, cords dest);

	void applySingleGate(cords target, std::function<void(Qureg, int)> gate);
	void applyRotation(cords target, Vector v, double angle);
	void hadamardGate(cords target){ applySingleGate(target, &hadamard); }
	void sqrtX(cords target){ applySingleGate(target, [](Qureg reg, int index){rotateX(reg, index, M_PI_2);}); }
	void sqrtY(cords target){ applySingleGate(target, [](Qureg reg, int index){rotateY(reg, index, M_PI_2);}); }
	void TGate(cords target){ applySingleGate(target, tGate); }
	void XGate(cords target){ applySingleGate(target, pauliX); }

	double calcBellFidelity(cords first, cords sec);
	double calcBellFidelityDirect(cords first, cords sec);
	int meas(cords target){ return measure(qubits, getIndex(target)); }
	double getSquaredAmp(int index);
	friend double fidelity(const QubitArray &arr1, const QubitArray &arr2){ return calcFidelity(arr1.qubits, arr2.qubits); }
	double calcProb(cords target){ return calcProbOfOutcome(qubits, getIndex(target), 1); }

	constexpr int getIndex(cords c){ return c.y * xSize + c.x; }

private:
	QuESTEnv env;
	Qureg qubits;

	int xSize;
	int ySize;
	double envCoupling;
	double singleGateCoupling;
	double multiGateCoupling;
	double singleErrRate;
	double multiErrRate;
	double singleGateTime;
	double multiGateTime;
	double loseTime;
	double dynamicNoise;

	void updateTotalTime();
	std::random_device rd{};
	std::mt19937 gen{rd()};

	void applyNoiseGate(int index, double coupling, double time);
	void applyNoise(int index);
	void applyNoise(cords target){ applyNoise(getIndex(target)); }
	void applyNoise();
	void applySingleGateErr(cords target);
	void applyMultiGateErr(cords target);
	void applyMultiGateErr(cords first, cords sec){ applyMultiGateErr(first); applyMultiGateErr(sec); }

	double totalTime;
	bool singleGateInCurLayer;
	bool multiGateInCurLayer;
	std::vector <int> lastNoiseTime;
	std::vector <bool> usedInCurLayer;
	std::vector <bool> isLost;
	void startNewLayer();
};

#endif // QUBITARRAY_H
