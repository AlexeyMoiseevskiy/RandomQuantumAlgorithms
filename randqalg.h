#ifndef RANDQALG_H
#define RANDQALG_H

#include <QuEST.h>
#include "qubitarray.h"
#include <array>

enum class gateName
{
	TGate,
	sqrtX,
	sqrtY,
	CZ,
	wait,
	init,
};

struct Gate
{
	constexpr Gate(gateName str, bool control = false, QubitArray::cords companionCords = {0, 0})
		: name(str), isControl(control), companion(companionCords){}

	gateName name;
	bool isControl;						//for multiqubit only
	QubitArray::cords companion;		//for multiqubit only
};

class RandQAlg : public QubitArray
{
public:
	RandQAlg(QuESTEnv externEnv, int a, int b, int depth);
	~RandQAlg();

	static constexpr std::array<Gate, 3> singleGates{Gate(gateName::TGate), Gate(gateName::sqrtX), Gate(gateName::sqrtY)};
	static constexpr auto gateInit = Gate(gateName::init);
	static constexpr auto gateWait = Gate(gateName::wait);

	Gate genSingleGate(cords target);
	void init();
	void generate();
	void evaluate();

private:
	int givenDepth;
	std::vector<std::vector<Gate>> layers;

	void applyGate(Gate gate, cords qubit);
	void addGate(Gate gate, cords qubit);
	Gate lastGateApplyed(cords qubit);
	Gate lastSingleGateApplyed(cords qubit);

	static constexpr int gapBeetwCZ = 2;
	static constexpr int lineSaplesInMask = 2;

	void fillLine(int lineNum, int begin);
	void fillCol(int colNum, int begin);
	void addLayer(const int beginPoints[], int numOfLines, void (RandQAlg::*fillerFunc)(int, int));

	/*struct Mask
	{
		int beginPoints[lineSaplesInMask];
		int endPoint;
		void (RandQAlg::*fillerFunc)(int, int);
	};
	static constexpr std::array<Mask, 1> masks{{0, 2}, 2, &(RandQAlg::fillLine)};*/
};

#endif // RANDQALG_H
