#ifndef RANDQALG_H
#define RANDQALG_H

#include <QuEST.h>
#include "qubitarray.h"
#include <vector>

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
	gateName name;
	bool isControl = false;				//for multiqubit only
	QubitArray::cords companion;		//for multiqubit only
};

class RandQAlg : public QubitArray
{
public:
	RandQAlg(QuESTEnv externEnv, int a, int b, int depth);
	~RandQAlg();

	const std::vector<gateName> singleGates{gateName::TGate, gateName::sqrtX, gateName::sqrtY};

	gateName randomSingleGate();
	gateName genSingleGate(cords target);
	void init();
	void generate();
	void evaluate();

private:
	int givenDepth;
	std::vector<std::vector<Gate>> layers;

	void applyGate(Gate gate, cords qubit);
	void applySingleGate(gateName name, cords qubit);
	void addGate(Gate gate, cords qubit);
	void addSingleGate(gateName name, cords qubit);
	gateName lastGateApplyed(cords qubit);
	gateName lastSingleGateApplyed(cords qubit);

	static constexpr int gapBtwCZ = 2;
	static constexpr int lineSaplesInMask = 2;
	void fillLine(int lineNum, int begin);
	void fillCol(int colNum, int begin);
	void addLayersBlock();
	bool newLayer();
};

#endif // RANDQALG_H
