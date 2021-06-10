#ifndef RANDQALG_H
#define RANDQALG_H

#include <QuEST.h>
#include "qubitarray.h"
#include <array>
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
	constexpr Gate(gateName newName, bool control = false, cords companionCords = {0, 0})
		: name(newName), isControl(control), companion(companionCords){}

	gateName name;
	bool isControl;				//for multiqubit only
	cords companion;			//for multiqubit only
};

class RandQAlg
{
public:
	RandQAlg(int xSize, int ySize, int depth);
	~RandQAlg();

	constexpr int getIndex(cords c){ return c.y * cols + c.x; }
	static constexpr std::array<Gate, 3> singleGates{Gate(gateName::TGate), Gate(gateName::sqrtX), Gate(gateName::sqrtY)};
	static constexpr auto gateInit = Gate(gateName::init);
	static constexpr auto gateWait = Gate(gateName::wait);

	Gate genSingleGate(cords target);
	void generate();
	void evaluate(QubitArray &qubits);
	void evaluateLayer(QubitArray &qubits, unsigned layerIndex);

private:
	int cols, lines;
	int givenDepth;
	int currentDepth;
	std::vector<std::vector<Gate>> layers;

	void applyGate(QubitArray &qubits, Gate gate, cords qubit);
	void addGate(Gate gate, cords qubit);
	Gate lastGateApplyed(cords qubit);
	Gate lastSingleGateApplyed(cords qubit);
	void init(QubitArray &qubits);

	static constexpr int gapBeetwCZ = 2;
	static constexpr int lineSamplesInMask = 2;

	void fillLine(unsigned lineNum, unsigned begin);
	void fillCol(unsigned colNum, unsigned begin);
	std::function<void(int,int)> lineFiller = [this](int lineNum, int begin) { this->fillLine(lineNum, begin); };
	std::function<void(int,int)> colFiller = [this](int lineNum, int begin) { this->fillCol(lineNum, begin); };
	void addLayer(const std::array<int, lineSamplesInMask> &beginPoints, int numOfLines, std::function<void(int,int)> fillerFunc);
};

#endif // RANDQALG_H
