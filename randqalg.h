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
	std::vector<std::vector<Gate>>& getAlgorithm();
	void setAlgorithm(std::vector<std::vector<Gate>> &algorithm);
	void generate();
	void evaluate(QubitArray &qubits);
	void evaluateLayer(QubitArray &qubits, int layerIndex);

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

	void fillLine(int lineNum, int begin);
	void fillCol(int colNum, int begin);
	void addLayer(const std::array<int, lineSamplesInMask> &beginPoints, int numOfLines, void (RandQAlg::*fillerFunc)(int, int));
};

#endif // RANDQALG_H
