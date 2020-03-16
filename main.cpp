#include <iostream>
#include <QuEST.h>
#include "qubitarray.h"
#include <unistd.h>
#include <vector>
#include <cmath>
#include <string>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

void startSimulation(int maxQubits, int avg,
					 double singleErrRate, double multiErrRate, double envCoupling,
					 double singleCoupling, double multiCoupling)
{
	QuESTEnv env = createQuESTEnv();
	QubitArray qubits(env, 1, 3);
	qubits.setEnvCoupling(envCoupling);
	qubits.setSingleErrRate(singleErrRate);
	qubits.setMultiErrRate(multiErrRate);
	qubits.setSingleGateCoupling(singleCoupling);
	qubits.setMultiGateCoupling(multiCoupling);

	cout << "Path \t fidelity \t fidelity std error\n";
	for(int path = 1; path <= maxQubits - 2; path++)
	{
		cerr << "Path " << path << "\n";
		qubits.setSize(1, path + 2);
		double avgFidelity = 0;
		vector <double> fidelity;
		for(int i = 0; i < avg; i++)
		{
			qubits.init();
			qubits.generateBell({0, 0}, {0, 1});
			qubits.move({0, 1}, {0, 1 + path});
			double curFidelity = qubits.calcBellFidelity({0, 0}, {0, 1 + path});
			avgFidelity += curFidelity;
			fidelity.push_back(curFidelity);
		}

		avgFidelity /= avg;
		double squaredErrSum = 0, standDiv = 0;
		if(avg > 1)
		{
			for(auto f: fidelity)
				squaredErrSum += (avgFidelity - f) * (avgFidelity - f);
			standDiv = sqrt( squaredErrSum / ((avg - 1) * avg) );
		}

		cout << path << " \t " << avgFidelity << " \t " << standDiv << "\n";
	}
}

int main()
{
	Config cfg;
	int qubits = 10, avg = 1000;
	double singleErrRate = 0.01, multiErrRate = 0.1, envCoupling = 0.2, singleCoupling = 0.2, multiCoupling = 0.2;
	try
	{
		cfg.readFile("simulatorConfig.cfg");
		cfg.lookupValue("regSize", qubits);
		cfg.lookupValue("averaging", avg);
		cfg.lookupValue("singleErrRate", singleErrRate);
		cfg.lookupValue("multiErrRate", multiErrRate);
		cfg.lookupValue("envCoupling", envCoupling);
		cfg.lookupValue("singleCoupling", singleCoupling);
		cfg.lookupValue("multiCoupling", multiCoupling);

		cerr << "Parameters found: "<< qubits << " " << avg << " " << singleErrRate << " " << multiErrRate
			 << " " << envCoupling << " " << singleCoupling << " " << multiCoupling << "\n";
	}
	catch(const FileIOException &fioex)
	{
		cerr << "I/O error while reading file." << std::endl;
	}
	catch(const ParseException &pex)
	{
		cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
				 << " - " << pex.getError() << std::endl;
	}

	string filename = "Simulation" + to_string(qubits) +  "Qubits"
			+ "ErrRateSingle" + to_string(singleErrRate) + "Multi" + to_string(multiErrRate)
			+ "Averaging" + to_string(avg) + ".txt";

	freopen(filename.c_str(), "w", stdout);
	startSimulation(qubits, avg, singleErrRate, multiErrRate, envCoupling, singleCoupling, multiCoupling);


	return 0;
}
