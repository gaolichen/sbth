#pragma once
#include <string>
#include <vector>
#include <fstream>
#include "StateType.h"
#include "Hamiltonian.h"
using namespace std;

class ScriptGenerator
{
private:
	string rootFolder;
	StateType type;
	void GroupStates(int bits, vector<vector<int> >& states);
	double NormMatrixScale(int bits, int index);
	const TraceState& GetState(int m, int index);
#ifdef HAM_PARAMETER
	Hamiltonian h0;
	Hamiltonian deltaH;
	Hamiltonian hPrime;
	void InitHamiltonians();
#endif
	void HamToMatlab(int bits, Hamiltonian& ham, ofstream& os, bool sparse = false);
	void VectorToMatlab(vector<int>& a, ofstream& os);
	void SparseMatrixToDataFile(vector<vector<Coefficient> > &mat, string file);
public:
	ScriptGenerator(string rootFolder, StateType stateType);

	void OutputHamToMatlab(int bits, Hamiltonian& ham, bool sparse = false);
	void OutputHamToMatlab(int bits, bool inverted);
	void OutputHamToDataFile(int bits, Hamiltonian& ham, bool isComplex = false);
	void OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename);
	void OutputNormToLaTeX(int minBits, int maxBits, string filename);
	void OutputNormToMatlab(int bits);
	void OutputNormToDataFile(int bits, int fNumber);
	void OutputStateToLaTeX(int minBit, int maxBit, string filename, bool compress = false);
	void OutputStateStructure(int maxBit = 11);
};
