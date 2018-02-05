#pragma once
#include <string>
#include <vector>
#include <map>
#include "BitUtility.h"
#include "EigenCommon.h"

using namespace std;

#define DefaultInvN EPS

// struct provides comparision function for sorting eigenstates.
struct EigenvalueLess
{
	Eigen::VectorXcd* values;
	EigenvalueLess(Eigen::VectorXcd& v)
	{
		this->values = &v;
	}

	bool operator() (int a, int b) const
	{
		return (*values)[a].real() <(*values)[b].real();
	}
};

class EigenEnergyLargeN
{
private:
	int M;
	int s;
	vector<i64> masks;

    // single trace eigen energies.
    // singleTraceEnergies[b][i]: the i-th eigenenergy of b-bit single trace eigenstate.
    vector<vector<double> > singleTraceEnergies;
    
    // states built by bosonic single trace states.
    // statesByBoson[k, b]: eigenenergies built out of k number of b-bit bosonic single trace state.
    vector<vector<vector<double> > > statesByBoson;

    // states built by fermionic single trace states.
    // statesByFermion[k, b]: eigenenergies built out of k number of b-bit fermionic single trace state.
    vector<vector<vector<double> > > statesByFermion;

    // statesByBoth[k, b]: eigenenergies built out of k number of b-bit single trace state.
    vector<vector<vector<double> > > statesByBoth;

    // all M-bit bosonic eigenenergies.
    vector<double> allStates;

	// list of eigen energies: for each element, the first component is the value of energy, 
	// and the second component is the degeneracy of the enerties.
	vector<pair<double, int> > energies;
	
	// calculate energy from two numbers.
	// nonZeroBits: the nonzero bit position of mask represents a nonzero mode, the mode can make 
	// either one unit position contribution or one unit negative contribution to the energy, 
	// where each unit mean sin(k * Pi/M);
	// bits: the nonzero bit position of bits represents a positive-contribution bit and zero 
	// bit position represents a negative-contribution bit.
	double CalcEnergy(vector<int>& nonZeroBits, int positiveBits);

	void Partition(int n, int k, i64 mask);

    // this function calculates all energies of single trace with bit number <= M and stores 
    // results in singleTraceEnergies.
    void CalcAllSingleTraceEnergies();

    // find all bosonic multi-trace energies built out of n number of b-bit single trace states.
    // store results in statesByBoson.
    void BuildBosonicMultiTraceEnergies(int b, int n);
    
    // the recursive function invoked by BuildBosonicMultiTraceEnergies.
    // index: pick n b-bit single traces from collection {0, 1, ... index}
    // energy: current value of energy.
    // res: the vector to stotre results.
    void BuildBosonicMultiTraceEnergiesRec(int b, int index, int n, double energy, vector<double>& res);

    // find all fermionic multi-trace energies built out of n number of b-bit single trace states.
    // store results in statesByFermion.
    void BuildFermionicMultiTraceEnergies(int b, int n);

    // the recursive function invoked by BuildFermionicMultiTraceEnergies.
    // index: pick n b-bit single traces from collection {0, 1, ... index}
    // energy: current value of energy.
    // res: the vector to stotre results.
    void BuildFermionicMultiTraceEnergiesRec(int b, int index, int n, double energy, vector<double>& res);
    
    // build statesByBoth
    void BuildStatesByBoth();
    
    // build allStates.
    void BuildAllStates();

    // The recursive function invoked by BuildAllStates().
    void BuildAllStatesRec(int bitRemain, int singleTraceBits, double energy, int deg);

    // merge two lists of eigenenergies to one.
    // res will be of size size(a)*size(b), and each element is of the form a[i] + b[j].
    static void MergeEnergy(vector<double>& a, vector<double>& b, vector<double>& res);

public:
	// constructor.
	// bits: number of bits of the stringbit system, i.e, the M in the paper
	// spin: number of spin d.o.f of the stringbit system, i.e, the s in the paper.
	EigenEnergyLargeN(int bits, int spin = 1)
	{
		this->M = bits;
		this->s = spin;
	}
	
	// calculate all eigen energies and eigenstates using EigenLibrary
	// and save the results to files. The file name for eigen energies data 
	// is of the form  EEs=1M=5.txt and the one for eigen states is of the form
	// ESs=1M=5.txt. 
	// invN: the value of 1/N, should be very small
	void CalculateByEigen(double invN = DefaultInvN, bool calcEigenvector = true);

	void Calculate();

	void CalculateByDynamics();

	//void LoadFromFile(string file);

	const vector<pair<double, int> > Energies() { return energies; }
};

