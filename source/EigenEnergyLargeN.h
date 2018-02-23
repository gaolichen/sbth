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

    // divide the energy eigenvalues into buckets. 
    // states: the energy eigenvalues to divide
    // buckets: number of buckets
    // range: all eigenvalues lie in the interval (-range, range).
    static vector<pair<double, int> > BucketEnergies(vector<double>& states, int buckets, double range);

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
    // buckets: number of buckets the energies are divided into. buckets = 0 means do not bucket.
    // calcEigenvector: true to calculate eigenvectors as well, otherwise only calculate eigenvalues.
	void CalculateByEigen(double invN = DefaultInvN, int buckets = 0, bool calcEigenvector = true);

    // this function calculates all energies of single trace with bit number <= M and stores 
    // results in singleTraceEnergies.
    void CalcAllSingleTraceEnergies();

	void Calculate();

	void CalculateByDynamics();

    // calculate thermodynamics of string bit.
    // T0: the string tension.
    // maxB: max value of beta, beta = 1/temperatur, the min value of beta is 0.0
    // steps: how many data sets to calculate. The interval between two connective data sets is maxB/steps.
    void CalculateThermo(double T0, double maxB, int steps = 100);

    // calculate thermodynamics of string bit for M=1,2,..,N.
    // the partition function Z is the sum of indidual Z of each M.
    // T0: the string tension.
    // maxB: max value of beta, beta = 1/temperatur, the min value of beta is 0.0
    // datafolder: the path to the foler that contains energy eigenvalue data.
    // steps: how many data sets to calculate. The interval between two connective data sets is maxB/steps.
    void CalculateThermoForN(double T0, double maxB, string datafolder, int steps = 100);

    void CalcFluctuation(double beta);

	// save eigen energies to file
	// buckets: number of buckets the energies are divided into. buckets = 0 means do not bucket.
	void SaveEnergies(int buckets = 0);

    // save single-trace eigen energies to file.
    // buckets: number of buckets the energies are divided into. buckets = 0 means do not bucket.
    void SaveSingleEnergies(int bit, int buckets = 0);

	const vector<pair<double, int> > Energies() { return energies; }
};

