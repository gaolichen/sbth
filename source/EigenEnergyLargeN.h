#pragma once
#include <string>
#include <vector>
#include "BitUtility.h"
#include "EigenCommon.h"

using namespace std;

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
	static const double DefaultInvN = EPS;

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
	void CalculateByEigen(double invN = DefaultInvN);

	void Calculate();

	//void LoadFromFile(string file);

	const vector<pair<double, int> > Energies() { return energies; }
};

