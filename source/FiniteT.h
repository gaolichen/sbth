#pragma once
#include <vector>
#include "Hamiltonian.h"

using namespace std;

class FiniteT
{
private:
	// temprature
	double temp;

	// number of bits
	int M;

	// all energies.
	vector<double> allenergies;

	vector<vector<double> > allstates;
public:
	// constructor
	FiniteT(int bits, double T) { this->M = bits; this->temp = T; }

	// function to calculate all eigen states.
	void CalcEigenStates();
}

