#pragma once
#include <complex>
#include <Eigen/Dense>
#include "HamOperator.h"
#include "StateCollection.h"
#include "SingleTrace.h"
#include "TraceState.h"
#include "StateType.h"

class Hamiltonian
{
protected:
	vector<HamOperator*> realOps;
	vector<HamOperator*> imaginaryOps;
	vector<int> rePrefactors;
	vector<int> imPrefactor;
	string filePrefix;
private:
	//bool inverted;
	//void Init(int xi);
public:
	Hamiltonian();
	//Hamiltonian(int xi, bool invert);
	~Hamiltonian();

	//inline bool Inverted() { return inverted; };
	HamOperator* RealOp(int);
	int RealOpSize();
	void Apply(const TraceState& state, MixState& real, MixState& imaginary);
	void Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm);
	
	// return numerical value of Hamiltonian matrix.
	// bits: number of bits of trace states for which the Hamiltonian apply to.
	// type: type of states, Boson or Fermion.
	// eps: value of 1/N
	// mat: the returned complex matrix.
	Eigen::MatrixXcd ToMatrixN(int bits, StateType type, double eps);
	void AddReadOp(HamOperator *op, int prefactor);
	void AddImaginaryOp(HamOperator *op, int prefactor);
	friend ostream& operator<<(ostream& os, const Hamiltonian& ham);
	string FilePrefix() { return this->filePrefix; }
};

class H0Hamiltonian : public Hamiltonian
{
public:
	H0Hamiltonian();
};

class ZeroHamiltonian : public Hamiltonian
{
public:
	ZeroHamiltonian();
};

class DeltaHamiltonian : public Hamiltonian
{
public:
	DeltaHamiltonian();
};
