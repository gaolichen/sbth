#pragma once
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"
#include "StringBitMatrices.h"
#include <Eigen/Dense>

#define USE_CACHE

using namespace std;

class VevCalculator
{
private:
	MatrixSB* matM;
	MatrixSB* omegaV;
	MatrixSB* omegaW;
	CDT gammaV;
	CDT gammaW;
	int M;
	map<int, CDT> cache;
#ifdef USE_CACHE
	map<int, CDT> cacheAll;
#endif
	CDT DoCalculate(int ops);
public:
	VevCalculator(MatrixSB &matM_, MatrixSB &omegaV_, MatrixSB &omegaW_, CDT gammaV_, CDT gammaW_)
	{
		matM = &matM_;
		omegaV = &omegaV_;
		omegaW = &omegaW_;
		gammaV = gammaV_;
		gammaW = gammaW_;
		M = matM_.innerSize();
	}

	CDT CalculateVev(int ops)
	{
		if (BitCount(ops) % 2 == 1) return 0.0;
		return DoCalculate(ops);
	}

	CDT VevV(int ops);
	CDT VevW(int ops);
	CDT VevAll(int ops);

	CDT CalculateVev(int ops, MatrixSB &omega);
};

class StateInfo
{
public:
	int Ops;
	DT Energy;
	StateInfo(int ops, DT energy):Ops(ops), Energy(energy)
	{
	}

};

class EnergyCalculator
{
private:
	int s;
	double xi;
	int totalStates;
	double calculateTime;
	Stopwatch watch;
	vector<DT> normalizedE;
	DT EnergyCorrection(int M, int L);
public:
	EnergyCalculator() : s(1) {}

	EnergyCalculator(int s_, double xi_ = .0) : s(s_), xi(xi_) {}

	//static CDT OperatorVev(int ops, int M, VevCalculator &calc, CDT &gamma, MatrixSB &omega);

	//static CDT OperatorVevAllZeros(int ops, int M, VevCalculator &calc, CDT &gammaW, MatrixSB &omegaW, CDT &gammaV, MatrixSB &omegaV);

	static vector<StateInfo> AllStates(int M);

	static DT Operators2Energy(vector<int>& ops, int M, int s);

	DT EnergyCorrection(int M, bool outputCorrections = false);

	int TotalStates()
	{
		return totalStates;
	}
	
	double CalculateTime()
	{
		return calculateTime;
	}
};
