#pragma once
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"
#include <Eigen/Dense>

#define SYMMETRIC_A

using namespace std;

#define esp 1e-10
typedef double DT;
typedef std::complex<DT> CDT;
typedef Eigen::MatrixXcd MatrixSB;
extern DT PI;
extern CDT I;
extern CDT TPI;

CDT ByPolar(DT r, DT angle);

CDT ByPolar(DT angle);

DT Chop(DT a);

CDT Chop(CDT& a);

class StringBitMatrices
{
private:
	int M;
	int L;
	CDT DotV1(int M, int L, int K, int m, int n);
	CDT DotV2(int M, int L, int K, int m, int n);
	CDT DotW(int M, int L, int K, int m);
	CDT Cm0(int M, int L, int K, int m);
	CDT Cmn1(int M, int L, int K, int m, int n);
	CDT Cmn2(int M, int L, int K, int m, int n);
	CDT CmM(int M, int L, int K, int m);
	CDT Smn1(int M, int L, int K, int m, int n);
	CDT Smn2(int M, int L, int K, int m, int n);
	CDT SmM(int M, int L, int K, int m);
	CDT ElementC(int M, int L, int K, int m, int n);
	CDT ElementS(int M, int L, int K, int m, int n);
	CDT ElementA(int k, int l, int M, int n, int m);
public:
	StringBitMatrices() {}

	MatrixSB MatrixC(int M, int L);
	MatrixSB MatrixS(int M, int L);
	MatrixSB MatrixM(int M, int L);
	MatrixSB MatrixA(int k, int l, int M);
	MatrixSB MatrixAV(int M, int L);
	MatrixSB MatrixAW(int M, int L);
	MatrixSB MatrixBV(int M, int L);
	MatrixSB MatrixBW(int M, int L);
	CDT MuPV(int M, int L, DT xi = .0);
	CDT MuPW(int M, int L, DT xi = .0);
};
