#pragma once
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
using namespace std;

#define EPS 1e-8

double chop(double v)
{
	if (abs(v) < EPS) return 0.0;
	return v;
}

complex<double> chop(complex<double>& v)
{
	return complex<double>(chop(v.real()), chop(v.imag()));
}

void template<T> chop(Eigen::MatrixBase<T>& mat)
{
	for (int i = 0; i < mat.rows(); i++)
	{
		for (int j = 0; j < mat.columns(); j++)
		{
			mat(i, j) = chop(mat(i, j));
		}
	}
}

