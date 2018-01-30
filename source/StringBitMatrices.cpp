#pragma warning(disable:4018)
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"
#include <Eigen/Dense>
#include "StringBitMatrices.h"
using namespace std;

DT PI = acos(-1.0);
CDT I(0,1);
CDT TPI = I * CDT(2 * PI);

CDT ByPolar(DT r, DT angle)
{
        return CDT(r*cos(angle) + I*r*sin(angle));
}

CDT ByPolar(DT angle)
{
        return CDT(cos(angle) + I*sin(angle));
}

DT Chop(DT a)
{
	if (abs(a) < esp) return 0;
	return a;
}

CDT Chop(CDT& a)
{
	return CDT(Chop(a.real()), Chop(a.imag()));
}


CDT StringBitMatrices::DotV1(int M, int L, int K, int m, int n)
{
        if((n * M - L * m) % (L * M) == 0)
        {
                return sqrt(L/(DT)M);
        }
        else
        {
                return -1/sqrt((DT)M * L) * ((DT)1 - ByPolar(-2 * PI * m * L / M)) /((DT)1 - ByPolar(-2 * PI * (n / (DT)L - m / (DT)M)));
        }
}

CDT StringBitMatrices::DotV2(int M, int L, int K, int m, int n)
{
        if((n * M - K * m) % (K * M) == 0)
        {
                return sqrt(K/(DT)M) * ByPolar(-2 * PI * n * L / K);
        }
        else
        {
                return 1/sqrt((DT)M * K) * ((DT)1 - ByPolar(-2 * PI * m * L / M) ) /((DT)1 - ByPolar(-2 * PI * (n / (DT)K - m / (DT)M)));
        }
}

CDT StringBitMatrices::DotW(int M, int L, int K, int m)
{
        if (m == 0) return 0;

        return -1/sqrt((DT)L * K) * ((DT)1.0 - ByPolar(- 2 * PI * m * L / M)) / ((DT)1.0 - ByPolar(2 * PI * m / M));
}

CDT StringBitMatrices::Cm0(int M, int L, int K, int m)
{
        if (m == 0) return 1;
        return 0;
}

CDT StringBitMatrices::Cmn1(int M, int L, int K, int m, int n)
{
        return DotV1(M, L, K, m, n) * cos(n * PI / (2 * L) - PI * m / (2 * M));
}

CDT StringBitMatrices::Cmn2(int M, int L, int K, int m, int n)
{
        return DotV2(M, L, K, m, n) * cos(n * PI / (2 * K) - PI * m / (2 * M));
}

CDT StringBitMatrices::CmM(int M, int L, int K, int m)
{
        return DotW(M, L, K, m) * cos(m * PI / (2 * M) - PI / 4);// * ByPolar(-PI / 4);
}

CDT StringBitMatrices::Smn1(int M, int L, int K, int m, int n)
{
        return DotV1(M, L, K, m, L - n) * cos(n * PI / (2 * L) + PI * m / (2 * M));
}

CDT StringBitMatrices::Smn2(int M, int L, int K, int m, int n)
{
        return DotV2(M, L, K, m, K - n) * cos(n * PI / (2 * K) + PI * m / (2 * M));
}

CDT StringBitMatrices::SmM(int M, int L, int K, int m)
{
        return DotW(M, L, K, m) * cos(m * PI / (2 * M) + PI / 4); //* ByPolar(PI / 4);
}

CDT StringBitMatrices::ElementC(int M, int L, int K, int m, int n)
{
        if (m == 0) return Cm0(M, L, K, n);
        if (n == 0) return 0;
        if (n < L) return Cmn1(M, L, K, m, n);
        if (n < M - 1) return Cmn2(M, L, K, m, n - L + 1);
        assert(n == M - 1);
        return CmM(M, L, K, m);
}

CDT StringBitMatrices::ElementS(int M, int L, int K, int m, int n)
{
        if (m == 0) return 0;
        if (n == 0) return 0;
        if (n < L) return Smn1(M, L, K, m, n);
        if (n < M - 1) return Smn2(M, L, K, m, n - L + 1);
        assert(n == M - 1);
        return SmM(M, L, K, m);
}

MatrixSB StringBitMatrices::MatrixC(int M, int L)
{
        MatrixSB ret(M, M);
        for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                        ret(i, j) = ElementC(M, L, M - L, i, j);

        return ret;
}

MatrixSB StringBitMatrices::MatrixS(int M, int L)
{
        MatrixSB ret(M, M);
        for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                        ret(i, j) = ElementS(M, L, M - L, i, j);

        return ret;
}

MatrixSB StringBitMatrices::MatrixM(int M, int L)
{
        return -MatrixC(M, L).inverse() * MatrixS(M, L);
}

CDT StringBitMatrices::ElementA(int k, int l, int M, int n, int m)
{
        return ByPolar(2 * PI * k * (m + n) / M) * sin((m - n) * PI / (2 * M))
                + ByPolar(2 * PI * (k + l) * (m + n) / (2 * M)) * sin((l - k - .5) * (m - n) * PI / M);
}

MatrixSB StringBitMatrices::MatrixA(int k, int l, int M)
{
        MatrixSB ret(M, M);
        for (int i = 0; i < M; i++)
        {
                for (int j = 0; j < M; j++)
                {
                        ret(i, j) = ElementA(k, l, M, i, j);
                }
        }

        return ret;
}

MatrixSB StringBitMatrices::MatrixAV(int M, int L)
{
	MatrixSB ret(M, M);
        for (int n = 0; n < M; n++)
        {
                for (int m = 0; m < M; m++)
                {
#ifdef SYMMETRIC_A
			//return (MatrixA(M, L + 1, M) + MatrixA(L, 1, M)) * .5;
			ret(n, m) = (1.0 + ByPolar(2 * PI * L * (m + n) / M)
				+ ByPolar(PI * (m + n + 2 * L * n) / M) 
				+ ByPolar(PI * (m + n + 2 * L * m) / M)) * 0.5 * sin((m-n) * PI / (2*M));
#else
                        ret(n, m) = sin((m - n) * PI / (2 * M))
                                + ByPolar(2 * PI * (L + 1) * (m + n) / (2 * M)) * sin((L + 0.5) * (m - n) * PI / M);
#endif
                }
        }

	return ret;
}

MatrixSB StringBitMatrices::MatrixAW(int M, int L)
{
        MatrixSB ret(M, M);
        for (int n = 0; n < M; n++)
        {
                for (int m = 0; m < M; m++)
                {
                        ret(n, m) = ((DT)1 + ByPolar(PI * (m + n) / M)) * ((DT)1
                                + ByPolar(2 * PI * L * (m + n) / M)) * sin((m - n) * PI / (2 * M));
                }
        }

        return ret;
}

MatrixSB StringBitMatrices::MatrixBV(int M, int L)
{
        MatrixSB invc = MatrixC(M, L).inverse();
        return invc * MatrixAV(M, L).adjoint() * invc.transpose();
}

MatrixSB StringBitMatrices::MatrixBW(int M, int L)
{
        MatrixSB invc = MatrixC(M, L).inverse();
        return invc * MatrixAW(M, L).adjoint() * invc.transpose();
}

CDT StringBitMatrices::MuPV(int M, int L, DT xi)
{
        CDT tr = (MatrixS(M, L).conjugate() * MatrixC(M, L).inverse() * MatrixAV(M, L).adjoint()).trace();
#ifdef SYMMETRIC_A
		return -1/tan(PI / (2 * M)) + (1 / tan((2 * L - 1) * PI / (2 * M)) - 1 / tan((2 * L + 1) * PI / (2 * M))) * 0.5 + M * xi - tr;
#else
        return -1/tan(PI / (2 * M)) - 1 / tan((2 * L + 1) * PI / (2 * M)) + M * xi - tr;
#endif
}

CDT StringBitMatrices::MuPW(int M, int L, DT xi)
{
        CDT tr = (MatrixS(M, L).conjugate() * MatrixC(M, L).inverse() * MatrixAW(M, L).adjoint()).trace();
        return -4/tan(PI / (2 * M)) + 2 * M * xi - tr;
}

