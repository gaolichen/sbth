
#include "EigenEnergyLargeN.h"
#include "Hamiltonian.h"
#include "BitUtility.h"
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;

void EigenEnergyLargeN::Partition(int n, int k, i64 mask)
{
	if (n <= 0)
	{
		// finish a partition, calculate the energies.
	        masks.push_back(mask);
		return;
	}

	if (n > k * (k + 1) / 2)
	{
		// n is too large to partition for the remaining k integers.
		return;
	};
	
	if (k == 0)
	{
		// should not reach here...
		return;
	}

	if (n >= k)
	{
		Partition(n - k, k - 1, mask | (1<<k));
	}
	
	Partition(n, k - 1, mask);
}

// TODO:
double EigenEnergyLargeN::CalcEnergy(vector<int>& modes, int positiveBits)
{
	double ret = .0;
	for (int i = 0; i < M/2; i++)
	{
		if (IsBitSet(positiveBits, i)) ret += sin((i+1) * PI / M);
		else ret -= sin((i+1) * PI / M);
	}

	return ret;
}

void EigenEnergyLargeN::Calculate()
{
	int M2 = M / 2;
	int half = 0;
	if (M % 2 == 0) half = M2;
	for (int n = 0; n * M + half <= M2 * (1 + M2) / 2; n++)
	{
		Partition(n * M + half, M2, (i64)0);
	}

	for (int i = 0; i < masks.size(); i++)
	{
		int zeroBits = BitCount(masks[i]);
		int nonZeroBits = (M2 - 1) / 2 - zeroBits;
		int nonZeroMask = ((1 << M2) - 1) ^ masks[i];
		vector<int> nonZeroDigits = Num2Digit(nonZeroMask, M2);

		for (int j = 0; j < (1 << nonZeroBits); j++)
		{
			if (M % 2 == 0 && (masks[i] & (1<<M2)) != 0)
			{
				energies.push_back(make_pair(CalcEnergy(nonZeroDigits, j), 1 << (zeroBits - 1)));
			}
			else 
			{
				energies.push_back(make_pair(CalcEnergy(nonZeroDigits, j), 1 << zeroBits));
			}
		}
	}
}

void EigenEnergyLargeN::CalculateByEigen(double invN)
{
	if (s != 1)
	{
		cout << "Error: s = " << s << ". CalculateByEigen only support for calculating the s=1 case.";
		return;
	}

	string file = "s=" + ToString(s) + "M=" + ToString(M);
	if (abs(invN) > EPS)
	{
		file += "N=" + ToString((int)round(1/invN + .5));
	}
	
	// file1 to save eigenvalues. 
	string file1 = "EE" + file + ".txt";
	// file2 to save eigenvectors.
	string file2 = "ES" + file + ".txt";
	
	// calculate eigenvalues and eigenvectors.
	H0Hamiltonian h0;
	Eigen::MatrixXcd mat = h0.ToMatrixN(M, Boson, invN);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(mat);

	Eigen::VectorXcd values = solver.eigenvalues();
	Eigen::MatrixXcd vectors = solver.eigenvectors();

	// sort them by eigenvalues.
	vector<int> pos(values.size(), 0);
	for (int i = 0; i < pos.size(); i++) pos[i] = i;
	EigenvalueLess comp(values);
	sort(pos.begin(), pos.end(), comp);

	ofstream ofs1(file1.c_str());

	double eps = 1e-6;
		
	for (int i = 0; i < values.rows(); i++)
	{
		ofs1 << Chop(values[pos[i]].real(), eps) << endl;
	}

	ofs1.close();

	ofstream ofs2(file2.c_str());

	for (int i = 0; i < vectors.cols(); i++)
	{
		ofs2 << i << endl;
		for (int j = 0; j < vectors.rows(); j++)
		{
			complex<double> res = Chop(vectors(j, pos[i]), eps);
			if (res.real() != 0.0 || res.imag() != 0.0)
			{
				ofs2 << j << ' ' << setprecision(3) << res.real() << ' ' << res.imag() << ' ';
			}	
			//ofs2 << setprecision(3) << Chop(vectors(j, pos[i]), eps) << ' ';
		}
		
		ofs2 << endl;
	}

	ofs2.close();
}
