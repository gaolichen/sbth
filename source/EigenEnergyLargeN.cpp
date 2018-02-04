
#include "EigenEnergyLargeN.h"
#include "Hamiltonian.h"
#include "BitUtility.h"
#include "StateCollection.h"
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

void EigenEnergyLargeN::CalcAllSingleTraceEnergies()
{
    for (int bit = 1; bit <= M; bit++)
    {
        vector<double> v;
        v.reserve(StateCollection::Inst()->SingleTraceStateNumber(bit));
        double e0 = -4 / tan(PI/(2 * bit));

	    for (int i = 0; i < (1<<bit); i+=2)
	    {
		    double deltE = 0.0;
		    int modes = 0;
		    for (int j = 1; (1<<j) <= i; j++)
		    {
			    if ((i & (1<<j)) == 0) continue;
			    deltE += 8 * sin(j * PI/bit);
			    modes += j;
		    }

		    if ((bit % 2 == 1 && modes % bit == 0) || (bit % 2 ==0 && modes % bit == bit / 2))
		    {
			    v.push_back(e0 + deltE);
		    }
	    }

        //cout << "bit = " << bit << ", state number = " << v.size() << ", expected = " << StateCollection::Inst()->SingleTraceStateNumber(bit) << endl;

        singleTraceEnergies.push_back(v);
    }
}

void EigenEnergyLargeN::BuildBosonicMultiTraceEnergies(int b, int n)
{
    vector<double>& res = statesByBoson[n - 1][b - 1];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b);
    res.reserve(BinomialCoefficient(size + n - 1, n));
    BuildBosonicMultiTraceEnergiesRec(b, singleTraceEnergies[b - 1].size() - 1, n, .0, res);
}

void EigenEnergyLargeN::BuildBosonicMultiTraceEnergiesRec(int b, int index, int n, double energy, vector<double>& res)
{
    if (n == 0)
    {
        res.push_back(Chop(energy));
        return;
    }

    if (index == 0)
    {
        res.push_back(Chop(energy + singleTraceEnergies[b - 1][0] * n));
        return;
    }

    for (int i = 0; i <= n; i++)
    {
        BuildBosonicMultiTraceEnergiesRec(b, index - 1, n - i, energy, res);
        energy += singleTraceEnergies[b - 1][index];
    }
}

void EigenEnergyLargeN::BuildFermionicMultiTraceEnergies(int b, int n)
{
    vector<double>& res = statesByFermion[n - 1][b - 1];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b);
    res.reserve(BinomialCoefficient(size, n));
    BuildFermionicMultiTraceEnergiesRec(b, singleTraceEnergies[b - 1].size() - 1, n, .0, res);
}

void EigenEnergyLargeN::BuildFermionicMultiTraceEnergiesRec(int b, int index, int n, double energy, vector<double>& res)
{
    if (n == 0)
    {
        res.push_back(Chop(energy));
        return;
    }

    if (index + 1 < n)
    {
        return;
    }
    else if (index + 1 == n)
    {
        for (int i = 0; i <= index; i++) 
        {
            energy += singleTraceEnergies[b - 1][i];
        }
        res.push_back(Chop(energy));
        return;
    }

    BuildFermionicMultiTraceEnergiesRec(b, index - 1, n - 1, energy + singleTraceEnergies[b - 1][index], res);
    BuildFermionicMultiTraceEnergiesRec(b, index - 1, n, energy, res);
}

void EigenEnergyLargeN::CalculateByDynamics()
{
        
    // first calculate all single trace state energies.
    CalcAllSingleTraceEnergies();

    // initialize statesByBoson and statesByFermion
    statesByBoson.reserve(M);
    statesByFermion.reserve(M);

    // the first element of statesByBoson is simply singleTraceEnergies
    statesByBoson.push_back(singleTraceEnergies);
    // for statesByFermion, we let its first element be empty, because we do
    // not need it in the following calculation.
    statesByFermion.push_back(vector<vector<double> >());

    for (int i = 2; i <= M; i++)
    {
        statesByBoson.push_back(vector<vector<double> >(M/i, vector<double>()));
        statesByFermion.push_back(vector<vector<double> >(M/i, vector<double>()));
        for (int bit = 1; bit * i <= M; bit++)
        {
            cout << "(" << i << ", " << bit << ") single trace states #:";
            cout << StateCollection::Inst()->SingleTraceStateNumber(bit);

            BuildBosonicMultiTraceEnergies(bit, i);
            cout << " statesByBoson: " << statesByBoson[i - 1][bit - 1];
            if (i <= StateCollection::Inst()->SingleTraceStateNumber(bit))
            {
                BuildFermionicMultiTraceEnergies(bit, i);
                cout << ", statesByFermion: " << statesByFermion[i - 1][bit - 1];
            }
            cout << endl;
        }
    }

    // build statesByBoth
    BuildStatesByBoth();

    BuildAllStates();
}

void EigenEnergyLargeN::BuildAllStates()
{
    allStates.reserve(StateCollection::Inst()->StateNumber(M));
    BuildAllStatesRec(M, M, .0, 0);
    sort(allStates.begin(), allStates.end());
    cout << "total states number: " << allStates.size() << endl;
    //cout << "eigenenergies: " << allStates << endl;
}

void EigenEnergyLargeN::BuildAllStatesRec(int bitRemain, int singleTraceBits, double energy, int deg)
{
    if (bitRemain == 0 || singleTraceBits == 1)
    {
        int n = (1 << (deg - 1));
        if (bitRemain > 0) n = (1 << deg);
        for (int i = 0; i < n; i++)        
            allStates.push_back(energy);
        return;
    }

    if (singleTraceBits > bitRemain)
    {
        BuildAllStatesRec(bitRemain, bitRemain, energy, deg);
        return;
    }

    // skip this single trace
    BuildAllStatesRec(bitRemain, singleTraceBits - 1, energy, deg);

    // pick one single trace.
    for (int i = 0; i < singleTraceEnergies[singleTraceBits - 1].size(); i++)
    {
        BuildAllStatesRec(bitRemain - singleTraceBits, singleTraceBits - 1, 
            energy + singleTraceEnergies[singleTraceBits - 1][i], deg + 1);
    }

    // pick two of more.
    for (int i = 2; i * singleTraceBits <= bitRemain; i++)
    {
        for (int j = 0; j < statesByBoth[i - 1][singleTraceBits - 1].size(); j++)
        {
            BuildAllStatesRec(bitRemain - singleTraceBits * i, singleTraceBits - 1, 
                energy + statesByBoth[i - 1][singleTraceBits - 1][j], deg + 1);
        }
    }
}

void EigenEnergyLargeN::BuildStatesByBoth()
{
    statesByBoth.reserve(M);
    statesByBoth.push_back(vector<vector<double> >());
    
    for (int i = 2; i <= M; i++)
    {
        statesByBoth.push_back(statesByBoson[i - 1]);
        for (int bit = 1; bit * i <= M; bit++)
        {
            vector<double>& v = statesByBoth[i - 1][bit - 1];
            for (int k = 2; k < i; k += 2)
            {
                if (statesByFermion[k - 1][bit - 1].size() > 0)
                {
                    vector<double> toAdd;
                    MergeEnergy(statesByFermion[k - 1][bit - 1], statesByBoson[i - k - 1][bit - 1], toAdd);
                    v.insert(v.end(), toAdd.begin(), toAdd.end());
                }
            }

            if (i % 2 == 0)
            {
                v.insert(v.end(), statesByFermion[i - 1][bit - 1].begin(), statesByFermion[i - 1][bit - 1].end());
            }

            cout << "(" << i << ", " << bit << "): " << statesByBoth[i - 1][bit - 1] << endl;
        }
    }
}

void EigenEnergyLargeN::MergeEnergy(vector<double>& a, vector<double>& b, vector<double>& res)
{
    res.reserve(a.size() * b.size());
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b.size(); j++)
        {
            res.push_back(Chop(a[i] + b[j]));
        }
    }
}

void EigenEnergyLargeN::CalculateByEigen(double invN, bool calcEigenvector)
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
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(mat, calcEigenvector);
	Eigen::VectorXcd values = solver.eigenvalues();

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

	if (calcEigenvector == false)
	{
		return;
	}

	Eigen::MatrixXcd vectors = solver.eigenvectors();
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
