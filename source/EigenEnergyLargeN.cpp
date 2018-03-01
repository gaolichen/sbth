#include "EigenEnergyLargeN.h"
#include "Hamiltonian.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
using namespace std;

#define HIGH_PRECISION

#ifdef HIGH_PRECISION
typedef long double Db;
#else
typedef double Db;
#endif

/*
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
}*/

void EigenEnergyLargeN::CalcAllSingleTraceEnergies()
{
    // if already calculated, return.
    if (singleTraceEnergies.size() > 0)
    {
        return;
    }

    singleTraceEnergies.reserve(M + 1);
    singleTraceEnergies.push_back(vector<double>(1, 0.0));
    i64 allOne = (i64)((1<<s) - 1);

    for (int bit = 1; bit <= M; bit++)
    {
        int Mhalf = 0;
        if ((s * (bit - 1)) % 2 == 1)
        {
            Mhalf = bit / 2;
        }

        vector<double> v;
        v.reserve(StateCollection::Inst()->SingleTraceStateNumber(bit, s));
        double e0 = -4 / tan(PI/(2 * bit)) * s;

	    for (i64 i = 0; i < ((i64)1<<(bit * s)); i += (1<<s))
	    {
		    double deltE = 0.0;
		    int modes = 0;
		    for (int j = 1; ((i64)1<<(j*s)) <= i; j++)
		    {
                i64 mask = i & (allOne << (j * s));
			    if (mask == 0) continue;
			    deltE += 8 * sin(j * PI/bit) * BitCount(mask);
			    modes += j * BitCount(mask);
		    }

		    if ((modes + Mhalf) % bit == 0)
		    {
			    v.push_back(e0 + deltE);
		    }
	    }

        //cout << "bit = " << bit << ", state number = " << v.size() << ", expected = " << StateCollection::Inst()->SingleTraceStateNumber(bit) << endl;

        sort(v.begin(), v.end());
        singleTraceEnergies.push_back(v);
    }
}

void EigenEnergyLargeN::BuildBosonicMultiTraceEnergies(int b, int n)
{
    vector<double>& res = statesByBoson[n][b];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b);
    res.reserve(BinomialCoefficient(size + n - 1, n));
    BuildBosonicMultiTraceEnergiesRec(b, singleTraceEnergies[b].size() - 1, n, .0, res);
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
        res.push_back(Chop(energy + singleTraceEnergies[b][0] * n));
        return;
    }

    for (int i = 0; i <= n; i++)
    {
        BuildBosonicMultiTraceEnergiesRec(b, index - 1, n - i, energy, res);
        energy += singleTraceEnergies[b][index];
    }
}

void EigenEnergyLargeN::BuildFermionicMultiTraceEnergies(int b, int n)
{
    vector<double>& res = statesByFermion[n][b];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b);
    res.reserve(BinomialCoefficient(size, n));
    BuildFermionicMultiTraceEnergiesRec(b, singleTraceEnergies[b].size() - 1, n, .0, res);
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
            energy += singleTraceEnergies[b][i];
        }
        res.push_back(Chop(energy));
        return;
    }

    BuildFermionicMultiTraceEnergiesRec(b, index - 1, n - 1, energy + singleTraceEnergies[b][index], res);
    BuildFermionicMultiTraceEnergiesRec(b, index - 1, n, energy, res);
}

void EigenEnergyLargeN::CalculateByDynamics()
{
        
    // first calculate all single trace state energies.
    CalcAllSingleTraceEnergies();

    // initialize statesByBoson and statesByFermion
    statesByBoson.reserve(M + 1);
    statesByFermion.reserve(M + 1);

    // the (0, i) elements are vacuum states.
    statesByBoson.push_back(vector<vector<double> >(M + 1, vector<double>(1, 0.0)));
    statesByFermion.push_back(vector<vector<double> >(M + 1, vector<double>(1, 0.0)));

    // the first element of statesByBoson is simply singleTraceEnergies
    statesByBoson.push_back(singleTraceEnergies);
    // for statesByFermion, we let its first element be empty, because we do
    // not need it in the following calculation.
    statesByFermion.push_back(vector<vector<double> >());

    for (int i = 2; i <= M; i++)
    {
        statesByBoson.push_back(vector<vector<double> >(M/i + 1, vector<double>()));
        statesByFermion.push_back(vector<vector<double> >(M/i + 1, vector<double>()));
        for (int bit = 1; bit * i <= M; bit++)
        {
            //cout << "(" << i << ", " << bit << ") single trace states #:";
            //cout << StateCollection::Inst()->SingleTraceStateNumber(bit);

            BuildBosonicMultiTraceEnergies(bit, i);
            //cout << " statesByBoson: " << statesByBoson[i][bit];
            if (i <= StateCollection::Inst()->SingleTraceStateNumber(bit))
            {
                BuildFermionicMultiTraceEnergies(bit, i);
                //cout << ", statesByFermion: " << statesByFermion[i][bit];
            }
            //cout << endl;
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
    cout << "total states number: " << allStates.size() << ", expected: " << StateCollection::Inst()->StateNumber(M) << endl;
    //cout << "eigenenergies: " << allStates << endl;
}

void EigenEnergyLargeN::CalculateThermo(double T0, double maxB, int steps)
{
    double minB = .0;
    double delta = (maxB - minB) / steps;

    string file = "THs=" + ToString(s) + "M=" + ToString(M) + "T0=" + ToString(T0) + ".txt";
	ofstream ofs(file.c_str());

    //ofs << "T Z Energy Entropy" << endl;

    for (Db beta = minB; beta <= maxB + 1e-8; beta += delta)
    {
        Db Z = .0;
        Db E = .0;

        for (int i = allStates.size() - 1; i >= 0; i--)
        {
            Db rho = exp(-beta * (T0 * allStates[i] + M)/sqrt(2));
            Z += rho;
            E += rho * (T0 * allStates[i] + M)/sqrt(2);
        }

        E /= Z;

        Db entropy = .0;

        for (int i = allStates.size() - 1; i >= 0; i--)
        {
            Db rho = exp(-beta * (T0 * allStates[i] + M)/sqrt(2)) / Z;
            entropy -= rho * log(rho);
        }

        ofs << beta << " " << log(Z) << " " << E << " " << entropy << endl;
    }

    ofs.close();
}

void EigenEnergyLargeN::CalculateThermoForN(double T0, double maxB, string datafolder, int steps)
{
    double minB = .0;
    double delta = (maxB - minB) / steps;

    // first, we need to find all the energies eigenvalues at finite N.
    vector<double> energies;
    energies.reserve(1<<(M+1));

    H0Hamiltonian h0;
    
    for (int bit = 1; bit <= M; bit++)
    {
        cout << "calculate energies for bit=" << bit << endl;
        if (bit > 8)
        {
            // load from file.
            string input = datafolder + "/EEs=" + ToString(s) + "M=" + ToString(bit) + "N=" + ToString(M) + ".txt";
            ifstream ifs(input.c_str());

            cout << "file: " << input << " " << (ifs.good() ? "good" : "fuck") << endl;

            if (!ifs.good())
            {
                return;
            }

            double energy;
            while (!ifs.eof())
            {
                ifs >> energy;
                // the energy is P_0 = (P_{+} + P_{-})/sqrt(2) and
                // to be consistent with Thorn's paper, we shift P_{-} by 8*T0/PI * M
                energies.push_back((T0 * (energy + 8 * bit /PI) + bit) / sqrt(2));
            }

            ifs.close();
        }
        else
        {
            // compute with eigen.
            Eigen::MatrixXcd mat = h0.ToMatrixN(bit, Boson, 1.0/M);
	        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
            solver.compute(mat, false);
	        Eigen::VectorXcd values = solver.eigenvalues();

            for (int i = 0; i < values.size(); i++)
            {
                // the energy is P_0 = (P_{+} + P_{-})/sqrt(2) and
                // to be consistent with Thorn's paper, we shift P_{-} by 8*T0/PI * M
                energies.push_back((T0 * (values[i].real() + 8 * bit /PI) + bit) / sqrt(2));
            }
        }
    }

    sort(energies.begin(), energies.end());

    cout << "Total number of energies states = " << energies.size() << endl;

    // then we calculate thermodynamics and save the results to file.
    string file = "THs=" + ToString(s) + "N=" + ToString(M) + "T0=" + ToString(T0) + ".txt";
	ofstream ofs(file.c_str());

    //ofs << "T Z Energy Entropy" << endl;

    for (Db beta = minB; beta <= maxB + 1e-8; beta += delta)
    {
        Db Z = .0;
        Db E = .0;

        for (int i = energies.size() - 1; i >= 0; i--)
        {
            Db rho = exp(-beta * energies[i]);
            Z += rho;
            E += rho * energies[i];
        }

        E /= Z;

        Db entropy = .0;

        for (int i = energies.size() - 1; i >= 0; i--)
        {
            Db rho = exp(-beta * energies[i]) / Z;
            entropy -= rho * log(rho);
        }

        ofs << beta << " " << log(Z) << " " << E << " " << entropy << endl;
    }

    ofs.close();
}

void EigenEnergyLargeN::CalcFluctuation(double beta)
{
    double maxT = 100.0;
    double delta = maxT / 100;

    string file = "FLs=" + ToString(s) + "M=" + ToString(M) + ".txt";
	ofstream ofs(file.c_str());

    double Z = .0;
    for (int i = allStates.size() - 1; i >= 0; i--)
    {
        double rho = exp(-beta * (allStates[i] + M)/sqrt(2));
        Z += rho;
    }

    ofs << beta << " " << Z << endl;
    for (double t = 0.0; t <= maxT + 1e-8; t += delta)
    {
        double ZtReal = .0;
        double ZtImag = .0;

        for (int i = allStates.size() - 1; i >= 0; i--)
        {
            double rho = exp(-beta * (allStates[i] + M) / sqrt(2));
            double angle = -(allStates[i] + M) / sqrt(2) * t;
            ZtReal += rho * cos(angle);
            ZtImag += rho * sin(angle);
        }

        ofs << t << " " << ZtReal << " " << ZtImag << endl;
    }

    ofs.close();
}

void EigenEnergyLargeN::SaveSingleEnergies(int bit, int buckets)
{
    string file = "EEs=" + ToString(s) + "M=" + ToString(bit) + "s";

    if (buckets > 0) file += "g.txt";
	else file += ".txt";

	ofstream ofs(file.c_str());

    if (buckets <= 0)
    {
        for (int i = 0; i < this->singleTraceEnergies[bit].size(); i++)
	    {
		    ofs << singleTraceEnergies[bit][i] << endl;
	    }
    }
    else
    {
        double range = 4.0 * s /tan(PI/(2 * bit));
		vector<pair<double, int> > res = BucketEnergies(this->singleTraceEnergies[bit], buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
			ofs << res[i].first << ' ' << res[i].second << endl;
		}
    }

    ofs.close();
}

vector<pair<double, int> > EigenEnergyLargeN::BucketEnergies(vector<double>& states, int buckets, double range)
{
    double delta = range * 2 / buckets;
	vector<pair<double, int> > res;
	for (int i = 0; i < buckets; i++)
	{
		res.push_back(make_pair(-range + delta * (i + .5), 0));
	}

	for (int i = 0; i < states.size(); i++)
	{
		int pos = (int)floor((states[i] + range) / delta);
		if (pos == -1) pos = 0;
		if (pos == buckets) pos = buckets - 1;
		res[pos].second++; 
	}

    return res;
}

void EigenEnergyLargeN::SaveEnergies(int buckets)
{
	if (s != 1)
	{
		cout << "Error: s = " << s << ". SaveEnergies only support for the s=1 case.";
		return;
	}

	string file = "EEs=" + ToString(s) + "M=" + ToString(M);
	if (buckets > 0) file += "g.txt";
	else file += ".txt";

	ofstream ofs(file.c_str());
	if (buckets <= 0)
	{
		for (int i = 0; i < this->allStates.size(); i++)
		{
			ofs << allStates[i] << endl;
		}
	}
	else
	{
		double range = 4.0 * s /tan(PI/(2 * M));
		vector<pair<double, int> > res = BucketEnergies(this->allStates, buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
			ofs << res[i].first << ' ' << res[i].second << endl;
		}
	}

	ofs.close();
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
    for (int i = 0; i < singleTraceEnergies[singleTraceBits].size(); i++)
    {
        BuildAllStatesRec(bitRemain - singleTraceBits, singleTraceBits - 1, 
            energy + singleTraceEnergies[singleTraceBits][i], deg + 1);
    }

    // pick two of more.
    for (int i = 2; i * singleTraceBits <= bitRemain; i++)
    {
        for (int j = 0; j < statesByBoth[i][singleTraceBits].size(); j++)
        {
            BuildAllStatesRec(bitRemain - singleTraceBits * i, singleTraceBits - 1, 
                energy + statesByBoth[i][singleTraceBits][j], deg + 1);
        }
    }
}

void EigenEnergyLargeN::BuildStatesByBoth()
{
    statesByBoth.reserve(M + 1);
    statesByBoth.push_back(vector<vector<double> >());
    statesByBoth.push_back(vector<vector<double> >());
    
    for (int i = 2; i <= M; i++)
    {
        statesByBoth.push_back(vector<vector<double> > (M/i + 1, vector<double>()));
        for (int bit = 1; bit * i <= M; bit++)
        {
            vector<double>& v = statesByBoth[i][bit];
			int size = 0;

			// go through first round to calculate number of states.
            for (int k = 0; k <= i; k += 2)
            {
				size += statesByFermion[k][bit].size() * statesByBoson[i - k][bit].size();
            }

			v.reserve(size);
			// second round to add states
			for (int k = 0; k <= i; k += 2)
            {
                if (statesByFermion[k][bit].size() > 0)
                {
                    MergeEnergy(statesByFermion[k][bit], statesByBoson[i - k][bit], v);
                }
            }

            //cout << "(" << i << ", " << bit << "): " << statesByBoth[i][bit] << endl;
        }
    }
}

void EigenEnergyLargeN::MergeEnergy(vector<double>& a, vector<double>& b, vector<double>& res)
{
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b.size(); j++)
        {
            res.push_back(Chop(a[i] + b[j]));
        }
    }
}

void EigenEnergyLargeN::CalculateByEigen(double invN, int buckets, bool calcEigenvector)
{
	if (s != 1)
	{
		cout << "Error: s = " << s << ". CalculateByEigen only support for calculating the s=1 case.";
		return;
	}

	string file = "s=" + ToString(s) + "M=" + ToString(M);
	if (abs(invN) > EPS)
	{
		file += "N=" + ToString((int)floor(1/invN + .5));
	}
	
	// file1 to save eigenvalues. 
	string file1 = "EE" + file;
    if (buckets > 0)
    {
        file1 += "g";
    }
    file1 += ".txt";

	// file2 to save eigenvectors.
	string file2 = "ES" + file + ".txt";
	
	// calculate eigenvalues and eigenvectors.
	H0Hamiltonian h0;
	Eigen::MatrixXcd mat = h0.ToMatrixN(M, Boson, invN);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(mat, calcEigenvector);
	Eigen::VectorXcd values = solver.eigenvalues();

	// sort them by eigenvalues.
	vector<int> pos(values.size(), 0);
	for (int i = 0; i < pos.size(); i++) pos[i] = i;
	EigenvalueLess comp(values);
	sort(pos.begin(), pos.end(), comp);

	ofstream ofs1(file1.c_str());

	double eps = 1e-6;

    if (buckets <= 0)
    {
	    for (int i = 0; i < values.rows(); i++)
	    {
		    ofs1 << Chop(values[pos[i]].real(), eps) << endl;
	    }
    }
    else
    {
        vector<double> energies(pos.size(), .0);
        for (int i = 0; i < pos.size(); i++)
        {
            energies[i] = Chop(values[pos[i]].real(), eps);
        }

        double range = max(abs(energies[1]), abs(energies.back()));
        vector<pair<double, int> > res = BucketEnergies(energies, buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
			ofs1 << res[i].first << ' ' << res[i].second << endl;
		}
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
