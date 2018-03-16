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

void EigenEnergyLargeN::CalcAllSingleTraceEnergies()
{
    // if already calculated, return.
    if (singleTraceEnergies.size() > 0)
    {
        return;
    }

    singleTraceEnergies.reserve(M + 1);
    singleTraceEnergies.push_back(vector<TE>(1, 0.0));
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

        sort(v.begin(), v.end());

#ifndef DEG_ENERGY_TYPE
        singleTraceEnergies.push_back(v);
        //cout << "bit = " << bit << ", state number = " << v.size() << ", expected = " << StateCollection::Inst()->SingleTraceStateNumber(bit, s) << endl;
#else
        // group energies.
        singleTraceEnergies.push_back(vector<TE>(0));
        singleTraceEnergies.back().reserve(v.size());
        double prev = -100000000.0;
        for (int i = 0; i < v.size(); i++)
        {
            if (v[i] - prev > EPS)
            {
                singleTraceEnergies.back().push_back(v[i]);
                prev = v[i];
            }
            else
            {
                singleTraceEnergies.back().back().Deg++;
            }
        }

        for (int i = 0; i < singleTraceEnergies.back().size(); i++)
        {
//            cout << "deg before shift: "  << singleTraceEnergies.back()[i].Deg;
            singleTraceEnergies.back()[i].Deg <<= (s-1);
//            cout << " after shift: "  << singleTraceEnergies.back()[i].Deg << endl;
        }

//        cout << "Single-trace: bit = " << bit << ", state number = " << TotalSize(singleTraceEnergies.back()) << ", expected = " << StateCollection::Inst()->SingleTraceStateNumber(bit, s) << endl;
#endif
    }
}

void EigenEnergyLargeN::BuildBosonicMultiTraceEnergies(int b, int n)
{
    vector<TE>& res = statesByBoson[n][b];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b, s);
    res.reserve(BinomialCoefficient(size + n - 1, n));
    BuildBosonicMultiTraceEnergiesRec(b, singleTraceEnergies[b].size() - 1, n, .0, res);
//    cout << "states built by " << n << " " << b << "-bit bosonic single trace: " << TotalSize(res) << endl;
}

void EigenEnergyLargeN::BuildBosonicMultiTraceEnergiesRec(int b, int index, int n, TE energy, vector<TE>& res)
{
    if (n == 0)
    {
        res.push_back(energy);
        return;
    }

    if (index == 0)
    {
        res.push_back(energy.Times(singleTraceEnergies[b][0].PowerBoson(n)));
        return;
    }

    BuildBosonicMultiTraceEnergiesRec(b, index - 1, n, energy, res);

    for (int i = 1; i <= n; i++)
    {
        BuildBosonicMultiTraceEnergiesRec(b, index - 1, n - i, energy.Times(singleTraceEnergies[b][index].PowerBoson(i)), res);
    }
}

void EigenEnergyLargeN::BuildFermionicMultiTraceEnergies(int b, int n)
{
    vector<TE>& res = statesByFermion[n][b];
    int size = StateCollection::Inst()->SingleTraceStateNumber(b, s);
    res.reserve(BinomialCoefficient(size, n));
    BuildFermionicMultiTraceEnergiesRec(b, singleTraceEnergies[b].size() - 1, n, .0, res);
//    cout << "states built by " << n << " " << b << "-bit fermionic single trace: " << TotalSize(res) << endl;
}

void EigenEnergyLargeN::BuildFermionicMultiTraceEnergiesRec(int b, int index, int n, TE energy, vector<TE>& res)
{
    if (n == 0)
    {
        res.push_back(energy);
        return;
    }

    if (index == 0)
    {
        if (singleTraceEnergies[b][0].Deg < n) return;
        res.push_back(energy.Times(singleTraceEnergies[b][0].PowerFermion(n)));
        return;
    }

    BuildFermionicMultiTraceEnergiesRec(b, index - 1, n, energy, res);

    for (int i = 1; i <= n && i <= singleTraceEnergies[b][index].Deg; i++)
    {
        BuildFermionicMultiTraceEnergiesRec(b, index - 1, n - i, energy.Times(singleTraceEnergies[b][index].PowerFermion(i)), res);
    }
}

void EigenEnergyLargeN::CalculateByDynamics()
{
        
    // first calculate all single trace state energies.
    CalcAllSingleTraceEnergies();

    // initialize statesByBoson and statesByFermion
    statesByBoson.reserve(M + 1);
    statesByFermion.reserve(M + 1);

    // the (0, i) elements are vacuum states.
    statesByBoson.push_back(vector<vector<TE> >(M + 1, vector<TE>(1, 0.0)));
    statesByFermion.push_back(vector<vector<TE> >(M + 1, vector<TE>(1, 0.0)));

    // the first element of statesByBoson is simply singleTraceEnergies
    statesByBoson.push_back(singleTraceEnergies);
    // for statesByFermion, we let its first element be empty, because we do
    // not need it in the following calculation.
    statesByFermion.push_back(vector<vector<TE> >());

    for (int i = 2; i <= M; i++)
    {
        statesByBoson.push_back(vector<vector<TE> >(M/i + 1, vector<TE>()));
        statesByFermion.push_back(vector<vector<TE> >(M/i + 1, vector<TE>()));
        for (int bit = 1; bit * i <= M; bit++)
        {
            //cout << "(" << i << ", " << bit << ") single trace states #:";
            //cout << StateCollection::Inst()->SingleTraceStateNumber(bit);

            BuildBosonicMultiTraceEnergies(bit, i);
            //cout << " statesByBoson: " << statesByBoson[i][bit];
            if (i <= StateCollection::Inst()->SingleTraceStateNumber(bit, s))
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

int EigenEnergyLargeN::EstimateStatesNumber(int s, int M)
{
    // The following fits are obtained by LibreOffice Calc
    switch (s)
    {
        case 1:
            return (int)floor(pow(2, 0.86 * M - 0.5) + 1);
        case 2:
            return (int)floor(pow(2, 1.26 * M - 1) + 1);
        case 3:
            return (int)floor(pow(2, 1.5071 * M - 1.49) + 1);
        case 4:
            return (int)floor(pow(2, 1.662 * M -1.431) + 1);
        case 5:
            return (int)floor(pow(2, 1.765 * M -1.521) + 1);
        default:
            return (int)floor(pow(2, 2.1 * M - 1.5 * M /s -1.5));
    }
}

void EigenEnergyLargeN::BuildAllStates()
{
    // To avoid allocating too much memory, we estimate the size
    // of the allStates vector. The estimation was obtained by fitting 
    // the run results data. 
    int cap = EstimateStatesNumber(s, M);
    allStates.reserve(cap);
    BuildAllStatesRec(M, M, .0, 0);
    sort(allStates.begin(), allStates.end());
    int size = allStates.size();
    Collect(allStates);

    if (TotalSize(allStates) != StateCollection::Inst()->StateNumber(M, s))
    {
        cout << "Error BuildAllStates: Total states=";
        cout << TotalSize(allStates) << ", expected=" << StateCollection::Inst()->StateNumber(M, s) << endl;
    }

    if (cap != allStates.capacity())
    {    
        cout << "expected capacity: " << cap << ", actual capacity: " << allStates.capacity() << endl; 
        cout << "allState vector size, before collect: " << size << ", after collect: " << allStates.size() << endl;
        cout << "Compress ratio: " << size / (double) TotalSize(allStates) << endl;
    }
}

void EigenEnergyLargeN::BuildAllStatesRec(int bitRemain, int singleTraceBits, TE energy, int deg)
{
    if (bitRemain == 0 || singleTraceBits == 1)
    {
        int n = (1 << (deg - 1));
        if (bitRemain > 0) n = (1 << deg);
        TE toAdd;
        if (bitRemain == 0)
        {
            toAdd = energy;
        }
        else if (bitRemain == 1)
        {
            toAdd = energy.Times(singleTraceEnergies[1][0]);
        }
        else
        {
            toAdd = energy.Times(statesByBoth[bitRemain][1][0]);
        }

        allStates.push_back(TE(toAdd.E, toAdd.Deg * n));
        return;
    }

    if (singleTraceBits == 0)
    {
        cout << "BuildAllStatesRec should not reach here!!!!!!!!!!!!!!!!!" << endl;
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
            energy.Times(singleTraceEnergies[singleTraceBits][i]), deg + 1);
    }

    // pick two of more.
    for (int i = 2; i * singleTraceBits <= bitRemain; i++)
    {
        for (int j = 0; j < statesByBoth[i][singleTraceBits].size(); j++)
        {
            BuildAllStatesRec(bitRemain - singleTraceBits * i, singleTraceBits - 1, 
                energy.Times(statesByBoth[i][singleTraceBits][j]), deg + 1);
        }
    }
}

void EigenEnergyLargeN::BuildStatesByBoth()
{
    statesByBoth.reserve(M + 1);
    statesByBoth.push_back(vector<vector<TE> >());
    statesByBoth.push_back(vector<vector<TE> >());
    
    for (int i = 2; i <= M; i++)
    {
        statesByBoth.push_back(vector<vector<TE> > (M/i + 1, vector<TE>()));
        for (int bit = 1; bit * i <= M; bit++)
        {
            vector<TE>& v = statesByBoth[i][bit];
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

            sort(v.begin(), v.end());
            Collect(v);
//            cout << "statesByBoth(" << i << ", " << bit << "): " << TotalSize(statesByBoth[i][bit]) << endl;
        }
    }
}

i64 TotalSize(vector<TE>& v)
{
    i64 tot = 0;
    for (int i = 0; i < v.size(); i++)
    {
        tot += v[i].Deg;
    }
    return tot;
}

void Collect(vector<TE>& v)
{
    int saved = 0;
    for (int i = 1; i < v.size(); i++)
    {
        if (v[i].E - v[saved].E < EPS)
        {
            v[saved].Deg += v[i].Deg;
        }
        else
        {
            v[++saved] = v[i];
        }
    }

    v.resize(saved + 1);
}

void MergeEnergy(vector<TE>& a, vector<TE>& b, vector<TE>& res)
{
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b.size(); j++)
        {
            res.push_back(a[i].Times(b[j]));
        }
    }
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
            Db rho = exp(-beta * (T0 * allStates[i].E + M)/sqrt(2)) * allStates[i].Deg;
            Z += rho;
            E += rho * (T0 * allStates[i].E + M)/sqrt(2);
        }

        E /= Z;

        Db entropy = .0;

        for (int i = allStates.size() - 1; i >= 0; i--)
        {
            Db rho = exp(-beta * (T0 * allStates[i].E + M)/sqrt(2)) / Z;
            entropy -= rho * log(rho) * allStates[i].Deg;
        }

        ofs << beta << " " << log(Z) << " " << E << " " << entropy << endl;
    }

    ofs.close();
}

void EigenEnergyFiniteN::CalculateThermoForN(double T0, double maxB, string datafolder, int steps)
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
            string input = datafolder + "/EEs=" + ToString(s) + "M=" + ToString(bit) + "N=" + ToString(N) + ".txt";
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
            Eigen::MatrixXcd mat = h0.ToMatrixN(bit, Boson, 1.0/N);
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
    string file = "THs=" + ToString(s) + "N=" + ToString(N) + "T0=" + ToString(T0) + ".txt";
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
        double rho = exp(-beta * (allStates[i].E + M)/sqrt(2));
        Z += rho * allStates[i].Deg;
    }

    ofs << beta << " " << Z << endl;
    for (double t = 0.0; t <= maxT + 1e-8; t += delta)
    {
        double ZtReal = .0;
        double ZtImag = .0;

        for (int i = allStates.size() - 1; i >= 0; i--)
        {
            double rho = exp(-beta * (allStates[i].E + M) / sqrt(2));
            double angle = -(allStates[i].E + M) / sqrt(2) * t;
            ZtReal += rho * cos(angle) * allStates[i].Deg;
            ZtImag += rho * sin(angle) * allStates[i].Deg;
        }

        ofs << t << " " << ZtReal << " " << ZtImag << endl;
    }

    ofs.close();
}

vector<DegEnergy> BucketEnergies(vector<TE>& states, int buckets, double range)
{
    double delta = range * 2 / buckets;
	vector<DegEnergy> res;
	for (int i = 0; i < buckets; i++)
	{
		res.push_back(DegEnergy(-range + delta * (i + .5), 0));
	}

	for (int i = 0; i < states.size(); i++)
	{
		int pos = (int)floor((states[i].E + range) / delta);
		if (pos == -1) pos = 0;
		if (pos == buckets) pos = buckets - 1;
		res[pos].Deg += states[i].Deg;
	}

    return res;
}

vector<DegEnergy> BucketEnergies(vector<double>& states, int buckets, double range)
{
    double delta = range * 2 / buckets;
	vector<DegEnergy> res;
	for (int i = 0; i < buckets; i++)
	{
		res.push_back(DegEnergy(-range + delta * (i + .5), 0));
	}

	for (int i = 0; i < states.size(); i++)
	{
		int pos = (int)floor((states[i] + range) / delta);
		if (pos == -1) pos = 0;
		if (pos == buckets) pos = buckets - 1;
		res[pos].Deg++;
	}

    return res;
}

void EigenEnergyLargeN::SaveSingleEnergies(int bit, int buckets)
{
    string file = "EEs=" + ToString(s) + "M=" + ToString(bit) + "s";

    if (buckets > 0) file += "g.txt";
	else file += ".txt";

    int tot = 0;

	ofstream ofs(file.c_str());

    if (buckets <= 0)
    {
        for (int i = 0; i < this->singleTraceEnergies[bit].size(); i++)
	    {
            tot += singleTraceEnergies[bit][i].Deg;
		    ofs << singleTraceEnergies[bit][i].E << ' ' << singleTraceEnergies[bit][i].Deg << endl;
	    }
    }
    else
    {
        double range = 4.0 * s /tan(PI/(2 * bit));
		vector<DegEnergy> res = BucketEnergies(this->singleTraceEnergies[bit], buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
            tot += res[i].Deg;
			ofs << res[i].E << ' ' << res[i].Deg << endl;
		}
    }

    cout << "Total sates: " << tot << ", expected: " << StateCollection::Inst()->SingleTraceStateNumber(bit, s) << endl;

    ofs.close();
}

void EigenEnergyLargeN::SaveEnergies(int buckets)
{
	string file = "EEs=" + ToString(s) + "M=" + ToString(M);
	if (buckets > 0) file += "g.txt";
	else file += ".txt";

	ofstream ofs(file.c_str());
	if (buckets <= 0)
	{
		for (int i = 0; i < this->allStates.size(); i++)
		{
			ofs << Chop(allStates[i].E) << ' ' << allStates[i].Deg << endl;
		}
	}
	else
	{
		double range = 4.0 * s /tan(PI/(2 * M));
		vector<DegEnergy> res = BucketEnergies(this->allStates, buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
			ofs << res[i].E << ' ' << res[i].Deg << endl;
		}
	}

	ofs.close();
}


void EigenEnergyFiniteN::CalculateByEigen(int buckets, bool calcEigenvector)
{
	if (s != 1)
	{
		cout << "Error: s = " << s << ". CalculateByEigen only support for calculating the s=1 case.";
		return;
	}

	string file = "s=" + ToString(s) + "M=" + ToString(M) + "N=" + ToString(N);;
	
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
	Eigen::MatrixXcd mat = h0.ToMatrixN(M, Boson, 1.0/N);
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
        vector<DegEnergy> res = BucketEnergies(energies, buckets, range);
		for (int i = 0; i < res.size(); i++)
		{
			ofs1 << res[i].E << ' ' << res[i].Deg << endl;
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

ostream& operator<<(ostream& os, const DegEnergy& de)
{
    os << '(' << de.E << ',' << de.Deg << ')';
    return os;
}
