#pragma warning(disable:4018)
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include "SingleTrace.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include "TraceState.h"
#include "StateGenerator.h"
#include "Hamiltonian.h"
#include "NormCalculator.h"
#include "ScriptGenerator.h"
#include "StateType.h"
#include "EigenSpectrum.h"
#include "StateCounter.h"
#include "Tests.h"

using namespace std;

#if WIN32
string scriptFolder = ".";
#else
string scriptFolder = ".";
#endif

void TestHamiltonian(string bits)
{
	TraceState state(bits);
	H0Hamiltonian ham;
	MixState re, im;
	ham.Apply(state, re, im);
	cout << "(" << ham << ")" << state << "|0> = ";
	cout << re.ToString();
	cout <<" + i{" << im.ToString() << "}" << endl;
}

void TestZeroHamiltonian()
{
	ZeroHamiltonian ham;
	MixState re, im;
	TraceState state("aaaaaa");
	ham.Apply(state, re, im);
	cout << ham  << state << "|0> = ";
	cout << re.ToString();
	cout << " + i {" << im.ToString() << "}" << endl;

	for (int i = 0; i < ham.RealOpSize(); i++)
	{
		HamOperator* op = ham.RealOp(i);
		MixState ms;
		op->ApplyOn(state, ms);
		cout << *op << state << "|0> = ";
		cout << ms.ToString() << endl; 
	}
}

void TestHamOperator()
{
	HamOperator *op;
	op = new BitNumberHamOperator();

	TraceState state("a,aabbb");
	MixState res;
	op->ApplyOn(state, res);
	cout << "H" << state << "|0> = ";
	cout << res.ToString() << endl;
	delete op;
}

void TestTraceState()
{
	SingleTrace single("aaaab");
	cout << single << endl;
	
	TraceState ts("aaa,ab");
	cout << ts << endl;
}

void GenerateStates(bool output)
{
	StateGenerator generator;
	/*for (int i = 1; i <= StateCounter::MAX_BIT_TO_COUNT; i++)
	{
		if (output)
		{
			cout << i << "\t" << generator.SingleStateNumber(i) << "\t" << generator.BosonNumber(i) << endl;
		}
	}*/

	generator.GenerateAllStates();
	generator.InitStateCollection(StateCollection::Inst());
}

void OutputHamiltonianMatrix(int bits)
{
	Hamiltonian ham;
	vector<vector<Coefficient> > rem, imm;
	ham.Matrix(bits, Boson, rem, imm);
	cout << bits << " bits Hamiltonian: " <<endl;
	for (int i = 0; i < rem.size(); i++)
	{
		for (int j = 0; j < rem[i].size(); j++)
		{
			if (!imm[i][j].IsZero())
			{
				cout << imm[i][j] << "i ";
			}
			else
			{
				cout << rem[i][j] << " ";
			}
			
		}
		cout << endl;
	}
}

void OutputHamToLatex()
{
	int arr[6] = {2, 3, 4, 5, 6, 7};
	double scale[6] = {1.0, 1.0, 1.0, 1, 0.57, 0.27};
	vector<int> bits(arr, arr + 6);
	StateType type = Boson;
	ScriptGenerator sc(scriptFolder, type);
	sc.OutputHamToLaTeX(bits, vector<double>(scale, scale + 6), "E:\\Dropbox\\proj\\stringbit\\hamiltonian.tex");
}

void CalculateNorm(int bits, bool output)
{
	BruteForceCalculator calc;
	StateCollection* sc = StateCollection::Inst();
	Stopwatch watch;
	
	watch.Start();
	cout << "Calculating norm matrix for " << bits << " bits states..." << endl;
	for (int i = 0; i < sc->StateNumber(bits); i++)
	{
		for (int j = 0; j < sc->StateNumber(bits); j++)
		{
			Polynomial res = calc.Calculate(sc->GetBosonState(bits, i), sc->GetBosonState(bits, j));
			if (output)
			{
				cout << res << "\t";
			}
		}

		if (output)
		{
			cout << endl;
		}
	}
	cout << "It takes " << watch.Stop() << " seconds." << endl;
}

void TestNormCalculator()
{
	int bits = 5;
	BruteForceCalculator calc;
	StateCollection* sc = StateCollection::Inst();

	const TraceState& a = sc->GetBosonState(bits, 7);
	const TraceState& b = sc->GetBosonState(bits, 8);

	cout << "(" << a << "," << b << ")=";
	cout << calc.Calculate(a, b) << endl;;
}

void GenerateHamDataFile(StateType type)
{
	ScriptGenerator sc(scriptFolder, type);
	
	H0Hamiltonian h0;
	DeltaHamiltonian delta;
	Stopwatch watch;
	watch.Start();
	for (int i = 3; i <= 9; i++)
	{
		sc.OutputHamToDataFile(i, h0, true);
		sc.OutputHamToDataFile(i, delta);
	}

	cout << "Time: " << watch.Stop() << " seconds." << endl;
}

void GenerateHamMatlab(StateType type)
{	
	ScriptGenerator sc(scriptFolder, type);
	H0Hamiltonian h0;
	DeltaHamiltonian delta;

	Stopwatch watch;
	watch.Start();
	//for (int i = 3; i <= StateGenerator::MAX_BIT_TO_GENERATE; i++)
	for (int i = 3; i <= 9; i++)
	{
		sc.OutputHamToMatlab(i, h0);
		sc.OutputHamToMatlab(i, delta);
	}

	cout << "Time: " << watch.Stop() << " seconds." << endl;
}

void GenerateNormMatlab(StateType type)
{
	Stopwatch watch;
	watch.Start();
	ScriptGenerator sc(scriptFolder, type);
	//for (int i = 3; i <= StateGenerator::MAX_BIT_TO_GENERATE; i++)
	for (int i = 3; i <= 11; i++)
	{
		//sc.OutputNormToMatlab(i);
	}

	cout << "Time: " << watch.Stop() << " seconds." << endl;
}

void GenerateNormDataFile(int minBits, int maxBits, StateType type)
{
	Stopwatch watch;
	watch.Start();
	ScriptGenerator sc(scriptFolder, type);
	int fNumber;
	for (int i = minBits; i <= maxBits; i++)
	{
		cout << "Processing bits=" << i << "..." << endl;
		if (type == Boson)
		{
			fNumber = 0;
		}
		else
		{
			fNumber = 1;
		}

		while (fNumber < i)
		{
			sc.OutputNormToDataFile(i, fNumber);
			fNumber += 2;
		}
	}

	cout << "Time: " << watch.Stop() << " seconds." << endl;
}

void GenerateAllHamDataFile()
{
	GenerateHamDataFile(Boson);
	GenerateHamDataFile(Fermion);
}

void GenerateAllHamMatlab()
{
	GenerateHamMatlab(Boson);
	GenerateHamMatlab(Fermion);
}

void GenerateLaTeX()
{
	ScriptGenerator sg1(scriptFolder, Boson);
	ScriptGenerator sg2(scriptFolder, Fermion);
	//string filename = "E:\\Dropbox\\proj\\stringbit\\boson-states-1-7.tex";
	string filename = "boson-states-1-7.tex";
	sg1.OutputStateToLaTeX(1, 7, filename, true);
	filename = "fermion-states-1-7.tex";
	sg2.OutputStateToLaTeX(1, 7, filename, true);
	OutputHamToLatex();
	sg1.OutputNormToLaTeX(2, 7, "E:\\Dropbox\\proj\\stringbit\\norms-2-7.tex");
}

void OutputStateStructure()
{
	ScriptGenerator sg1(scriptFolder, Boson);
	ScriptGenerator sg2(scriptFolder, Fermion);
	sg1.OutputStateStructure();
	sg2.OutputStateStructure(); 
}

void CalculateAllEigenvalues(int spin, int bits, int buckets)
{
    cout << "Calculating large N energies for spin=" << spin << ", bits=" << bits << " ..." << endl;
	EigenEnergyLargeN ee(bits, spin);
    ee.CalculateByDynamics();
    ee.SaveEnergies(buckets);
    //ee.CalcFluctuation(1.5);
}

void SaveSingleEnergies(int s, int M, int buckets = 0)
{
    cout << "Calculating large N single-trace energies for s= " << s << " bits=" << M << " ..." << endl;
	EigenEnergyLargeN ee(M, s);

    ee.CalcAllSingleTraceEnergies();
    ee.SaveSingleEnergies(M, buckets);
}

void CalculateThermodynamics(int bits, double T0, double maxB, int steps = 100)
{
    cout << "Calculating thermodynamics for bits=" << bits << " ..." << endl;
	EigenEnergyLargeN ee(bits);
    ee.CalculateByDynamics();
    ee.CalculateThermo(T0, maxB, steps);
}

void CalculateAllEigenvaluesByEigen(int bits, int N, int buckets)
{
    cout << "Calculating large N energies for bits=" << bits << " ..." << endl;
	EigenEnergyFiniteN ee(bits, (double)N);
    ee.SaveEnergies(buckets, false);
}

void CalculateThermodynamicsForN(int N, double T0, double maxB, string datafolder, int steps = 1000)
{
    cout << "Calculating thermodynamics for N=" << N << " ..." << endl;
	EigenEnergyFiniteN ee(N, (double)N);
    ee.SaveThermoData(T0, maxB, datafolder, steps);
}

int main(int argc, char* argv[])
{
	GenerateStates(false);    
	//TestTraceState();
	//TestHamOperator();
	//TestHamiltonian("bbb");
	//TestHamiltonian("aab");
	//TestZeroHamiltonian();
	//CalculateNorm(9, false);
	//TestNormCalculator();
	//OutputHamiltonianMatrix(3);
	//GenerateLaTeX();

	//GenerateAllHamMatlab();
	//OutputStateStructure();
	//GenerateAllHamDataFile();
	//GenerateNormDataFile(8, 10, Boson);

    Stopwatch watch;
	watch.Start();

    cout << "Start computation at " << watch.Now() << endl;
	string cmd = "";
	if (argc > 1) cmd = argv[1];

    if (cmd == "-t")
    {
        Tests::Run();
    }
    else if (cmd == "-sn")
    {
        if (argc < 3)
		{
			cout << "need one more integer type  parameter." << endl;
            cout << "Usage: sbth -sn s" << endl;
            cout << "Example: sbth -sn 3" << endl;
			return -1;
		}

        Tests test(atoi(argv[2]));
        test.DemoStateNumber();
    }
    else if (cmd == "-e")
	{
		if (argc < 4)
		{
			cout << "need one more integer type  parameter." << endl;
            cout << "Usage: sbth -e spin bits [buckets=101]" << endl;
            cout << "Example: sbth -e 1 19 50" << endl;
			return -1;
		}

        int bucket = 101;

        if (argc >= 5)
        {
            bucket = atoi(argv[4]);
        }

	    CalculateAllEigenvalues(atoi(argv[2]), atoi(argv[3]), bucket);
	}
    else if (cmd == "-es")
    {
        if (argc < 4)
		{
			cout << "need more integer  parameters: s M buckets" << endl;
			return -1;
		}

        SaveSingleEnergies(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    }
    else if (cmd == "-ebe")
    {
        // calculate all eigen energies using the Eigen library.
        if (argc < 5)
        {
            cout << "please input 3 integer parameters: M, N, and buckets." << endl;
            return -1;
        }

        CalculateAllEigenvaluesByEigen(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    }
    else if (cmd == "-th")
    {
        // calculate thermodunamics.
        if (argc < 5)
        {
            cout << "please input 3 integer parameters: M, T0, and maxBeta." << endl;
            return -1;
        }
        CalculateThermodynamics(atoi(argv[2]), atof(argv[3]), atof(argv[4]));
    }
    else if (cmd == "-thn")
    {
        // calculate thermodunamics for finite N.
        if (argc < 6)
        {
            cout << "please input 3 integer parameters: N, T0, maxBeta, data folder" << endl;
            return -1;
        }

        CalculateThermodynamicsForN(atoi(argv[2]), atof(argv[3]), atof(argv[4]), argv[5]);
    }
    
    cout << "Finish computation at " << watch.Now() << endl;
    cout << "Elapsed time: " << watch.Stop() << " seconds." << endl;
	return 0;
}
