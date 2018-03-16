#include "Tests.h"
#include <iomanip>
#include <fstream>

Tests::Tests(int s)
{
    this->s = s;
    this->datapath = "/users/gaolichen/gitroot/sbth/data/";
}

void Tests::DemoStateNumber()
{
    StateCounter* inst = StateCounter::Inst();
    cout << "Number of states of spin-" << s << " system:" << endl;
    cout.width(6);
    cout << right << "M" << " ";
    cout.width(22);
    cout << right << "SingleTrace" << " ";
    cout.width(15);
    cout << right << "ration1" << " ";
    cout.width(22);
    cout << right << "MultiTrace" << " ";
    cout.width(15);
    cout << right << "ratio2" << endl;
    for (int i = 1; i * s <= StateCounter::MAX_BIT_TO_COUNT; i++)
	{
        cout.width(6);
		cout << right << i << " ";
        cout.width(22);
        cout << right << inst->SingleTrace(i, s) << " ";
        cout.width(15);
        cout << setprecision(10) << right << inst->SingleTrace(i, s) / (long double)((i64)1<<(s * i)) * i << " ";
        cout.width(22);
        cout << right << inst->MultiTrace(i, s) << " ";
        cout.width(15);
        cout << setprecision(10) << right << inst->MultiTrace(i, s) / (long double)((i64)1<<(s * i)) << endl;
	}
}

void Tests::SingleTraceEnergies_StateNumber(int M)
{
    int maxSpin = 4;

    StateCounter* inst = StateCounter::Inst();
    for (int i = 1; i <= maxSpin; i++)
    {
        EigenEnergyLargeN calc(M, i);
        calc.CalcAllSingleTraceEnergies();

        for (int bit = 1; bit <= M; bit++)
        {
            if (calc.SingleEnergySize(bit) != inst->SingleTrace(bit, i))
            {
                cout << "SingleTraceEnergies failed for s=" << i << ", M=" << bit << ": ";
                cout << "\t returned " << calc.SingleEnergySize(bit) << ", expected " << inst->SingleTrace(bit, i) << endl;
                return;
            }
        }
    }

    cout << "Test::SingleTraceEnergies passed." << endl;
}

void Tests::MultiTraceEigenEnergies_StateNumber(int maxM)
{
    StateCounter* inst = StateCounter::Inst();
    for (int bit = 1; bit <= maxM; bit++)
    {
        EigenEnergyLargeN calc(bit, s);
        calc.CalculateByDynamics();

        if (calc.MultiEnergySize() != inst->MultiTrace(bit, s))
        {
            cout << "MultiTraceEigenEnergies failed for s=" << s << ", M=" << bit << ": ";
            cout << "\t returned " << calc.MultiEnergySize() << ", expected " << inst->MultiTrace(bit, s) << endl;
            return;
        }
    }

    cout << "Test::MultiTraceEigenEnergies passed for s=" << s << ", maxM=" << maxM << endl;
}

vector<TE> Tests::ReadEnergyFile(string file)
{
    ifstream ifs(file.c_str());
    if (!ifs.good())
    {
        cout << "Tests::ReadEnergyFile: cannot open " << file << endl;
        ifs.close();
        return vector<TE>();
    }

    vector<TE> ret;

    double e1, e2;
    ifs >> e1 >> e2;

    int deg = (int)floor(e2 + 0.5);

    // if e2 is integer, then the format of each line is [E degenerate]
    if (abs(deg-e2) < EPS && deg > 0)
    {
        ret.push_back(TE(e1, deg));
        while (ifs >> e1 >> deg)
        {
            ret.push_back(TE(e1,deg));
        }
    }
    else
    {
        // e2 is also an energy eigenvalue
        ret.push_back(e1);
        if (e2 - e1 < EPS)
        {
            ret.back().Deg++;
        }
        else
        {
            ret.push_back(e2);
        }

        while (ifs >> e1)
        {
            if (e1 - ret.back().E < EPS)
            {
                ret.back().Deg++;
            }
            else
            {
                ret.push_back(e1);
            }
        }
    }

    ifs.close();
    return ret;
}

void Tests::MultiTraceEigenEnergies_States(int maxM)
{
    for (int bit = 1; bit <= maxM; bit++)
    {
        string file = this->datapath + "EEs=" + ToString(s) + "M=" + ToString(bit) + ".txt";

        vector<TE> expected = ReadEnergyFile(file);

        // no existing results, cannot compare, so ignore it.
        if (expected.size() == 0) continue;

        EigenEnergyLargeN calc(bit, s);
        calc.CalculateByDynamics();

        const vector<TE>& ret = calc.AllStates();
        Assert(EnergyListEqual(expected, ret), "MultiTraceEigenEnergies_States, bit=" + ToString(bit) + " expected != ret");
    }

    cout << "Test::MultiTraceEigenEnergies_States passed for s=" << s << ", maxM=" << maxM << endl;
}

void Tests::SingleTraceEnergies_States(int maxM)
{
    EigenEnergyLargeN calc(maxM, s);
    calc.CalcAllSingleTraceEnergies();

    for (int bit = 1; bit <= maxM; bit++)
    {
        string file = this->datapath + "EEs=" + ToString(s) + "M=" + ToString(bit) + "s.txt";
        vector<TE> expected = ReadEnergyFile(file);
        // no existing results, cannot compare, so ignore it.
        if (expected.size() == 0) continue;

        const vector<TE>& ret = calc.SingleEnergies(bit);
        //cout << "ret=" << ret << "\n expected=" << expected << endl;
        //if (EnergyListEqual(expected, ret)) cout << "equal" << endl;
        //else cout << "not equal" << endl;

        Assert(EnergyListEqual(expected, ret), "SingleTraceEnergies_States: bit=" + ToString(bit) + " expected != ret");
    }

    cout << "Test::SingleTraceEnergies_States passed for s=" << s << ", maxM=" << maxM << endl;
}

void Tests::DemoAverageEnergy()
{
    StateCounter* inst = StateCounter::Inst();

    for (int i = 2; i <= StateCounter::MAX_BIT_TO_COUNT; i+=2)
    {
        cout << i << ": " << inst->NoHalfMode(i) << endl;
    }

    cout << "average energy of each bit:" << endl;    
    for (int i = 1; i <= StateCounter::MAX_BIT_TO_COUNT; i++)
    {
        cout << i << ": " << setprecision(10) << inst->AverageEnergy(i) << endl;
    }
}

void Tests::Run()
{
    cout << ToString(2.5) << endl;
    Tests test1(1);
    test1.SingleTraceEnergies_StateNumber(7);
    test1.MultiTraceEigenEnergies_StateNumber(11);

    Tests test2(2);
    test2.MultiTraceEigenEnergies_StateNumber(8);

    Tests test3(3);
    test3.MultiTraceEigenEnergies_StateNumber(6);

    test1.SingleTraceEnergies_States(13);
    test1.MultiTraceEigenEnergies_States(13);
}
