#include "Tests.h"
#include <iomanip>

Tests::Tests(int s)
{
    this->s = s;
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

void Tests::SingleTraceEnergies(int M)
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
    Tests test1(1);
    test1.SingleTraceEnergies(7);
}
