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
    cout << "M\tSingle-Trace\t\tMulti-Trace" << endl;
    for (int i = 1; i * s <= StateCounter::MAX_BIT_TO_COUNT; i++)
	{
		cout << i << "\t" << inst->SingleTrace(i, s) << "\t\t" << inst->MultiTrace(i, s) << endl;
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
    Tests test2(2);
    Tests test3(3);

    test1.DemoStateNumber();
    test2.DemoStateNumber();
    test3.DemoStateNumber();

    test1.SingleTraceEnergies(7);
}
