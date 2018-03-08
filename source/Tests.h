// This is the class for test code. 
//
#pragma once
#include "BitUtility.h"
#include "StateCounter.h"
#include "EigenEnergyLargeN.h"

class Tests
{
private:
    int s;

public:
	Tests(int s = 1);

    // test functions
    void SingleTraceEnergies(int M);

    void MultiTraceEigenEnergies(int maxM);

    // demo functions
    void DemoStateNumber();

    void DemoAverageEnergy();

    static void Run();
};

