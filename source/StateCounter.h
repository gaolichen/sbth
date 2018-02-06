#pragma once
#include <vector>
#include "BitUtility.h"

using namespace std;

// class implements all counting functions.
class StateCounter
{
private:
    vector<vector<vector<snum> > > stateNumbers;
    vector<vector<snum> > oddOnlyStateNumbers;
    vector<vector<snum> > evenTraceNumbers;
	vector<snum> singleTraceNumbers;
    vector<double> aveESingle;
    vector<vector<long double> > aveE;
	vector<snum> noHalfMode;

	StateCounter(void);
	static StateCounter* inst;

	void InitNoHalfMode();

    snum StateNumbers(int bit, int remain, int b);
    snum OddOnlyStateNumbers(int bit, int remain, int parity);
    snum EvenTraceNumber(int bit, int remain, int parity);
    long double AverageEnergy(int bit, int remain, int parity);

    void InitSingleTraceNumber();
    void InitAverageEnergy();
public:

    const static int MAX_BIT_TO_COUNT = 62;
	static StateCounter* Inst();

	// function to initialize all data of the class.
	void Init();

	// Returns the number of partitions of M without M/2.
	// A valid partition is M = n_1 + n_2 + ..., where 1 <= n_1 < n_2 < ... < M and n_k != M/2 for all k.
	snum NoHalfMode(int M);

    snum SingleTrace(int M);

    snum MultiTrace(int M);

    double AverageEnergy(int bit);

	~StateCounter(void);
};
