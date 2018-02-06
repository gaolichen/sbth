#pragma once
#include <vector>
#include "BitUtility.h"
#include "StateGenerator.h"

using namespace std;

// class implements all counting functions.
// TODO: move all counting functions to the class.
class StateCounter
{
private:
	vector<snum> noHalfMode;

	StateCounter(void);
	static StateCounter* inst;

	void StateCounter::InitNoHalfMode();
public:
	static StateCounter* Inst();

	// function to initialize all data of the class.
	void Init();

	// Returns the number of partitions of M without M/2.
	// A valid partition is M = n_1 + n_2 + ..., where 1 <= n_1 < n_2 < ... < M and n_k != M/2 for all k.
	snum NoHalfMode(int M);

	~StateCounter(void);
};

