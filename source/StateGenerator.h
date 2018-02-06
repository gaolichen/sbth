#pragma once
#include "SingleTrace.h"
#include "TraceState.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include "StateCounter.h"
#include <vector>

using namespace std;

class StateGenerator
{
private:
    StateCounter* counter;
	bool* myFlags;
	bool visited[30][30][2];

	vector<vector<SingleTrace> > singleFermions;
	vector<vector<SingleTrace> > singleBosons;
	vector<vector<vector<TraceState> > > fermions;
	vector<vector<vector<TraceState> > > bosons;

	void FindSingleStates(int n);
	void GeneratSingleStates();
	static void DoPickFermion(vector<SingleTrace>& allstates, int index, int remain, vector<SingleTrace>& curr, vector<TraceState>& ret);
	static void DoPickBoson(vector<SingleTrace>& allstates, int index, int remain, vector<SingleTrace>& curr, vector<TraceState>& ret);
	static vector<TraceState> PickFermionFromSingleState(vector<SingleTrace> &states, int number);
	static vector<TraceState> PickBosonFromSingleState(vector<SingleTrace> &states, int number);
	static TraceState CombineStates(TraceState& a, TraceState& b, TraceState& c);

	void BuildSingleOperatorStates(int remBit, int currBits, vector<int>& res);
public:
	const static int MAX_BIT_TO_GENERATE = 11;

	StateGenerator();
	~StateGenerator();
	
	void GenerateAllStates();
    
	snum BosonNumber(int n);
	snum FermionNumber(int n);
	snum SingleStateNumber(int n);
	TraceState BosonState(int n, int index);
	TraceState FermionState(int n, int index);
	void InitStateCollection(StateCollection* collection);
	void GenerateSingleOperatorStates(int bits, vector<TraceState> &res);
};
