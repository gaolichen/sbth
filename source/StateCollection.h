#pragma once
#include <map>
#include <vector>
#include "TraceState.h"
#include "StateId.h"
#include "StateType.h"
using namespace std;

class StateCollection
{
private:
	map<TraceState, StateId> state2Id;
	vector<vector<TraceState> > stateList;
	StateCollection();
	static StateCollection* inst;
public:
	static StateCollection* Inst();
	void Init(int bits, vector<TraceState>& bosons, vector<TraceState> fermions);
	StateId GetId(const TraceState& state) const;
	const TraceState& GetState(const StateId& id) const;
	int StateNumber(int bits) const;
	const TraceState& GetBosonState(int bits, int index) const;
	const TraceState& GetFermionState(int bits, int index) const;
	const TraceState& GetState(int bits, int index, StateType type) const;
};