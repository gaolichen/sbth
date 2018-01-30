#pragma warning(disable:4018)
#include "StateCollection.h"
#include "StateGenerator.h"

StateCollection* StateCollection::inst = NULL;

StateCollection::StateCollection()
{
	// .. do nothing.
	stateList.resize(StateGenerator::MAX_BIT_TO_GENERATE + 1, vector<TraceState>());
}

StateCollection* StateCollection::Inst()
{
	if (inst == NULL)
	{
		inst = new StateCollection();
	}

	return inst;
}

void StateCollection::Init(int bits, vector<TraceState>& bosons, vector<TraceState> fermions)
{
	vector<TraceState>& vt = stateList[bits];
	for (int i = 0; i < bosons.size(); i++)
	{
		vt.push_back(bosons[i]);
		vt.push_back(fermions[i]);
		this->state2Id[bosons[i]] = StateId(bits, 2 * i);
		this->state2Id[fermions[i]] = StateId(bits, 2 * i + 1);
	}
}

StateId StateCollection::GetId(const TraceState& state) const
{
	map<TraceState, StateId>::const_iterator it;
	it = state2Id.find(state); 
	if (it == state2Id.end())
	{
		return StateId(-1, -1);
	}

	return it->second;
}

const TraceState& StateCollection::GetState(const StateId& id) const
{
	return stateList[id.BitNumber][id.Index];
}

int StateCollection::StateNumber(int bits) const
{
	return stateList[bits].size() / 2;
}

const TraceState& StateCollection::GetBosonState(int bits, int index) const
{
	return stateList[bits][2 * index];
}

const TraceState& StateCollection::GetFermionState(int bits, int index) const
{
	return stateList[bits][2 * index + 1];
}

const TraceState& StateCollection::GetState(int bits, int index, StateType type) const
{
	if (type == Boson)
	{
		return GetBosonState(bits, index);
	}
	else
	{
		return GetFermionState(bits, index);
	}
}