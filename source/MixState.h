#pragma once
#include <map>
#include <string>
#include <iostream>
#include "Coefficient.h"
#include "TraceState.h"
#include "StateCollection.h"
#include "StateId.h"

using namespace std;

class MixState
{
private:
	map<StateId, Coefficient> coefficients;
public:
	MixState();
	void AddState(const StateId& stateId, const Coefficient& coefficient);
	void Merge(const MixState&, int prefactor);
	void Merge(const MixState& state);
	void Extend(TraceState& state, int parity);
	friend ostream& operator << (ostream& os, const MixState& ms);
	MixState& operator*= (int n);
	string ToString();
	map<StateId, Coefficient>::const_iterator Begin() const;
	map<StateId, Coefficient>::const_iterator End() const;
};