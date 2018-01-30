#pragma once
#include <vector>
#include <iostream>
#include <string>
#include "SingleTrace.h"
#include "Coefficient.h"
using namespace std;

class TraceState
{
private:
	vector<SingleTrace> vs;
public:
	TraceState();
	TraceState(string);
	TraceState(vector<SingleTrace> &);
	int TraceNumber() const;
	int FermionNumber() const;
	int TotalBits() const;
	void AddTrace(SingleTrace trace);
	void QuickNormalize();
	Coefficient Normalize();
	vector<SingleTrace>& Traces();
	const SingleTrace& Trace(int n) const;
	string ToString() const;
	string ToLaTeX() const;
	void CopyTo(TraceState& state, int exclude) const;
	void CopyTo(TraceState& state, int exclude1, int exclude2) const;
	bool operator< (const TraceState&) const;
	friend ostream& operator<<(ostream& os, const TraceState& ms);
	static void Merge(const TraceState& a, const TraceState& b, TraceState& res);
};
