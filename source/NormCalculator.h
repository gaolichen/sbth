#pragma once
#include "TraceState.h"
#include "Polynomial.h"
#include <vector>
#include <map>
using namespace std;

class NormCalculator
{
private:
protected:
	Polynomial zero;
	virtual Polynomial& DoCalculate(const TraceState&, const TraceState&) = 0;
public:
	NormCalculator();
	virtual Polynomial& Calculate(const TraceState&, const TraceState&);
};

class BruteForceCalculator : public NormCalculator
{
private:
	bool visited[1<<5];
	map< pair<TraceState, TraceState>, Polynomial>::iterator cIt;
	map< pair<TraceState, TraceState>, Polynomial> res;
	static void Positions(const TraceState&, vector<int>& bPos, vector<int>& fPos);
	static void IndexPosMap(const TraceState&, vector<int>& index2Pos, vector<int>& pos2Index);
protected:
	virtual Polynomial& DoCalculate(const TraceState&, const TraceState&);
public:
	BruteForceCalculator();
};

class RecursiveCalculator : public NormCalculator
{
private:

protected:
	virtual Polynomial& DoCalculate(const TraceState&, const TraceState&);
public:
	RecursiveCalculator();
};

