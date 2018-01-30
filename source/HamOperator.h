#pragma once
#include <map>
#include "MixState.h"
#include "SingleTrace.h"
#include "TraceState.h"
#include "StateCollection.h"
using namespace std;

class HamOperator
{
private:
	map<SingleTrace, MixState> cache1;
	map<pair<SingleTrace, SingleTrace>, MixState> cache2;
	MixState& ApplyOn(const SingleTrace& single);
	MixState& ApplyOn(const SingleTrace& single1, const SingleTrace& single2);
protected:
	void AddState(TraceState& state, int parity, MixState& res, int order2Change = 0);
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res) = 0;
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res) = 0;
public:
	const static int AA = 0;
	const static int AB = 1;
	const static int BA = 2;
	const static int BB = 3;
	HamOperator();
	virtual void ApplyOn(const TraceState& state, MixState& res);
	virtual string ToString() const = 0;
	friend ostream& operator<<(ostream& os, const HamOperator& op);
	friend bool operator == (const HamOperator& op1, const HamOperator& op2);
};

class BitNumberHamOperator : public HamOperator
{
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
	BitNumberHamOperator();
	virtual void ApplyOn(const TraceState& state, MixState& res);
	virtual string ToString() const;
};

class HamOperatorA : public HamOperator
{
private:
	int annihilator;
	int creator;
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
	HamOperatorA(int creator, int annihilator);
	virtual string ToString() const;
};

class HamOperatorB : public HamOperator
{
private:
	int annihilator;
	int creator;
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
	HamOperatorB(int creator, int annihilator);
	virtual string ToString() const;
};
