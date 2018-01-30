#pragma once
#include <string>
#include <vector>
#include "BitUtility.h"
using namespace std;


class SingleTrace
{
private:
	int mask;
public:
	SingleTrace();
	SingleTrace(string trace);
	SingleTrace(int trace, int bitNumber);
	int BitNumber() const;
	int Bit(int index) const;
	int Trace() const;
	int FermiNumber() const;
	bool IsFermionic() const;

	// normalize the trace to it's minimum form
	// return the change of parity
	int Normalize();
	void Split(int pos, SingleTrace& a, SingleTrace& b) const;
	void Split(int pos1, int pos2, SingleTrace& a, SingleTrace& b, SingleTrace& c) const;
	string ToString() const;
	string ToLaTeX() const;

	bool operator< (const SingleTrace& other) const;
	friend ostream& operator<<(ostream& os, const SingleTrace& single);

	static SingleTrace Merge(SingleTrace& a, SingleTrace& b);
	static SingleTrace Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c);
};
