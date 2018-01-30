#pragma warning(disable:4018)
#include <string>
#include <iostream>
#include "BitUtility.h"
#include "SingleTrace.h"

using namespace std;

SingleTrace::SingleTrace()
{
	SingleTrace(0, 0);
}

SingleTrace::SingleTrace(int trace, int bitNumber)
{
	mask = BuildMask(trace, bitNumber);
}

SingleTrace::SingleTrace(string trace)
{
	int n = 0;
	for (int i = 0; i < trace.length(); i++)
	{
		if (trace[i] == 'b')
		{
			n |= (1 << (trace.length() - i - 1));
		}
	}

	mask = BuildMask(n, trace.length());
}

int SingleTrace::BitNumber() const
{
	return mask >> MAX_TRACE_BITS;
}

int SingleTrace::Bit(int index) const
{
	if ((mask & (1 << index)) != 0)
	{
		return 1;
	}

	return 0;
}

int SingleTrace::Trace() const
{
	return mask & FULL_TRACE_BITS;
}

int SingleTrace::FermiNumber() const
{
	return BitCount(Trace());
}

bool SingleTrace::IsFermionic() const
{
	return (FermiNumber() & 1) == 1;
}

string SingleTrace::ToString() const
{
	string ret;
	int trace = Trace();
	for (int i = 0; i < BitNumber(); i++)
	{
		if ((trace & (1<<i)) != 0)
		{
			ret = "b" + ret;
		}
		else
		{
			ret = "a" + ret;
		}
	}

	return ret;
}

ostream& operator<<(ostream& os, const SingleTrace& single)
{
	os << "Tr(" << single.ToString() << ")";

	return os;
}

string SingleTrace::ToLaTeX() const
{
	string ret;
	int trace = Trace();
	for (int i = 0; i < BitNumber(); i++)
	{
		if ((trace & (1<<i)) != 0)
		{
			ret = "\\bar{b}" + ret;
		}
		else
		{
			ret = "\\bar{a}" + ret;
		}
	}

	return "\\mathrm{Tr}" + ret;
}

bool SingleTrace::operator< (const SingleTrace& other) const
{
	return mask < other.mask;
}

void SingleTrace::Split(int pos, SingleTrace& a, SingleTrace& b) const
{
	b.mask = BuildMask(PickBits(mask, pos), pos);
	int remain = (mask >> (pos + 1));
	a.mask = BuildMask(PickBits(remain, BitNumber() - pos - 1), BitNumber() - pos - 1);
}

void SingleTrace::Split(int pos1, int pos2, SingleTrace& a, SingleTrace& b, SingleTrace& c) const
{
	SingleTrace d(0, 0);
	Split(pos1, a, d);
	d.Split(pos2, b, c);
}

SingleTrace SingleTrace::Merge(SingleTrace& a, SingleTrace& b)
{
	int trace = (a.Trace() << b.BitNumber()) | b.Trace();
	return SingleTrace(trace, a.BitNumber() + b.BitNumber());	
}

SingleTrace SingleTrace::Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c)
{
	int trace = (a.Trace() << b.BitNumber()) | b.Trace();
	trace = trace << c.BitNumber() | c.Trace();
	return SingleTrace(trace, a.BitNumber() + b.BitNumber() + c.BitNumber());	
}

int SingleTrace::Normalize()
{
	int trace = this->Trace();
	int n = this->BitNumber();
	int min = trace;
	int pos = 0;
	for (int i = 1; i < n; i++)
	{
		trace = CyclicRotation(trace, n);
		if (trace < min)
		{
			min = trace; 
			pos = i;
		}
	}

	int b = BitCount(PickBits(this->Trace(), pos));

	int ret = b * (this->FermiNumber() - b);
	this->mask = BuildMask(min, BitNumber());
	return ret;
}
