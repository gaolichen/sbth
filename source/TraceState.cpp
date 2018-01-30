#pragma warning(disable:4018)
#include "TraceState.h"
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

TraceState::TraceState()
{
}

TraceState::TraceState(vector<SingleTrace>& traces)
{
	vs = traces;
}

TraceState::TraceState(string traces)
{
	int pos, next = -1;
	do
	{
		pos = next + 1;
		next = traces.find_first_of(',', pos);
		if (next != string::npos)
			vs.push_back(SingleTrace(traces.substr(pos, next - pos)));
		else vs.push_back(SingleTrace(traces.substr(pos, next)));
	} while (next != string::npos);
}

int TraceState::FermionNumber() const
{
	int ret = 0;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].FermiNumber();
	}

	return ret;
}

int TraceState::TotalBits() const
{
	int ret = 0;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].BitNumber();
	}

	return ret;
}

void TraceState::AddTrace(SingleTrace trace)
{
	vs.push_back(trace);
}

vector<SingleTrace>& TraceState::Traces()
{
	return vs;
}

const SingleTrace& TraceState::Trace(int n) const
{
	return vs[n];
}

void TraceState::QuickNormalize()
{
	if (vs.size() <= 1) return;
	sort(vs.begin(), vs.end());
}

bool TraceState::operator< (const TraceState& other) const
{
	int f1 = FermionNumber();
	int f2 = other.FermionNumber();
	if (f1 != f2)
	{
		return f1 < f2;
	}

	if (vs.size() != other.vs.size())
	{
		return vs.size() < other.vs.size();
	}

	for (int i = 0; i < vs.size(); i++)
	{
		if (vs[i] < other.vs[i]) return true;
		else if (other.vs[i] < vs[i]) return false;
	}

	return false;
}

int TraceState::TraceNumber() const
{
	return vs.size();
}

string TraceState::ToString() const
{
	string ret;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += "(" + vs[i].ToString() + ")";
	}

	return ret;
}

ostream& operator<<(ostream& os, const TraceState& ts)
{
	for (int i = 0; i < ts.vs.size(); i++)
	{
		os << ts.vs[i];
	}

	return os;
}

string TraceState::ToLaTeX() const
{
	string ret;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].ToLaTeX();
	}

	return ret + "\\left|0\\right\\rangle";
}

Coefficient TraceState::Normalize()
{
	int totF = 0;
	Coefficient ret(0, 0);
	for (int i = vs.size() - 1; i >= 0; i--)
	{
		if (vs[i].BitNumber() == 0)
		{
			ret.N++;
			vs.erase(vs.begin() + i);
		}
		else
		{
			totF += vs[i].Normalize();
		}
	}

	if (ret.N == 0)
	{
		ret.One = 1;
	}

	for (int i = 0; i < vs.size(); i++)
	{
		int m = i;
		for (int j = i + 1; j < vs.size(); j++)
		{
			if (vs[j] < vs[m])
			{
				m = j;
			}
		}

		if (m != i)
		{
			// swap vs[i], vs[m]
			int f = 0;
			for (int j = i; j < m; j++)
			{
				f += vs[j].FermiNumber();
			}
			// move vs[m] to pos i.
			totF += f * vs[m].FermiNumber();
			// move vs[i] to pos m.
			totF += vs[i].FermiNumber() * (f - vs[i].FermiNumber());

			swap(vs[i], vs[m]);
		}
	}

	if (totF % 2 == 1)
	{
		ret.Opposite();
	}
	return ret;
}

void TraceState::Merge(const TraceState& a, const TraceState& b, TraceState& res)
{
	res.vs = a.vs;
	res.vs.insert(res.vs.end(), b.vs.begin(), b.vs.end());
}

void TraceState::CopyTo(TraceState& state, int exclude) const
{
	state.vs.insert(state.vs.end(), vs.begin(), vs.begin() + exclude);
	state.vs.insert(state.vs.end(), vs.begin() + exclude + 1, vs.end());
}

void TraceState::CopyTo(TraceState& state, int exclude1, int exclude2) const
{
	state.vs.insert(state.vs.end(), vs.begin(), vs.begin() + exclude1);
	state.vs.insert(state.vs.end(), vs.begin() + exclude1 + 1, vs.begin() + exclude2);
	state.vs.insert(state.vs.end(), vs.begin() + exclude2 + 1, vs.end());
}
