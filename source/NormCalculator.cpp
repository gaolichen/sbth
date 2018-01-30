#pragma warning(disable:4018)
#include "NormCalculator.h"
#include <vector>
#include <algorithm>
#include <cstring>
using namespace std;

NormCalculator::NormCalculator()
{
	// ..
}

Polynomial& NormCalculator::Calculate(const TraceState& a, const TraceState& b)
{
	if (a.FermionNumber() != b.FermionNumber() || a.TotalBits() != b.TotalBits())
	{
		return zero;
	}

	if (a < b)
	{
		return DoCalculate(a, b);
	}
	else
	{
		return DoCalculate(b, a);
	}
}

BruteForceCalculator::BruteForceCalculator()
{
	// ..
}

void BruteForceCalculator::Positions(const TraceState& state, vector<int>& bPos, vector<int>& fPos)
{
	int p = 0;
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		const SingleTrace& single = state.Trace(i);
		for (int j = 0; j < single.BitNumber(); j++)
		{
			if (single.Bit(single.BitNumber() - j - 1) == 0)
			{
				bPos.push_back(p + j);
			}
			else
			{
				fPos.push_back(p + j);
			}
		}

		p += single.BitNumber();
	}
}

void BruteForceCalculator::IndexPosMap(const TraceState& state, vector<int>& index2Pos, vector<int>& pos2Index)
{
	int p = 0;
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		int bits = state.Trace(i).BitNumber();
		for (int j = 0; j < bits; j++)
		{
			index2Pos[p + j] = p + (j + 1) % bits;
			pos2Index[p + (j + 1) % bits] = j + p;
		}

		p += bits;
	}
}

/// Calculate <0| A B |0> where A, B are two tracestates.
/// It bruteforcely go through all the possible contractions of A and B.
/// For each contraction, find how many pow of N it produce and parity of it (even or odd)
Polynomial& BruteForceCalculator::DoCalculate(const TraceState& a, const TraceState& b)
{
	vector<int> bPos1, bPos2, fPos1, fPos2;
	Positions(a, bPos1, fPos1);
	Positions(b, bPos2, fPos2);
	pair<TraceState, TraceState> key = make_pair(a, b);
	cIt = res.find(key);
	if (cIt != res.end())
	{
		return cIt->second;
	}

	Polynomial ret;
	int bits = a.TotalBits();
	vector<int> pos2Tindex(bits);
	vector<int> dIndex2Pos(bits);
	vector<int> temp(bits);
	IndexPosMap(a, dIndex2Pos, temp);
	IndexPosMap(b, temp, pos2Tindex);
	vector<int> creator2Annihilator(bits);
	vector<int> annihilator2Creator(bits);

	bool odd = false;

	sort(fPos1.begin(), fPos1.end());
	do
	{
		odd = (InverseNumber(fPos1) % 2 == 1);
		for (int i = 0; i < fPos1.size(); i++)
		{
			creator2Annihilator[fPos2[i]] = fPos1[i];
			annihilator2Creator[fPos1[i]] = fPos2[i];
		}

		sort(bPos1.begin(), bPos1.end());
		do
		{
			for (int i = 0; i < bPos1.size(); i++)
			{
				creator2Annihilator[bPos2[i]] = bPos1[i];
				annihilator2Creator[bPos1[i]] = bPos2[i];
			}

			memset(visited, 0x0, sizeof(visited));

			// find how many loops in this contraction configuration
			int loop = 0;
			for (int dIndex = 0; dIndex < bits; dIndex++)
			{
				if (visited[dIndex]) continue;
				loop++;
				int i = dIndex;
				do{
					visited[i] = true;
					int j = creator2Annihilator[i];
					int k = dIndex2Pos[j];
					int x = annihilator2Creator[k];
					i = pos2Tindex[x];
				} while (!visited[i]);
			}

			if (odd)
			{
				ret.Increase(loop, -1);
			}
			else
			{
				ret.Increase(loop, 1);
			}
		} while(next_permutation(bPos1.begin(), bPos1.end()));
	} while(next_permutation(fPos1.begin(), fPos1.end()));

	res[key] = ret;

	return res[key];
}

Polynomial& RecursiveCalculator::DoCalculate(const TraceState&, const TraceState&)
{
	return zero;
}
