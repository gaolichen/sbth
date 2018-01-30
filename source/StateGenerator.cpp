#pragma warning(disable:4018)
#include "StateGenerator.h"
#include "BitUtility.h"
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;
//#define TEST_STATENUMBER

StateGenerator::StateGenerator()
{
	vector<vector<snum> > numbers = vector<vector<snum> >(MAX_BIT_TO_COUNT + 1, vector<snum>(2, (snum)(-1)));
	stateNumbers.resize(MAX_BIT_TO_COUNT + 1, numbers);
	myFlags = new bool[1 << MAX_BIT_TO_GENERATE];
	InitSingleTraceNumber();
}

StateGenerator::~StateGenerator()
{
	delete[] myFlags;
}

void StateGenerator::InitSingleTraceNumber()
{
	singleTraceNumbers.resize(MAX_BIT_TO_COUNT + 1);
	for (int m = 1; m <= MAX_BIT_TO_COUNT; m++)
	{
#ifdef TEST_STATENUMBER
		singleTraceNumbers[m] = ((((i64)1) << (m - 1)) + m - 1) / m;
#else
		int p = m;
		while (p % 2 == 0)
		{
			p /= 2;
		}

		snum res = 0;
		for (int i = 1; i <= p; i++)
		{
			int toShift = m / p * Gcd(i, p);
			res += ((i64)1 << toShift);
		}

		singleTraceNumbers[m] = res / (2 * m);
#endif
	}
}

void StateGenerator::GeneratSingleStates()
{
	singleFermions.resize(MAX_BIT_TO_GENERATE + 1);
	singleBosons.resize(MAX_BIT_TO_GENERATE + 1);
	for (int i = 1; i <= MAX_BIT_TO_GENERATE; i++)
	{
		FindSingleStates(i);
	}
}

void StateGenerator::FindSingleStates(int bit)
{
	memset(myFlags, 0, sizeof(bool) * (1<<bit));

	for (int i = 0; i < (1<<bit); i++)
	{
		if (myFlags[i]) continue;
		int n = i;
		for (int j = 1; j < bit; j++)
		{
			n = CyclicRotation(n, bit);
			if (n == i)
			{
				// .. check if the states is zero.
				if (BitCount(i) % 2 == 0 && BitCount(i & ((1<<j)-1)) % 2 == 1)
				{
					myFlags[i] = true;
				}
			}
			else myFlags[n] = true;
		}

		if (myFlags[i]) continue;
		if ((BitCount(i) & 1) == 0)
		{
			singleBosons[bit].push_back(SingleTrace(i, bit));
		}
		else
		{
			singleFermions[bit].push_back(SingleTrace(i, bit));
		}
	}
}

snum StateGenerator::SingleStateNumber(int n)
{
	return singleTraceNumbers[n];
}

snum StateGenerator::BosonNumber(int n)
{
	return StateNumbers(n, n, 0);
}

snum StateGenerator::FermionNumber(int n)
{
	return StateNumbers(n, n, 1);
}

snum StateGenerator::StateNumbers(int bit, int remain, int b)
{
	if (bit == 0)
	{
		if (remain == 0 && b == 0) return 1;
		else return 0;
	}
	if (remain == 0)
	{
		if (b == 0) return 1;
		return 0;
	}

	snum ret = stateNumbers[bit][remain][b];
	if (ret >= 0) return ret;
	int a = remain / bit;

	ret = 0;
	for (int i = 0; i <= a; i++)
		for (int j = 0; i + j <= a; j++)
		{
			//int n1 = BinomialCoefficient(singleFermions[bit].size(), i);
			//int n2 = BinomialCoefficient(singleBosons[bit].size() - 1 + j, j);
			snum n1 = BinomialCoefficient(singleTraceNumbers[bit], (i64)i);
			snum n2 = BinomialCoefficient(singleTraceNumbers[bit] - 1 + j, (i64)j);
			ret += n1 * n2 * StateNumbers(bit - 1, remain - (i+j) * bit, (b+i)%2);
		}
		stateNumbers[bit][remain][b] = ret;

	return ret;
}

void StateGenerator::GenerateAllStates()
{
	GeneratSingleStates();

	memset(visited, 0, sizeof(visited));
	fermions.resize(MAX_BIT_TO_GENERATE + 1, vector<vector<TraceState> >(MAX_BIT_TO_GENERATE + 1, vector<TraceState>()));
	bosons.resize(MAX_BIT_TO_GENERATE + 1, vector<vector<TraceState> >(MAX_BIT_TO_GENERATE + 1, vector<TraceState>()));
	bosons[0][0] = vector<TraceState>(1, TraceState());

	vector<vector<TraceState> > fv;
	vector<vector<TraceState> > bv;

	for (int i = 1; i <= MAX_BIT_TO_GENERATE; i++)
	{
		// find TraceStates for i-bit.
		fv.clear();
		bv.clear();

		// maximum number of i-bit states can be in one state.
		int n = MAX_BIT_TO_GENERATE / i;
		// fv[k]: all states consisting of k single i-bit fermions
		// bv[k]: all states cosnsiting of k single i-bit bosons.
		for (int k = 0; k <= n; k++)
		{
			fv.push_back(PickFermionFromSingleState(singleFermions[i], k));
			bv.push_back(PickBosonFromSingleState(singleBosons[i], k));
		}

		for (int j = 0; j <= MAX_BIT_TO_GENERATE; j++)
		{
			int a = j / i;
			for (int k = 0; k <= a; k++)
			{
				// pick k different single states fermions from i bit states.
				for (int h = 0; h + k <= a; h++)
				{
					// pick h single state bosons from i bit states.
					for (int x = 0; x < fv[k].size(); x++)
						for (int y = 0; y < bv[h].size(); y++)
						{
							// even number of i-bit fermions
							if (k % 2 == 0)
							{
								vector<TraceState>& vts1 = fermions[i - 1][j - (k + h) * i];
								for (int z = 0; z < vts1.size(); z++)
								{
									fermions[i][j].push_back(CombineStates(vts1[z], fv[k][x], bv[h][y]));
								}
								vector<TraceState>& vts2 = bosons[i - 1][j - (k + h) * i];
								for (int z = 0; z < vts2.size(); z++)
								{
									bosons[i][j].push_back(CombineStates(vts2[z], fv[k][x], bv[h][y]));
								}
							}
							else
							{
								// odd number of i-bit fermions.
								vector<TraceState>& vts3 = fermions[i - 1][j - (k + h) * i];
								for (int z = 0; z < vts3.size(); z++)
								{
									bosons[i][j].push_back(CombineStates(vts3[z], fv[k][x], bv[h][y]));
								}
								vector<TraceState>& vts4 = bosons[i - 1][j - (k + h) * i];
								for (int z = 0; z < vts4.size(); z++)
								{
									fermions[i][j].push_back(CombineStates(vts4[z], fv[k][x], bv[h][y]));
								}
							}
						}
				}
			}
		}
	}

	for (int i = 1; i <= MAX_BIT_TO_GENERATE; i++)
	{
		for (int j = 0; j < bosons[i][i].size(); j++)
		{
			bosons[i][i][j].QuickNormalize();
			fermions[i][i][j].QuickNormalize();
		}

		sort(bosons[i][i].begin(), bosons[i][i].end());
		sort(fermions[i][i].begin(), fermions[i][i].end());
	}
}

void StateGenerator::DoPickFermion(vector<SingleTrace>& allstates, int index, int remain, 
	vector<SingleTrace>& curr, vector<TraceState>& ret)
{
	if (remain == 0)
	{
		ret.push_back(TraceState(curr));
		return;
	}

	if (allstates.size() - index < remain)
	{
		return;
	}

	curr.push_back(allstates[index]);
	DoPickFermion(allstates, index + 1, remain - 1, curr, ret);
	curr.pop_back();

	DoPickFermion(allstates, index + 1, remain, curr, ret);
}

void StateGenerator::DoPickBoson(vector<SingleTrace>& allstates, int index, int remain,
	vector<SingleTrace>& curr, vector<TraceState>& ret)
{
	if (remain == 0)
	{
		ret.push_back(TraceState(curr));
		return;
	}

	if (index == allstates.size())
	{
		if (remain == 0)
		{
			ret.push_back(TraceState(curr));
		}
		else
		{
			//.. do nothing.
		}

		return;
	}

	for (int i = 0; i <= remain; i++)
	{
		DoPickBoson(allstates, index + 1, remain - i, curr, ret);
		curr.push_back(allstates[index]);
	}

	for (int i = 0; i <= remain; i++)
	{
		curr.pop_back();
	}
}

vector<TraceState> StateGenerator::PickFermionFromSingleState(vector<SingleTrace> &states, int number)
{
	vector<TraceState> ret;
	vector<SingleTrace> curr;
	DoPickFermion(states, 0, number, curr, ret);

	return ret;
}

vector<TraceState> StateGenerator::PickBosonFromSingleState(vector<SingleTrace> &states, int number)
{
	vector<TraceState> ret;
	vector<SingleTrace> curr;
	DoPickBoson(states, 0, number, curr, ret);

	return ret;
}

TraceState StateGenerator::CombineStates(TraceState& a, TraceState& b, TraceState& c)
{
	TraceState ret;
	for (int i = 0; i < a.TraceNumber(); i++)
	{
		ret.AddTrace(a.Traces()[i]);
	}

	for (int i = 0; i < b.TraceNumber(); i++)
	{
		ret.AddTrace(b.Traces()[i]);
	}

	for (int i = 0; i < c.TraceNumber(); i++)
	{
		ret.AddTrace(c.Traces()[i]);
	}

	return ret;
}

TraceState StateGenerator::BosonState(int n, int index)
{
	return bosons[n][n][index];
}

TraceState StateGenerator::FermionState(int n, int index)
{
	return fermions[n][n][index];
}

void StateGenerator::InitStateCollection(StateCollection* collection)
{
	for (int i = 1; i <= MAX_BIT_TO_GENERATE; i++)
	{
		collection->Init(i, bosons[i][i], fermions[i][i]);
	}
}

void StateGenerator::BuildSingleOperatorStates(int remBit, int currBits, vector<int>& res)
{
	if (remBit == 0)
	{
		// todo 
		return;
	}

	if (currBits > remBit)
	{
		return;
	}

	for (int i = 0; i * currBits <= remBit; i++)
	{
		res[currBits] = i;
		BuildSingleOperatorStates(remBit - i * currBits, currBits + 1, res);
	}
}

void StateGenerator::GenerateSingleOperatorStates(int bits, vector<TraceState>& res)
{
	vector<int> states(bits + 1, 0);
	BuildSingleOperatorStates(bits, 1, states);
}