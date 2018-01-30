#pragma warning(disable:4018)
#include <iostream>
#include <vector>
#include <cassert>
#include "BitUtility.h"
#include "Partition.h"
#include "EnergyCalculator.h"
using namespace std;

void BitManager::Init(int maxBit)
{
	for (int i = 0; i <= maxBit; i++)
	{
		bitCollection.push_back(vector<vector<int> >(i + 1));
		bitchecks.push_back(vector<vector<int> >(i + 1));
	}

	for (int i = 0; i < (1 << maxBit); i++)
	{
		int check = 0;
		for (int j = 0; j < maxBit; j++)
		{
			if (IsBitSet(i, j) && !IsBitSet(i, j + 1))
			{
				check |= (1 << (j + 1));
			}
		}

		int bits = BitCount(i);
		for (int j = maxBit; (1<<j) > i; j--)
		{
			bitCollection[j][bits].push_back(i);
			bitchecks[j][bits].push_back(PickBits(check, j));
		}
	}
}

Partitioner::Partitioner(int M_, int s_): M(M_), s(s_)
{
	bitUnit = 1;
	while ((1 << bitUnit) <= s) bitUnit++;
	allOnes = (1 << bitUnit) - 1;
	bm.Init(s + 1);
}

vector<i64>& Partitioner::AllPartitions()
{
	i64 n = 1;
	for (int i = 0; i < M - 1; i++) n *= s + 1;
	//cout << "bitUnit=" << bitUnit << endl << " n=" << n << endl;

	for (int i = 0; i < n; i++)
	{
		i64 toadd = 0;
		int ii = i;
		int count = 0;
		int k = 0;
		//cout << "i=" << i << endl;
		while (ii > 0)
		{
			int j = ii % (s + 1);
			ii /= s + 1;
			toadd |= (j<<(k * bitUnit));
			k++;
			count += j * k;
		}

		if (((M & 1) == 1) || ((s & 1) == 0))
		{
			if (count % M == 0) 
			{ 
				modes.push_back(toadd);
			}
		}
		else if (count % M == M / 2)
		{
			modes.push_back(toadd);
		}
	}

	return modes;	
}

vector<int> Partitioner::ToArray(i64 mode)
{
	vector<int> ret(M);
	for (int i = 0; i < M - 1; i++)
	{
		ret[i] = (mode & allOnes);
		mode >>= bitUnit;
	}

	return ret;
}

vector<vector<int> >& Partitioner::AllPartitionsBruteForce()
{
	vector<int> partition(M);
	AllPartitionsRecursive(partition, 0, M - 1);
	return allPartitionsTest;
}

void Partitioner::AllPartitionsRecursive(vector<int>& partition, int total, int bit)
{
	if (bit == 0)
	{
		if (s * (M - 1) % 2 == 0 && total % M == 0)
		{
			allPartitionsTest.push_back(partition);
		}
		else if (s * (M - 1) % 2 == 1 && total % M == M / 2)
		{
			allPartitionsTest.push_back(partition);
		}

		return;
	}
	
	for (int i = 0; i <= s; i++)
	{
		partition[bit - 1] = i;
		AllPartitionsRecursive(partition, total + i * bit, bit - 1);
	}
}

vector<vector<i64> >& DividePartition::Divide(vector<int>& partition, int offset)
{
	vector<i64> division(s);
	res.clear();
	DoDivide(partition, partition.size() - 1, offset, division);

	return res;
}

vector<vector<i64> >& DividePartition::Divide2(vector<int>& partition, vector<i64>& division)
{
	res.clear();
	DoDivide(partition, partition.size() - 1, 0, division);
	return res;
}

void DividePartition::DoDivide(vector<int>& partition, int index, int offset, vector<i64>& division)
{
	if (index < 0)
	{
		res.push_back(division);
		return;
	}

	int n = partition[index];
	vector<int> positions = bm.GetNumbers(s, n);
	for (int i = 0; i < positions.size(); i++)
	{
		int pos = positions[i];
		bool ok = true;
		for (int j = 0; j < s; j++)
		{
			if (!IsBitSet(pos, j)) continue;
			if (j > 0 && !IsBitSet(pos, j-1) && division[j] == division[j-1]) 
			{
				ok = false; 
				break;
			}
		}

		if (!ok) continue;
		for (int j = 0; j < s; j++)
		{
			if (IsBitSet(pos, j)) division[j] |= ((i64)1<<(index + 1 + offset));
		}

		DoDivide(partition, index - 1, offset, division);

		for (int j = 0; j < s; j++)
		{
			if (IsBitSet(pos, j)) division[j] ^= ((i64)1<<(index + 1 + offset));
		}
	}
}

void DividePartitionForTest::DoDivide(vector<int>& partition, int index, int spin, int offset, vector<i64>& division)
{
	if (index < 0)
	{
		res.push_back(division);
		return;
	}

	if (partition[index] == 0)
	{
		DoDivide(partition, index - 1, 0, offset, division);
		return;
	}

	if (spin == s)
	{
		return;
	}

	i64 n = (i64)1<<(index + 1 + offset);
	assert((division[spin] & n) == 0);
	
	//if ((division[spin] & n) != 0)
	//{
	//	cout << "error! DoDivide()! " << endl;
	//}

	if (spin == 0 || division[spin - 1] >= (division[spin] | n))
	{
		division[spin] |= n;
		partition[index]--;
		DoDivide(partition, index, spin + 1, offset, division);
		division[spin] ^= n;
		partition[index]++;
	}

	DoDivide(partition, index, spin + 1, offset, division);
}

vector<vector<i64> >& DividePartitionForTest::Divide2(vector<int>& partition, vector<i64>& division)
{
	res.clear();
	DoDivide(partition, partition.size() - 1, 0, 0, division);
	return res;
}

vector<vector<i64> >& DividePartitionForTest::Divide(vector<int> &partition, int offset)
{
	vector<i64> division(s);
	res.clear();
	DoDivide(partition, partition.size() - 1, 0, offset, division);
	return res;
}

ModesGenerator::ModesGenerator(int M_, int L_, int s_) : M(M_), L(L_), s(s_)
{
	this->partitioner1 = new Partitioner(M - L, s);
	this->partitioner2 = new Partitioner(L, s);
	this->divide1 = new DividePartition(s);
	this->divide2 = new DividePartition(s);
}

ModesGenerator::~ModesGenerator()
{
	delete this->partitioner1;
	delete this->partitioner2;
	delete this->divide1;
	delete this->divide2;
}

vector<vector<i64> >& ModesGenerator::Generate()
{
	vector<vector<int> >& partition1 = partitioner1->AllPartitionsBruteForce();
	vector<vector<int> >& partition2 = partitioner2->AllPartitionsBruteForce();
	//cout << "partition1 size: " << partition1.size() << endl;
	//cout << "partition2 size: " << partition2.size() << endl;
	vector<DT> energies1(partition1.size());
	vector<DT> energies2(partition2.size());

	for (int i = 0; i < partition1.size(); i++) 
		energies1[i] = EnergyCalculator::Operators2Energy(partition1[i], M - L, s);
	for (int i = 0; i < partition2.size(); i++)
		energies2[i] = EnergyCalculator::Operators2Energy(partition2[i], L, s);

	for (int i = 0; i < partition1.size(); i++)
	{
		vector<vector<i64> >& division1 = divide1->Divide(partition1[i], L - 1);
		for (int k = 0; k < division1.size(); k++)
		{
			for (int j = 0; j < partition2.size(); j++)
			{
				vector<vector<i64> >& division2 = divide2->Divide2(partition2[j], division1[k]);
				allModes.insert(allModes.end(), division2.begin(), division2.end());
				allEnergies.insert(allEnergies.end(), division2.size(), energies1[i] + energies2[j]);
			}
		}
	}

	return allModes;
}
