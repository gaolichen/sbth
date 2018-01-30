#pragma once
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"
#include "StringBitMatrices.h"

class BitManager
{
private:
	vector<vector<vector<int> > > bitCollection;
	vector<vector<vector<int> > > bitchecks;
public:
	void Init(int maxBit);

	vector<int>& GetNumbers(int maxBit, int bits)
	{
		return bitCollection[maxBit][bits];
	}

	vector<int>& GetChecks(int maxBit, int bits)
	{
		return bitchecks[maxBit][bits];
	}
};

class Partitioner
{
private:
	int M;
	int s;
	int bitUnit;
	int allOnes;
	vector<i64> modes;
	vector<vector<int> > allPartitionsTest;
	vector<vector<int> > res;
	BitManager bm; 
	
	void AllPartitionsRecursive(vector<int>& partition, int total, int bit);
public:
	Partitioner(int M_, int s_);
	
	// return all possible partitions: a_1 * 1 + a_2 * 2 + ... + a_M * M = n * M if s*(M-1) is even 
	// or a_1 * 1 + a_2 * 2 + ... + a_M * M = (n + 1/2)M when s*(M-1) is odd.
	// where 0 <= a_i <= s.  
	// Each element of the returned in array is an integer mode = a_1 + a_2 * S + a_3 * S^2 + .. + a_M * S^(M-1), 
	// where S is smallest number of the form 2^n grater than s. 
	vector<i64>& AllPartitions();

	// covert the mode in to array. For each mode returned by AllPartitions function,
	// it returns {a_1, a_2, a_3, ..., a_M }. 
	vector<int> ToArray(i64 mode);

	// return all partitions bruteforcely. this function is for test purpose. 
	// The returned values are equivalent to: v = AllPartitions(); foreach (e in v) ret.push_back(ToArray(e));
	vector<vector<int> >& AllPartitionsBruteForce();
};

class DividePartition
{
private:
	int s;
	BitManager bm; 
	vector<vector<i64> > res;

	void DoDivide(vector<int>& partition, int index, int offset, vector<i64>& division);
public:
	DividePartition(int s_) : s(s_)
	{
		bm.Init(s + 1);
	}

	vector<vector<i64> >& Divide(vector<int>& partition, int offset);

	vector<vector<i64> >& Divide2(vector<int>& partition, vector<i64>& division);
};

class DividePartitionForTest
{
private:
	int s;
	vector<vector<i64> > res;
	void DoDivide(vector<int>& partition, int index, int spin, int offset, vector<i64>& division);
public:
	DividePartitionForTest(int s_) : s(s_) {}

	vector<vector<i64> >& Divide(vector<int>& partiton, int offset);

	vector<vector<i64> >& Divide2(vector<int>& partition, vector<i64>& division);
};

class ModesGenerator
{
private:
	int M;
	int L;
	int s;
	Partitioner* partitioner1;
	Partitioner* partitioner2;
	DividePartition* divide1;
	DividePartition* divide2;
	vector<vector<i64> > allModes;
	vector<DT> allEnergies;
	
public:
	ModesGenerator(int M_, int L_, int s_);

	~ModesGenerator();

	vector<vector<i64> >& Generate();
	vector<DT>& AllEnergies() { return allEnergies; }
};

