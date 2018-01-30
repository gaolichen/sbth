#pragma warning(disable:4018)
#include "BitUtility.h"
#include <time.h>
using namespace std;

int BitCount(int i)
{
	 i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

int CyclicRotation(int n, int bitNumber)
{
	return ((n & 1) << (bitNumber -1)) | (n >> 1);
}

snum BinomialCoefficient(snum a,i64 b)
{
	if (b > a) return 0;
	if (b + b > a) b = a - b;
	snum ret = 1;
	for (int i = 0; i < b; i++)
	{
		ret *= (a - i);
		ret /= (i+1);
	}

	return ret;
}

int PickBits(int n, int bitNumber)
{
	return n & ((1 << bitNumber) - 1);
}

bool IsBitSet(int n, int bit)
{
	return (n & (1 << bit)) != 0;
}

vector<int> Num2Digit(i64 n, int maxBit)
{
	vector<int> ret;
	ret.reserve(maxBit);
	for (int i = 0; i < maxBit; i++)
	{
		if (IsBitSet(n, i)) ret.push_back(i);
	}

	return ret;
}

i64 Digit2Num(int* v, int size)
{
	i64 ret = 0;
	for (int i = 0; i < size; i++)
	{
		ret |= (1<< v[i]);
	}

	return ret;
}

i64 Factorial(int n)
{
	i64 ret = 1;
	while (n > 0) { ret *= n; n--; }
	return ret;
}

int SymmetryFactor(vector<i64>& mode, int s)
{
	i64 ret = Factorial(s);
	int r = 1;
	for (int j = 1; j < mode.size(); j++)
	{
		if (mode[j] == mode[j-1]) r++;
		else
		{
			ret /= Factorial(r);
			r = 1;
		}
	}

	ret /= Factorial(r);

	return (int)ret;
}

int InverseNumber(const vector<int>& v)
{
	int ret = 0;
	for (int i = 1; i < v.size(); i++)
	{
		for (int j = 0; j < i; j++)
			if (v[i] < v[j]) ret++;
	}

	return ret;
}

int BuildMask(int trace, int bitNumber)
{
	return (bitNumber << MAX_TRACE_BITS) | trace;
}

int Gcd(int a, int b)
{
	if (b == 0) return a;
	return Gcd(b, a % b);
}

string CombinePath(string path1, string path2)
{
	return path1 + "/" + path2;
}

string Bits2String(int bits, int bitNumber)
{
	string ret(bitNumber, 'a');
	for (int i = 0; i < bitNumber; i++)
	{
		if ((bits & (1 << i)) != 0)
		{
			ret[bitNumber - 1 - i] = char('b');
		}
	}

	return ret;
}

string ToUpper(string s)
{
	string ret = s;
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] >= 'a' && s[i] <= 'z')
		{
			ret[i] = char((s[i] - 'a') + 'A');
		}
	}

	return ret;
}

Stopwatch::Stopwatch()
{
	start = 0;
}

void Stopwatch::Start()
{
	start = clock();
}

/// return how much time in second since last Start() function call
double Stopwatch::Stop()
{
	double ret = (clock() - start) / (double)CLOCKS_PER_SEC;
	start = 0;

	return ret;
}
