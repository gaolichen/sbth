#include "StateCounter.h"
#include <iostream>

StateCounter* StateCounter::inst = 0;

StateCounter::StateCounter()
{
}

StateCounter* StateCounter::Inst()
{
	if (inst == 0)
	{
		inst = new StateCounter();
		inst->Init();
	}

	return inst;
}


StateCounter::~StateCounter(void)
{
}

snum StateCounter::SingleTrace(int M, int s)
{
    if (s * M > MAX_BIT_TO_COUNT)
    {
        cout << "StateCounter::SingleTrace: Exceed max bit number to count" << endl;
        return -1;
    }

    return singleTraceNumbers[s][M];
}

snum StateCounter::MultiTrace(int M, int s)
{
    if (s * M > MAX_BIT_TO_COUNT)
    {
        cout << "StateCounter::MultiTrace: Exceed max bit number to count" << endl;
        return -1;
    }

    return StateNumbers(s, M, M, 0) >> (s - 1);
}

void StateCounter::InitSingleTraceNumber()
{
    singleTraceNumbers.resize(MAX_SPIN_TO_COUNT + 1, vector<snum>());

    for (int s = 1; s <= MAX_SPIN_TO_COUNT; s++)
    {
	    singleTraceNumbers[s].resize(MAX_BIT_TO_COUNT / s + 1);
	    for (int m = 1; m * s <= MAX_BIT_TO_COUNT; m++)
	    {
            // the formula for number of M-bit single trace states is:
            // 1/M * sum_{k=1}^{M} 2^[s * gcd(M,k) - s] where the sum is run over k such that M/gcd(M,k) is odd.

            // p is largest odd factor of m.
		    int p = m;
		    while (p % 2 == 0)
		    {
			    p /= 2;
		    }

            // since in the formula we sum over m/gcd(m,k) is odd,
            // it implies that we only need to enumurate all intergers divide p.
            // it also implies that m/p is a factor of gcd(m,k).
		    snum res = 0;
		    for (int i = 1; i <= p; i++)
		    {
			    int toShift = (m / p * Gcd(i, p) - 1) * s;
			    res += ((i64)1 << toShift);
		    }

		    singleTraceNumbers[s][m] = res / m;
	    }
    }
}

snum StateCounter::NoHalfMode(int M)
{
	if (M % 2 == 1)
	{
		cout << "Error: M should be even." << endl;
		return -1;
	}
	return noHalfMode[M / 2];
}

void StateCounter::Init()
{
    // initialize all state numbers to be zero.
    stateNumbers.resize(MAX_SPIN_TO_COUNT + 1);
    for (int s = 1; s <= MAX_SPIN_TO_COUNT; s++)
    {
        vector<vector<snum> > numbers = vector<vector<snum> >(MAX_BIT_TO_COUNT / s + 1, vector<snum>(2, (snum)(-1)));
	    stateNumbers[s].resize(MAX_BIT_TO_COUNT / s + 1, numbers);
    }

    oddOnlyStateNumbers.resize(MAX_BIT_TO_COUNT / 2 + 1, vector<snum>(MAX_BIT_TO_COUNT + 1, (snum)(-1)));
    evenTraceNumbers.resize(MAX_BIT_TO_COUNT + 1, vector<snum>(MAX_BIT_TO_COUNT + 1, (snum)(-1)));

    InitSingleTraceNumber();
	InitNoHalfMode();
    InitAverageEnergy();
}

void StateCounter::InitNoHalfMode()
{
	noHalfMode.resize(MAX_BIT_TO_COUNT / 2 + 1, 0);
	for (int M = 2; M <= MAX_BIT_TO_COUNT; M += 2)
	{
		vector<vector<snum> > num(2, vector<snum>(M, 0));
		int curr = 0;
		int prev = 1;
		num[prev][0] = 1;

		for (int i = 1; i < M; i++)
		{
			if (i * 2 == M) continue;

			for (int j = 0; j < M; j++)
			{
				num[curr][j] = num[prev][j] + num[prev][(j + M - i) % M];
			}

			swap(curr, prev);
		}

		noHalfMode[M / 2] = num[prev][M / 2];
	}
}

void StateCounter::InitAverageEnergy()
{
    aveESingle.resize(MAX_BIT_TO_COUNT / 2 + 1, 0.0);
    aveE.resize(MAX_BIT_TO_COUNT + 1, vector<long double>(MAX_BIT_TO_COUNT + 1, (long double)-1.0));

    for (int i = 2; i <= MAX_BIT_TO_COUNT; i += 2)
    {
        aveESingle[i / 2] = 4 - 8 * NoHalfMode(i) / (long double)SingleTrace(i);
    }
}

double StateCounter::AverageEnergy(int bit)
{
    return AverageEnergy(bit, bit, 0) / MultiTrace(bit);
}

long double StateCounter::AverageEnergy(int bit, int remain, int parity)
{
    if (remain == 0) return 0.0;
    if (bit <= 1) return 0.0;

    long double ret = aveE[bit][remain];
	if (ret >= 0.0) return ret;
	int a = remain / bit;

	ret = 0.0;
	for (int i = 0; i <= a; i++)
    {
		for (int j = 0; i + j <= a; j++)
		{
			snum n1 = BinomialCoefficient(SingleTrace(bit), (i64)i);
			snum n2 = BinomialCoefficient(SingleTrace(bit) - 1 + j, (i64)j);
            if (bit % 2 == 0)
            {
                ret += (long double)aveESingle[bit/2] * n1 * n2 * (i + j) * StateNumbers(1, bit - 1, remain - (i+j) * bit, (parity + i) % 2);
            }

			ret += n1 * n2 * AverageEnergy(bit - 1, remain - (i+j) * bit, (parity + i) % 2);
		}
    }

    aveE[bit][remain] = ret;

    return ret;
}

snum StateCounter::EvenTraceNumber(int bit, int remain, int parity)
{
    if (remain == 0) return 0;
    if (bit <= 1) return 0;

    snum ret = evenTraceNumbers[bit][remain];
	if (ret >= 0) return ret;
	int a = remain / bit;

	ret = 0;
	for (int i = 0; i <= a; i++)
    {
		for (int j = 0; i + j <= a; j++)
		{
			snum n1 = BinomialCoefficient(SingleTrace(bit), (i64)i);
			snum n2 = BinomialCoefficient(SingleTrace(bit) - 1 + j, (i64)j);
            if (bit % 2 == 0)
            {
                ret += n1 * n2 * (i + j) * StateNumbers(1, bit - 1, remain - (i+j) * bit, (parity + i) % 2);
            }

			ret += n1 * n2 * EvenTraceNumber(bit - 1, remain - (i+j) * bit, (parity + i) % 2);
		}
    }

    evenTraceNumbers[bit][remain] = ret;

    return ret;
}

snum StateCounter::OddOnlyStateNumbers(int bit, int remain, int parity)
{
	if (remain == 0)
	{
        if (parity == 0) return 1;
        return 0;
	}

    if (bit == 1)
	{
		return 1;
	}

    if (bit < 1)
    {
        cout << "OddOnlyStateNumbers: should not reach here..." << endl;
        return -1;
    }

	snum ret = oddOnlyStateNumbers[bit / 2][remain];
	if (ret >= 0) return ret;
	int a = remain / bit;

	ret = 0;
	for (int i = 0; i <= a; i++)
    {
		for (int j = 0; i + j <= a; j++)
		{
			snum n1 = BinomialCoefficient(SingleTrace(bit), (i64)i);
			snum n2 = BinomialCoefficient(SingleTrace(bit) - 1 + j, (i64)j);
			ret += n1 * n2 * OddOnlyStateNumbers(bit - 2, remain - (i+j) * bit, (parity + i) % 2);
		}
    }

	oddOnlyStateNumbers[bit / 2][remain] = ret;
	return ret;
}

snum StateCounter::StateNumbers(int s, int bit, int remain, int b)
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

	snum ret = stateNumbers[s][bit][remain][b];
	if (ret >= 0) return ret;
	int a = remain / bit;

	ret = 0;
	for (int i = 0; i <= a; i++)
    {
		for (int j = 0; i + j <= a; j++)
		{
			snum n1 = BinomialCoefficient(SingleTrace(bit, s) << (s - 1), (i64)i);
			snum n2 = BinomialCoefficient((SingleTrace(bit, s) << (s - 1)) - 1 + j, (i64)j);
			ret += n1 * n2 * StateNumbers(s, bit - 1, remain - (i+j) * bit, (b+i)%2);
		}
    }

	stateNumbers[s][bit][remain][b] = ret;
	return ret;
}

void StateCounter::Test(int s)
{
    StateCounter* inst = StateCounter::Inst();

    cout << "Number of states of spin-" << s << " system:" << endl;
    cout << "M\tSingle-Trace\t\tMulti-Trace" << endl;
    for (int i = 1; i * s <= StateCounter::MAX_BIT_TO_COUNT; i++)
	{
		cout << i << "\t" << inst->SingleTrace(i, s) << "\t\t" << inst->MultiTrace(i, s) << endl;
	}

    //cout << "statecounter" << endl;
    for (int i = 2; i <= StateCounter::MAX_BIT_TO_COUNT; i+=2)
    {
        //cout << i << ": " << inst->NoHalfMode(i) << endl;
    }

    //cout << "average energy of each bit:" << endl;    
    for (int i = 1; i <= StateCounter::MAX_BIT_TO_COUNT; i++)
    {
        //cout << i << ": " << setprecision(10) << inst->AverageEnergy(i) << endl;
    }
}
