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
	InitNoHalfMode();
}

void StateCounter::InitNoHalfMode()
{
	noHalfMode.resize(StateGenerator::MAX_BIT_TO_COUNT / 2 + 1, 0);
	for (int M = 2; M <= StateGenerator::MAX_BIT_TO_COUNT; M += 2)
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
