
#include "EigenEnergyLargeN.h"

void EigenEnergyLargeN::Partition(int n, int k, i64 mask)
{
	if (n <= 0)
	{
		// finish a partition, calculate the energies.
	        masks.push_back(mask);
		return;
	}

	if (n > k * (k + 1) / 2)
	{
		// n is too large to partition for the remaining k integers.
		return;
	};
	
	if (k == 0)
	{
		// should not reach here...
		return;
	}

	if (n >= k)
	{
		Partition(n - k, k - 1, mask | (1<<k));
	}
	
	Partition(n, k - 1, mask);
}

// TODO:
double EigenEnergyLargeN::CalcEnergy(i64 mask, int bits)
{
	return 0.0;
}

void EigenEnergyLargeN::Calculate()
{
	int M2 = bits / 2;
	int half = 0;
	if (bits % 2 == 0) half = M2;
	for (int n = 0; n * bits + half <= M2 * (1 + M2) / 2; n++)
	{
		Partition(n * bits + half, M2, (i64)0);
	}

	for (int i = 0; i <= masks.size(); i++)
	{
		int zeroBits = BitCount(masks[i]);
		int count = M2 - zeroBits;
		if (bits % 2 == 1) count++;

		for (int j = 0; j < (1 << count); j++)
		{
			energies.push_back(make_pair(CalcEnergy(masks[i], j), 1 << zeroBits));
		}
	}
}
