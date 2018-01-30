
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
double EigenEnergyLargeN::CalcEnergy(vector<int>& modes, int positiveBits)
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
		int nonZeroBits = M2 - zeroBits;
		int nonZeroMask = ((1 << M2) - 1) ^ masks[i];
		vector<int> nonZeroDigits = Num2Digit(nonZeroMask, M2);

		for (int j = 0; j < (1 << nonZeroBits); j++)
		{
			energies.push_back(make_pair(CalcEnergy(nonZeroDigits, j), 1 << zeroBits));
		}
	}
}
