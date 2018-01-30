#include "BitUtility.h"

using namespace std;

class EigenEnergyLargeN
{
private:
	int bits;
	vector<i64> masks;

	// list of eigen energies: for each element, the first component is the value of energy, 
	// and the second component is the degeneracy of the enerties.
	vector<pair<double, int> > energies;
	
	// calculate energy from two numbers.
	// nonZeroBits: the nonzero bit position of mask represents a nonzero mode, the mode can make 
	// either one unit position contribution or one unit negative contribution to the energy, 
	// where each unit mean sin(k * Pi/M);
	// bits: the nonzero bit position of bits represents a positive-contribution bit and zero 
	// bit position represents a negative-contribution bit.
	double CalcEnergy(vector<int>& nonZeroBits, int positiveBits);

	void Partition(int n, int k, i64 mask);

public:
	EigenEnergyLargeN(int bits_)
	{
		this->bits = bits_;
	}
	
	void Calculate();

	const vector<pair<double, int> > Energies() { return energies; }
};

