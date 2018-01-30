#include"StateId.h"

StateId::StateId()
	: BitNumber(-1), Index(-1)
{
}

StateId::StateId(int bitNumber, int index)
	: BitNumber(bitNumber), Index(index)
{
	// ...
}

bool StateId::IsValid()
{
	return BitNumber >= 0 && Index >= 0;
}

bool operator< (const StateId& a, const StateId& b)
{
	if (a.BitNumber != b.BitNumber)
	{
		return a.BitNumber < b.BitNumber;
	}

	return a.Index < b.Index;
}

ostream& operator<<(ostream& os, const StateId& id)
{
	os << (id.Index / 2) + 1;
	return os;
}