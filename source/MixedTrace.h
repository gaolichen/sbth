#pragma once
#include "BitUtility.h"
#include "SingleTrace.h"

class MixedTrace
{
private:
	int left;
	int right;
public:
	MixedTrace();
	inline int LeftBits() const;
	inline int RightBits() const;
	inline int LeftTrace() const;
	inline int RightTrace() const;
	MixedTrace(int left, int leftBits, int right, int rightBits);
	void Contract(int, SingleTrace&, MixedTrace&, SingleTrace&);
};