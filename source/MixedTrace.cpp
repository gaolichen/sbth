#include "MixedTrace.h"

MixedTrace::MixedTrace()
{
	left = 0;
	right = 0;
	//..
}

MixedTrace::MixedTrace(int left, int leftBits, int right, int rightBits)
{
	this->left = BuildMask(left, leftBits);
	this->right = BuildMask(right, rightBits);
}

int MixedTrace::LeftBits() const
{
	return left >> MAX_TRACE_BITS;
}

int MixedTrace::RightBits() const
{
	return right >> MAX_TRACE_BITS;
}

int MixedTrace::LeftTrace() const
{
	return left & FULL_TRACE_BITS;
}

int MixedTrace::RightTrace() const
{
	return right & FULL_TRACE_BITS;
}

void MixedTrace::Contract(int index, SingleTrace& left, MixedTrace& mid, SingleTrace& right)
{
	;
}
