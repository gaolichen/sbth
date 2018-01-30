#pragma warning(disable:4018)
#include <sstream>
#include <cmath>
#include <algorithm>
#include "Polynomial.h"
using namespace std;

Polynomial::Polynomial()
{
	//..
}

bool Polynomial::IsZero() const
{
	if (this->coef.size() == 0)
	{
		return true;
	}

	for (int i = 0; i < coef.size(); i++)
	{
		if (coef[i] != 0)
		{
			return false;
		}
	}

	return true;
}

int Polynomial::GetCoef(int order) const
{
	if (order >= this->coef.size()) return 0;
	return coef[order];
}

int Polynomial::MaxPow() const
{
	for (int i = this->coef.size() - 1; i >= 0; i--)
	{
		if (this->coef[i] != 0) return i;
	}

	return 0;
}

Polynomial& Polynomial::operator+=(const Polynomial& poly)
{
	if (this->coef.size() < poly.coef.size())
	{
		coef.resize(poly.coef.size());
	}

	for (int i = 0; i < this->coef.size(); i++)
	{
		this->coef[i] += poly.coef[i];
	}

	return *this;
}

Polynomial& Polynomial::operator+=(int value)
{
	if (this->coef.size() == 0)
	{
		this->coef.push_back(value);
	}
	else
	{
		this->coef[0] += value;
	}

	return *this;
}

Polynomial& Polynomial::operator*= (const Polynomial& poly)
{
	if (this->IsZero()) return *this;
	if (poly.IsZero())
	{
		this->coef.resize(0);
		return *this;
	}

	int increasePow = poly.MaxPow();
	if (increasePow > 0)
	{
		this->coef.resize(this->coef.size() + increasePow);
	}

	for (int i = increasePow + this->MaxPow(); i >= 0; i--)
	{
		int res = 0;
		for (int j = max(i - increasePow, 0); j < this->coef.size() && j <= i; j++)
		{
			res += this->coef[j] * poly.coef[i - j];
		}

		this->coef[i] = res;
	}

	return *this;
}

Polynomial& Polynomial::operator*= (int value)
{
	if (value == 0)
	{
		this->coef.resize(0);
	}
	else
	{
		for (int i = 0; i < this->coef.size(); i++)
		{
			this->coef[i] *= value;
		}
	}

	return *this;
}

Polynomial operator+ (const Polynomial& a, const Polynomial& b)
{
	Polynomial ret = a;
	ret += b;

	return ret;
}

Polynomial operator* (const Polynomial& a, const Polynomial& b)
{
	Polynomial ret = a;
	ret *= b;

	return ret;
}

void Polynomial::IncreaseOrder(int pow)
{
	if (pow <= 0) return;

	this->coef.resize(this->MaxPow() + pow + 1);
	for (int i = this->MaxPow() + pow; i >= pow; i--)
	{
		this->coef[i] = this->coef[i - pow];
	}

	for (int i = 0; i < pow; i++)
	{
		this->coef[i] = 0;
	}
}

void Polynomial::Increase(int order, int value)
{
	if (this->coef.size() <= order)
	{
		this->coef.resize(order + 1);
	}

	this->coef[order] += value;
}

ostream& operator << (ostream& os, const Polynomial& poly)
{
	if (poly.IsZero())
	{
		os << 0;
		return os;
	}

	int maxP = poly.MaxPow();
	for (int i = maxP; i >= 0; i--)
	{
		if (poly.GetCoef(i) == 0) continue;
		else if (i != maxP && poly.GetCoef(i) > 0)
		{
			os << "+";
		}

		os << poly.GetCoef(i);
		if (i > 0)
		{
			os << "N^" << i;
		}
	}

	return os;
}

string Polynomial::ToMatlab(int orderToDeduct) const
{
	if (IsZero())
	{
		return "0";
	}

	ostringstream oss;
	int maxP = this->MaxPow();
	for (int i = maxP; i >= 0; i--)
	{
		if (GetCoef(i) == 0) continue;
		else if (GetCoef(i) < 0)
		{
			oss << '-';
		}
		else if (i != maxP)
		{
			oss << '+';
		}

		int c = abs(GetCoef(i));
		if (i - orderToDeduct == 0)
		{
			oss << c;
		}
		else if (i - orderToDeduct > 0)
		{
			oss << c << "N";
			if (i - orderToDeduct > 1)
			{
				oss << '^'<< (i - orderToDeduct);
			}
		}
		else
		{
			
			oss << c << "/N"; 
			if (orderToDeduct - i > 1)
			{
				oss << '^'<< (orderToDeduct - i);
			}
		}
	}

	return oss.str();
}

string Polynomial::ToLaTeX(int orderToDeduct) const
{
	if (IsZero())
	{
		return "0";
	}

	ostringstream oss;
	int maxP = this->MaxPow();
	for (int i = maxP; i >= 0; i--)
	{
		if (GetCoef(i) == 0) continue;
		else if (GetCoef(i) < 0)
		{
			oss << '-';
		}
		else if (i != maxP)
		{
			oss << '+';
		}

		int c = abs(GetCoef(i));
		if (i - orderToDeduct == 0)
		{
			oss << c;
		}
		else if (i - orderToDeduct > 0)
		{
			oss << c << "N";
			if (i - orderToDeduct > 1)
			{
				oss << '^'<< (i - orderToDeduct);
			}
		}
		else
		{
			
			oss << "\\frac{" << c << "}{N"; 
			if (orderToDeduct - i > 1)
			{
				oss << "^{"<< (orderToDeduct - i) << "}";
			}
			oss << "}";
		}
	}

	return oss.str();
}
