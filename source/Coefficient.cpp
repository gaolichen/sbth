#include "Coefficient.h"
#include <sstream>
#include <cmath>
using namespace std;

Coefficient::Coefficient()
	:N(0), One(0), InvN(0)
{
}

Coefficient::Coefficient(int n, int one, int invN)
	:N(n), One(one), InvN(invN)
{
}

Coefficient Coefficient::operator+ (const Coefficient& a) const
{
	return Coefficient(N + a.N, One + a.One, InvN + a.InvN);
}

Coefficient& Coefficient::operator+= (const Coefficient& a)
{
	this->N += a.N;
	this->One += a.One;
	this->InvN += a.InvN;

	return *this;
}

void Coefficient::Opposite()
{
	N = -N;
	One = -One;
	InvN = -InvN;
}

void Coefficient::ChangeOrder(int order2Change)
{
	if (order2Change == 0) return;
	if (order2Change > 1 || order2Change < -1)
	{
		cout << "Only 1 or -1 order change is allowd." << endl;
		return;
	}
	
	if ((InvN != 0 && order2Change == -1) || (N != 0 && order2Change == 1))
	{
		cout << "ChangeOrder(): Unexpected: InvN or N should be 0.";
	}

	if (order2Change == -1)
	{
		InvN = One;
		One = N;
		N = 0;
	}
	else
	{
		N = One;
		One = InvN;
		InvN = 0;
	}
}

bool Coefficient::IsZero() const
{
	return N == 0 && One == 0 && InvN == 0;
}

Coefficient& Coefficient::operator*= (int n)
{
	N *= n;
	One *= n;
	InvN *= n;

	return *this;
}

Coefficient operator* (const Coefficient& a, int n)
{
	return Coefficient(a.N * n, a.One * n, a.InvN * n);
}

ostream& operator<<(ostream& os, const Coefficient& coef)
{
	if (coef.IsZero())
	{
		os << "0";
		return os;
	}

	int v[3] = {coef.N, coef.One, coef.InvN};
	bool isFirst = true;
	for (int i = 0; i < 3; i++)
	{
		if (v[i] == 0) continue;
		if (isFirst)
		{
			os << v[i];
			isFirst = false;
		}
		else
		{
			if (v[i] > 0)
			{
				os << "+";
			}
			os << v[i];
		}

		if (i == 1)
		{
			os << "/N";
		}
		else if (i == 2)
		{
			os << "/N^2";
		}
	}

    return os;
}

string Coefficient::ToLatex() const
{
	if (IsZero())
	{
		return "0";
	}

	int v[3] = {N, One, InvN};
	bool isFirst = false;

	ostringstream oss;
	for (int i = 0; i < 3; i++)
	{
		if (v[i] == 0) continue;
		if (v[i] < 0)
		{
			oss << "-";
		}
		else if (!isFirst)
		{
			oss << "+";
		}

		isFirst = false;
		if (i == 0)
		{
			oss << abs(v[i]);
		}
		else
		{
			oss << "\\frac{" << abs(v[i]) << "}{N";
			if (i == 2) oss << "^2}";
			else oss << "}";
		}
	}

	return oss.str();
}
