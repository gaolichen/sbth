#pragma once
#include <iostream>
using namespace std;

struct Coefficient
{
	int N;
	int One;
	int InvN;
	Coefficient();
	Coefficient(int n, int one, int invN = 0);
	Coefficient& operator+= (const Coefficient& a);
	Coefficient operator+ (const Coefficient& a) const;
	void Opposite();
	void ChangeOrder(int order2Change);
	bool IsZero() const;
	string ToLatex() const;
	friend ostream& operator<<(ostream& os, const Coefficient& coef);
	Coefficient& operator*= (int n);
};

Coefficient operator* (const Coefficient& a, int n);

