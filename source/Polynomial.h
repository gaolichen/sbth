#pragma once
#include <vector>
#include <iostream>

using namespace std;

class Polynomial
{
private:
	vector<int> coef;
public:
	Polynomial();

	bool IsZero() const;
	int MaxPow() const;
	int GetCoef(int order) const;
	void IncreaseOrder(int pow);
	void Increase(int order, int value);
	Polynomial& operator+= (const Polynomial& poly);
	Polynomial& operator*= (const Polynomial& poly);
	Polynomial& operator+= (int value);
	Polynomial& operator*= (int value);
	string ToLaTeX(int orderToDeduct) const;
	string ToMatlab(int orderToDeduct) const;

	friend Polynomial operator+ (const Polynomial& a, const Polynomial& b);
	friend Polynomial operator* (const Polynomial& a, const Polynomial& b);
	friend ostream& operator<<(ostream& os, const Polynomial& poly);
};

