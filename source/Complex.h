#pragma once

template<class T> class Complex
{
private:
	T re;
	T im;
public:
	Complex(T real);
	Complex(T real, T im);
	T Re();
	T Im();

	Complex<T>& operator+= (const Complex<T>&);
	Complex<T>& operator-= (const Complex<T>&);
	Complex<T>& operator*= (const Complex<T>&);

	friend Complex<T> operator+ (const Complex<T>&, const Complex<T>&);
	friend Complex<T> operator- (const Complex<T>&, const Complex<T>&);
	friend Complex<T> operator* (const Complex<T>&, const Complex<T>&);
};