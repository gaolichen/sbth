#include "Complex.h"

template<class T> Complex<T>::Complex(T real)
{
	Complex(real, (T)0);
}

template<class T> Complex<T>::Complex(T real, T imagine)
	: re(real), im(imagine)
{
}

template<class T> T Complex<T>::Re()
{
	return re;
}

template<class T> T Complex<T>::Im()
{
	return im;
}

template<class T> Complex<T>& Complex<T>::operator+= (const Complex<T>& a)
{
	re += a.re;
	im += a.im;

	return *this;
}

template<class T> Complex<T>& Complex<T>::operator-= (const Complex<T>& a)
{
	re -= a.re;
	im -= a.im;
	return *this;
}

template<class T> Complex<T>& Complex<T>::operator*= (const Complex<T>& a)
{
	int real = a.re * re - a.im * im;
	int imagine = a.im * re + a.re * im;
	re = real;
	im = image;
	return *this;
}

template<class T> Complex<T> operator+ (const Complex<T>& a, const Complex<T>& b)
{
	return Complex<T>(a.re + b.re, a.im + b.im);
}

template<class T> Complex<T> operator- (const Complex<T>& a, const Complex<T>& b)
{
	return Complex<T>(a.re - b.re, a.im - b.im);
}

template<class T> Complex<T> operator* (const Complex<T>& a, const Complex<T>& b)
{
	return Complex<T>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
