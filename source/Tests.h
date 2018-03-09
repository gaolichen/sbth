// This is the class for test code. 
//
#pragma once
#include <exception>
#include "BitUtility.h"
#include "StateCounter.h"
#include "EigenEnergyLargeN.h"


class TestError : public exception
{
private:
    string message;

public:
    TestError(string message)
    {
        this->message = message;
    }

    const char* what() const throw()
    {
        return message.c_str();
    }

    ~TestError() throw() {}
};

template<class T> bool VectorEqual(const vector<T>& a, const vector<T>& b)
{
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++)
    {
        if (!(a[i] == b[i]))
        {
            return false;
        }
    }

    return true;
};

/*template<class T> void Assert(const T& a, const T& b, string message)
{
    if (a == b) return;
    throw TestError(message + "\n\t" + "returned " + ToString(a) + ". Expected: " + ToString(b));
}*/


class Tests
{
private:
    int s;

    string datapath;

    static vector<TE> ReadEnergyFile(string file);

public:
	Tests(int s = 1);

    static void Assert(bool condition, string message)
    {
        if (!condition)
        {
            throw TestError(message);
        }
    }

    static bool EnergyListEqual(const vector<TE>& a, const vector<TE>& b)
    {
        if (a.size() != b.size()) return false;
        for (int i = 0; i < a.size(); i++)
        {
            if ((abs(a[i].E-b[i].E) < EPS || abs(a[i].E-b[i].E) <= 1e-5 * abs(a[i].E)) && a[i].Deg == b[i].Deg)
            {
                continue;
            }

            cout <<a[i].E << ' ' << b[i].E << ' ' << a[i].E - b[i].E << endl;

            return false;
        }

        return true;
    };

    // test functions
    void SingleTraceEnergies_StateNumber(int M);
    void MultiTraceEigenEnergies_StateNumber(int maxM);

    void SingleTraceEnergies_States(int maxM);
    void MultiTraceEigenEnergies_States(int maxM);
    

    // demo functions
    void DemoStateNumber();

    void DemoAverageEnergy();

    static void Run();
};

