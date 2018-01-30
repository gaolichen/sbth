#pragma warning(disable:4018)
#include "ScriptGenerator.h"
#include "Hamiltonian.h"
#include "BitUtility.h"
#include "NormCalculator.h"
#include <fstream>
using namespace std;

ScriptGenerator::ScriptGenerator(string rootFolder, StateType stateType)
{
	if (stateType == Boson)
	{
		this->rootFolder = CombinePath(rootFolder, "boson");
	}
	else
	{
		this->rootFolder = CombinePath(rootFolder, "fermion");
	}

	this->type = stateType;
#ifdef HAM_PARAMETER
	InitHamiltonians();
#endif
}

#ifdef HAM_PARAMETER
void ScriptGenerator::InitHamiltonians()
{
	h0.AddReadOp(new HamOperatorA(HamOperator::AA, HamOperator::AA), 2);
	h0.AddReadOp(new HamOperatorA(HamOperator::BB, HamOperator::BB), -2);
	h0.AddReadOp(new HamOperatorA(HamOperator::AB, HamOperator::BA), 2);
	h0.AddReadOp(new HamOperatorA(HamOperator::BA, HamOperator::BA), 2);
	h0.AddReadOp(new HamOperatorA(HamOperator::AB, HamOperator::AB), 2);
	h0.AddReadOp(new HamOperatorA(HamOperator::BA, HamOperator::AB), -2);

	h0.AddImaginaryOp(new HamOperatorA(HamOperator::BB, HamOperator::AA), -2);
	h0.AddImaginaryOp(new HamOperatorA(HamOperator::AA, HamOperator::BB), 2);

	hPrime.AddReadOp(new HamOperatorB(HamOperator::AA, HamOperator::AA), 1);
	hPrime.AddReadOp(new HamOperatorB(HamOperator::BA, HamOperator::BA), 1);
	hPrime.AddReadOp(new HamOperatorB(HamOperator::AB, HamOperator::BA), -1);

	deltaH.AddReadOp(new HamOperatorA(HamOperator::AA, HamOperator::AA), 2);
	deltaH.AddReadOp(new HamOperatorA(HamOperator::BB, HamOperator::BB), 2);
	deltaH.AddReadOp(new HamOperatorA(HamOperator::AB, HamOperator::BA), 2);
	deltaH.AddReadOp(new HamOperatorA(HamOperator::BA, HamOperator::AB), 2);
	deltaH.AddReadOp(new BitNumberHamOperator(), -2);
}

void ScriptGenerator::OutputHamToMatlab(int bits, bool inverted)
{
	string function;
	string sign;
	if (!inverted)
	{
		sign = "";
		function = "tham" + ToString(bits);
	}
	else
	{
		sign = "-";
		function = "tpham" + ToString(bits);
	}

	string file = CombinePath(rootFolder, function + ".m");
	cout << file << endl;
	ofstream ofs(file.c_str());
	ofs << "function f = " << function << "(N, xi)" << endl;
	ofs << "A = ";
	HamToMatlab(bits, h0, ofs);
	ofs << ";" << endl << endl;
	
	ofs << "B = xi * ";
	HamToMatlab(bits, deltaH, ofs);
	ofs << ";" << endl;

	ofs << "if " << sign;
	ofs << "(2 * xi + 2) >= 0" << endl;
	ofs << "  f = " << sign << "(A + B);" << endl;
	ofs << "else" << endl;
	ofs << "  f = " << sign << "(A + B)";
	if (inverted)
	{
		ofs << " + (2 * xi + 2) *";
	}
	else
	{
		ofs << " - (2 * xi + 2) *";
	}

	HamToMatlab(bits, hPrime, ofs);
	ofs << ";" << endl;
	ofs << "end" << endl;
	
	ofs << "end" << endl;
	ofs.close();
}

#endif

void ScriptGenerator::VectorToMatlab(vector<int>& a, ofstream& ofs)
{
	ofs << "[";
	for (int i = 0; i < a.size(); i++)
	{
		if (i > 0 && i % 20 == 0)
		{
			ofs << "..." << endl;
			ofs << "\t";
		}
		
		ofs << a[i] + 1 << " ";
	}
	
	ofs << "];" << endl;
}

void ScriptGenerator::HamToMatlab(int bits, Hamiltonian& ham, ofstream& os, bool sparse)
{
        vector<vector<Coefficient> > rem, imm;
        ham.Matrix(bits, type, rem, imm);
	if (sparse)
	{
		vector<int> a, b;
		for (int j = 0; j < rem.size(); j++)
		{
			for (int k = 0; k < rem[j].size(); k++)
			{
				if (imm[j][k].IsZero() && rem[j][k].IsZero()) continue;
				a.push_back(j);
				b.push_back(k);
			}
		}
		
		os << "a = ";
		VectorToMatlab(a, os);
		os << "b = ";
		VectorToMatlab(b, os);
		os << endl;
		os << "v = [";
		for (int i = 0; i < a.size(); i++)
		{
			if (i != 0 && i % 20 == 0)
			{
				os << "..." << endl;
				os << "\t";
			}

			int x = a[i];
			int y = b[i];
			if (!imm[x][y].IsZero())
			{
				os << imm[x][y] << "*1i ";
			}
			else
			{
				os << rem[x][y] << ' ';
			}
		}
		os << "];" << endl;
		
	}
	else
	{
	        os << "[" << endl;
	        for (int j = 0; j < rem.size(); j++)
	        {
	                os << "   ";
	                for (int k = 0; k < rem[j].size(); k++)
	                {
	                        os << ' ';
	                        if (!imm[j][k].IsZero())
	                        {
	                                os << imm[j][k] << "*1i";
	                        }
	                        else
	                        {
	                                os << rem[j][k];
	                        }
        	        }

	                if (j != rem.size() - 1)
	                {
	                        os << ";" << endl;
                	}
        	}

        	os << "]";
	}
}


void ScriptGenerator::OutputHamToMatlab(int bits, Hamiltonian& ham, bool sparse)
{
	string function;
	function = ham.FilePrefix() + ToString(bits);

	string file = CombinePath(rootFolder, function + ".m");
	cout << file << endl;
	ofstream ofs(file.c_str());
	ofs << "function f = " << function << "(N)" << endl;
	if (sparse)
	{
		HamToMatlab(bits, ham, ofs, sparse);
		ofs << "f = sparse(a, b, v);" << endl;
		ofs << "end" << endl;
	}
	else
	{
		ofs << "f=";
		HamToMatlab(bits, ham, ofs);
		ofs << ";" << endl;
		ofs << "end" << endl;
	}

	ofs.close();
}

void ScriptGenerator::SparseMatrixToDataFile(vector<vector<Coefficient> > &mat, string file)
{
	string file1 = CombinePath(rootFolder, file + "-0.dat");
	string file2 = CombinePath(rootFolder, file + "-1.dat");

	cout << file1 << endl << file2 << endl;
        ofstream ofs1(file1.c_str());
	ofstream ofs2(file2.c_str());

        for (int j = 0; j < mat.size(); j++)
        {
                for (int k = 0; k < mat[j].size(); k++)
                {
			if (mat[j][k].N != 0)
			{
				ofs1 << j + 1 << ' ' << k + 1 << ' ' << mat[j][k].N << endl;
			}

                        if (mat[j][k].One != 0)
                        {
                                ofs2 << j + 1 << ' ' << k + 1 << ' ' << mat[j][k].One << endl;
                        }
               }

        }
        ofs1 << mat.size() <<' ' << mat.size() << " 0" << endl;
        ofs2 << mat.size() <<' ' << mat.size() << " 0" << endl;
        ofs1.close();
	ofs2.close();

}

void ScriptGenerator::OutputHamToDataFile(int bits, Hamiltonian& ham, bool isComplex)
{
	string filePrefix = ToString(bits) + ham.FilePrefix();
	
        vector<vector<Coefficient> > rem, imm;
        ham.Matrix(bits, type, rem, imm);
	
	SparseMatrixToDataFile(rem, filePrefix + "-re");

	if (isComplex)
	{
		SparseMatrixToDataFile(imm, filePrefix + "-im");
	}
}

void ScriptGenerator::OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename)
{
	ofstream ofs(filename.c_str());
	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;
	ofs << "\\usepackage{graphicx}" << endl;
	ofs << "\\usepackage{geometry}" << endl;
	ofs << "\\geometry{verbose,tmargin=0.1cm,bmargin=0.1cm,lmargin=0.1cm,rmargin=0.1cm}" << endl;
	ofs << "\\usepackage{pdflscape}" << endl;
	//ofs << "\\usepackage{lscape}" << endl;
	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs <<"\\begin{document}" << endl << endl;
	ofs << "\\newcommand\\scalemath[2]{\\scalebox{#1}{\\mbox{\\ensuremath{\\displaystyle #2}}}}";
	ofs << endl << endl;
	ofs << "\\begin{landscape}" << endl;

	Hamiltonian ham;
	
	for (int i = 0; i < bits.size(); i++)
	{
		ofs << "\\subsection*{"<< bits[i]<< " bits:}" << endl << endl;
		ofs << "\\[" << endl;

		vector<vector<Coefficient> > rem, imm;
		ham.Matrix(bits[i], type, rem, imm);
		
		string s(rem.size(), 'c');
		ofs << "\\left(\\scalemath{" << scales[i] <<"}" << endl;
		ofs << "{\\begin{array}{"<< string(rem.size(), 'c') << "}" << endl;
		for (int j = 0; j < rem.size(); j++)
		{
			for (int k = 0; k < rem[j].size(); k++)
			{
				if (k > 0)
				{
					ofs << "& ";
				}
				if (!imm[j][k].IsZero())
				{
					ofs << imm[j][k].ToLatex() << "i ";
				}
				else
				{
					ofs << rem[j][k].ToLatex() << " ";
				}
			}

			if (j != rem.size() - 1)
			{
				ofs << "\\\\";
			}
			ofs << endl;
		}

		ofs << "\\end{array}}\\right)" << endl;
		ofs << "\\]" << endl << endl << endl;
	}

	ofs << "\\end{landscape}" << endl;

	ofs << "\\end{document}" << endl;
	ofs.close();
}

void ScriptGenerator::GroupStates(int bits, vector<vector<int> >& states)
{
	StateCollection* sc = StateCollection::Inst();
	int n = sc->StateNumber(bits);
	int b = 0;
	if (type == Fermion)
	{
		b = 1;
	}

	vector<int> v;
	if (bits < 5)
	{
		for (int i = 0; i < n; i++)
		{
			v.push_back(i);
		}
		states.push_back(v);
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			if (sc->GetState(bits, j, type).FermionNumber() == b)
			{
				v.push_back(j);
			}
			else
			{
				b = sc->GetState(bits, j, type).FermionNumber();
				states.push_back(v);
				v.clear();
				v.push_back(j);
			}
		}

		states.push_back(v);
	}
}

double ScriptGenerator::NormMatrixScale(int bits, int index)
{
	map<int, vector<double> > scales;
	double a[4] = {1.0, 0.6, 0.79, 1.0};
	double b[4] = {0.51, 0.25, 0.26, 1.0};
	scales[6] = vector<double>(a, a + 4);
	scales[7] = vector<double>(b, b + 4);

	if (scales.find(bits) == scales.end())
	{
		return 1.0;
	}
	else
	{
		return scales[bits][index];
	}
}

void ScriptGenerator::OutputNormToLaTeX(int minBits, int maxBits, string filename)
{
	StateCollection* sc = StateCollection::Inst();
	BruteForceCalculator calc;

	ofstream ofs(filename.c_str());

	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;
	ofs << "\\usepackage{graphicx}" << endl;
	ofs << "\\usepackage{geometry}" << endl;
	ofs << "\\geometry{verbose,tmargin=0.1cm,bmargin=0.1cm,lmargin=0.1cm,rmargin=0.1cm}" << endl;
	ofs << "\\usepackage{pdflscape}" << endl;
	//ofs << "\\usepackage{lscape}" << endl;
	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs <<"\\begin{document}" << endl << endl;
	ofs << "\\newcommand\\scalemath[2]{\\scalebox{#1}{\\mbox{\\ensuremath{\\displaystyle #2}}}}";
	ofs << endl << endl;
	ofs << "\\begin{landscape}" << endl;

	ofs << "\\title{States norm matrices}" << endl << endl;
	ofs << "\\maketitle" << endl << endl;

	for (int i = minBits; i <= maxBits; i++)
	{
		ofs << "\\subsection*{" << i << " bits:}" << endl;
		vector<vector<int> > states;
		GroupStates(i, states);
		if (states.size() > 1)
		{
			ofs << "\\begin{itemize}" << endl;
		}

		for (int j = 0; j < states.size(); j++)
		{
			if (states.size() > 1)
			{
				ofs << "\\item State " << states[j][0] + 1 << " to State " << states[j].back() + 1 << ": " << endl; 
			}
			
			ofs << "\\[" << endl;
			ofs << "\\left(\\scalemath{" << NormMatrixScale(i, j) <<"}" << endl;
			ofs << "{\\begin{array}{" << string(states[j].size(), 'c') << "}" << endl;

			for (int k = 0; k < states[j].size(); k++)
			{
				for (int h = 0; h < states[j].size(); h++)
				{
					if (h > 0)
					{
						ofs << " & "; 
					}
					ofs << calc.Calculate(sc->GetState(i, states[j][k], type), sc->GetState(i, states[j][h], type)).ToLaTeX(i);
				}

				if (k + 1 != states[j].size())
				{
					ofs << "\\\\";
				}

				ofs << endl;
			}

			ofs << "\\end{array}}\\right)" << endl;
			ofs << "\\]" << endl;
		}

		if (states.size() > 1)
		{
			ofs << "\\end{itemize}" << endl;
		}
	}

	ofs << "\\end{landscape}" << endl;
	ofs << "\\end{document}" << endl;
	ofs.close();
}

const TraceState& ScriptGenerator::GetState(int m, int index)
{
	if (type == Boson)
	{
		return StateCollection::Inst()->GetBosonState(m, index);
	}
	else
	{
		return StateCollection::Inst()->GetFermionState(m, index);
	}
}

void ScriptGenerator::OutputStateToLaTeX(int minBit, int maxBit, string filename, bool compress)
{
	ofstream ofs(filename.c_str());
	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;

	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs << "\\usepackage{amsmath}" << endl;
	ofs << "\\allowdisplaybreaks" << endl;
	ofs <<"\\begin{document}" << endl;

	for (int i = minBit; i <= maxBit; i++)
	{
		ofs << "\\subsection*{" << i << " bits}" << endl << endl;
		
		if (type == Boson)
		{
			ofs << StateCollection::Inst()->StateNumber(i) << " bosonic states:" << endl;
		}
		else
		{
			ofs << StateCollection::Inst()->StateNumber(i) << " fermionic states:" << endl;
		}

		if (compress)
		{
			int columns = 3;
			if (StateCollection::Inst()->StateNumber(i) < columns)
			{
				columns = StateCollection::Inst()->StateNumber(i);
			}

			string cc = string(columns, 'l');
			ofs << "\\begin{align}" << endl;
			
			for (int j = 0; j < StateCollection::Inst()->StateNumber(i); j++)
			{
				if (j % columns > 0)
				{
					ofs << " & ";
				}
				if (j > 0 && j % columns == 0)
				{
					ofs << " \\nonumber\\\\" << endl;
				}
				
				ofs << "\\left|" <<  j + 1 << "\\right\\rangle &=";	
				ofs << GetState(i, j).ToLaTeX();
			}
			
			ofs << " \\nonumber" << endl;
			ofs << "\\end{align}" << endl;
		}
		else
		{
			ofs << "\\begin{itemize}" << endl;
			for (int j = 0; j < StateCollection::Inst()->StateNumber(i); j++)
			{
				ofs << "\\item $\\left|" << j + 1 << "\\right\\rangle=";
				ofs << GetState(i, j).ToLaTeX() <<"$" <<  endl;
			}
			ofs << "\\end{itemize}" << endl;
		}
	}
	ofs << "\\end{document}" << endl;

	ofs.close();
}

void ScriptGenerator::OutputNormToMatlab(int bits)
{
	StateCollection* sc = StateCollection::Inst();
	BruteForceCalculator calc;
	vector<vector<int> > states;
	GroupStates(bits, states);

	bool output = (bits >= 10);
	string function = "norm" + ToString(bits);
	string file = CombinePath(rootFolder, function + ".m");
	ofstream ofs(file.c_str());
	ofs << "function f = " << function << "(N)" << endl;
	for (int i = 0; i < states.size(); i++)
	{
		if (output)
		{
			cout << "Calculating state " << i + 1 << "..." << endl;
		}

		ofs << (char)((int)'A' + i) << "=[" << endl;
		for (int k = 0; k < states[i].size(); k++)
		{
			ofs << "    ";
			for (int h = 0; h < states[i].size(); h++)
			{
				if (h > 0)
				{
					ofs << ", ";
				}
				ofs << calc.Calculate(sc->GetState(bits, states[i][k], type), sc->GetState(bits, states[i][h], type)).ToMatlab(bits);
			}

			if (k + 1 != states[i].size())
			{
				ofs << ";" << endl;
			}
		}

		ofs << endl << "];" << endl << endl;
	}

	if (states.size() == 1)
	{
		ofs << "f=A;" << endl;
	}
	else
	{
		ofs << "f=blkdiag(";
		for (int i = 0; i < states.size(); i++)
		{
			if (i > 0)
			{
				ofs << ", ";
			}

			ofs << (char)((int)'A' + i);
		}

		ofs << ");" << endl;
	}
	
	ofs << "end" << endl;
}

void ScriptGenerator::OutputNormToDataFile(int bits, int fNumber)
{
	StateCollection* sc = StateCollection::Inst();
        BruteForceCalculator calc;
        vector<int > stateIndex;
	int tot = sc->StateNumber(bits);
	
	for (int i = 0; i < tot; i++)
	{
		TraceState state = sc->GetState(bits, i, type);
		if (state.FermionNumber() == fNumber)
		{
			stateIndex.push_back(i);
		}
	}

	if (stateIndex.size() == 0)
	{
		cout << "no states found for bits= " << bits << " fermion number=" << fNumber << endl;
		return;
	}

	string name = "norm" + ToString(bits) + "_" + ToString(fNumber);
	string file = CombinePath(rootFolder, name + ".dat");
	//cout << "file =" << file << endl;
	ofstream ofs(file.c_str());

	ofs << stateIndex.size() << endl << endl;
	
	for (int i = 0; i < stateIndex.size(); i++)
	{
		TraceState state = GetState(bits, stateIndex[i]);
		for (int j = i; j < stateIndex.size(); j++)
		{
			Polynomial poly = calc.Calculate(state, GetState(bits, stateIndex[j]));
			if (poly.IsZero())
			{
				ofs << 0 << endl;
			}
			else
			{
				int maxPow = poly.MaxPow();
				int offset = (bits - maxPow) % 2;
				if (offset == 0) { ofs << 2; }
				else { ofs << 1; }

				for (int k = bits - offset; k >=0; k -= 2)
				{
					ofs << ' ' << poly.GetCoef(k);
				}
				ofs << endl;
			}
		}

		ofs << endl;
	}

}

void ScriptGenerator::OutputStateStructure(int maxBit)
{	
	StateCollection* sc = StateCollection::Inst();
	string filePrefix = "state_structure_";
	
	for (int i = 3; i<= maxBit; i++)
	{
		string file = CombinePath(rootFolder, filePrefix + ToString(i) + ".txt");
		cout << "output state structure to " << file << endl;
		ofstream os(file.c_str());
		os << "index\tfermion\ttrace" << endl;
		
		for (int j = 0; j < sc->StateNumber(i); j++)
		{
			TraceState state = sc->GetState(i, j, type);
			os << j + 1 << "\t" << state.FermionNumber() << "\t" << state.TraceNumber() << endl;
		}

		os.close();
	}
}
