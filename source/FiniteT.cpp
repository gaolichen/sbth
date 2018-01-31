#include "FiniteT.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

void FiniteT::CalcEigenStates(double invN)
{
	H0Hamiltonian h0;
	Eigen::MatrixXcd mat = h0.ToMatrixN(bits, Boson, invN);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(mat);
}

