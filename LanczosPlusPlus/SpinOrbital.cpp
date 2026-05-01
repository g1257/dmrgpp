#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

template <typename ComplexOrRealType> class Hamiltonian {

public:

	typedef std::vector<ComplexOrRealType>           VectorType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType>    MatrixType;

	Hamiltonian(int twiceJ, int nsites)
	    : twiceJ_(twiceJ)
	    , nsites_(nsites)
	    , sPerSite_(nsites)
	    , lPerSite_(nsites)
	{
		statesS_ = pow(twiceJ + 1, nsites);
		statesL_ = statesS_;
	}

	void fill(SparseMatrixType& sparse, MatrixType* dense = nullptr)
	{
		const int  total = statesS_ * statesL_;
		VectorType rowVector(total);
		sparse.resize(total, total);

		int counter = 0;
		if (dense)
			dense->resize(total, total);

		for (int idL = 0; idL < statesL_; ++idL) {
			for (int idS = 0; idS < statesS_; ++idS) {
				int row = packSandL(idS, idL);
				sparse.setRow(row, counter);
				fillMatrixRow(rowVector, idS, idL);
				for (int col = 0; col < total; ++col) {
					const ComplexOrRealType value = rowVector[col];
					if (value == 0)
						continue;

					sparse.pushCol(col);
					sparse.pushValue(value);
					if (dense)
						dense->operator()(row, col) = value;
					rowVector[col] = 0;
					++counter;
				}
			}
		}

		sparse.setRow(total, counter);
		sparse.checkValidity();
		if (dense && !isHermitian(*dense, true)) {
			throw std::runtime_error("Dense is not Hermitian\n");
		}

		if (!isHermitian(sparse, true)) {
			throw std::runtime_error("H is not Hermitian\n");
		}
	}

private:

	void fillMatrixRow(VectorType& rowVector, const int idS, const int idL)
	{
		indexToVector(sPerSite_, idS);
		indexToVector(lPerSite_, idL);

		for (int i = 0; i < nsites_; ++i) {
			for (int x = 0; x < 1; ++x) {
				const int j = (x == 0) ? i + 1 : i - 1;
				if (j < 0 || j >= nsites_)
					continue;
				for (int which0 = 0; which0 < 3; ++which0) { // +- -+ or zz
					for (int which1 = 0; which1 < 3; ++which1) { // +- -+ or zz
						ComplexOrRealType valueS = 0;
						int               braS
						    = createOneTerm(valueS, idS, i, j, which0, 0);
						if (braS < 0)
							continue;

						ComplexOrRealType valueL = 0;
						int               braL
						    = createOneTerm(valueL, idL, i, j, which1, 1);
						if (braL < 0)
							continue;
						int col = packSandL(braS, braL);

						rowVector[col] += valueS * valueL;
					}
				}
			}
		}
	}

	int createOneTerm(ComplexOrRealType& value, int id, int i, int j, int which, int what)
	{
		const double      jValue = twiceJ_ * 0.5;
		std::vector<int>& v      = (what == 0) ? sPerSite_ : lPerSite_;
		int               ret    = -1;
		switch (which) {
		case 0:
			value = 0.5 * factorSpSm(jValue, v[i] - jValue);
			ret   = createSpSm(i, v[i], j, v[j], v);
			// assert(tret < 0 || value == 1);
			break;
		case 1:
			value = 0.5 * factorSpSm(jValue, v[j] - jValue);
			ret   = createSmSp(i, v[i], j, v[j], v);
			// assert(ret < 0 || value == 1);
			break;
		case 2:
			value = (v[i] - jValue) * (v[j] - jValue);
			ret   = id;
			break;
		default:
			throw std::runtime_error("createOneTerm\n");
		}

		return ret;
	}

	double factorSpSm(double j, double m) const
	{
		const double value = j * (j + 1) - m * (m + 1);
		return value;
	}

	int createSpSm(int i, int mi, int j, int mj, std::vector<int>& v) const
	{
		if (mi == twiceJ_)
			return -1;
		if (mj == 0)
			return -1;
		return createOneState(i, mi + 1, j, mj - 1, v);
	}

	int createSmSp(int i, int mi, int j, int mj, std::vector<int>& v) const
	{
		if (mj == twiceJ_)
			return -1;
		if (mi == 0)
			return -1;
		return createOneState(i, mi - 1, j, mj + 1, v);
	}

	int createOneState(int ind, int mind, int jnd, int mjnd, std::vector<int>& v) const
	{
		const int savedInd = v[ind];
		v[ind]             = mind;
		const int savedJnd = v[jnd];
		v[jnd]             = mjnd;
		const int braIndex = vectorToIndex(v);
		v[ind]             = savedInd;
		v[jnd]             = savedJnd;
		return braIndex;
	}

	int packSandL(int idS, int idL) const
	{
		assert(idS < statesS_);
		return idS + idL * statesS_;
	}

	void indexToVector(std::vector<int>& v, int ind) const
	{
		assert(ind < statesS_);
		const int three = twiceJ_ + 1;
		int       n     = v.size();
		int       tmp   = ind;
		for (int i = 0; i < n; ++i) {
			v[i] = tmp % three;
			tmp -= v[i];
			tmp /= three;
		}
	}

	int vectorToIndex(const std::vector<int>& v) const
	{
		const int three  = twiceJ_ + 1;
		int       n      = v.size();
		int       tmp    = 0;
		int       factor = 1;
		for (int i = 0; i < n; ++i) {
			assert(v[i] < three);
			tmp += v[i] * factor;
			factor *= three;
		}

		return tmp;
	}

	int              twiceJ_;
	int              statesS_;
	int              statesL_;
	int              nsites_;
	std::vector<int> sPerSite_;
	std::vector<int> lPerSite_;
};

template <typename ComplexOrRealType> void solveMatrix(PsimagLite::Matrix<ComplexOrRealType>& dense)
{
	const int           n = dense.rows();
	std::vector<double> e(n);
	diag(dense, e, 'N');
	std::cout << "LAPACK energy=" << e[0] << "\n";
}

template <typename ComplexOrRealType>
void solveMatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& sparse)
{
	typedef PsimagLite::ParametersForSolver<ComplexOrRealType>   SolverParametersType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	SolverParametersType params;
	params.lotaMemory = true;
	params.options    = "reortho";
	params.tolerance  = -1;

	PsimagLite::LanczosSolver<SolverParametersType,
	                          PsimagLite::CrsMatrix<ComplexOrRealType>,
	                          VectorType>
	    lanczosSolver(sparse, params);

	const int  n = sparse.rows();
	double     e = 0;
	VectorType z(n, 0.0);
	VectorType initial(n);
	PsimagLite::fillRandom(initial);
	lanczosSolver.computeOneState(e, z, initial, 0);

	std::cout << "Lanczos energy=" << e << "\n";
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "USAGE: " << argv[0] << " nsites [twiceJ]\n";
		return 1;
	}

	int       nsites = atoi(argv[1]);
	const int twiceJ = (argc > 2) ? atoi(argv[2]) : 2;

	typedef Hamiltonian<double>               HamiltonianType;
	typedef HamiltonianType::SparseMatrixType SparseMatrixType;
	typedef HamiltonianType::MatrixType       MatrixType;

	HamiltonianType  h(twiceJ, nsites);
	SparseMatrixType sparse;
	MatrixType       dense;

	h.fill(sparse, &dense);

	solveMatrix(sparse);

	solveMatrix(dense);
}
