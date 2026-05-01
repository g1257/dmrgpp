#ifndef EXACTDIAG_H
#define EXACTDIAG_H
#include "InputCheck.h"
#include "InputNg.h"
#include "Matrix.h"
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename InternalProductType> class ExactDiag {

public:

	typedef typename InternalProductType::ModelType              ModelType;
	typedef typename ModelType::ComplexOrRealType                ComplexOrRealType;
	typedef typename ModelType::SparseMatrixType                 SparseMatrixType;
	typedef typename ModelType::RealType                         RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type          VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::InputNg<InputCheck>                      InputNgType;

	enum class ParamType
	{
		TEMPERATURE,
		BETA
	};

	ExactDiag(InputNgType::Readable& io, const ModelType&, InternalProductType& h)
	    : tbWhat_(ParamType::TEMPERATURE)
	    , tbTotal_(0)
	    , tbStart_(0)
	    , tbStep_(0)
	{
		PsimagLite::String tmp;
		io.readline(tmp, "TemperatureOrBeta=");
		if (tmp != "beta" && tmp != "temperature")
			err("TemperatureOrBeta= can only be beta or temperature\n");
		tbWhat_ = (tmp == "beta") ? ParamType::BETA : ParamType::TEMPERATURE;

		io.readline(tbStart_, "TemperatureOrBetaStart=");
		io.readline(tbTotal_, "TemperatureOrBetaTotal=");
		io.readline(tbStep_, "TemperatureOrBetaStep=");

		eigenvectors_.resize(h.rows(), h.rows());
		eigenvalues_.resize(h.rows());
		h.fullDiag(eigenvalues_, eigenvectors_);
	}

	void printEnergiesVsTemperature(std::ostream& os) const
	{
		os << "#tb=" << tbWhat_ << "\n";
		os << "#Parameter Energy\n";
		for (SizeType i = 0; i < tbTotal_; ++i) {
			RealType tb = tbStart_ + i * tbStep_;
			RealType e  = calculateOneParam(tb);
			os << tb << " " << e << "\n";
		}
	}

private:

	friend std::ostream& operator<<(std::ostream& os, const ParamType& x)
	{
		PsimagLite::String str = (x == ParamType::BETA) ? "beta" : "temperature";
		os << str;
		return os;
	}

	RealType calculateOneParam(RealType tb) const
	{
		RealType       numer = 0;
		RealType       denom = 0;
		const SizeType total = eigenvalues_.size();
		for (SizeType n = 0; n < total; ++n) {
			RealType weight = weightAt(n, tb);
			numer += eigenvalues_[n] * weight;
			denom += weight;
		}

		return numer / denom;
	}

	RealType weightAt(SizeType n, RealType tb) const
	{
		return (tbWhat_ == ParamType::BETA) ? exp(-tb * eigenvalues_[n])
		                                    : exp(-eigenvalues_[n] / tb);
	}

	ParamType      tbWhat_;
	SizeType       tbTotal_;
	RealType       tbStart_;
	RealType       tbStep_;
	MatrixType     eigenvectors_;
	VectorRealType eigenvalues_;
};
}
#endif // EXACTDIAG_H
