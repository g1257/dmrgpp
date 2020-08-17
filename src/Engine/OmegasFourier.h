#ifndef OMEGASFOURIER_H
#define OMEGASFOURIER_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ComplexOrRealType, typename InputNgType>
class OmegasFourier {

	class OmegasGeometry {

	public:

		OmegasGeometry(typename InputNgType::Readable& io) : subname_("NONE")
		{
			io.readline(name_, "GeometryName=");
			try {
			io.readline(subname_, "GeometrySubname=");
			} catch (std::exception&) {}
		}

		PsimagLite::String name() const
		{
			return name_;
		}

		PsimagLite::String subname() const
		{
			return subname_;
		}

	private:

		PsimagLite::String name_;
		PsimagLite::String subname_;
	};

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef ProgramGlobals ProgramGlobalsType;
	//typedef PsimagLite::Geometry<ComplexOrRealType, InputNgType, ProgramGlobalsType> GeometryType;

	static const SizeType M_MAX = 0;

	OmegasFourier(bool skipFourier, typename InputNgType::Readable& io)
	    : skipFourier_(skipFourier), geometry_(io)
	{
		if (skipFourier_) return;

		io.readline(numberOfSites_, "TotalNumberOfSites=");
		int tmp = 0;
		io.readline(tmp, "IsPeriodicX=");
		isPeriodicX_ = (tmp > 0);
		VectorSizeType v;
		io.read(v, "TSPSites");
		if (v.size() != 1)
			err("TSPSites must be a vector of exactly one entry\n");
		centralSite_ = v[0];
		io.readline(orbitals_, "Orbitals=");

		qValues_.resize(numberOfSites_);
	}

	void fourier(const VectorRealType& values1, const VectorRealType& values2)
	{
		if (skipFourier_) return;

		//use qvalues_;
		PsimagLite::String subname = geometry_.subname();
		if (subname=="average") {
			return fourierLadderAverage(values1, values2);
		}

		PsimagLite::String name = geometry_.name();
		if (name == "chain") {
			return fourierChain(values1, values2);
		}

		if (name == "ladder" && subname == "GrandCanonical") {
			return fourierChainGC(values1, values2);
		}

		if (name == "ladder") {
			return fourierLadder(values1, values2);
		}

		if (name == "LongRange" || name == "General") {

			if (subname == "chain") {
				if (orbitals_ == 1) {
					return fourierChain(values1, values2);
				} else if (orbitals_ == 2) {
					return fourierChain2orb(values1, values2);
				}
			}

			if (subname == "ladder") {
				if (orbitals_ == 1) {
					return fourierLadder(values1, values2);
				} else if (orbitals_ == 2) {
					return fourierLadder2orb(values1, values2);
				}
			}

			const PsimagLite::String honey = "HoneyComb";
			if (subname.find(honey) == 0) {
				PsimagLite::String type = subname.substr(honey.length(),
				                                         subname.size() - honey.length());
				return fourierHoneycomb(values1, values2, type);
			}
		}

		err("OmegasFourier: undefined geometry " + name  + "\n");
	}

private:

	void fourierChain(const VectorRealType& values1, const VectorRealType& values2)
	{
		const SizeType nOverTwo = static_cast<SizeType>(numberOfSites_/2);
		if (!isPeriodicX_) {
			bool b = (centralSite_ == nOverTwo);
			if (!b && (centralSite_ != nOverTwo - 1)) {
				err("Chain of " + ttos(numberOfSites_) + "sites, but central site is " +
				    ttos(centralSite_) + ", makes no sense!?\n");
			}
		}

		const SizeType numberOfQs = (M_MAX > 0) ? M_MAX : numberOfSites_;
		if (qValues_.size() < numberOfQs)
			err("INTERNAL ERROR at fourierChain\n");

		for (SizeType m = 0; m < numberOfQs; ++m) {
			ComplexType sum = 0;
			RealType q = getQ(m, numberOfQs, isPeriodicX_);
			for (SizeType i = 0; i < numberOfSites_; ++i) {
				RealType arg = q*(i - centralSite_);
				RealType carg = cos(arg);
				RealType sarg = sin(q*(i + 1))*sin(q*(centralSite_ + 1));
				RealType cOrSarg = (isPeriodicX_) ? carg : sarg;
				sum += ComplexType(values1[i]*cOrSarg, values2[i]*cOrSarg);
			}

			assert(m < qValues_.size());
			qValues_[m] = sum;
		}
	}

	void fourierLadderAverage(const VectorRealType& values1, const VectorRealType& values2)
	{
		err("unimplemented fourierLadderAverage\n");
	}

	void fourierChainGC(const VectorRealType& values1, const VectorRealType& values2)
	{
		err("unimplemented fourierLadderAverage\n");
	}

	void fourierLadder(const VectorRealType& values1, const VectorRealType& values2)
	{
		err("unimplemented fourierLadderAverage\n");
	}

	void fourierChain2orb(const VectorRealType& values1, const VectorRealType& values2)
	{
		err("unimplemented fourierLadderAverage\n");
	}

	void fourierLadder2orb(const VectorRealType& values1, const VectorRealType& values2)
	{
		err("unimplemented fourierLadder2orb\n");
	}

	void fourierHoneycomb(const VectorRealType& values1,
	                      const VectorRealType& values2,
	                      PsimagLite::String type)
	{
		err("unimplemented fourierLadderAverage\n");
	}

	static RealType getQ(SizeType m, SizeType n, bool isPeriodic)
	{
		return (isPeriodic) ? 2.0*M_PI*m/n : (m + 1.0)*M_PI/(n+1.0);
	}

	bool skipFourier_;
	OmegasGeometry geometry_;
	SizeType numberOfSites_;
	SizeType centralSite_;
	SizeType orbitals_;
	bool isPeriodicX_;
	VectorComplexType qValues_;
};
}
#endif // OMEGASFOURIER_H
