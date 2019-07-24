#ifndef PRINTERINDETAIL_H
#define PRINTERINDETAIL_H
#include <iostream>

namespace Dmrg {

template<typename LeftRightSuperType>
class PrinterInDetail {

	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::QnType QnType;

public:

	PrinterInDetail(const LeftRightSuperType& lrs, bool extended)
	    : lrs_(lrs), extended_(extended)
	{}

	void print(std::ostream& os, PsimagLite::String msg) const
	{
		if (!extended_) return;
		printOneSide(os, "left", lrs_.left());
		printOneSide(os, "right", lrs_.right());
	}

private:

	void printOneSide(std::ostream& os,
	                  PsimagLite::String msg,
	                  const BasisWithOperatorsType& basis) const
	{
		SizeType sites = basis.block().size();
		os<<"Side="<<msg<<"\n";
		os<<"SitesOnThisSide ";
		for (SizeType i = 0; i < sites; ++i) {
			os<<basis.block()[i]<<" ";
		}

		os<<"\n";

		SizeType n = basis.partition();
		os<<"Partitions "<<n<<"\n";
		for (SizeType i = 0; i < n - 1; ++i) {
			SizeType s = basis.partition(i + 1) - basis.partition(i);
			const typename BasisWithOperatorsType::QnType& j = basis.qnEx(i);
			os<<j<<" "<<s<<"\n";
		}

		assert(sites > 0);
		SizeType site = basis.block()[sites - 1];
		SizeType end = basis.operatorsPerSite(0);
		SizeType siteC = site;
		if (msg == "right") {
			assert(site >= basis.block()[0]);
			siteC = site - basis.block()[0];
		}

		os<<"Operators at site "<<site<<" ("<<siteC<<")\n";
		for (SizeType sigma = 0; sigma < end; ++sigma) {
			typename BasisWithOperatorsType::PairType p = basis.getOperatorIndices(siteC, sigma);
			os<<sigma<<" non-zeroes="<<basis.getOperatorByIndex(p.first).getStorage().nonZeros();
			os<<" rows="<<basis.getOperatorByIndex(p.first).getStorage().rows()<<"\n";
		}
	}

	const LeftRightSuperType& lrs_;
	bool extended_;
}; // class PrinterInDetail
}
#endif // PRINTERINDETAIL_H
