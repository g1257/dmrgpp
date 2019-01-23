#ifndef OBSERVEDRIVER1_H
#define OBSERVEDRIVER1_H

#include "ObserveDriver.h"

namespace Dmrg {
template<typename VectorWithOffsetType,
         typename ModelType>
bool observeOneFullSweep(IoInputType& io,
                         const ModelType& model,
                         const PsimagLite::String& list,
                         SizeType orbitals)
{
	typedef typename ModelType::GeometryType GeometryType;
	typedef Observer<VectorWithOffsetType,ModelType,IoInputType> ObserverType;
	typedef ObservableLibrary<ObserverType> ObservableLibraryType;

	static SizeType start = 0;

	const GeometryType& geometry = model.geometry();
	bool verbose = false;
	SizeType n  = geometry.numberOfSites();
	SizeType rows = n; // could be n/2 if there's enough symmetry
	SizeType cols = n;
	SizeType nf = n - 2;
	SizeType trail = 0;
	SizeType end = start + nf;

	PsimagLite::Vector<PsimagLite::String>::Type vecOptions;
	PsimagLite::split(vecOptions, list, ",");
	bool hasTrail = false;

	for (SizeType i = 0; i < vecOptions.size(); ++i) {
		PsimagLite::String item = vecOptions[i];

		PsimagLite::String label = "%nf=";
		std::size_t labelIndex = item.find(label);
		if (labelIndex == 0) {
			nf = atoi(item.substr(label.length()).c_str());
			rows = nf;
			cols = nf;
			std::cerr<<"observe: Found "<<label<<" = "<<nf;
			std::cerr<<" (rows and cols also set to it)\n";
		}

		label = "%trail=";
		labelIndex = item.find(label);
		if (labelIndex == 0) {
			trail = atoi(item.substr(label.length()).c_str());
			std::cerr<<"observe: Found "<<label<<" = "<<trail<<"\n";
			hasTrail = true;
		}

		label = "%rows=";
		labelIndex = item.find(label);
		if (labelIndex == 0) {
			std::cerr<<"observe: Found %rows= with index "<<labelIndex<<"\n";
			rows = atoi(item.substr(label.length()).c_str());
		}

		label = "%cols=";
		labelIndex = item.find(label);
		if (labelIndex == 0) {
			std::cerr<<"observe: Found %cols= with index "<<labelIndex<<"\n";
			cols = atoi(item.substr(label.length()).c_str());
		}
	}

	if (!hasTrail)
		trail = n - 2 - nf;

	ObservableLibraryType observerLib(io,
	                                  n,
	                                  model,
	                                  start,
	                                  nf,
	                                  trail,
	                                  verbose);

	for (SizeType i = 0; i < vecOptions.size(); ++i) {
		PsimagLite::String item = vecOptions[i];

		if (item.find("%") == 0) continue;

		if (item.length() > 0 && item[0] != '<')
			observerLib.measure(item, rows, cols, orbitals);
		else
			observerLib.interpret(item, rows, cols);
	}

	start = end;
	return observerLib.endOfData();
}
}
#endif // OBSERVEDRIVER1_H
