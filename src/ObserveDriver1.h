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
	typedef typename ObservableLibraryType::ManyPointActionType ManyPointActionType;
	typedef typename PsimagLite::OneOperatorSpec::SiteSplit SiteSplitType;

	static SizeType start = 0;

	const GeometryType& geometry = model.geometry();
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
	                                  trail);

	ManyPointActionType* manyPointAction = new ManyPointActionType(false, "");
	for (SizeType i = 0; i < vecOptions.size(); ++i) {
		PsimagLite::String item = vecOptions[i];

		if (item.find("%") == 0) continue;

		SiteSplitType braceContent = PsimagLite::OneOperatorSpec::extractSiteIfAny(item,
		                                                                           '{',
		                                                                           '}');

		PsimagLite::String actionString;
		if (braceContent.hasSiteString) {
			PsimagLite::String actionContent = braceContent.siteString;
			size_t index = actionContent.find("action=");

			if (index == PsimagLite::String::npos)
				err("Only action=something accepted for brace option, not "
				    + braceContent.siteString + "\n");

			actionString = actionContent.substr(index + 7, actionContent.length() - 7);
			delete manyPointAction;
			manyPointAction = new ManyPointActionType(braceContent.hasSiteString, actionString);
			if (braceContent.root != "")
				err("Garbage trailing after brace close\n");
			continue;
		}

		if (item.length() > 0 && item[0] != '<')
			observerLib.measure(item, rows, cols, *manyPointAction, orbitals);
		else
			observerLib.interpret(item, rows, cols, *manyPointAction);
	}

	start = end;
	return observerLib.endOfData();
}
}
#endif // OBSERVEDRIVER1_H
