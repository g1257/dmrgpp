#ifndef LANCZOSDRIVER_1_H
#define LANCZOSDRIVER_1_H
#include "Io/IoSimple.h"
#include "LabeledOperator.h"
#include "LanczosDriver.h"
#include "LanczosGlobals.h"

template <typename ModelType> SizeType maxOrbitals(const ModelType& model)
{
	SizeType res = 0;
	for (SizeType i = 0; i < model.geometry().numberOfSites(); i++) {
		if (res < model.orbitals(i))
			res = model.orbitals(i);
	}
	return res;
}

template <typename EngineType>
void extendedStatic(PsimagLite::String                   manypoint,
                    const EngineType&                    engine,
                    const typename EngineType::PairType& braAndKet)
{
	typedef typename EngineType::VectorSizeType  VectorSizeType;
	PsimagLite::Vector<PsimagLite::String>::Type str;
	PsimagLite::split(str, manypoint, ";");

	VectorSizeType                                             sites;
	VectorSizeType                                             spins;
	VectorSizeType                                             orbs;
	PsimagLite::Vector<LanczosPlusPlus::LabeledOperator>::Type whats;
	for (SizeType i = 0; i < str.size(); ++i) {
		PsimagLite::Vector<PsimagLite::String>::Type str2;
		PsimagLite::split(str2, str[i], "?");
		if (str2.size() < 3)
			throw PsimagLite::RuntimeError("-M option malformed\n");
		whats.push_back(LanczosPlusPlus::LabeledOperator(str2[0]));
		sites.push_back(atoi(str2[1].c_str()));
		spins.push_back(atoi(str2[2].c_str()));
		if (str2.size() == 4)
			orbs.push_back(atoi(str2[3].c_str()));
		else
			orbs.push_back(0);
	}

	std::cout << "<gs|" << manypoint << "|gs>=";
	std::cout << engine.manyPoint(sites, whats, spins, orbs, braAndKet);
	std::cout << "\n";
}

template <typename ModelType,
          typename SpecialSymmetryType,
          template <typename, typename> class InternalProductTemplate>
void mainLoop3(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType         GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef LanczosPlusPlus::Engine<ModelType, InternalProductTemplate, SpecialSymmetryType>
	    EngineType;
	using ContinuedFractionCollectionType = PsimagLite::ContinuedFractionCollection<RealType>;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	const GeometryType& geometry = model.geometry();
	EngineType          engine(model, io);

	//! get the g.s.:
	RealType Eg = engine.energies(0);
	std::cout.precision(8);
	std::cout << "Energy=" << Eg << "\n";
	PsimagLite::String filename = PsimagLite::basenameOf(io.filename());

	const SizeType nmeas = lanczosOptions.measure.size();
	for (SizeType i = 0; i < nmeas; ++i) {
		VectorStringType tokens;
		PsimagLite::split(tokens, lanczosOptions.measure[i], ",");
		const SizeType ntokens = tokens.size();
		for (SizeType j = 0; j < ntokens; ++j) {
			VectorStringType braOpKet;
			PsimagLite::split(braOpKet, tokens[j], "|");
			engine.measure(braOpKet);
		}
	}

	bool needsDos = false;
	try {
		int tmp = 0;
		io.readline(tmp, "ComputeDensityOfStates=");
		needsDos = (tmp > 0);
	} catch (std::exception&) { }

	typedef std::pair<SizeType, SizeType>  PairSizeType;
	PsimagLite::Vector<PairSizeType>::Type pairOfSites;
	const SizeType                         n = geometry.numberOfSites();

	if (needsDos) {
		lanczosOptions.gf.push_back(LanczosPlusPlus::LabeledOperator("c"));
		for (SizeType i = 0; i < n; ++i)
			pairOfSites.push_back(PairSizeType(i, i));
	}

	try {
		io.read(lanczosOptions.sites, "TSPSites");

		if (lanczosOptions.sites.size() == 0)
			err("TSPSites must have at least one site\n");

		if (lanczosOptions.sites.size() == 1)
			lanczosOptions.sites.push_back(lanczosOptions.sites[0]);

		pairOfSites.push_back(
		    PairSizeType(lanczosOptions.sites[0], lanczosOptions.sites[1]));
	} catch (std::exception&) { }

	bool     hasCenter  = false;
	SizeType centerSite = 0;
	try {
		io.readline(centerSite, "TSPCenter=");
		std::cout << "TSPCenter=" << centerSite << "\n";

		for (SizeType i = 0; i < n; ++i)
			pairOfSites.push_back(PairSizeType(centerSite, i));
		hasCenter = true;
	} catch (std::exception&) { }

	bool doAllPairs = false;
	try {
		int tmp = 0;
		io.readline(tmp, "DoAllPairs=");
		doAllPairs = (tmp > 0);
	} catch (std::exception&) { }

	if (doAllPairs && hasCenter)
		err("You cannot have both TSPCenter and DoAllPairs\n");

	if (doAllPairs) {
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = 0; j < n; ++j)
				pairOfSites.push_back(PairSizeType(i, j));
	}

	for (SizeType gfi = 0; gfi < lanczosOptions.gf.size(); ++gfi) {
		SizeType       counter  = 0;
		const SizeType nIndices = pairOfSites.size();
		for (SizeType sIndex = 0; sIndex < nIndices; ++sIndex) {
			const SizeType site0 = pairOfSites[sIndex].first;
			const SizeType site1 = pairOfSites[sIndex].second;

			std::cout << "#gf(i=" << site0 << ", j=" << site1 << ")\n";

			typename EngineType::VectorStringType vstr;
			PsimagLite::IoSimple::Out ioOut(filename + ttos(counter) + ".comb");

			ioOut.write(site0, "Site0");
			ioOut.write(site1, "Site1");

			if (hasCenter)
				ioOut.write(centerSite, "TSPCenter");

			ContinuedFractionCollectionType cfCollection(PsimagLite::FreqEnum::REAL);
			SizeType                        norbitals = maxOrbitals(model);
			for (SizeType orb1 = 0; orb1 < norbitals; orb1++) {
				for (SizeType orb2 = orb1; orb2 < norbitals; orb2++) {
					engine.spectralFunction(
					    cfCollection,
					    vstr,
					    lanczosOptions.gf[gfi],
					    site0,
					    site1,
					    lanczosOptions.spins,
					    std::pair<SizeType, SizeType>(orb1, orb2));
				}
			}

			ioOut << "#INDEXTOCF ";
			for (SizeType i = 0; i < vstr.size(); ++i)
				ioOut << vstr[i] << " ";
			ioOut << "\n";
			cfCollection.write(ioOut);
			std::cerr << "LanczosDriver1.h: Written to " << ioOut.filename() << "\n";
			++counter;
		}
	}

	for (SizeType cicji = 0; cicji < lanczosOptions.cicj.size(); cicji++) {
		SizeType                              total = geometry.numberOfSites();
		PsimagLite::Matrix<ComplexOrRealType> cicjMatrix(total, total);
		SizeType                              norbitals = maxOrbitals(model);
		for (SizeType orb1 = 0; orb1 < norbitals; orb1++) {
			for (SizeType orb2 = 0; orb2 < norbitals; orb2++) {
				engine.twoPoint(cicjMatrix,
				                lanczosOptions.cicj[cicji],
				                lanczosOptions.spins,
				                std::pair<SizeType, SizeType>(orb1, orb2),
				                std::pair<SizeType, SizeType>(0, 0));
				std::cout << cicjMatrix;
			}
		}
	}

	if (lanczosOptions.split >= 0) {
		LanczosPlusPlus::ReducedDensityMatrix<ModelType> reducedDensityMatrix(
		    model, engine.eigenvector(0), lanczosOptions.split);
		reducedDensityMatrix.printAll(std::cout);
	}

	if (lanczosOptions.extendedStatic != "") {
		PsimagLite::Vector<PsimagLite::String>::Type str;
		PsimagLite::split(str, lanczosOptions.extendedStatic, ",");
		for (SizeType i = 0; i < str.size(); ++i)
			extendedStatic(str[i], engine, typename EngineType::PairType(0, 0));
	}
}

template <typename ModelType, typename SpecialSymmetryType>
void mainLoop2(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	PsimagLite::String tmp;
	io.readline(tmp, "SolverOptions=");
	bool onthefly = (tmp.find("InternalProductOnTheFly") != PsimagLite::String::npos);

	if (onthefly) {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductOnTheFly>(
		    model, io, lanczosOptions);
	} else {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductStored>(
		    model, io, lanczosOptions);
	}
}

#endif // LANCZOSDRIVER_1_H
