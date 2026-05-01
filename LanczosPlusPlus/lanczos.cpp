#include "LanczosDriver.h"

PsimagLite::String license = "Copyright (c) 2009-2012, UT-Battelle, LLC\n"
                             "All rights reserved\n"
                             "\n"
                             "[Lanczos++, Version 1.0]\n"
                             "\n"
                             "-------------------------------------------------------------\n"
                             "THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
                             "CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
                             "WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
                             "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
                             "PARTICULAR PURPOSE ARE DISCLAIMED. \n"
                             "\n"
                             "Please see full open source license included in file LICENSE.\n"
                             "-------------------------------------------------------------\n"
                             "\n";

using namespace LanczosPlusPlus;

void fillOrbsOrSpin(PsimagLite::Vector<LanczosPlusPlus::LanczosOptions::PairSizeType>::Type& spinV,
                    const PsimagLite::Vector<PsimagLite::String>::Type&                      strV)
{
	for (SizeType i = 0; i < strV.size(); i++) {
		PsimagLite::Vector<PsimagLite::String>::Type strV2;
		PsimagLite::split(strV2, strV[i], ",");
		if (strV2.size() != 2)
			throw std::runtime_error("-o needs pairs\n");
		LanczosPlusPlus::LanczosOptions::PairSizeType spins;
		spins.first  = atoi(strV2[0].c_str());
		spins.second = atoi(strV2[1].c_str());
		spinV.push_back(spins);
	}
}

template <typename ModelType>
void mainLoop(InputNgType::Readable&           io,
              const ModelType&                 model,
              LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType  GeometryType;
	typedef typename ModelType::BasisBaseType BasisBaseType;

	int tmp = 0;
	try {
		io.readline(tmp, "UseTranslationSymmetry=");
	} catch (std::exception& e) { }

	bool useTranslationSymmetry = (tmp == 1) ? true : false;

	try {
		io.readline(tmp, "UseReflectionSymmetry=");
	} catch (std::exception& e) { }

	bool useReflectionSymmetry = (tmp == 1) ? true : false;

	if (useTranslationSymmetry) {
		mainLoop2<ModelType,
		          LanczosPlusPlus::TranslationSymmetry<GeometryType, BasisBaseType>>(
		    model, io, lanczosOptions);
	} else if (useReflectionSymmetry) {
		mainLoop2<ModelType,
		          LanczosPlusPlus::ReflectionSymmetry<GeometryType, BasisBaseType>>(
		    model, io, lanczosOptions);
	} else {
		mainLoop2<ModelType, LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>>(
		    model, io, lanczosOptions);
	}
}

template <typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io, LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef PsimagLite::
	    Geometry<ComplexOrRealType, InputNgType::Readable, LanczosPlusPlus::LanczosGlobals>
	        GeometryType;
	typedef LanczosPlusPlus::
	    ModelSelector<ComplexOrRealType, GeometryType, InputNgType::Readable>
	                                                  ModelSelectorType;
	typedef typename ModelSelectorType::ModelBaseType ModelBaseType;

	GeometryType geometry(io);

	std::cout << geometry;

	ModelSelectorType    modelSelector(io, geometry);
	const ModelBaseType& modelPtr = modelSelector();

	std::cout << modelPtr;
	mainLoop(io, modelPtr, lanczosOptions);
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp                           application("lanczos++", &argc, &argv, 1);
	int                                          opt = 0;
	LanczosOptions                               lanczosOptions;
	PsimagLite::String                           file = "";
	PsimagLite::Vector<PsimagLite::String>::Type str;
	InputCheck                                   inputCheck;
	int                                          precision        = 6;
	bool                                         versionOnly      = false;
	SizeType                                     threadsInCmdLine = 0;

	/* PSIDOC LanczosDriver
	\begin{itemize}
	\item[-f file] Input file to use. DMRG++ inputs can be used.
	\item[-c label] Computes the two-point correlation for label.
	\item[-M string] Computes many-point static correlations.
	string = opsec0;opesec1;... is a semicolon-separated list of operator specifications.
	An opsec ( equal to id?site?spin?orb ) is a question-mark separated list of
	operator properties.
	Id is the id of the operator; see operator ids in Lanczos++ below in this manual;
	site, spin, and orbital are 0-based specifications with their respective meanings.
	Orbital is optional, so that opsec = id?site?spin is also valid.
	\item[-S number] Override Threads= in input line if preset, and set threads
	to number given here.
	\item[-g label] Computes the spectral function (continued fraction) for label.
	\item[-s ``s1,s2''] computes correlations or spectral functions for spin s1,s2.
	Only s1==s2 is supported for now.
	\item[-r siteForSplit] Calculates the reduced density matrix with a lattice
	split at the siteForSplit.
	\item[-p precision] precision in decimals to use.
	\item[-V] prints version and exits.
	\end{itemize}
	*/
	while ((opt = getopt(argc, argv, "g:c:m:f:s:r:p:M:S:V")) != -1) {
		switch (opt) {
		case 'g':
			lanczosOptions.gf.push_back(LabeledOperator(optarg));
			break;
		case 'f':
			file = optarg;
			break;
		case 'c':
			lanczosOptions.cicj.push_back(LabeledOperator(optarg));
			break;
		case 'm':
			lanczosOptions.measure.push_back(optarg);
			break;
		case 's':
			lanczosOptions.spins.clear();
			PsimagLite::split(str, optarg, ";");
			fillOrbsOrSpin(lanczosOptions.spins, str);
			str.clear();
			break;
		case 'r':
			lanczosOptions.split = atoi(optarg);
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'M':
			lanczosOptions.extendedStatic = optarg;
			break;
		case 'S':
			threadsInCmdLine = atoi(optarg);
			break;
		case 'V':
			versionOnly = true;
			break;
		default: /* '?' */
			inputCheck.usage(argv[0]);
			return 1;
		}
	}

	if (file == "" && !versionOnly) {
		inputCheck.usage(argv[0]);
		return 1;
	}

	// print license
	if (ConcurrencyType::root()) {
		std::cerr << license;
		std::cerr << "Lanczos++ Version " << LANCZOSPP_VERSION << "\n";
		std::cerr << "PsimagLite version " << PSIMAGLITE_VERSION << "\n";
	}

	if (versionOnly)
		return 0;

	// Setup the Geometry
	InputNgType::Writeable ioWriteable(file, inputCheck);
	InputNgType::Readable  io(ioWriteable);

	bool isComplex     = false;
	bool setAffinities = false;

	PsimagLite::String solverOptions;
	io.readline(solverOptions, "SolverOptions=");
	try {
		int fermionSign = -1;
		io.readline(fermionSign, "FermionSign=");
		std::cerr << "WARNING= FermionSign=" << fermionSign << "\n";
		std::cout << "WARNING= FermionSign=" << fermionSign << "\n";
		LanczosPlusPlus::LanczosGlobals::FERMION_SIGN = fermionSign;
	} catch (std::exception&) { }

	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::split(tokens, solverOptions, ",");
	for (SizeType i = 0; i < tokens.size(); ++i) {
		if (tokens[i] == "useComplex") {
			isComplex = true;
		} else if (tokens[i] == "setAffinities") {
			setAffinities = true;
		}
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	try {
		io.readline(npthreads, "Threads=");
	} catch (std::exception&) { }

	if (threadsInCmdLine > 0)
		npthreads = threadsInCmdLine;

	PsimagLite::CodeSectionParams codeSectionParams(npthreads, 1, setAffinities, 0);
	ConcurrencyType::setOptions(codeSectionParams);

	typedef std::complex<RealType> ComplexType;

	if (isComplex)
		mainLoop0<ComplexType>(io, lanczosOptions);
	else
		mainLoop0<RealType>(io, lanczosOptions);
}
