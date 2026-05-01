#include "ExactDiag.h"
#include "LanczosDriver.h"

PsimagLite::String license = "Copyright (c) 2009-2020, UT-Battelle, LLC\n"
                             "All rights reserved\n"
                             "\n"
                             "[Lanczos++, Version 1.]\n"
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

template <typename ModelType> void mainLoop(InputNgType::Readable& io, const ModelType& model)
{
	typedef typename ModelType::GeometryType                              GeometryType;
	typedef typename ModelType::BasisBaseType                             BasisBaseType;
	typedef LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType> DefaultSymmetryType;
	typedef InternalProductStored<ModelType, DefaultSymmetryType>         InternalProductType;
	typedef ExactDiag<InternalProductType>                                ExactDiagType;

	DefaultSymmetryType rs(model.basis(), model.geometry(), "");
	InternalProductType hamiltonian(model, rs);
	ExactDiagType       exactDiag(io, model, hamiltonian);

	exactDiag.printEnergiesVsTemperature(std::cout);
}

template <typename ComplexOrRealType> void mainLoop0(InputNgType::Readable& io)
{
	typedef PsimagLite::
	    Geometry<ComplexOrRealType, InputNgType::Readable, LanczosPlusPlus::ProgramGlobals>
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
	mainLoop(io, modelPtr);
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp                           application("lanczos++ED", &argc, &argv, 1);
	int                                          opt  = 0;
	PsimagLite::String                           file = "";
	PsimagLite::Vector<PsimagLite::String>::Type str;
	InputCheck                                   inputCheck;
	int                                          precision        = 10;
	bool                                         versionOnly      = false;
	SizeType                                     threadsInCmdLine = 0;

	/* PSIDOC EdDriver
	\begin{itemize}
	\item[-f file] Input file to use. DMRG++ inputs can be used.
	\item[-S number] Override Threads= in input line if preset, and set threads
	to number given here.
	\item[-p precision] precision in decimals to use.
	\item[-V] prints version and exits.
	\end{itemize}
	*/
	while ((opt = getopt(argc, argv, "f:p:S:V")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
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
		std::cerr << "Lanczos++ ED Version " << LANCZOSPP_VERSION << "\n";
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
		mainLoop0<ComplexType>(io);
	else
		mainLoop0<RealType>(io);
}
