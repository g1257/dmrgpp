\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

%\title{The DynamicTargetting Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o DynamicTargetting.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Dynamic Targetting To get Frequency Dependent Observables}
This class implements the dynamic DMRG algorithm as described in \cite{re:jeckelmann02}.
Note that this class implements the \verb=Targetting= interface also
used by GroundStateTargetting (static DMRG) and TimeStepTargetting (for
time dependent DMRG)

Following paper reference~\cite{re:jeckelmann02}, the name \emph{dyn-vectors} will be used for the four
vector: (i) the ground state $|\psi_{gs}\rangle$,
(ii) the vector $A|\psi_{gs}\rangle$, (iii)
the ``correction vector'' $|Y_A\rangle$, and (iv)
the ``correction vector'' $|X_A\rangle$, as defined in the paper.
The last 3 vectors will be stored in the private member \verb|targetVectors_|,
whereas the first one will be stored in the private member \verb|psi_|.

@o DynamicTargetting.h -t
@{
#ifndef DYNAMICTARGETTING_H
#define DYNAMICTARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "TimeSerializer.h"
#include "DynamicDmrgParams.h"
#include "DynamicFunctional.h"
#include "Minimizer.h"

namespace Dmrg {
	@<theClassHere@>
} // namespace
#endif // DYNAMICTARGETTING_H

@}

This class is templated on 7 templates, which are:
\begin{enumerate}
\item \verb|LanczosSolverTemplate|, being usually the \verb|LanczosSolver| class.
\item \verb|InternalProductTemplate|, being usually the \verb|InternalProductOnTheFly| class.
This is a very short class that allows to compute the superblock matrix either on-the-fly or to store it.
(by using \verb|InternalProductStored|). The latter option is limited to small systems due to memory constraints.
\item \verb|WaveFunctionTransformationType| is usually the \verb|WaveFunctionTransformation| class. 
The wave function transformation is too long to explain here but it is a standard computational trick in DMRG, and
was introduce in 1996 by S. White (need to write here the corresponding PRB article FIXME).
\item \verb|ModelType| is the model in question. These are classes under the directory Models.
\item \verb|ConcurrenyType| is the type to deal with parallelization or lack thereof.
\item \verb|IoType| is usually the \verb|IoSimple| class, and deals with writing to disk the dyn-vectors produced by this class.
\item \verb|VectorWithOffsetTemplate| is usually the \verb|VectorWithOffsets| class that encapsulates
the functionality of a vector that is mostly zero except for chunks of non-zero numbers at certain offsets.
Note that there is (for efficiency reasons) a \verb|VectorWithOffset| class that encapsulates the functionality
of a vector with a single chunk. That class is used in \verb|GroundStateTargetting| but not here.
Why do vectors in chunks appear here (you might be wondering)? Well, because of symmetries the vectors are zero mostly
everywhere except on the (targetted) symmetry sector(s).
\end{enumerate}

@d theClassHere
@{
template<
	template<typename,typename,typename> class LanczosSolverTemplate,
	template<typename,typename> class InternalProductTemplate,
	typename WaveFunctionTransformationType_,
	typename ModelType_,
	typename ConcurrencyType_,
	typename IoType_,
	template<typename> class VectorWithOffsetTemplate>
class DynamicTargetting  {
public:
	@<publicTypedefs@>
	@<enumsAndConstants@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class DynamicTargetting
@}

A long series of typedefs follow. Need to explain these maybe (FIXME).
@d publicTypedefs
@{
typedef WaveFunctionTransformationType_ WaveFunctionTransformationType;
typedef ModelType_ ModelType;
typedef ConcurrencyType_ ConcurrencyType;
typedef IoType_ IoType;
typedef typename ModelType::RealType RealType;
typedef std::complex<RealType> ComplexType;
typedef InternalProductTemplate<ComplexType,ModelType> InternalProductType;
typedef typename ModelType::OperatorsType OperatorsType;
//typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
typedef typename ModelType::MyBasisWithOperators BasisWithOperatorsType;
typedef std::vector<ComplexType> ComplexVectorType;
typedef LanczosSolverTemplate<RealType,InternalProductType,ComplexVectorType> LanczosSolverType;
//typedef std::vector<RealType> VectorType;
typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
//typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
typedef typename BasisWithOperatorsType::OperatorType OperatorType;
typedef typename BasisWithOperatorsType::BasisType BasisType;
typedef DynamicDmrgParams<ModelType> TargettingParamsType;
typedef typename BasisType::BlockType BlockType;
typedef VectorWithOffsetTemplate<ComplexType> VectorWithOffsetType;
typedef typename VectorWithOffsetType::VectorType VectorType;
typedef ComplexVectorType TargetVectorType;
typedef BlockMatrix<ComplexType,ComplexMatrixType> ComplexBlockMatrixType;
typedef ApplyOperatorLocal<BasisWithOperatorsType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
typedef TimeSerializer<RealType,VectorWithOffsetType> TimeSerializerType;
@}

And now a few enums and other constants. The first refers to the 3 stages in which the dynamic algorithm can be.
\begin{enumerate}
\item DISABLED
Dynamic targetting is disabled if we are not computing any dyn-dependent operator yet (like when we're in the infinite algorithm)%'
or if the user specified TSTLoops with numbers greater than zero, those numbers indicate the loops that must pass
before dyn-dependent observables are computed.
\item OPERATOR
In this stage we're applying an operator%'
\item WFT\CONVERGING
In this stage we're converging .%'
\end{enumerate}

@d enumsAndConstants
@{
enum {DISABLED,OPERATOR,CONVERGING};
enum {	EXPAND_ENVIRON=WaveFunctionTransformationType::EXPAND_ENVIRON,
		EXPAND_SYSTEM=WaveFunctionTransformationType::EXPAND_SYSTEM,
		INFINITE=WaveFunctionTransformationType::INFINITE};

static const size_t parallelRank_ = 0; // DYNT needs to support concurrency FIXME
@}

Now comes the constructor which takes 6 arguments.
The first 3 arguments are the system (left-block), environment (right-block), and superblock (system + environment).
As usual, the first 2 are heavy objects---with operators---, and the superblock is light.
The 4th argument is the model object. The 5th argument is a \verb|TargettingStructureType| object 
which is a \verb|TargettingStructureParms| object 
A structure is just a bunch of data bundled together, and you can see this in the file \verb|TargetStructureParams.h|.
The last argument is a \verb|WaveFunctionTransformation| object. More info about this class is
in \verb|WaveFunctionTransformation.h|.
@d constructor
@{
DynamicTargetting(
		const BasisWithOperatorsType& basisS,
		const BasisWithOperatorsType& basisE,
		const BasisType& basisSE,
		const ModelType& model,
		const TargettingParamsType& tstStruct,
		const WaveFunctionTransformationType& wft)
:	@<stackInitialization@>
@<constructorBody@>
@}

Now let us look at the private data of this class:
@d privateData
@{
std::vector<size_t> stage_;
VectorWithOffsetType psi_;
const BasisWithOperatorsType& basisS_;
const BasisWithOperatorsType& basisE_;
const BasisType& basisSE_;
const ModelType& model_;
const TargettingParamsType& tstStruct_;
const WaveFunctionTransformationType& waveFunctionTransformation_;
PsimagLite::ProgressIndicator progress_;
RealType currentOmega_;
std::vector<VectorWithOffsetType> targetVectors_;
std::vector<RealType> weight_;
RealType gsWeight_;
typename IoType::Out io_;
ApplyOperatorType applyOpLocal_;
@}

Now we get to the stack initialization of this object.
We said before that the algorithm could be in 4 stages. In reality, there is not a stage for the full
algorithm but there's a stage for each operator to be applied (like holon and then doublon).%'
These operators are specified by the user in the input file in TSPSites. 
All stages are set to DISABLED at the beginning.
We make reference copies to the bases for system (\verb|basisS|), environment (\verb|basisE|), and
superblock (\verb|basisSE|). We also make a reference copy for the model and the tst t(ime)s(tep)t(argetting)Struct(ure),
and of the \verb|waveFunctionTransformation| object.
We initialize the \verb|progress| object that helps with printing progress to the terminal.
The frequency that we are calculating here needs to be described FIXME.
We will think this to do one frequency at a time, as suggested by the reference above.
Multiple frequencies should be parallelized, and we need to provide a frequency range and stepping
in the input file (FIXME)
The \verb|weight|, which is a vector of weights, for each target state (except possibly the ground state)
is set to 4 vectors.
Next an \verb|io| or input/output object is constructed. This is needed to dump the dyn vectors to disk
since we don't do computations \emph{in-situ} here.%'
All right, we may do something \emph{in situ} just to check.
The \verb|applyLocal| operator described before is also initalized on the stack.
@d stackInitialization
@{
stage_(tstStruct.sites.size(),DISABLED),
basisS_(basisS),
basisE_(basisE),
basisSE_(basisSE),
model_(model),
tstStruct_(tstStruct),
waveFunctionTransformation_(wft),
progress_("DynamicTargetting",0),
currentOmega_(tstStruct_.omega),
targetVectors_(3),
weight_(targetVectors_.size()),
io_(tstStruct_.filename,parallelRank_),
applyOpLocal_(basisS,basisE,basisSE)
@}

The body of the constructor follows:
@d constructorBody
@{
{
	if (!wft.isEnabled()) throw std::runtime_error(" DynamicTargetting "
			"needs an enabled wft\n");
	RealType sum = 0;
	size_t n = weight_.size();
	for (size_t i=0;i<n;i++) {
		weight_[i] = 1.0/(n+1);
		sum += weight_[i];
	}

	gsWeight_=1.0-sum;
	sum += gsWeight_;
	if (fabs(sum-1.0)>1e-5) throw std::runtime_error("Weights don't amount to one\n");
	printHeader();
}
@}

@d publicFunctions
@{
@<weight@>
@<gsWeight@>
@<normSquared@>
@<setGs@>
@<operatorBracket@>
@<gs@>
@<includeGroundStage@>
@<size@>
@<operatorParens@>
@<evolve@>
@<initialGuess@>
@<BasisGetFunctions@>
@}

The public member function \verb|weight| returns the weight of target state $i$. This is needed for
the \verb|DensityMatrix| class to be able weight the states properly.
Note that it throws if you ask for weights of dyn-vectors when all stages are disabled, since
this would be an error.%
@d weight
@{
RealType weight(size_t i) const
{
	if (allStages(DISABLED)) throw std::runtime_error("TST: What are you doing here?\n");
	return weight_[i];
	//return 1.0;
}
@}

The public member function \verb|gsWeight| returns the weight of the ground state.
During the disabled stages it is $1$ since there are no other vectors to target. 
@d gsWeight
@{
RealType gsWeight() const
{
	if (allStages(DISABLED)) return 1.0;
	return gsWeight_;
}
@}

This member function returns the squared norm of dynamic vector number $i$.
(Do we really need this function?? (FIXME))
@d normSquared
@{
RealType normSquared(size_t i) const
{
	// call to mult will conjugate one of the vector
	return real(multiply(targetVectors_[i],targetVectors_[i]));
}
@}

The function below sets the ground state to whatever is passed in \verb|v|.
The basis to which this state belongs need be passed in \verb|someBasis| because of
the chunking of this vector.
@d setGs
@{
template<typename SomeBasisType>
void setGs(const std::vector<TargetVectorType>& v,
		const SomeBasisType& someBasis)
{
	psi_.set(v,someBasis);
}
@}

The functions below returns the $i-$th element of the ground state \verb|psi|.
@d operatorBracket
@{
const ComplexType& operator[](size_t i) const { return psi_[i]; }
			
ComplexType& operator[](size_t i) { return psi_[i]; }
@}

The function below returns the full ground state vector as a vector with offset:
@d gs
@{
const VectorWithOffsetType& gs() const { return psi_; }
@}

The function below tells if the ground state will be included in the density matrix.
If using this class it
will always be.
@d includeGroundStage
@{
bool includeGroundStage() const {return true; }
@}

How many dyn vectors does the \verb|DensityMatrix| need to include, excepting
the ground state? The function below tells.
Note that when all stages are disabled no dyn vectors are
included.
@d size
@{
size_t size() const
{
	if (allStages(DISABLED)) return 0;
	return targetVectors_.size();
}
@}

The function below returns the full dyn vector number $i$ as a vector with offsets:
@d operatorParens
@{
const VectorWithOffsetType& operator()(size_t i) const
{
	return targetVectors_[i];
}
@}

This function provides a hook to (possibly) start the computation of
dynamic observables. Five arguments are passed. First $Eg$, the ground state energy,
then the \verb|direction| of expansion (system or environment), then the \verb|block|
being currently grown or shrunk, then the \verb|loopNumber| of the finite algorithm,
and finally a flag \verb|needsPrinting| 
that indicates if dyn-vectors need to be printed to disk for post-processing or not.

Here the main work is done two functions,  a different function \verb|evolve|
is called to either WFT transform the vector or to apply the operators to the ground state.

This function call other functions. We'll continue linearly describing each one %'
in order of appearance.
@d evolve
@{
void evolve(RealType Eg,size_t direction,const BlockType& block,
		size_t loopNumber, bool needsPrinting)
{
	size_t count =0;
	VectorWithOffsetType phiOld = psi_;
	VectorWithOffsetType phiNew;
	size_t max = tstStruct_.sites.size();

	if (noStageIs(DISABLED)) max = 1;

	// Loop over each operator that needs to be applied
	// in turn to the g.s.
	for (size_t i=0;i<max;i++) {
		count += evolve(i,phiNew,phiOld,Eg,direction,block,loopNumber,max-1);
		phiOld = phiNew;
	}

	if (count==0) {
		// always print to keep observer driver in sync
		if (needsPrinting) {
			zeroOutVectors();
			printVectors(block);
		}
		return;
	}

	ComplexType val = calcDynVectors(Eg,phiNew,direction);

	cocoon(val,direction,block); // in-situ

	if (needsPrinting) printVectors(block); // for post-processing
}
@}

The function below provides an initial guess for the Lanczos vector.
Traditionally, when DMRG is only targetting the ground state this is a standard procedure
(see file \verb|GroundStateTargetting|). Here, when \verb|DynamicTargetting|, we need be
concerned with all target states and we the stages of the application of operators.
When all stages are disabled then the initial guess is just delegated to one call of the WFT's%'
\verb|setInitialVector| function.
When stages are advancing we need to weight each target wave-function-transformed  state
with the appropriate weights:
@d initialGuess
@{
void initialGuess(VectorWithOffsetType& v) const
{
	waveFunctionTransformation_.setInitialVector(v,psi_,basisS_,basisE_,basisSE_);
	if (!allStages(CONVERGING)) return;
	std::vector<VectorWithOffsetType> vv(targetVectors_.size());
	for (size_t i=0;i<targetVectors_.size();i++) {
		waveFunctionTransformation_.setInitialVector(vv[i],
				targetVectors_[i],basisS_,basisE_,basisSE_);
		if (norm(vv[i])<1e-6) continue;
		VectorWithOffsetType w= weight_[i]*vv[i];
		v += w;
	}
}
@}

Finally, the following 3 member public functions return the superblock object, the system (left-block)
or the environment objects, or rather, the references held by this class.
@d BasisGetFunctions
@{
const BasisType& basisSE() const { return basisSE_; }

const BasisWithOperatorsType& basisS() const { return basisS_; }

const BasisWithOperatorsType& basisE() const { return basisE_; }
@}

This completes the list of public functions.
What remains are private (i.e. non-exported) code used only by this class.
We'll visit one function at a time. %'

@d privateFunctions
@{
@<evolvePrivate@>
@<computephi@>
@<cocoon@>
@<checkOrderOfSites@>
@<allStages@>
@<noStageIs@>
@<getStage@>
@<calcDynVectors@>
@<minimizeFunctional@>
@<minimizeFunctional2@>
@<obtainXA@>
@<obtainXA2@>
@<guessPhiSectors@>
@<zeroOutVectors@>
@<printVectors@>
@<printHeader@>
@<test@>
@<areAllTargetsSensible@>
@<isThisTargetSensible@>
@}

The below function is called from the \verb|evolve| above and, if appropriate, applies operator $i$ to
\verb|phiOld| storing the result in \verb|phiNew|. In some cases it just advances, through the WFT,
state \verb|phiOld| into \verb|phiNew|.
Let's look at the algorithm in detail.%'
@d evolvePrivate
@{
size_t evolve(
		size_t i,
		VectorWithOffsetType& phiNew,
		VectorWithOffsetType& phiOld,
		RealType Eg,
		size_t direction,
		const BlockType& block,
		size_t loopNumber,
		size_t lastI)
{
	@<checkIfWeAreInTheRightLoop@>
	@<checkIfAddedBlockIsSizeOne@>
	@<checkStage@>
	@<checkOperator@>
	@<computeAtimesPsi@>
	return 1;
}
@}

If we have not yet reached the finite loop that the user specified as a starting loop,
or if we are in the infinite algorithm phase, then we do nothing:
@d checkIfWeAreInTheRightLoop
@{
if (tstStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return 0;
@}
Currently this class can only deal with a DMRG algorithm that uses single site blocks for growth and 
shrinkage:
@d checkIfAddedBlockIsSizeOne
@{
if (block.size()!=1) throw
std::runtime_error("DynamicTargetting::evolve(...):"
	" blocks of size != 1 are unsupported (sorry)\n");
size_t site = block[0];
@}

If the stage is disabled and this is not the site on which the user specified, through \verb|tstStruct_.sites|,
to apply the operator, then do nothing:
@d checkStage
@{
if (site != tstStruct_.sites[i] && stage_[i]==DISABLED) return 0;
@}

If we are on the site specified by the user to apply the operator,
and we were disabled before, change the stage to \verb|OPERATOR|.
Otherwise, do not apply the operator, just advance in space one site using the WFT.
We also check the order in which sites were specified by the user against the order in which sites
are presented to us by the DMRG sweeping. This is explained below under function \verb|checkOrder|
@d checkOperator
@{
if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
else stage_[i]=CONVERGING;
if (stage_[i] == OPERATOR) checkOrder(i);
@}

We now print some progress.
Up to now, we simply set the stage but we are ready to apply the operator to the state \verb|phiOld|.
We delegate that to function \verb|computePhi| that will be explain below.
@d computeAtimesPsi
@{
std::ostringstream msg;
msg<<"Evolving, stage="<<getStage(i)<<" site="<<site<<" loopNumber="<<loopNumber;
msg<<" Eg="<<Eg;
progress_.printline(msg,std::cout);
				
// phi = A|psi>
computePhi(i,phiNew,phiOld,direction);
@}

Let us look at \verb|computephi|, the next private function.
@d computephi
@{
void computePhi(size_t i,VectorWithOffsetType& phiNew,
	VectorWithOffsetType& phiOld,size_t systemOrEnviron)
{
	if (stage_[i]==OPERATOR) {
		@<computePhiOperator@>
	} else if (stage_[i]== CONVERGING) {
		@<computePhiAdvance@>
	} else {
		throw std::runtime_error("It's 5 am, do you know what line "
			" your code is exec-ing?\n");
	}
}
@}
If we're in the stage of applying operator $i$, then we %'
call \verb|applyLocal| (see function \verb|operator()| in file \verb|ApplyLocalOperator.h|)
to apply this operator to state \verb|phiOld| and store the result in \verb|phiNew|.
@d computePhiOperator
@{
std::ostringstream msg;
msg<<"I'm applying a local operator now";
progress_.printline(msg,std::cout);
FermionSign fs(basisS_,tstStruct_.electrons);
applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],fs,systemOrEnviron);
RealType norma = norm(phiNew);
if (norma==0) throw std::runtime_error("Norm of phi is zero\n");
//std::cerr<<"Norm of phi="<<norma<<" when i="<<i<<"\n";
@}
Else we need to advance in space with the WFT. In principle, to do this we just call
function \verb|setInitialVector| in file \verb|WaveFunctionTransformation.h| as you can see below.
There is, however, a slight complication, in that the \verb|WaveFunctionTransformation| class expects
to know which sectors in the resulting vector (\verb|phiNew|) will turn out to be non-zero.
So, we need to either guess which sectors will be non-zero by calling \verb| guessPhiSectors| as described below,
or just populate all sectors with and then ``collapse'' the non-zero sectors for efficiency.
@d computePhiAdvance
@{
std::ostringstream msg;
msg<<"I'm calling the WFT now";
progress_.printline(msg,std::cout);

if (tstStruct_.aOperators.size()==1) guessPhiSectors(phiNew,i,systemOrEnviron);
else phiNew.populateSectors(basisSE_);

// OK, now that we got the partition number right, let's wft:
waveFunctionTransformation_.setInitialVector(phiNew,targetVectors_[0],
	basisS_,basisE_,basisSE_); // generalize for su(2)
phiNew.collapseSectors();
@}

The \verb|cocoon| function measures the density of all time vectors \emph{in situ}.
This is done only for debugging purposes, and uses the function \verb|test|.
@d cocoon
@{
void cocoon(ComplexType& val,size_t direction,const BlockType& block) const
{
	size_t site = block[0];
	std::cerr<<"-------------&*&*&* Cocoon output starts\n";
	test(psi_,psi_,direction,"<PSI|A|PSI>",site);
	std::cerr<<"OMEGA "<<currentOmega_<<" "<<imag(val)<<" "<<real(val)<<" "<<site<<"\n";
	for (size_t j=0;j<targetVectors_.size();j++) {
		std::string s = "<P"+utils::ttos(j)+"|A|P"+utils::ttos(j)+">";
		test(targetVectors_[j],targetVectors_[0],direction,s,site);
	}
	std::cerr<<"-------------&*&*&* Cocoon output ends\n";
}
@}

If we see $site[i]$ then we need to make sure we've seen all sites $site[j]$ for $j\le i$.%'
In other words, the order in which the user specifies the affected sites for the application of operators
needs to be the same as the order in which the DMRG sweeping process encounters those sites. Else we throw.
@d checkOrderOfSites
@{
void checkOrder(size_t i) const
{
	if (i==0) return;
	for (size_t j=0;j<i;j++) {
		if (stage_[j] == DISABLED) {
			std::string s ="TST:: Seeing dynamic site "+utils::ttos(tstStruct_.sites[i]);
			s =s + " before having seen";
			s = s + " site "+utils::ttos(j);
			s = s +". Please order your dynamic sites in order of appearance.\n";
			throw std::runtime_error(s);
		}
	}
}
@}

The little function below returns true if the stages of all the operators to be applied 
(or of all the sites on which those operators are to be applied) is equal to $x$.
Else it returns false.
Valid stages were noted before (cross reference here FIXME).
@d allStages
@{
bool allStages(size_t x) const
{
	for (size_t i=0;i<stage_.size();i++)
		if (stage_[i]!=x) return false;
	return true;
}
@}

The function below returns true if no stage is $x$, else false.
@d noStageIs
@{
bool noStageIs(size_t x) const
{
	for (size_t i=0;i<stage_.size();i++)
		if (stage_[i]==x) return false;
	return true;
}
@}

This function returns a string (human-readable) representation of the stage given by $i$.
@d getStage
@{
std::string getStage(size_t i) const
{
	switch (stage_[i]) {
	case DISABLED:
		return "Disabled";
		break;
	case OPERATOR:
		return "Applying operator for the first time";
		break;
	case CONVERGING:
		return "Converging DDMRG";
		break;
	}
	return "undefined";
}
@}

The below function computes steps 3 and 4 of
the algorithm described in page 3 of reference~\cite{re:jeckelmann02}.
The incoming arguments are the ground state energy \verb=Eg=,
the vector \verb|phi| or $|\phi\rangle$ which is $|\phi\rangle\equiv A|\psi_{gs}\rangle$, and
the direction of growth specified in \verb|systemOrEnviron|.
Note that \verb|phi| will be stored in \verb|targetVectors_[0]|,
$|Y_A\rangle$ in \verb|targetVectors_[1]|,
$|X_A\rangle$ in \verb|targetVectors_[2]|.
@d calcDynVectors
@{
ComplexType calcDynVectors(
		RealType Eg,
		const VectorWithOffsetType& phi,
		size_t systemOrEnviron)
{
	RealType retIm = minimizeFunctional(targetVectors_[1],Eg,phi,systemOrEnviron);
	obtainXA(targetVectors_[2],targetVectors_[1],Eg);
	RealType retRe = -real(targetVectors_[2]*phi)/M_PI; // Eq.~(12a)
	targetVectors_[0] = phi;
	areAllTargetsSensible();
	return ComplexType(retRe,retIm);
}
@}

After \verb=calcDynVectors= is called, \verb=psiMin= contains $|Y_A\rangle$,
		as explained in Eq.~(15).
From Eq.~(16), \verb=iaw= contains $I_A(\omega)$.

Below we minimize Eq.~(14) of reference~\cite{re:jeckelmann02}, and obtain $\psi_{min}$ which is stored in psiMin.
@d minimizeFunctional
@{
RealType minimizeFunctional(
		VectorWithOffsetType& psiMin,
		RealType Eg,
		const VectorWithOffsetType&phi,
		size_t systemOrEnviron)
{
	VectorWithOffsetType phiCopy = phi;
	psiMin = phi;
	RealType ret = 0;
	for (size_t i=0;i<phiCopy.sectors();i++) {
		VectorType sv;
		size_t ii = phiCopy.sector(i);
		psiMin.extract(sv,ii);
		if (sv.size()==0) throw std::runtime_error("Non-zero sector is zero!\n");
		ret += minimizeFunctional(sv,Eg,phi,ii);
		psiMin.setDataInSector(sv,ii);
	}
	return ret;
}
@}

The function computes the minimum of the $W$ functional and returns the complex number $Im[G(\omega+i\eta)]$.
Note that the return values use (16).
@d minimizeFunctional2
@{
RealType minimizeFunctional(VectorType& sv,RealType Eg,const VectorWithOffsetType& phi,size_t ind)
{
	size_t p = basisSE_.findPartitionNumber(phi.offset(ind));
	typename ModelType::ModelHelperType modelHelper(p,basisSE_,basisS_,basisE_,model_.orbitals());
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	LanczosMatrixType h(&model_,&modelHelper);
	typedef DynamicFunctional<RealType,LanczosMatrixType,VectorType> DynamicFunctionalType;
	VectorType aVector;
	phi.extract(aVector,ind);
	DynamicFunctionalType wFunctional(h,aVector,currentOmega_,Eg,tstStruct_.eta);
	size_t maxIter = 1000;

	PsimagLite::Minimizer<RealType,DynamicFunctionalType> min(wFunctional,maxIter);
	std::vector<RealType> svReal(2*sv.size());
	//wFunctional.packComplexToReal(svReal,sv);
	for (size_t i=0;i<svReal.size();i++) svReal[i]=drand48();
	RealType norma = std::norm(svReal);
	for (size_t i=0;i<svReal.size();i++) svReal[i]/=norma;

	int iter = -1;
	RealType delta = 1e-3;
	RealType tolerance = 1e-3;
	size_t counter = 0;
	while (iter<0 && counter<100) {
		iter = min.simplex(svReal,delta,tolerance);
		if (iter>=0) {
			std::cerr<<"delta="<<delta<<" tolerance="<<tolerance<<"\n";
		}
		delta /= 2;
		tolerance *= 1.2;
		counter++;
	}
	if (iter<0) {
		std::cerr<<"delta="<<delta<<" tol="<<tolerance<<"\n";
		throw std::runtime_error
			("DynTargetting::minimizeFunctional(...):No minimum found\n");
	}
	wFunctional.packRealToComplex(sv,svReal);
	return  -wFunctional(svReal)/(M_PI*tstStruct_.eta); // Eq.~(16)
}
@}

The function below implements Eq.~(11) of reference \cite{re:jeckelmann02}.

@d obtainXA
@{
void obtainXA(
		VectorWithOffsetType& xa,
		const VectorWithOffsetType& ya,
		RealType Eg)
{
	xa = ya;
	for (size_t i=0;i<ya.sectors();i++) {
		size_t ii = ya.sector(i);
		obtainXA(xa,Eg,ya,ii);
	}
}
@}

@d obtainXA2
@{
void obtainXA(VectorWithOffsetType& xa,RealType Eg,const VectorWithOffsetType& ya,size_t i)
{
	size_t p = basisSE_.findPartitionNumber(ya.offset(i));
	typename ModelType::ModelHelperType modelHelper(p,basisSE_,basisS_,basisE_,model_.orbitals());
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	LanczosMatrixType h(&model_,&modelHelper);
	VectorType yaThisSector;
	ya.extract(yaThisSector,i);
	VectorType sv(yaThisSector.size(),0.0);
	h.matrixVectorProduct(sv,yaThisSector); // sv = H * yaThisSector
	RealType factor =  (Eg+currentOmega_);
	sv -= (yaThisSector * factor);
	sv *= (1/tstStruct_.eta);
	xa.setDataInSector(sv,i);
}
@}

void areAllTargetsSensible

As explained above (cross reference here), we need to know before applying an operator
were the non-zero sectors are going to be. The operator (think $c^\dagger$)
does not necessarily have the symmetry of the Hamiltonian, so non-zero sectors of the original
vector are not---in general---going to coincide with the non-zero sectors of the result vector, neither
will the empty sectors be the same.
Note that using simpy
\begin{verbatim}:
size_t partition = targetVectors_[0].findPartition(basisSE_);
\end{verbatim}
doesn't work, since \verb|targetVectors_[0]| is stale at this point%'
This function should not be called when more than one operator will be applied.
@d guessPhiSectors
@{
void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
{
	FermionSign fs(basisS_,tstStruct_.electrons);
	if (allStages(CONVERGING)) {
		VectorWithOffsetType tmpVector = psi_;
		for (size_t j=0;j<tstStruct_.aOperators.size();j++) {
			applyOpLocal_(phi,tmpVector,tstStruct_.aOperators[j],fs,
					systemOrEnviron);
			tmpVector = phi;
		}
		return;
	}
	applyOpLocal_(phi,psi_,tstStruct_.aOperators[i],fs,
			systemOrEnviron);
}
@}

The function below makes all target vectors empty:
@d zeroOutVectors
@{
void zeroOutVectors()
{
	for (size_t i=0;i<targetVectors_.size();i++)
		targetVectors_[i].resize(basisSE_.size());
}
@}
The function below prints all target vectors to disk, using the \verb|TimeSerializer| class.
@d printVectors
@{
void printVectors(const std::vector<size_t>& block)
{
	if (block.size()!=1) throw std::runtime_error(
			"DynamicTargetting only supports blocks of size 1\n");

	TimeSerializerType ts(currentOmega_,block[0],targetVectors_);
	ts.save(io_);
}
@}

Print header to disk to index the dyn vectors. This indexing wil lbe used at postprocessing.
@d printHeader
@{
void printHeader()
{
	io_.print(tstStruct_);
	std::string label = "omega";
	std::string s = "Omega=" + utils::ttos(currentOmega_);
	io_.printline(s);
	label = "weights";
	io_.printVector(weight_,label);
	s = "GsWeight="+utils::ttos(gsWeight_);
	io_.printline(s);
}
@}

The \verb|test| function below performs a measurement \emph{in situ}.
This is mainly for testing purposes, since measurements are better done, post-processing.
@d test
@{
void test(
		const VectorWithOffsetType& src1,
		const VectorWithOffsetType& src2,
		size_t systemOrEnviron,
		const std::string& label,
		size_t site) const
{
	VectorWithOffsetType dest;
	OperatorType A = tstStruct_.aOperators[0];
	CrsMatrix<ComplexType> tmpC(model_.getOperator("c",0,0));
	/*CrsMatrix<ComplexType> tmpCt;
				transposeConjugate(tmpCt,tmpC);
				multiply(A.data,tmpCt,tmpC);*/
	A.fermionSign = 1;
	A.data.makeDiagonal(tmpC.rank(),1.0);
	FermionSign fs(basisS_,tstStruct_.electrons);
	applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

	ComplexType sum = 0;
	for (size_t ii=0;ii<dest.sectors();ii++) {
		size_t i = dest.sector(ii);
		size_t offset1 = dest.offset(i);
		for (size_t jj=0;jj<src2.sectors();jj++) {
			size_t j = src2.sector(jj);
			size_t offset2 = src2.offset(j);
			if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
			for (size_t k=0;k<dest.effectiveSize(i);k++)
				sum+= dest[k+offset1] * conj(src2[k+offset2]);
		}
	}
	std::cerr<<site<<" "<<sum<<" "<<" "<<currentOmega_;
	std::cerr<<" "<<label<<std::norm(src1)<<" "<<std::norm(src2)<<" "<<std::norm(dest)<<"\n";
}
@}

The function below is just for checking:
@d areAllTargetsSensible
@{
void areAllTargetsSensible() const
{
	for (size_t i=0;i<targetVectors_.size();i++)
		isThisTargetSensible(i);
}
@}

@d isThisTargetSensible
@{
void isThisTargetSensible(size_t i) const
{
	RealType norma = std::norm(targetVectors_[i]);
	if (norma<1e-6) throw std::runtime_error("Norma is zero\n");
}
@}
\bibliographystyle{plain}
\bibliography{thesis}

\end{document}
