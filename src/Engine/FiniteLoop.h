#ifndef DMRG_FINITELOOP_H
#define DMRG_FINITELOOP_H

#include "Io/IoSerializerStub.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "TruncationControl.h"
#include "Vector.h"
#include <iostream>
#include <map>

namespace Dmrg {
/* PSIDOC FiniteLoop
\subsection{Enabling finite loops}

To enable finite loops make sure that the option `nofiniteloops` is \emph{not}
present under `SolverOptions=`. Remember that the entry \verb!FiniteLoops!
in the input file is a series of space-separated numbers. More than one space is allowed.
The first
number is the number of finite algorithm ``movements,'' followed by series of three
numbers for
each movement. Of the three numbers, the first is the number of sites to go forward
if positive
or backward if negative. The second number is the $m$ for this movement and the
last number
is a bitwise option described in \S\ref{sec:thirdnumber}.
The first movement starts from where the infinite loop left off, at the
 middle of the
lattice.

\subsection{Example of a Finite loops line in the input file}

\begin{verbatim}
        FiniteLoops 4 7 200 0 -7 200 0 7 200 1 7 200 1
\end{verbatim}

The number 4 implies 4 finite loops. The first fine loop is \verb!7 200 0!, meaning
go forward 7 steps, use \emph{m=200} for this finite sweep, and 0: do not store
transformation in disk.
The next is \verb!-7 200 0!, which goes backwards 7 sites, etc.
Remember that the finite loops start at the middle of the lattice, where the
infinite loop left off.
\todo{ADD FIGURE SHOWING WHAT THIS DOES.}

\subsection{The third number in the triplet}\label{sec:thirdnumber}
The third number in the triplet is a bitwise option where the
first bit means save or don't save, the second bit
compute the g.s. or WFT it without eigenvalue or eigenvector updates,
and the third bit compute the g.s. of WFT it updating eigenvalue and eigenvector.
It is a fatal error to have both bits 1 and 2 set.
If bit 5 is set then so must bit 3.
\begin{table}
\begin{tabular}{ll}\toprule
Bit & What Happens if Set\\
0       & Saves state for the observe code\\
1       & WFTs the ground state in a fast way instead of computing it\\
2       & WFTs the ground state slowly  instead of computing it\\
3       & Forces random guess for ground state\\
4       & MultiSitePush (whatever that means)\\
5       & OneSiteTruncation\\
\end{tabular}
\caption{Meaning of each bit of the third number in the
finite loop triplet. It is a fatal error to have both bits 1 and 2 set.}
\end{table}
\subsection{Caveats and Troubleshooting}

If `nofiniteloops` is an option in the options line of the input file then
the \verb!FiniteLoops! line in the input file is ignored, and no finite loops are done.
In this case, DMRG++ stops when the infinite algorithm has finished.

Make sure the first number is the number of triplets that follow.

Make sure you don't fall off the lattice, by going forward or backwards too much.
Remember that at least one site must remain for the system part of the lattice.
So on a 16 site chain, when you start the finite loops you're at the middle, you
can go forward at most 7 sites, and backwards at most 7 sites.

There is some checking done to the finite loops input,
but you might find that it's not comprehensive.
*/
template <typename RealType>
class FiniteLoop {

public:

	using VectorStringType = PsimagLite::Vector<PsimagLite::String>::Type;
	using TruncationControlType = TruncationControl<RealType>;

	FiniteLoop(int sl,
	           SizeType ks,
	           PsimagLite::String str,
	           const TruncationControlType& truncationControl)
	    : stepLength_(sl)
	    , keptStates_(ks)
	    , bitField_(0)
	    , truncationControl_(truncationControl)
	{
		str = removeCharsAtStart(str, "\"");
		str = removeCharsAtEnd(str, "\"");

		setMap();

		VectorStringType tokens;
		PsimagLite::split(tokens, str, ":");
		bool hasBitField = false;
		bool hasAtField = false;
		for (SizeType i = 0; i < tokens.size(); ++i)
			procToken(hasBitField, hasAtField, tokens[i]);

		warnIfFiniteMlessThanMin();
		checkBitField();
	}

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer,
	           PsimagLite::IoNgSerializer::WriteMode mode = PsimagLite::IoNgSerializer::NO_OVERWRITE) const
	{
		PsimagLite::String root = label;

		if (mode != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			ioSerializer.createGroup(root);
		ioSerializer.write(root + "/stepLength", stepLength_, mode);
		ioSerializer.write(root + "/keptStates", keptStates_, mode);
		ioSerializer.write(root + "/saveOption", bitField_, mode);
	}

	int stepLength() const { return stepLength_; }

	SizeType keptStates() const { return keptStates_; }

	const TruncationControlType& truncationControl() const { return truncationControl_; }

	bool wants(PsimagLite::String what) const
	{
		PsimagLite::String whatLower = ProgramGlobals::toLower(what);
		assert(mapNameBit_.find(whatLower) != mapNameBit_.end());
		return (bitField_ & mapNameBit_.at(whatLower));
	}

	void print(std::ostream& os) const
	{
		os << stepLength_ << " " << keptStates_ << " " << bitField_;
	}

private:

	void setMap()
	{
		const VectorStringType str = { "save",
			                       "onlyfastwft",
			                       "onlyslowwft",
			                       "randomguess",
			                       "multisitepush",
			                       "onesitetruncation" };

		SizeType bit = 1;
		for (SizeType i = 0; i < str.size(); ++i) {
			mapNameBit_[str[i]] = bit;
			bit <<= 1;
		}
	}

	// either @something or
	// setting=value
	void procToken(bool& hasBitField, bool& hasAtField, PsimagLite::String str)
	{
		if (isAdigit(str)) {
			dieIfTrue(hasBitField, "Bit field already set for this finite loop\n");
			dieIfTrue(hasAtField, "Bit field already set for this finite loop\n");
			hasBitField = true;
			bitField_ = PsimagLite::atoi(str);
			return;
		}

		if (str.length() < 2)
			err("FiniteLoop with " + str + " : syntax error (1)\n");

		if (str[0] == '@') {
			dieIfTrue(hasBitField, "Bit field already set for this finite loop\n");
			procAtValue(str.substr(1, str.length() - 1));
			hasAtField = true;
			return;
		}

		VectorStringType tokens;
		PsimagLite::split(tokens, str, "=");
		if (tokens.size() != 2)
			err("FiniteLoop with " + str + " : syntax error (2)\n");

		truncationControl_.setNameValue(tokens[0], tokens[1]);
	}

	void procAtValue(PsimagLite::String what)
	{
		assert(mapNameBit_.find(what) != mapNameBit_.end());
		bitField_ |= mapNameBit_.at(what);
	}

	void warnIfFiniteMlessThanMin() const
	{
		if (keptStates_ >= truncationControl().mMin())
			return;

		std::cout << "WARNING: Triplet has m= " << keptStates_;
		std::cout << " which is less than its minimum m = " << truncationControl().mMin();
		std::cout << " as found in TruncationTolerance or set only for this loop\n";
	}

	void checkBitField() const
	{
		// Only 1 bit of bit 1, bit 2, and bit 3 may be set
		SizeType saveOption = bitField_;
		bool flag = false;
		for (SizeType i = 0; i < 3; ++i) {
			saveOption >>= 1; // ATTENTION: Discards bit 0 when i === 0
			if (!(saveOption & 1))
				continue;
			if (flag) {
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg << "Triplet 3rd: Only one bit of bit 1, bit 2, and bit 3 may be set";
				msg << " VALUE= " << bitField_ << "\n";
				err(msg.str());
			} else {
				flag = true;
			}
		}

		bool wantsOneSiteTrunc = (bitField_ & 32);
		if (wantsOneSiteTrunc && !(bitField_ & 8))
			err(PsimagLite::String("A finite loop that wants one site truncation must also want")
			    + " random guess.\nIn other words, if bit 5 is set so must bit 3.\n");
	}

	static void dieIfTrue(bool b, PsimagLite::String msg)
	{
		if (!b)
			return;
		std::cout << msg;
	}

	static bool isAdigit(PsimagLite::String str)
	{
		for (SizeType i = 0; i < str.length(); ++i)
			if (!std::isdigit(str[i]))
				return false;
		return true;
	}

	static PsimagLite::String removeCharsAtStart(PsimagLite::String str,
	                                             PsimagLite::String chars)
	{
		size_t start = str.find_first_not_of(chars);
		return (start == PsimagLite::String::npos) ? "" : str.substr(start);
	}

	static PsimagLite::String removeCharsAtEnd(PsimagLite::String str,
	                                           PsimagLite::String chars)
	{
		size_t end = str.find_last_not_of(chars);
		return (end == PsimagLite::String::npos) ? "" : str.substr(0, end + 1);
	}

	static std::map<PsimagLite::String, SizeType> mapNameBit_;
	const int stepLength_; // how much to go right (+) or left (-)
	const SizeType keptStates_; // kept states
	SizeType bitField_;
	TruncationControlType truncationControl_;
};

template <typename T>
std::map<PsimagLite::String, SizeType> FiniteLoop<T>::mapNameBit_;
} // namespace Dmrg

#endif // DMRG_FINITELOOP_H
