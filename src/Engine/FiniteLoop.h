#ifndef DMRG_FINITELOOP_H
#define DMRG_FINITELOOP_H

#include "Vector.h"
#include <iostream>
#include "Io/IoSerializerStub.h"
#include "PsimagLite.h"

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
class FiniteLoop {

public:

	FiniteLoop(int sl, SizeType ks, SizeType so)
	    : stepLength_(sl), keptStates_(ks), bitField_(so)
	{
		checkBitField();
	}

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer,
	           PsimagLite::IoNgSerializer::WriteMode mode =
	        PsimagLite::IoNgSerializer::NO_OVERWRITE) const
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

	bool wantsSave() const { return (bitField_ & 1); }

	bool wantsOnlyFastWft() const { return (bitField_ & 2); }

	bool wantsOnlySlowWft() const { return (bitField_ & 4); }

	bool wantsRandomGuess() const { return (bitField_ & 8); }

	bool wantsMultiSitePush() const { return (bitField_ & 16); }

	bool wantsOneSiteTruncation() const { return (bitField_ & 32); }

private:

	void checkBitField() const
	{
		// Only 1 bit of bit 1, bit 2, and bit 3 may be set
		SizeType saveOption = bitField_;
		bool flag = false;
		for (SizeType i = 0; i < 3; ++i) {
			saveOption >>= 1; // ATTENTION: Discards bit 0 when i === 0
			if (!(saveOption & 1)) continue;
			if (flag) {
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"Triplet 3rd: Only one bit of bit 1, bit 2, and bit 3 may be set";
				msg<<" VALUE= "<<bitField_<<"\n";
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

	int stepLength_; // how much to go right (+) or left (-)
	SizeType keptStates_; // kept states
	SizeType bitField_;
};

} // namespace Dmrg

#endif // DMRG_FINITELOOP_H
