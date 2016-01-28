#ifndef DMRG_FINITELOOP_H
#define DMRG_FINITELOOP_H

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
\begin{table}
\begin{tabular}{ll}\toprule
Bit & Description\\
0       & save or don't save state for the observe code\\
1       & compute the g.s. or fast WFT it\\
2       & compute the g.s. or slowly WFT it\\
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

There is some checking done to the finite loops input, see PTEXREF{139},
but you might find that it's not comprehensive.
*/
struct FiniteLoop {
	int stepLength; // how much to go right (+) or left (-)
	unsigned int keptStates; // kept states
	int saveOption; // to save or not to save

	FiniteLoop() : stepLength(0),keptStates(0),saveOption(0)
	{}

	FiniteLoop(int sl,unsigned int ks,int so)
	    : stepLength(sl),keptStates(ks),saveOption(so)
	{}

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & stepLength;
		ar & keptStates;
		ar & saveOption;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "FiniteLoop";
		const char* start = (const char *)&stepLength;
		const char* end = (const char*)&keptStates;
		SizeType total = mres.memResolv(&stepLength,end-start,str + " stepLength");

		start = end;
		end = (const char*)&saveOption;
		total += mres.memResolv(&keptStates,end-start,str + " keptStates");
		total += mres.memResolv(&saveOption,sizeof(*this)-total,str + " saveOption");

		return total;
	}
};

} // namespace Dmrg

#endif // DMRG_FINITELOOP_H

