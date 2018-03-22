/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file RegisterSignals.h
 *
 *
 */

#ifndef DMRG_REGISTER_SIGNALS_H
#define DMRG_REGISTER_SIGNALS_H

#ifdef USE_SIGNALS
#include <signal.h>
#include "ProgressIndicator.h"
#include <iostream>
#endif

namespace Dmrg {

/* PSIDOC RegisterSignals
   PLEASE NOTE: This is an experimental (CITATION NEEDED FIXME) feature. To use it
                you must add \verb!-DUSE_SIGNALS! to
                \verb!CPPFLAGS! in the Makefile.

   \subsection{SIGUSR1}

   Rationale: When running a process in a queue batching system the standard output
   and standard error might be buffered, and, thus, might not be seen until program
   completion.
   DMRG++ allows the user to store (a fragment of) the stdout and stderr buffers into
   a temporary file to monitor program process in situations where stdout and stderr
   would not normally
   be accessible.

   Sending the signal SIGUSR1 to the DMRG++ process will result in switching the state
   of the ProgressIndicator buffer: if the state was inactive it will become active, and
   viceversa. Only when the state of the ProgressIndicator buffer is active
   does ProgressIndicator store its
   stream in memory. This stream contains the standard output and standard error printed
   by DMRG++. When the state of the ProgressIndicator is switched back from active to
   inactive,
   DMRG++ dumps the buffer into a temporary file, and closes the buffer.
   The temporary file is named bufferN.txt where N is the PID of the DMRG++ process.

   HINT: qsig might be used to send a signal if the DMRG++ process is running in
         PBS or torque.

   CAVEATS: Leaving the buffer on for long periods of time might cause high memory
   consumption. The temporary buffer file is overwritten if the buffer is used more
   than once by the same process. The temporary buffer file is not deleted at the end of
   program execution.
 */
void registerSignals()
{
#ifdef USE_SIGNALS
	int signum = SIGUSR1;
	sigset_t *maskset = new sigset_t;

	int ret = sigemptyset(maskset);
	if (ret != 0) {
		PsimagLite::String str("RegisterSignals:");
		throw PsimagLite::RuntimeError(str + " sigemptyset returned non zero\n");
	}

	struct sigaction act;
	act.sa_handler = PsimagLite::ProgressIndicator::updateBuffer;
	act.sa_flags = 0;
	act.sa_mask = *maskset;

	ret = sigaction(signum,&act,0);
	delete maskset;
	if (ret != 0) {
		PsimagLite::String msg("FATAL: sigaction failed\n");
		throw PsimagLite::RuntimeError(msg);
	}

	std::cerr<<"EXPERIMENTAL: Signal support enabled\n";

#endif
}

} // namespace Dmrg
/*@}*/
#endif

