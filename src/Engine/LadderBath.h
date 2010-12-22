// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file LadderBath.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDER_BATH_H
#define LADDER_BATH_H

#include "Utils.h"

namespace Dmrg {
	
	class LadderBath {
		public:
			enum {DIRECTION_X,DIRECTION_Y,DIRECTION_BATH};

			LadderBath(size_t linSize) : linSize_(linSize)
			{
			}

			template<typename IoInputter>
			void init(IoInputter& io)
			{
				int x = 0;
				io.readline(x,"BathSitesPerSite=");
				if (x<0) throw std::runtime_error("BathSitesPerSite<0 is an error\n");
				bathSitesPerSite_ = x;
				clusterSize_ = linSize_/(1+bathSitesPerSite_);
			}

			size_t getVectorSize(size_t leg,size_t dirId) const
			{
				if (dirId==DIRECTION_X) return clusterSize_-leg;
				else if (dirId==DIRECTION_Y) return clusterSize_ - clusterSize_/leg;
				else return bathSitesPerSite_*clusterSize_;

			}

			bool connected(bool& bypass,size_t& i1,size_t& i2) const
			{
				size_t j1 = i1;
				size_t j2 = i2;
				bypass = false;
				if (connectedBathLadder(i1,i2)) return true;
				// only chance of connection now is that both are cluster sites:
				if (getClusterSite(j1)>=0 || getClusterSite(j2)>=0) return false;
				bypass = true;
				return true; // bogus
			}

			size_t calcDir(bool& bypass,size_t& i1,size_t& i2) const
			{
				bypass = false;
				if (connectedBathLadder(i1,i2)) return DIRECTION_BATH;
				bypass = true;
				return 0; // bogus
			}

			bool fringe(bool& bypass,size_t& i,size_t& smax,size_t& emin) const
			{
				bypass = false;
				if (getClusterSite(i)>=0) return false; // no bath site is ever fringe
				bypass = true;
				fringeBathLadder(i,smax,emin);
				return true; //bogus
			}

			// are sites i1 and site i2 connected thru the bath?
			bool connectedBathLadder(size_t& i1,size_t& i2) const
			{
				size_t middle = linSize_/2;
				size_t clustersize = linSize_/(1+bathSitesPerSite_);
				size_t totalBathSites = clustersize*bathSitesPerSite_;
				size_t ismall = (i1<i2) ? i1 : i2;
				size_t ibig = (i1<i2) ? i2 : i1;
				size_t ii1 = (i1<middle) ? i1 : i1-totalBathSites;
				size_t ii2 = (i2<middle) ? i2 : i2-totalBathSites;
				i1 = ii1; i2=ii2;
				if (ismall<middle && ibig>=middle) return false;
				if (ismall<middle) {
					// only chance is that ismall is in the cluster and ibig in in the bath:
					int x = getClusterSite(ibig);
					if (x<0) return false; // both in the cluster
					if (size_t(x)==ismall) return true;
					return false;
				}
				// ismall>=middle
				// only chance is that ismall is in the bath and ibig in in the cluster:
				int x = getClusterSite(ismall);
				if (x<0) return false;
				if (size_t(x)==ibig) return true;
				return false;
			}

			void fringeBathLadder(size_t& i,size_t& smax,size_t& emin) const
			{
				size_t middle = linSize_/2;
				size_t clustersize = linSize_/(1+bathSitesPerSite_);
				size_t totalBathSites = clustersize*bathSitesPerSite_;
				size_t ii = (i<middle) ? i : i-totalBathSites;
				i = ii;
				// find the maximum *cluster* site that is less or equal smax:
				size_t newsmax = 0;
				for (size_t i=0;i<=smax;i++) {
					int x = getClusterSite(i);
					if (x<0) { // i is in the cluster
						if (i>newsmax) newsmax=i;
					} else { // i is in the bath with cluster site x:
						if (size_t(x)>newsmax) newsmax=x;
					}
				}
				smax = (newsmax<middle) ? newsmax : newsmax-totalBathSites;

				// find the minimum *cluster* site that is greater or equal emin:
				size_t newemin = linSize_;
				for (size_t i=0;i<=emin;i++) {
					int x = getClusterSite(i);
					if (x<0) { // i is in the cluster
						if (i<newemin) newemin=i;
					} else { // i is in the bath with cluster site x:
						if (size_t(x)<newemin) newemin=x;
					}
				}
				emin = (newemin<middle) ? newemin : newemin-totalBathSites;

			}

			// if i is in the cluster return -1
			// else return the corresponding cluster site
			int getClusterSite(size_t i) const
			{
				size_t middle = linSize_/2;
				size_t c = linSize_/(1+bathSitesPerSite_); // clustersize
				if (i<c/2 || i>=linSize_-c/2) return -1;
				// ok, i is in the bath:
				int s = -1;
				if (i<middle) {
					//find s such that i = s + c/2*x for some x=1,..bathSitesPerSite
					for (size_t x=1;x<=bathSitesPerSite_;x++) {
						if (i<c/2*x) continue;
						s = i-c/2*x;
						if (size_t(s)>=c/2) continue;
						else break;
					}
					return s;
				}
				// i>=middle
				//find s such that s = i + c/2*x for some x=1,..bathSitesPerSite
				for (size_t x=1;x<=bathSitesPerSite_;x++) {
					s = i + c/2*x;
					if (size_t(s)<linSize_-c/2) continue;
					else break;
				}
				return s;
			}

		private:
			size_t linSize_;
			size_t bathSitesPerSite_;
			size_t clusterSize_;
	}; // class LadderBath
} // namespace Dmrg 

/*@}*/
#endif // GEOMETRY_H
