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

#include "Ladder.h"

namespace Dmrg {
	
	class LadderBath  {
			typedef std::pair<int,int> PairType;
			typedef Ladder LadderType;
		public:
			enum {DIRECTION_X=LadderType::DIRECTION_X,DIRECTION_Y=LadderType::DIRECTION_Y,DIRECTION_BATH};

			LadderBath(size_t linSize,size_t leg,size_t bathSitesPerSite) :
				linSize_(linSize),bathSitesPerSite_(bathSitesPerSite),
				clusterSize_(linSize_/(1+bathSitesPerSite_)),ladder_(clusterSize_,leg)
			{
			}

			size_t getVectorSize(size_t dirId) const
			{
				try {
					ladder_.getVectorSize(dirId);
				} catch (std::exception& e) {
				}
				return bathSitesPerSite_*clusterSize_;
			}

			bool connected(size_t i1,size_t i2) const
			{
				if (i1==i2) return false;
				int c1 = getClusterSite(i1).first;
				int c2 = getClusterSite(i2).first;
				//4 possibilites
				// 1. both in the cluster
				if (c1<0 && c2<0) return connectedInCluster(i1,i2);
				// both in the bath:
				if (c1>=0 && c2>=0) return false;
				// cluster - bath
				if (c1<0) {
					return (size_t(c2)==i1) ? true : false;
				}
				// bath - cluster
				return (size_t(c1)==i2) ? true : false;
			}

			// assumes i1 and i2 are connected
			size_t calcDir(size_t i1,size_t i2) const
			{
				int c1 = getClusterSite(i1).first;
				int c2 = getClusterSite(i2).first;
				// two possibilities
				// 1. both in the cluster
				if (c1<0 && c2<0) return calcDirInCluster(i1,i2);
				// cluster - bath or bath cluster:
				return DIRECTION_BATH;
			}

			bool fringe(size_t i,size_t smax,size_t emin) const
			{
				int c = getClusterSite(i).first;
				if (c>=0) return false; // no bath site is ever fringe
				return fringeInCluster(i,smax,emin);
			}

			// assumes i1 and i2 are connected
			size_t handle(size_t i1,size_t i2) const
			{
				PairType c1 = getClusterSite(i1);
				PairType c2 = getClusterSite(i2);
				// two possibilities
				// 1. both in the cluster
				if (c1.first<0 && c2.first<0) return handleInCluster(i1,i2);
				// cluster - bath or bath cluster
				PairType x =  (c1.first<0) ? c2 : c1;
				if (x.first<0 || x.second<0) throw std::runtime_error("Internal error in handle\n");
				size_t firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
				x.first -= firstClusterSite;

				return x.first*bathSitesPerSite_+x.second;
			}

			// siteNew2 is fringe in the environment
			size_t getSubstituteSite(size_t smax,size_t emin,size_t siteNew2) const
			{
				throw std::runtime_error("Umhph, ouch, ayyayyayy, what?\n");
			}

			std::string label() const
			{
				return "ladderbath";
			}

		private:

			// if i is in the cluster return the pair (-1,-1)
			// else return the corresponding cluster site c and the number
			// of this bath site as a pair (c,b)
			PairType getClusterSite(size_t i) const
			{
				size_t firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
				size_t lastP1ClusterSite = firstClusterSite + clusterSize_;
				if (i>=firstClusterSite && i<lastP1ClusterSite) return PairType(-1,-1);

				size_t middle = linSize_/2;
				size_t cs = clusterSize_/2;
				// now i is in the bath:
				if (i<middle) { // i is in the system
					return PairType(i%cs + firstClusterSite,i/cs);
				}
				// is in the bath and in the environ:
				size_t iprime = i - lastP1ClusterSite;
				size_t offset = lastP1ClusterSite - cs;
				return PairType(iprime%cs + offset,iprime/cs);
			}

			// assumes i1 and i2 are in the cluster
			// if connected return true, else false
			bool connectedInCluster(size_t i1,size_t i2) const
			{
				ladderize(i1,i2);
				return ladder_.connected(i1,i2);
			}

			// assumes i1 and i2 are connected and in the cluster
			size_t calcDirInCluster(size_t i1,size_t i2) const
			{
				ladderize(i1,i2);
				return ladder_.calcDir(i1,i2);
			}

			// assumes i1 and i2 are in the cluster
			bool fringeInCluster(size_t i,size_t smax,size_t emin) const
			{
				size_t firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
				i -= firstClusterSite;
				smax -= firstClusterSite;
				emin -= firstClusterSite;
				return ladder_.fringe(i,smax,emin);
			}

			// assumes i1 and i2 are connected and in the cluster
			size_t handleInCluster(size_t i1,size_t i2) const
			{
				ladderize(i1,i2);
				return ladder_.handle(i1,i2);
			}

			void ladderize(size_t& i1,size_t& i2) const
			{
				size_t firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
				i1 -= firstClusterSite;
				i2 -= firstClusterSite;
			}

			size_t linSize_;
			size_t bathSitesPerSite_;
			size_t clusterSize_;
			LadderType ladder_;
	}; // class LadderBath
} // namespace Dmrg 

/*@}*/
#endif // GEOMETRY_H
