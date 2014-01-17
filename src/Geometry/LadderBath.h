/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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

/** \ingroup PsimagLite */
/*@{*/

/*! \file LadderBath.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDER_BATH_H
#define LADDER_BATH_H

#include "Ladder.h"
#include "String.h"

namespace PsimagLite {

template<typename InputType>
class LadderBath : public GeometryBase<InputType> {

	typedef std::pair<int,int> PairType;
	typedef Ladder<InputType> LadderType;

public:

	enum {DIRECTION_X=LadderType::DIRECTION_X,
		  DIRECTION_Y=LadderType::DIRECTION_Y,
		  DIRECTION_BATH};

	LadderBath(SizeType linSize,InputType& io)
	    : linSize_(linSize),ladder_(0)
	{
		io.readline(bathSitesPerSite_,"BathSitesPerSite=");
		if (bathSitesPerSite_ < 0)
			throw RuntimeError("BathSitesPerSite<0 is an error\n");

		clusterSize_ = linSize_/(1+bathSitesPerSite_);

		ladder_ = new LadderType(clusterSize_,io);
	}

	~LadderBath()
	{
		if (ladder_) delete ladder_;
	}

	virtual SizeType dirs() const { return 3; }

	virtual SizeType length(SizeType i) const
	{
		return this->unimplemented("length");
	}

	virtual SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		return this->unimplemented("translate");
	}

	SizeType getVectorSize(SizeType dirId) const
	{
		if (dirId==DIRECTION_BATH) return bathSitesPerSite_*clusterSize_;
		return ladder_->getVectorSize(dirId);
	}

	bool connected(SizeType i1,SizeType i2) const
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
			return (SizeType(c2)==i1) ? true : false;
		}
		// bath - cluster
		return (SizeType(c1)==i2) ? true : false;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		int c1 = getClusterSite(i1).first;
		int c2 = getClusterSite(i2).first;
		// two possibilities
		// 1. both in the cluster
		if (c1<0 && c2<0) return calcDirInCluster(i1,i2);
		// cluster - bath or bath cluster:
		return DIRECTION_BATH;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		int c = getClusterSite(i).first;
		if (c>=0) return false; // no bath site is ever fringe
		return fringeInCluster(i,smax,emin);
	}

	// assumes i1 and i2 are connected
	SizeType handle(SizeType i1,SizeType i2) const
	{
		PairType c1 = getClusterSite(i1);
		PairType c2 = getClusterSite(i2);
		// two possibilities
		// 1. both in the cluster
		if (c1.first<0 && c2.first<0) return handleInCluster(i1,i2);
		// cluster - bath or bath cluster
		PairType x =  (c1.first<0) ? c2 : c1;
		if (x.first<0 || x.second<0) throw RuntimeError("Internal error in handle\n");
		SizeType firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
		x.first -= firstClusterSite;

		return x.first*bathSitesPerSite_+x.second;
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		throw RuntimeError("Umhph, ouch, ayyayyayy, what?\n");
	}

	String label() const
	{
		return "ladderbath";
	}

	SizeType maxConnections() const
	{
		return clusterSize_+1;
	}

	SizeType findReflection(SizeType site) const
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

private:

	// if i is in the cluster return the pair (-1,-1)
	// else return the corresponding cluster site c and the number
	// of this bath site as a pair (c,b)
	PairType getClusterSite(SizeType i) const
	{
		SizeType firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
		SizeType lastP1ClusterSite = firstClusterSite + clusterSize_;
		if (i>=firstClusterSite && i<lastP1ClusterSite) return PairType(-1,-1);

		SizeType middle = linSize_/2;
		SizeType cs = clusterSize_/2;
		// now i is in the bath:
		if (i<middle) { // i is in the system
			return PairType(i%cs + firstClusterSite,i/cs);
		}
		// is in the bath and in the environ:
		SizeType iprime = i - lastP1ClusterSite;
		SizeType offset = lastP1ClusterSite - cs;
		return PairType(iprime%cs + offset,iprime/cs);
	}

	// assumes i1 and i2 are in the cluster
	// if connected return true, else false
	bool connectedInCluster(SizeType i1,SizeType i2) const
	{
		ladderize(i1,i2);
		return ladder_->connected(i1,i2);
	}

	// assumes i1 and i2 are connected and in the cluster
	SizeType calcDirInCluster(SizeType i1,SizeType i2) const
	{
		ladderize(i1,i2);
		return ladder_->calcDir(i1,i2);
	}

	// assumes i1 and i2 are in the cluster
	bool fringeInCluster(SizeType i,SizeType smax,SizeType emin) const
	{
		SizeType firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
		i -= firstClusterSite;
		smax -= firstClusterSite;
		emin -= firstClusterSite;
		return ladder_->fringe(i,smax,emin);
	}

	// assumes i1 and i2 are connected and in the cluster
	SizeType handleInCluster(SizeType i1,SizeType i2) const
	{
		ladderize(i1,i2);
		return ladder_->handle(i1,i2);
	}

	void ladderize(SizeType& i1,SizeType& i2) const
	{
		SizeType firstClusterSite = (clusterSize_/2)*bathSitesPerSite_;
		i1 -= firstClusterSite;
		i2 -= firstClusterSite;
	}

	SizeType linSize_;
	SizeType bathSitesPerSite_;
	SizeType clusterSize_;
	LadderType* ladder_;
}; // class LadderBath
} // namespace PsimagLite 

/*@}*/
#endif // GEOMETRY_H

