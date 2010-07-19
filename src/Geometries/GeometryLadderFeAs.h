// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file GeometryLadderFeAs.h
 *
 *  This class implements GeometryBase for n-leg ladders with multiple bands
 */
#ifndef GEOMETRYLADDERFEAS
#define GEOMETRYLADDERFEAS

#include "Matrix.h" // in psimag
#include "GeometryBase.h"
#include "ProgressIndicator.h"

namespace Dmrg {
	//! Class to handle the connection between System and Environ
	//! in the dmrg algorithm: case of a ladder
	//! Numbering is
	//! 0---- leg -- 2*leg ----...
	//! 1 ---leg+1 --2*leg+1---...
	//! 2 ---leg+2 --2*leg+2 ---...	
	//! ...
	//! RESTRICTION: leg must be even
	template<typename Field,typename ConnectorsType>
	class GeometryLadderFeAs : public GeometryBase<Field,ConnectorsType> {
		static const size_t SystemSystem=ProgramGlobals::SYSTEM_SYSTEM;
		static const size_t SystemEnviron=ProgramGlobals::SYSTEM_ENVIRON;
		static const size_t EnvironSystem=ProgramGlobals::ENVIRON_SYSTEM;
		static const size_t EnvironEnviron=ProgramGlobals::ENVIRON_ENVIRON;
	public:
		//static const int SystemSystem=0,SystemEnviron=1,EnvironSystem=2,EnvironEnviron=3;
		typedef  typename GeometryBase<Field,ConnectorsType>::BlockType BlockType;
		
		GeometryLadderFeAs(const ConnectorsType& connectors,int sizeOfInitialBlock=1,int leg=2)
			 : connectors_(connectors), leg_(leg),
			dof_(connectors.dof()),maxSites_(connectors.linSize()),progress_("GeometryLadderFeAs",0)
		{
			split(sizeOfInitialBlock);
		}
		
		const BlockType& systemBlock() const { return systemBlock_; }
		
		//! split into S X Y and E
		void setBlocksOfSites(BlockType &S,std::vector<BlockType> &X,std::vector<BlockType> &Y,BlockType &E) const 
		{
			S=S_;
			X=X_;
			Y=Y_;
			E=E_;
		}
		
		// FIXME: spread the what!!
		Field calcConnectorValue(int type,int ind,int jnd,int smax,int emin,size_t what) const 
		{
			//! There are four cases:
			//! 1. (ind,jnd) in SUX --> use input connectors
			//! 2. (ind,jnd) in YUE --> use reflected connectors
			//! 3. ind in SUX , jnd in YUE --> delegate
			//! 4. ind in YUE, jnd in SUX --> delegate
			
			//int total = S_.size()+X_.size()+Y_.size()+E_.size(); 
			Field x=0;
			size_t indR=0,jndR=0;
			switch (type) {
				case EnvironEnviron:
					indR = findReflection(ind);
					jndR = findReflection(jnd);
					// x= connectors_(dof1+indR*dof_,dof2+jndR*dof_);
					// hack to deal with the "reflection changes x+y into x-y problem":
					x= connectors_(dof1+findXReflection(indR)*dof_,dof2+findXReflection(jndR)*dof_);
					break;
				
				case EnvironSystem:
					//x=calcHoppingAux(ind,jnd,smax,emin,0);
					x=connector(ind,dof1,jnd,dof2,smax,emin,0);
					break;
				
				case SystemEnviron:
					//x=calcHoppingAux(ind,jnd,smax,emin,1);
					x=connector(ind,dof1,jnd,dof2,smax,emin,1);
					break;
				case SystemSystem:
					x= connectors_(dof1+ind*dof_,dof2+jnd*dof_);
					break;
				
			}
			return x;
			
		}

		int calcConnectorType(int ind,int jnd) const 
		{
			return  calcConnectorType(ind,jnd,systemBlock_);
		}

		//! given i in the environment returns the site symmetric to i (in the system)
		int findReflection(int i) const
		{
			int r=(i%leg_)-int(leg_/2);
			return 2*maxSites_-i+2*r;
		}

		//!given i in the system, reflect by an horizontal axis that cuts the ladder in two
		int findXReflection(int i) const
		{
			if (i%2==0) return i+1;
			return i-1;
		}

		int calcConnectorType(int ind,int jnd,BlockType const &sBlock) const 
		{
			//! There are four cases:
			//! 1. (ind,jnd) in SUX --> use input connectors
			//! 2. (ind,jnd) in YUE --> use reflected connectors
			//! 3. ind in SUX , jnd in YUE --> delegate
			//! 4. ind in YUE, jnd in SUX --> delegate
			int i = utils::isInVector(sBlock,ind);
			int j = utils::isInVector(sBlock,jnd);
			
			if (i<0) { // i is in env
				if (j<0) { // j is in env
					//! reflect
					return EnvironEnviron;
				} else { // j is in the system
					return EnvironSystem;
				}
			} else { // i is in the system
				if (j<0) { // j is in the env
					return SystemEnviron;
				} else { // j is in the system
					return SystemSystem;
				}
			}
			
		}

		size_t connectorValues() const
		{
			return connectors_.size();
		}

	private:
		const ConnectorsType& connectors_;
		BlockType S_,E_;
		std::vector<BlockType> X_,Y_;
		BlockType systemBlock_; //,envBlock;
		int leg_;
		size_t dof_;
		size_t maxSites_;
		ProgressIndicator progress_;
		
		Field connector(int ind,size_t dof1,int jnd,size_t dof2,int smax,int emin,int type) const
		{
			if (type==0) {// ind in env, jnd in the system
				int n = smax+1;
				if (n%leg_==0) return connectorComplete(jnd,dof2,ind,dof1,smax,emin,type);
				return connectorIncomplete(jnd,dof2,ind,dof1,smax,emin,type);
			} else { // ind in system, jnd in env
				return connector(jnd,dof2,ind,dof1,smax,emin,0);
			}
			return static_cast<Field>(0);
		}

		Field connectorIncomplete(int ind,size_t dof1,int jnd,size_t dof2,int smax,int emin,int type) const
		{
			// use default connector...
			size_t direction=1; // y-direction
			Field t=connectors_.defaultValue(direction,dof1,dof2);
			// ...unless it is possible to take connectors as if current environment were spatially next to jnd
			/*if (ind<maxSites_) {
				t2=connectors_(dof1+ind*dof_,dof2+jnd*dof_);
				if (t2!=2) throw std::logic_error("connectorIncomplete direction==1\n");
			}*/
			if (ind==smax && jnd==emin) return t;
			
			if (leg_!=2) std::cerr<<"EEEEEEEEEEEEEERRRRRRRRRRRRRRROOOOOOOOOOOOOOORRRRRRRRRRRRRR\n";
			
			for (int j=0;j<leg_;j++) {
				// use default connector...
				direction=0; // x-direction
				t=connectors_.defaultValue(direction,dof1,dof2);
				// ...unless it is possible to take connectors as if current environment were spatially next to jnd
				//if (jnd+leg_<maxSites_) t=connectors_(dof2+jnd*dof_,dof1+(jnd+leg_)*dof_);
				if (ind==smax-j && jnd==emin+leg_-1-j) return t;
			}
			// smin with emin + 2 (x-y)
			direction=2;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax && jnd==emin+2) return t;

			//smin -2 with emin (x-y)
			direction=2;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax-2 && jnd==emin) return t;

			return static_cast<Field>(0);
		}

		Field connectorComplete(int ind,size_t dof1,int jnd,size_t dof2,int smax,int emin,int type) const
		{
			size_t direction=0;
			Field t=connectors_.defaultValue(direction,dof1,dof2);

			// smax with emin +1 (x)
			direction=0;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax && jnd==emin+1) return t;

			// smax -1 with emin (x)
			direction=0;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax-1 && jnd==emin) return t;

			// smax-1 with emin+1 (x-y)
			direction=2;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax-1 && jnd==emin+1) return t;

			// smax with emin (x+y)
			direction=3;
			t=connectors_.defaultValue(direction,dof1,dof2);
			if (ind==smax && jnd==emin) return t;

			return static_cast<Field>(0);
		}

		void split(int sizeOfInitialBlock) 
		{
			int i,j;
			int legOver2=int(leg_/2);
			int n=maxSites_;
			BlockType tmpBlock;
			if (sizeOfInitialBlock>n) 
				throw std::runtime_error("GeometryBase::split(): sizeofinitialblock="+
						utils::ttos(sizeOfInitialBlock)+" but n="+utils::ttos(n)+"\n");

			tmpBlock.resize(legOver2);

			// S U X = system, and Y U E is the reflexion
			for (i=0;i<sizeOfInitialBlock;i++) {
				S_.push_back(i);
				systemBlock_.push_back(i);
				std::ostringstream msg;
				msg<<"systemBlock_["<<i<<"]="<<systemBlock_[i];
				progress_.printline(msg,std::cerr);
			}
			for (i=sizeOfInitialBlock;i<n;i+=legOver2) {
				for (j=0;j<legOver2;j++) {
					tmpBlock[j]=i+j;
					systemBlock_.push_back(tmpBlock[j]);
				}
				X_.push_back(tmpBlock);
				std::ostringstream msg;
				msg<<"X added block of size "<<tmpBlock.size();
				progress_.printline(msg,std::cerr);
			}
			std::ostringstream msg0;
			msg0<<"X total size is "<<X_.size();
			progress_.printline(msg0,std::cerr);
			
			for (i=2*n-1-sizeOfInitialBlock;i>=n+legOver2-1;i-=legOver2) {
				for (j=0;j<legOver2;j++) tmpBlock[j]=i-j;
				Y_.push_back(tmpBlock);
				std::ostringstream msg;
				msg<<"Y added block of size "<<tmpBlock.size();
				progress_.printline(msg,std::cerr);
			}
			std::ostringstream msg;
			msg<<"Y_ total size is "<<Y_.size();
			progress_.printline(msg,std::cerr);
			
			for (i=0;i<sizeOfInitialBlock;i++) E_.push_back(2*n-1-i);
			std::ostringstream msg2;
			msg2<<"E total size is "<<E_.size();
			progress_.printline(msg2,std::cerr);
		}

		void createHoppings(
					psimag::Matrix<Field>& connectors,
					const std::vector<psimag::Matrix<Field> >& hoppingsOneSite,
					size_t linSize,
					size_t numberOfOrbitals,
					size_t dof)
		{
			size_t legOfLadder=leg_;
			connectors.resize(linSize*dof,linSize*dof);
			for (size_t x=0;x<connectors.n_row();x++) 
				for (size_t y=0;y<connectors.n_col();y++) 
					connectors(x,y)=0;
			for (size_t x=0;x<size_t(linSize/legOfLadder);x++) {
				for (size_t y=0;y<legOfLadder;y++) {
					size_t i = y +x*legOfLadder;
					for (size_t orb1=0;orb1<numberOfOrbitals;orb1++) {
						for (size_t orb2=0;orb2<numberOfOrbitals;orb2++) {
							if (y+1<legOfLadder) {
								size_t j=y+1 + x*legOfLadder;
								size_t direction=1;
								connectors(orb1+i*dof,orb2+j*dof)=
								connectors(orb2+j*dof,orb1+i*dof)=hoppingsOneSite[direction](orb1,orb2);
							}
							if (x+1<size_t(linSize/legOfLadder)) {
								size_t j=y + (x+1)*legOfLadder;
								size_t direction=0;
								connectors(orb1+i*dof,orb2+j*dof)=
								connectors(orb2+j*dof,orb1+i*dof)=hoppingsOneSite[direction](orb1,orb2);
							}
							if (x+1<size_t(linSize/legOfLadder) && y+1<legOfLadder) {
								size_t j=y+1 + (x+1)*legOfLadder;
								size_t direction=2;
								connectors(orb1+i*dof,orb2+j*dof)=
								connectors(orb2+j*dof,orb1+i*dof)=hoppingsOneSite[direction](orb1,orb2);
							}
							if (x+1<size_t(linSize/legOfLadder) && y>0) {
								size_t j=y-1 + (x+1)*legOfLadder;
								size_t direction=3;
								connectors(orb1+i*dof,orb2+j*dof)=
								connectors(orb2+j*dof,orb1+i*dof)=hoppingsOneSite[direction](orb1,orb2);
							}
							
						}
					}
				}
			}
			// add spin down (=spin up and hoppings is diagonal in spin)
			for (size_t i=0;i<linSize;i++) for (size_t j=0;j<linSize;j++) 
				for (size_t orb1=0;orb1<numberOfOrbitals;orb1++) 
						for (size_t orb2=0;orb2<numberOfOrbitals;orb2++) 
							connectors(orb1+1*numberOfOrbitals+i*dof,orb2+1*numberOfOrbitals+j*dof)=
								connectors(orb1+0*numberOfOrbitals+i*dof,orb2+0*numberOfOrbitals+j*dof);
		}

		
		//! Make sure no one tries to copy construct or assign objects of this class
		GeometryLadderFeAs(const GeometryLadderFeAs&);

		GeometryLadderFeAs& operator=(const GeometryLadderFeAs&);		
	}; // class GeometryLadderFeAs
} // namespace Dmrg

/*@}*/
#endif
