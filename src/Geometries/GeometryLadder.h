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

/*! \file GeometryLadder.h
 *
 *  This class implements DmrgGeometryBase for n-leg ladders
 */
#ifndef GEOMETRYLADDER
#define GEOMETRYLADDER

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
	//! RESTRICTION: leg must be even (tested only for leg==2)
	//! RESTIRCTION: homogeneous system only!!
	template<typename Field,typename ConnectorsType_>
	class GeometryLadder : public GeometryBase<Field,ConnectorsType_> {
		static const size_t SystemSystem=ProgramGlobals::SYSTEM_SYSTEM;
		static const size_t SystemEnviron=ProgramGlobals::SYSTEM_ENVIRON;
		static const size_t EnvironSystem=ProgramGlobals::ENVIRON_SYSTEM;
		static const size_t EnvironEnviron=ProgramGlobals::ENVIRON_ENVIRON;
	public:
		typedef ConnectorsType_ ConnectorsType;
		//static const int SystemSystem=0,SystemEnviron=1,EnvironSystem=2,EnvironEnviron=3;
		typedef  typename GeometryBase<Field,ConnectorsType>::BlockType BlockType;
		
				
		GeometryLadder(ConnectorsType& connectors,int sizeOfInitialBlock,int leg):
			connectors_(connectors), leg_(leg),progress_("GeometryLadder",0)
		{
			split(sizeOfInitialBlock);
// 			std::ostringstream msg;
// 			msg<<"Connectors:\n";
// 			size_t n = connectors_.n_row();
// 			for (size_t i=0;i<n;i++) {
// 				for (size_t j=0;j<n;j++) {
// 					msg<<connectors_(i,j);
// 				}
// 				msg<<"\n";
// 			}
// 			msg<<"End of Connectors\n";
// 			progress_.printline(msg,std::cerr);
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

		Field calcConnectorValue(int type,int ind,int jnd,int smax,int emin,size_t what) const 
		{
			//! There are four cases:
			//! 1. (ind,jnd) in SUX --> use input connectors
			//! 2. (ind,jnd) in YUE --> use reflected connectors
			//! 3. ind in SUX , jnd in YUE --> delegate
			//! 4. ind in YUE, jnd in SUX --> delegate

			Field x=0;
			switch (type) {
				case EnvironEnviron:
					x= connectors_(findReflection(ind),findReflection(jnd));
					break;

				case EnvironSystem:
					x=connector(ind,jnd,smax,emin,EnvironSystem);
					break;

				case SystemEnviron:
					x=connector(ind,jnd,smax,emin,SystemEnviron);
					break;
				case SystemSystem:
					x= connectors_(ind,jnd);
					break;
			}
			return x;
		}

		int calcConnectorType(int ind,int jnd) const { return  calcConnectorType(ind,jnd,systemBlock_); }

		//! given i in the environment returns the site symmetric to i (in the system)
		int findReflection(int i) const
		{
			int r=(i%leg_)-int(leg_/2);
			return 2*connectors_.n_row()-i+2*r;
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
		int leg_;
		BlockType S_,E_;
		std::vector<BlockType> X_,Y_;
		BlockType systemBlock_;
		ProgressIndicator progress_;
		
		Field connector(int ind,int jnd,int smax,int emin,int type) const
		{
			int n = smax+1;
			if (type==EnvironSystem) {// ind in env, jnd in the system
				if (n%leg_==0) return connectorComplete(ind,jnd,smax,emin,type);
				return connectorIncomplete(ind,jnd,smax,emin,type);
				
			} else { // ind in system, jnd in env
				return connector(jnd,ind,smax,emin,EnvironSystem);
			}
		}
		
		// ind in env, jnd in the system
		Field connectorIncomplete(int ind,int jnd,int smax,int emin,int type) const
		{
			Field t=connectors_.defaultValue(1); // y-direction
			if (size_t(smax+1)<connectors_.n_row()) t = connectors_(smax,smax+1);
			if (jnd==smax && ind==emin) return t;

			t=connectors_.defaultValue(0); // x-direction
			/*int indR = findReflection(ind); 
			if (jnd<int(connectors_.n_row()) && indR<int(connectors_.n_row()))
				t=connectors_(jnd,indR);
			*/
			for (int j=0;j<leg_;j++) 
				if (jnd==smax-j && ind==emin+leg_-1-j) return t;

			return static_cast<Field>(0);
		}
		
		// ind in env, jnd in the system
		Field connectorComplete(int ind,int jnd,int smax,int emin,int type) const
		{
			Field t=connectors_.defaultValue(0);
			//int indR = findReflection(ind);
			int jndR = findReflection(jnd);
			//if (jnd<int(connectors_.n_row()) && indR<int(connectors_.n_row()))
			//	t=connectors_(jnd,indR);

			for (int j=smax;j>smax-leg_;j--) {
				if (j<0) break;
				if (jnd==j && ind==jndR) return t;
			}
			return static_cast<Field>(0);
		}
		
		void split(int sizeOfInitialBlock) 
		{
			
			int i,j;
			int legOver2=int(leg_/2);
			int n=connectors_.linSize();
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


		void splitOld(int sizeOfInitialBlock) 
		{
			int i;
			int n=connectors_.n_row();
			BlockType tmpBlock;
			if (sizeOfInitialBlock>n) 
				throw std::runtime_error("GeometryBase::split(): sizeofinitialblock="+
						utils::ttos(sizeOfInitialBlock)+" but n="+utils::ttos(n)+"\n");

			tmpBlock.resize(1);

			// S U X = system, and Y U E is the reflexion
			for (i=0;i<sizeOfInitialBlock;i++) {
				S_.push_back(i);
				systemBlock_.push_back(i);
			}
			for (i=sizeOfInitialBlock;i<n;i++) {
				tmpBlock[0]=i;
				X_.push_back(tmpBlock);
				systemBlock_.push_back(i);
			}

			for (i=2*n-1-sizeOfInitialBlock;i>=n;i--) {
				tmpBlock[0]=i;
				Y_.push_back(tmpBlock);
			}
			for (i=0;i<sizeOfInitialBlock;i++) E_.push_back(2*n-1-i);
		}

		//! Make sure no one tries to copy construct or assign objects of this class
		GeometryLadder(const GeometryLadder&);

		GeometryLadder& operator=(const GeometryLadder&);	
	};
} // namespace Dmrg

/*@}*/
#endif
