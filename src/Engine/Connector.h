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

/*! \file Connector.h
 *
 *  FIXME
 *
 */
#ifndef CONNECTOR_H
#define CONNECTOR_H

namespace Dmrg {
	template<typename FieldType>
	class Connector {
	public:
		Connector(const psimag::Matrix<FieldType>& connectorMatrix,FieldType defaultConnector) 
		{
			connectorMatrix_=connectorMatrix;
			std::vector<psimag::Matrix<FieldType> > connectorsOneSite(1);
			psimag::Matrix<FieldType> tmp(1,1);
			tmp(0,0) = defaultConnector;
			connectorsOneSite[0] = tmp;
			createDefaultHoppings(defaultConnector_,connectorsOneSite,1);
			//std::cerr<<"DEF="<<defaultConnector_[0](0,0)<<"\n";
		}
		
		Connector(const psimag::Matrix<FieldType>& connectorMatrix,const std::vector<FieldType>& defaultConnector) 
		{
			connectorMatrix_=connectorMatrix;
			std::vector<psimag::Matrix<FieldType> > connectorsOneSite(defaultConnector.size());
			psimag::Matrix<FieldType> tmp(1,1);
			for (size_t i=0;i<defaultConnector.size();i++) {
				tmp(0,0) = defaultConnector[i];
				connectorsOneSite[i] = tmp;
			}
			createDefaultHoppings(defaultConnector_,connectorsOneSite,1);
			//std::cerr<<"DEF="<<defaultConnector_[0](0,0)<<"\n";
		}
		
		
		//For ladderFeAS
		Connector(	const std::vector<psimag::Matrix<FieldType> >& connectorsOneSite,
      				size_t linSize,
	  			size_t numberOfOrbitals,
				size_t dof,size_t leg) 
		{
		
			createHoppings(connectorMatrix_,connectorsOneSite,linSize,numberOfOrbitals,dof,leg);
			createDefaultHoppings(defaultConnector_,connectorsOneSite,numberOfOrbitals);
			
		}
		
		FieldType operator()(size_t i,size_t j) const
		{
			return connectorMatrix_(i,j);
		}
		
		FieldType defaultValue(size_t dir = 0,size_t a = 0,size_t b = 0) const
		{
			return defaultConnector_[dir](a,b);
		}
		
		size_t n_row() const { return connectorMatrix_.n_row(); }

	private:
		std::vector<psimag::Matrix<FieldType> > defaultConnector_;
		psimag::Matrix<FieldType>  connectorMatrix_;
		
		void createHoppings(
					psimag::Matrix<FieldType>& connectors,
					const std::vector<psimag::Matrix<FieldType> >& hoppingsOneSite,
					size_t linSize,
					size_t numberOfOrbitals,
					size_t dof,
				   	size_t legOfLadder)
		{
			size_t max = hoppingsOneSite.size();
			connectors.resize(linSize*dof,linSize*dof);
			for (size_t x=0;x<connectors.n_row();x++) 
				for (size_t y=0;y<connectors.n_col();y++) 
					connectors(x,y)=0;
			for (size_t x=0;x<size_t(linSize/legOfLadder);x++) {
				for (size_t y=0;y<legOfLadder;y++) {
					size_t i = y +x*legOfLadder;
					for (size_t orb1=0;orb1<numberOfOrbitals;orb1++) {
						for (size_t orb2=0;orb2<numberOfOrbitals;orb2++) {
							if (y+1<legOfLadder && max>1) {
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
							if (x+1<size_t(linSize/legOfLadder) && y+1<legOfLadder && max>1) {
								size_t j=y+1 + (x+1)*legOfLadder;
								size_t direction=2;
								connectors(orb1+i*dof,orb2+j*dof)=
								connectors(orb2+j*dof,orb1+i*dof)=hoppingsOneSite[direction](orb1,orb2);
							}
							if (x+1<size_t(linSize/legOfLadder) && y>0 && max>1) {
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
		
		void createDefaultHoppings(psimag::Matrix<FieldType>& matrixTmp,size_t direction,
					   const std::vector<psimag::Matrix<FieldType> >& hoppingsOneSite,
					  size_t numberOfOrbitals)
		{
			matrixTmp.resize(numberOfOrbitals*hoppingsOneSite[direction].n_row(),
					 numberOfOrbitals*hoppingsOneSite[direction].n_col());
			for (size_t i=0;i<matrixTmp.n_row();i++) for (size_t j=0;j<matrixTmp.n_col();j++) matrixTmp(i,j)=0;
			size_t row = hoppingsOneSite[direction].n_row();
			size_t col = hoppingsOneSite[direction].n_col();
			for (size_t i=0;i<row;i++) 
				for (size_t j=0;j<col;j++) 
					for (size_t spin=0;spin<numberOfOrbitals;spin++)
						matrixTmp(i+spin*row,j+spin*col)=hoppingsOneSite[direction](i,j);
		}

		void createDefaultHoppings(std::vector<psimag::Matrix<FieldType> >& dh,
					   const std::vector<psimag::Matrix<FieldType> >& hoppingsOneSite,
					  size_t numberOfOrbitals)
		{
			dh.clear();
			psimag::Matrix<FieldType> matrixTmp;
			size_t numberOfDirections  = hoppingsOneSite.size();
			for (size_t direction=0;direction<numberOfDirections;direction++) {
				createDefaultHoppings(matrixTmp,direction,hoppingsOneSite,numberOfOrbitals);
				dh.push_back(matrixTmp);
			}
		}
			
	}; //Connector
	
// 	template<typename FieldType>
// 	std::ostream& operator<<(std::ostream& os,Connector<FieldType>& connector)
// 	{
// 		os<<connector.connectorMatrix_;
// 		return os;
// 	}
} // namespace Dmrg
/*@}*/
#endif //TJ_ONEORBITAL_HEADER_H
