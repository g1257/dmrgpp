/*2:*/
#line 11 "DensityMatrix.w"

#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

/*:2*//*3:*/
#line 18 "DensityMatrix.w"

#include "Utils.h"

/*:3*//*4:*/
#line 27 "DensityMatrix.w"

#include "BlockMatrix.h"

/*:4*//*5:*/
#line 32 "DensityMatrix.w"

#include "DensityMatrixLocal.h"
#include "DensityMatrixSu2.h"

namespace Dmrg{

/*:5*//*6:*/
#line 44 "DensityMatrix.w"

template<
typename RealType,
typename DmrgBasisType,
typename DmrgBasisWithOperatorsType,
typename TargettingType
> 
class DensityMatrix{

/*:6*//*7:*/
#line 55 "DensityMatrix.w"

typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
typedef BlockMatrix<DensityMatrixElementType,psimag::Matrix<DensityMatrixElementType> > 
BlockMatrixType;
typedef typename DmrgBasisType::FactorsType FactorsType;
enum{EXPAND_SYSTEM= TargettingType::EXPAND_SYSTEM};
typedef DensityMatrixLocal<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
TargettingType> DensityMatrixLocalType;
typedef DensityMatrixSu2<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
TargettingType> DensityMatrixSu2Type;
typedef DensityMatrixBase<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
TargettingType> DensityMatrixBaseType;

/*:7*//*8:*/
#line 71 "DensityMatrix.w"

public:
typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;

/*:8*//*9:*/
#line 87 "DensityMatrix.w"

DensityMatrix(
const TargettingType&target,
const DmrgBasisWithOperatorsType&pBasis,
const DmrgBasisWithOperatorsType&pBasisSummed,
const DmrgBasisType&pSE,
size_t direction,
bool debug= false,
bool verbose= false)

/*:9*//*10:*/
#line 102 "DensityMatrix.w"

:densityMatrixLocal_(target,pBasis,pBasisSummed,pSE,direction,debug,verbose),
densityMatrixSu2_(target,pBasis,pBasisSummed,pSE,direction,debug,verbose)
{

/*:10*//*11:*/
#line 112 "DensityMatrix.w"

if(DmrgBasisType::useSu2Symmetry()){
densityMatrixImpl_= &densityMatrixSu2_;
}else{
densityMatrixImpl_= &densityMatrixLocal_;
}
/*:11*//*12:*/
#line 123 "DensityMatrix.w"

densityMatrixImpl_->init(target,pBasis,pBasisSummed,pSE,direction);
}

/*:12*//*13:*/
#line 135 "DensityMatrix.w"

BlockMatrixType&operator()()
{
return densityMatrixImpl_->operator()();
}

/*:13*//*14:*/
#line 145 "DensityMatrix.w"

size_t rank(){return densityMatrixImpl_->rank();}

/*:14*//*15:*/
#line 149 "DensityMatrix.w"

void check(int direction)
{
return densityMatrixImpl_->check(direction);
}

/*:15*//*16:*/
#line 156 "DensityMatrix.w"

void check2(int direction)
{
densityMatrixImpl_->check2(direction);
}

/*:16*//*17:*/
#line 165 "DensityMatrix.w"

template<typename ConcurrencyType> 
void diag(std::vector<RealType> &eigs,char jobz,ConcurrencyType&concurrency)
{
if(!DmrgBasisType::useSu2Symmetry()){
densityMatrixLocal_.diag(eigs,jobz,concurrency);
}else{
densityMatrixSu2_.diag(eigs,jobz,concurrency);
}

}

/*:17*//*18:*/
#line 179 "DensityMatrix.w"

template<
typename RealType_,
typename DmrgBasisType_,
typename DmrgBasisWithOperatorsType_,
typename TargettingType_
> 
friend std::ostream&operator<<(std::ostream&os,
const DensityMatrix<RealType_,
DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_> &dm);

/*:18*//*19:*/
#line 195 "DensityMatrix.w"

private:
DensityMatrixLocalType densityMatrixLocal_;
DensityMatrixSu2Type densityMatrixSu2_;
DensityMatrixBaseType*densityMatrixImpl_;


};

/*:19*//*20:*/
#line 207 "DensityMatrix.w"

template<
typename RealType,
typename DmrgBasisType,
typename DmrgBasisWithOperatorsType,
typename TargettingType
> 
std::ostream&operator<<(std::ostream&os,
const DensityMatrix<RealType,DmrgBasisType,
DmrgBasisWithOperatorsType,TargettingType> &dm)
{
os<<(*dm.densityMatrixImpl_);
return os;
}
}

#endif

/*:20*/
