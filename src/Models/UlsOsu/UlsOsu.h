/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/** \ingroup DMRG */
/*@{*/

/*! \file UlsOsu.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_ULSOSU_H
#define DMRG_ULSOSU_H
#include "../Models/UlsOsu/ParametersModelUlsOsu.h"
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"
#include "Utils.h"
#include "Complex.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class UlsOsu : public ModelBaseType {

    enum InternalDir {DIR_SX, DIR_SY, DIR_SZ, DIR_LX, DIR_LY, DIR_LZ};

public:

    typedef typename ModelBaseType::VectorRealType VectorRealType;
    typedef ModelHubbard<ModelBaseType> ModelHubbardType;
    typedef typename ModelBaseType::ModelHelperType ModelHelperType;
    typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
    typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
    typedef typename ModelBaseType::LinkType LinkType;
    typedef typename ModelHelperType::OperatorsType OperatorsType;
    typedef typename OperatorsType::OperatorType OperatorType;
    typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
    typedef typename ModelHelperType::RealType RealType;
    typedef typename ModelBaseType::QnType QnType;
    typedef typename QnType::VectorQnType VectorQnType;
    typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
    typedef typename SparseMatrixType::value_type ComplexOrRealType;
    typedef typename ModelBaseType::HilbertBasisType HilbertBasisFeAsType;
    typedef typename HilbertBasisFeAsType::value_type HilbertStateFeAs;
    typedef HilbertSpaceFeAs<HilbertStateFeAs> HilbertSpaceFeAsType;
    typedef	typename ModelBaseType::MyBasis BasisType;
    typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
    typedef typename ModelHubbardType::HilbertState HilbertStateType;
    typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
    typedef typename ModelHelperType::BlockType BlockType;
    typedef typename ModelBaseType::SolverParamsType SolverParamsType;
    typedef typename ModelBaseType::VectorType VectorType;
    typedef typename ModelBaseType::InputValidatorType InputValidatorType;
    typedef typename OperatorType::PairType PairType;
    typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
    typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
    typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
    typedef typename ModelBaseType::OpsLabelType OpsLabelType;
    typedef typename ModelBaseType::OpForLinkType OpForLinkType;
    typedef typename ModelBaseType::ModelTermType ModelTermType;

    static const int FERMION_SIGN = -1;
    enum {STATE_EMPTY = 0, STATE_UP_A = 1, STATE_DOWN_A = 4};
    enum {SPIN_UP, SPIN_DOWN};

    UlsOsu(const SolverParamsType& solverParams,
               InputValidatorType& io,
               const SuperGeometryType& geometry)
        : ModelBaseType(solverParams, geometry, io),
          modelParameters_(io),
          superGeometry_(geometry),
          offset_(6), // Sx, Sy, Sz, Lx, Ly, Lz
          spinSquared_(spinSquaredHelper_,modelParameters_.orbitals,2*modelParameters_.orbitals)
    {
        if (modelParameters_.orbitals != 1)
            throw PsimagLite::RuntimeError("UlsOsu: must use Orbital=1 \n");

        SizeType n = superGeometry_.numberOfSites();
        SizeType mx = modelParameters_.magneticFieldX.size();
        SizeType my = modelParameters_.magneticFieldY.size();
        SizeType mz = modelParameters_.magneticFieldZ.size();
        SizeType m = mz;

        if (mx != my || my != mz || mz != mx ) {
            PsimagLite::String msg("tJKitaev: If provided, ");
            msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
            msg += " MagneticFieldX, MagneticFieldY, MagneticFieldZ must be ";
            msg += "provided in all 3 (x,y,z) directions.\n";
            throw PsimagLite::RuntimeError(msg);
        }

        if (m > 0 && m != n) {
            PsimagLite::String msg("Kitaev: If provided, ");
            msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
            throw PsimagLite::RuntimeError(msg);
        }

        if (BasisType::useSu2Symmetry())
            err("UlsOsu does not have SU(2) implimentation \n");

        // fill caches
        ProgramGlobals::init(superGeometry_.numberOfSites() + 1);
        BlockType block(1,0);
        setNaturalBasis(basis_, block, true);
    }

    //! Find c^\dagger_isigma in the natural basis natBasis
    SparseMatrixType findCreationMatrices(int,
                                          SizeType sigma,
                                          const VectorHilbertStateType&) const
    {
        assert(sigma < creationMatrix_.size());
        return creationMatrix_[sigma].getCRS();
    }

    SizeType maxElectronsOneSpin() const
    {
        return 1*superGeometry_.numberOfSites() + 1;
    }

    void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
    {
        if (!io.doesGroupExist(label1))
            io.createGroup(label1);

        PsimagLite::String label = label1 + "/" + this->params().model;
        io.createGroup(label);
        modelParameters_.write(label, io);
        io.write(label + "/offset_", offset_);
        spinSquaredHelper_.write(label, io);
        spinSquared_.write(label, io);
        io.write(label + "/basis_", basis_);
        io.write(label + "/qq_", qq_);
        io.write(label + "/creationMatrix_", creationMatrix_);
    }

protected:

    void fillLabeledOperators(VectorQnType& qns)
    {
        SizeType site = 0;
        VectorSizeType block(1, site);
        VectorHilbertStateType natBasis;
        SparseMatrixType tmpMatrix;
        typename MatrixType::value_type dummy = 0.0;
        setNaturalBasis(natBasis, block, false);
        OpsLabelType& sx = this->createOpsLabel("sx");
        OpsLabelType& sy = this->createOpsLabel("sy");
        OpsLabelType& sz = this->createOpsLabel("sz");
        OpsLabelType& lx = this->createOpsLabel("lx");
        OpsLabelType& ly = this->createOpsLabel("ly");
        OpsLabelType& lz = this->createOpsLabel("lz");

        this->makeTrackable("sx");
        this->makeTrackable("sy");
        this->makeTrackable("sz");
        this->makeTrackable("lx");
        this->makeTrackable("ly");
        this->makeTrackable("lz");

        // Set local spin and charge matrices
        for (SizeType i=0;i<block.size();i++) {

            typename OperatorType::Su2RelatedType su2related;

            // Sx
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_SX, dummy);
            OperatorType myOp1(tmpMatrix,
                              ProgramGlobals::FermionOrBosonEnum::BOSON,
                              PairType(0, 0),
                              1.0,
                              su2related);
            sx.push(myOp1);

            // Sy
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_SY, dummy);
            OperatorType myOp2(tmpMatrix,
                               ProgramGlobals::FermionOrBosonEnum::BOSON,
                               PairType(0, 0),
                               1.0,
                               su2related);
            sy.push(myOp2);

            // Sz
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_SZ, dummy);
            OperatorType myOp3(tmpMatrix,
                               ProgramGlobals::FermionOrBosonEnum::BOSON,
                               PairType(0, 0),
                               1.0,
                               su2related);
            sz.push(myOp3);


            // Lx
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_LX, dummy);
            OperatorType myOp4(tmpMatrix,
                              ProgramGlobals::FermionOrBosonEnum::BOSON,
                              PairType(0, 0),
                              1.0,
                              su2related);
            lx.push(myOp4);

            // Ly
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_LY, dummy);
            OperatorType myOp5(tmpMatrix,
                               ProgramGlobals::FermionOrBosonEnum::BOSON,
                               PairType(0, 0),
                               1.0,
                               su2related);
            ly.push(myOp5);

            // Lz
            tmpMatrix = findSdirMatrices(i, natBasis, DIR_LZ, dummy);
            OperatorType myOp6(tmpMatrix,
                               ProgramGlobals::FermionOrBosonEnum::BOSON,
                               PairType(0, 0),
                               1.0,
                               su2related);
            lz.push(myOp6);
        }
    }

    void fillModelLinks()
    {
        OpForLinkType sx("sx");
        OpForLinkType sy("sy");
        OpForLinkType sz("sz");
        OpForLinkType lx("lx");
        OpForLinkType ly("ly");
        OpForLinkType lz("lz");

        ModelBaseType::createTerm("sxsx").push(sx, 'N', sx, 'N');
        ModelBaseType::createTerm("sysy").push(sy, 'N', sy, 'N');
        ModelBaseType::createTerm("szsz").push(sz, 'N', sz, 'N');

        ModelBaseType::createTerm("lxlx").push(lx, 'N', lx, 'N');
        ModelBaseType::createTerm("lyly").push(ly, 'N', ly, 'N');
        ModelBaseType::createTerm("lzlz").push(lz, 'N', lz, 'N');
    }

private:

    SparseMatrixType findSdirMatrices(SizeType,
                                      const HilbertBasisType&,
                                      InternalDir,
                                      RealType) const
    {
        err("Kitaev needs useComplex in SolverOptions in the input file\n");
        throw PsimagLite::RuntimeError("FATAL\n");
    }

    SparseMatrixType findSdirMatrices(SizeType,// site,
                                      const HilbertBasisType& natBasis,
                                      InternalDir dir,
                                      std::complex<RealType>) const
    {
        SizeType total = natBasis.size();
        double val = sqrt(2.0);
        MatrixType Ax(3,3), Ay(3,3), Az(3,3), Iden(3,3);
        Ax.setTo(0.0); Ay.setTo(0.0); Az.setTo(0.0); Iden.setTo(0.0);

        // -- Spin-1 matricies --
        Ax(0,1) = 1.0/val; Ax(1,0)= 1.0/val;
        Ax(1,2) = 1.0/val; Ax(2,1) = 1.0/val;

        Ay(0,1) = 1.0/val * std::complex<RealType>(0.0, -1.0);
        Ay(1,2) = 1.0/val * std::complex<RealType>(0.0, -1.0);
        Ay(1,0) = 1.0/val * std::complex<RealType>(0.0,  1.0);
        Ay(2,1) = 1.0/val * std::complex<RealType>(0.0,  1.0);

        Az(0,0) =  1.0;
        Az(2,2) = -1.0;

        Iden(0,0) = 1.0;
        Iden(1,1) = 1.0;
        Iden(2,2) = 1.0;


        // -- Spin-1 matricies --
        MatrixType cm(total,total);
        cm.setTo(0.0);
        assert(total == 9);

        if (dir == DIR_SX) {
            outerProduct(cm,Ax,Iden);
        } else if (dir == DIR_SY) {
            outerProduct(cm,Ay,Iden);
        } else if (dir == DIR_SZ) {
            outerProduct(cm,Az,Iden);
        } else if (dir == DIR_LX) {
            outerProduct(cm,Iden,Ax);
        } else if (dir == DIR_LY) {
            outerProduct(cm,Iden,Ay);
        } else if (dir == DIR_LZ) {
            outerProduct(cm,Iden,Az);
        } else {
            assert(false);
        }

        SparseMatrixType operatorMatrix(cm);
        return operatorMatrix;
    }

    void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
                                    const BlockType& block,
                                    RealType) const
    {
        SizeType n=block.size();
        SizeType linSize = superGeometry_.numberOfSites();
        if (modelParameters_.magneticFieldX.size() != linSize)
            return; // <<---- PLEASE NOTE EARLY EXIT HERE
        if (modelParameters_.magneticFieldY.size() != linSize)
            return; // <<---- PLEASE NOTE EARLY EXIT HERE
        if (modelParameters_.magneticFieldZ.size() != linSize)
            return; // <<---- PLEASE NOTE EARLY EXIT HERE

        for (SizeType i = 0; i < n; ++i) {
            SizeType site = block[i];

            // magnetic field x
            const OperatorType& sx = ModelBaseType::naturalOperator("sx", site, 0);
            RealType tmp = modelParameters_.magneticFieldX[block[0]];
            hmatrix += tmp*sx.getCRS();

            // magnetic field y
            const OperatorType& sy = ModelBaseType::naturalOperator("sy", site, 0);
            tmp = modelParameters_.magneticFieldY[block[0]];
            hmatrix += tmp*sy.getCRS();

            // magnetic field z
            const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
            tmp = modelParameters_.magneticFieldZ[block[0]];
            hmatrix += tmp*sz.getCRS();
        }
    }

    void setNaturalBasis(HilbertBasisType& basis,
                         const VectorSizeType& block,
                         bool truncated) const
    {
        assert(block.size()==1);
        assert(modelParameters_.orbitals == 1);
        HilbertStateType total = 9;

        basis.resize(total);
        //for (SizeType a = 0; a< total; ++a) basis[a] = a;
        if (modelParameters_.orbitals == 1 && basis.size() == 9) {
            basis[0] = 0; // e - e
            basis[1] = 1; // e - u
            basis[2] = 2; // e - d

            basis[3] = 3; // u - e
            basis[4] = 4; // u - u
            basis[5] = 5; // u - d

            basis[6] = 6; // d - e
            basis[7] = 7; // d - u
            basis[8] = 8; // d - d
        }
    }


    ParametersModelUlsOsu<RealType, QnType>  modelParameters_;
    const SuperGeometryType& superGeometry_;
    SizeType offset_;
    SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
    SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;
    HilbertBasisType basis_;
    VectorQnType qq_;
    VectorOperatorType creationMatrix_;
};	//class UlsOsu

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_Anisotropic_H

