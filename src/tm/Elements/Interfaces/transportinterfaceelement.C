/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "tm/Elements/Interfaces/transportinterfaceelement.h"
#include "tm/Materials/InterfaceMaterials/transportinterfacematerial.h"

#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"



namespace oofem {
TransportInterfaceElement :: TransportInterfaceElement(int n, Domain *aDomain) : Element(n, aDomain)
{
}

int TransportInterfaceElement :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.resize(this->giveDofManager(1)->giveCoordinates().giveSize());
    answer.zero();

    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= numNodes/2; i++ ) {
        const auto &nodeCoord = this->giveDofManager(i)->giveCoordinates();
        answer.add(N.at(i), nodeCoord );
    }

    return true;
}


  
    
void
TransportInterfaceElement :: computeConductivityMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
  {
    FloatArray T, dT, F, dF;
    this->computeVectorOf(VM_Total, tStep, T);
    this->giveInternalForcesVector_num(F, tStep, T);
    double eps = 1.e-7;
    for(int i = 1; i <= 8; i++) {
      dT = T;
      dT.at(i) += eps;
      this->giveInternalForcesVector_num(dF, tStep, dT);
      for(int j = 1; j <= 8; j++) {
	answer.at(i,j) = dF.at(j) - F.at(j);
      }
      
    }
    answer.times(1./eps);
    
    
  }

  


void
TransportInterfaceElement :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // Computes the stiffness matrix of the receiver K_cohesive = int_A ( N^t * dT/dj * N ) dA
    FloatMatrix N, Nf, D, DN;
    //  bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    bool matStiffSymmFlag = true;
    answer.clear();

    FloatMatrix rotationMatGtoL, Dn;
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
      D = static_cast< TransportInterfaceMaterial * >( this->giveMaterial() )->giveConductivityMatrix(rMode, ip, tStep);
      
        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
	Dn.beProductOf(rotationMatGtoL,D);
        //D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system
        this->computeNmatrixAt(ip, N);
	this->computeFullNmatrixAt(ip, Nf);
        DN.beProductOf(Dn, N);
        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(Nf, DN, dA);
        } else {
            answer.plusProductUnsym(Nf, DN, dA);
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }

    /*
    FloatMatrix A(8,8);
    this->computeConductivityMatrix_num(A, rMode, tStep);
    */
}


void
TransportInterfaceElement :: computeTemperatureJump(FloatArray &answer, IntegrationPoint *ip, TimeStep *tStep)
{
    // Computes the spatial jump vector at the Gauss point (ip) of
    // the receiver, at time step (tStep). jump = N*u
    FloatMatrix N;
    FloatArray T;

    if ( !this->isActivated(tStep) ) {
        //match dimensions
        //answer.resize(3);
        this->computeNmatrixAt(ip, N);
        answer.resize(N.giveNumberOfRows());
        answer.zero();
        return;
    }

    this->computeNmatrixAt(ip, N);
    this->computeVectorOf(VM_Total, tStep, T);

    answer.beProductOf(N, T);
}


void
TransportInterfaceElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces
    // if useGpRecord == 1 then data stored in ip->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

  FloatMatrix N, Nf;
    FloatArray T, flux, jump;

    this->computeVectorOf(VM_TotalIntrinsic, tStep, T); //u in GCS (LCS only if defined computeGtoLRotationMatrix() )
   
    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        this->computeNmatrixAt(ip, N);
	this->computeFullNmatrixAt(ip, Nf);
	jump.beProductOf(N, T);
	this->computeFluxVector(flux, ip, jump, tStep);
        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(Nf, flux, dA);
    }

}


void
TransportInterfaceElement :: giveInternalForcesVector_num(FloatArray &answer, TimeStep *tStep, const FloatArray &T)
{
  // Computes internal forces
    // if useGpRecord == 1 then data stored in ip->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

    FloatMatrix N, Nf;
    FloatArray flux, jump;
   
    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        this->computeNmatrixAt(ip, N);
	this->computeFullNmatrixAt(ip, Nf);
	jump.beProductOf(N, T);
	this->computeFluxVector(flux, ip, jump, tStep);
        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(Nf, flux, dA);
    }

}




  
void
TransportInterfaceElement :: computeFluxVector(FloatArray &flux, IntegrationPoint *ip, const FloatArray &jump, TimeStep *tStep)
{
    TransportInterfaceMaterial *mat = static_cast< TransportInterfaceMaterial* >( this->giveMaterial() );
    // Returns the traction in global coordinate system
    FloatMatrix rotationMatGtoL;
    this->computeTransformationMatrixAt(ip, rotationMatGtoL);

    FloatArray jumpRot = jump;
    jumpRot.rotatedWith(rotationMatGtoL, 'n');      // transform jump to local coord system
    flux = mat->giveFluxVector(jump, ip, tStep);
    flux.rotatedWith(rotationMatGtoL,'n');     // transform traction to global coord system
}


void
TransportInterfaceElement :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{

  if ( mtrx == TangentStiffnessMatrix ) {
    FloatMatrix tmp;
    this->computeConductivityMatrix(answer, Conductivity, tStep);
    //   this->computeBCMtrxAt(tmp, tStep, VM_TotalIntrinsic);
    //answer.add(tmp);
  } else if ( mtrx == ConductivityMatrix ) {
    this->computeConductivityMatrix(answer, Conductivity, tStep);
  } else if ( mtrx == CapacityMatrix || mtrx == MassMatrix ) {
    //do nothing this->computeCapacityMatrix(answer, tStep);
  } else if ( mtrx == LumpedMassMatrix ) {
    //do nothing
  } else {
    OOFEM_ERROR("Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
  }
}


void
TransportInterfaceElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                                       TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
  if ( mtrx == InternalForcesVector ) {
    this->giveInternalForcesVector(answer, tStep);
  } else if ( mtrx == InertiaForcesVector ) {
    //do noting this->computeInertiaForcesVector(answer, tStep);
    answer.zero();
    } else if ( mtrx == LumpedMassMatrix ) {
    answer.zero();
    //do nothung this->computeLumpedCapacityVector(answer, tStep);
  } else if ( mtrx == ExternalForcesVector ) {
    answer.zero();    
  } else {
    OOFEM_ERROR( "Unknown Type of characteristic mtrx (%s)",
		 __CharTypeToString(mtrx) );
  }
}

void
TransportInterfaceElement :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}


void
TransportInterfaceElement :: updateInternalState(TimeStep *tStep)
{
    // Updates the receiver at end of step.
    FloatArray tractionG, jumpL;

    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            this->computeTemperatureJump(jumpL, gp, tStep);
            this->computeFluxVector(tractionG, gp, jumpL, tStep);
        }
    }
}



int
TransportInterfaceElement :: giveIPValue(FloatArray &answer, IntegrationPoint *aIntegrationPoint, InternalStateType type, TimeStep *tStep)
{
    return Element :: giveIPValue(answer, aIntegrationPoint, type, tStep);
}


void
TransportInterfaceElement :: initializeFrom(InputRecord &ir)
{
    Element :: initializeFrom(ir);
}

void TransportInterfaceElement :: giveInputRecord(DynamicInputRecord &input)
{
    Element :: giveInputRecord(input);
}





} // end namespace oofem

