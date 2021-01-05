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

#include "tm/Elements/Interfaces/intelsurfquadlin_ht.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(IntElSurfQuadLin_ht);

  FEI2dQuadLin IntElSurfQuadLin_ht :: interpolation(1,2);

IntElSurfQuadLin_ht :: IntElSurfQuadLin_ht(int n, Domain *aDomain) :
    TransportInterfaceElement(n, aDomain)
{
    numberOfDofMans = 8;
}


void
IntElSurfQuadLin_ht :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    this->interpolation.evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 8);
    answer.zero();

    answer.at(1, 1) = -N.at(1);
    answer.at(1, 2) = -N.at(2);
    answer.at(1, 3) = -N.at(3);
    answer.at(1, 4) = -N.at(4);
    answer.at(1, 5) =  N.at(1);
    answer.at(1, 6) =  N.at(2);
    answer.at(1, 7) =  N.at(3);
    answer.at(1, 8) =  N.at(4);
}



void
IntElSurfQuadLin_ht :: computeFullNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    this->interpolation.evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 8);
    answer.zero();

    answer.at(1, 1) = -N.at(1);
    answer.at(1, 2) = -N.at(2);
    answer.at(1, 3) = -N.at(3);
    answer.at(1, 4) = -N.at(4);
    answer.at(1, 5) =  N.at(1);
    answer.at(1, 6) =  N.at(2);
    answer.at(1, 7) =  N.at(3);
    answer.at(1, 8) =  N.at(4);


    answer.at(2, 1) = -N.at(1);
    answer.at(2, 2) = -N.at(2);
    answer.at(2, 3) = -N.at(3);
    answer.at(2, 4) = -N.at(4);
    answer.at(2, 5) =  N.at(1);
    answer.at(2, 6) =  N.at(2);
    answer.at(2, 7) =  N.at(3);
    answer.at(2, 8) =  N.at(4);



    answer.at(3, 1) = -N.at(1);
    answer.at(3, 2) = -N.at(2);
    answer.at(3, 3) = -N.at(3);
    answer.at(3, 4) = -N.at(4);
    answer.at(3, 5) =  N.at(1);
    answer.at(3, 6) =  N.at(2);
    answer.at(3, 7) =  N.at(3);
    answer.at(3, 8) =  N.at(4);

}


void
IntElSurfQuadLin_ht :: computeGaussPoints()
{
    // Sets up the array of Gauss Points of the receiver.
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        //integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>(1,domain, 1, 2);
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3); 
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Square, 4, _3dHeInterface); 
    }
}


void
IntElSurfQuadLin_ht :: computeCovarBaseVectorsAt(IntegrationPoint *ip, FloatArray &G1, FloatArray &G2)
{
    FloatMatrix dNdxi;
    this->interpolation.evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    G1.resize(3);
    G2.resize(3);
    G1.zero();
    G2.zero();

    FloatArray meanNode;
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        meanNode = 0.5 * ( this->giveNode(i)->giveCoordinates() + this->giveNode(i + numNodes / 2)->giveCoordinates() );
        G1 += dNdxi.at(i, 1) * meanNode;
        G2 += dNdxi.at(i, 2) * meanNode;
    }
}


double
IntElSurfQuadLin_ht :: computeAreaAround(IntegrationPoint *ip)
{
    FloatArray G1, G2, G3;
    this->computeCovarBaseVectorsAt(ip, G1, G2);
    double weight  = ip->giveWeight();
    G3.beVectorProductOf(G1, G2);
    return G3.computeNorm() * weight;
}


void
IntElSurfQuadLin_ht :: initializeFrom(InputRecord &ir)
{
    TransportInterfaceElement :: initializeFrom(ir);
}


void
IntElSurfQuadLin_ht :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { T_f };
}

#if 0
bool
IntElSurfQuadLin_ht :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    computeTransformationMatrixAt(this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), lcs);
    answer.resize(18, 18);
    for ( int i = 0; i < 6; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i * 3 + 1, i * 3 + j) = lcs.at(3, j);
            answer.at(i * 3 + 2, i * 3 + j) = lcs.at(1, j);
            answer.at(i * 3 + 3, i * 3 + j) = lcs.at(2, j);
        }
    }
    
    
    return 1;
}
#endif

void
IntElSurfQuadLin_ht :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    FloatArray G1, G2, Normal;
    this->computeCovarBaseVectorsAt(gp, G1, G2);
    Normal.beVectorProductOf(G1, G2);
    Normal.normalize();
    answer.initFromVector(Normal,false);
    //xLocalCoordSys(Normal);
}



int
IntElSurfQuadLin_ht :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatArray N, meanNode;
    this->interpolation.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        meanNode = 0.5 * ( this->giveNode(i)->giveCoordinates() + this->giveNode(i + 3)->giveCoordinates() );
        answer += N.at(i) * meanNode;
    }

    return 1;
}


bool
IntElSurfQuadLin_ht :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    OOFEM_ERROR("Not implemented");
    return false;
}

} // end namespace oofem
