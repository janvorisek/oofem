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

#ifndef intelsurfquadlin_ht_h
#define intelsurfquadlin_ht_h

#include "tm/Elements/Interfaces/transportinterfaceelement.h"
#include "fei2dquadlin.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

#define _IFT_IntElSurfQuadLin_ht_Name "intelsurfquadlin_ht"

namespace oofem {
/**
 * This class implements 3d triangular surface interface element with linear interpolation.
 */
class IntElSurfQuadLin_ht : public TransportInterfaceElement
{
protected:
    static FEI2dQuadLin interpolation;

public:
    IntElSurfQuadLin_ht(int n, Domain *d);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;
    virtual void computeCovarBaseVectorsAt(IntegrationPoint *ip, FloatArray &G1, FloatArray &G2);
    //bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 8; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeAreaAround(IntegrationPoint *ip) override;

   
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_IntElSurfQuadLin_ht_Name; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_hexa_1; }


protected:
    void computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer) override;
    void computeFullNmatrixAt(GaussPoint *ip, FloatMatrix &answer) override;
    void computeGaussPoints() override;

};
} // end namespace oofem
#endif // intelsurfquadlin_ht_h
