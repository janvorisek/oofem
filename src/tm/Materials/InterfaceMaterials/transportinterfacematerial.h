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

#ifndef transportinterfacematerial_h
#define transportinterfacematerial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "matconst.h"
#include "matstatus.h"
#include "floatarrayf.h"

///@name Input fields for StructuralInterfaceMaterial
//@{

#define _IFT_TransportInterfaceMaterial_h "h"
#define _IFT_TransportInterfaceMaterial_Name "transportintermat"

//@}

namespace oofem {
class GaussPoint;


class TransportInterfaceMaterial : public Material
{
private:
  double h;
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    TransportInterfaceMaterial(int n, Domain * d);

    /**
     * Computes the first Piola-Kirchoff traction vector for given total jump/gap and integration point.
     * The total gap is computed from the displacement field at the given time step.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reducedF Deformation gradient in in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
  FloatMatrixF<1,1> giveConductivityMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
  FloatArrayF<1> giveFluxVector(FloatArrayF<1> jump, GaussPoint *gp, TimeStep *tStep) const;

    // identification and auxiliary functions
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    //virtual int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_TransportInterfaceMaterial_Name; }
    const char *giveClassName() const override { return "TransportInterfaceMaterial"; }

  
};
} // end namespace oofem
#endif // structuralinterfacematerial_h
