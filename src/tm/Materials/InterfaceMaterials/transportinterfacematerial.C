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

#include "transportinterfacematerial.h"
//#include "termalinterfacematerialstatus.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_Material(TransportInterfaceMaterial);

  
TransportInterfaceMaterial :: TransportInterfaceMaterial(int n, Domain *d) : Material(n, d)
{
}


int
TransportInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  //TermalInterfaceMaterialStatus *status = static_cast< TermalInterfaceMaterialStatus * >( this->giveStatus(gp) );
    /*if ( type == IST_InterfaceJump ) {
        answer = status->giveJump();
        return 1;
    } else if ( type == IST_InterfaceTraction ) {
        answer = status->giveTraction();
        return 1;
    } else {
        return Material :: giveIPValue(answer, gp, type, tStep);
	}*/
  return 1;
}


void
TransportInterfaceMaterial :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->h, _IFT_TransportInterfaceMaterial_h);
}


void
TransportInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
}


FloatMatrixF<1,1>
TransportInterfaceMaterial :: giveConductivityMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
  FloatMatrixF<1,1> d;
  d.at(1,1) = h;
  return d;
}



FloatArrayF<1>
TransportInterfaceMaterial :: giveFluxVector(FloatArrayF<1> jump, GaussPoint *gp, TimeStep *tStep) const
{
  return h * jump;
}



} // end namespace oofem
