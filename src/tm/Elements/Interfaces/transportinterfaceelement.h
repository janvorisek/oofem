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

#ifndef transportinterfaceelement_h
#define transportinterfaceelement_h

#include "element.h"
#include "floatmatrix.h"
#include "function.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "integrationdomain.h"
#include "dofmantransftype.h"

namespace oofem {
class TimeStep;
class Node;
class TransportInterfaceMaterial;
class GaussPoint;
class FloatArray;
class IntArray;
class FEInterpolation;


/**
 * Abstract base class for all structural interface elements. It declares a common interface available
 * to all derived elements. The implementation of these services is partly left to the derived classes,
 * some general services are implemented here (but they can be overload by more efficient
 * element implementations).
 * The general implementation provided here is intended for both linear and nonlinear computations.
 */
class TransportInterfaceElement : public Element
{
protected:

public:
    /**
     * Constructor. Creates structural element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    TransportInterfaceElement(int n, Domain * d);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType, TimeStep *tStep) override;
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override;


  virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
  virtual void computeFluxVector(FloatArray &flux, IntegrationPoint *ip, const FloatArray &temperaturejump, TimeStep *tStep);
    virtual void computeTemperatureJump(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType) override { return nullptr; }

    //@}

    // Overloaded methods.
    void updateInternalState(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "TransportInterfaceElement"; }
   virtual double computeAreaAround(GaussPoint *gp) = 0;

    Element_Geometry_Type giveGeometryType() const override { return EGT_unknown; }

    //virtual methods that should be overloaded by the elements
  void computeConductivityMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
  void giveInternalForcesVector_num(FloatArray &answer, TimeStep *tStep, const FloatArray &T);



  void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
  

protected:
    /**
     * Computes modified interpolation matrix (N) for the element which multiplied
     * with the unknowns vector (u) produces the spatial jump.
     * @param gp Integration point for which answer is assembled.
     * @param answer Interpolation matrix evaluated at gp.
     */
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;
    virtual void computeFullNmatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;

    // transformation matrix from local to global
    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;

};
} // end namespace oofem
#endif // structuralinterfaceelement_h
