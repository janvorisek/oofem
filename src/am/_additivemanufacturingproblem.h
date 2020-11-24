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

#ifndef AdditiveManufacturingProblem_h
#define AdditiveManufacturingProblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "sparsenonlinsystemnm.h"

#include <memory>
#include <string>

///@name Input fields for AdditiveManufacturingProblem
//@{
#define _IFT_AdditiveManufacturingProblem_Name "additivemanufacturing"
#define _IFT_AdditiveManufacturingProblem_gcode "gcode"
#define _IFT_AdditiveManufacturingProblem_alpha "alpha" ///< Defines the time discretization of the value: @f$ u = u_n (1-\alpha) + u_{n+1} \alpha @f$.
#define _IFT_AdditiveManufacturingProblem_deltaT "deltat" ///< Fixed timestep.
#define _IFT_AdditiveManufacturingProblem_initt "initt" ///< Initial time
#define _IFT_AdditiveManufacturingProblem_dtFunction "dtfunction" ///< Function that determines size of time step.
#define _IFT_AdditiveManufacturingProblem_prescribedTimes "prescribedtimes" ///< Discrete times for each time step.
#define _IFT_AdditiveManufacturingProblem_keepTangent "keeptangent" ///< Fixes the tangent to be reused on each step.
#define _IFT_AdditiveManufacturingProblem_lumped "lumped" ///< Use of lumped "mass" matrix
#define _IFT_AdditiveManufacturingProblem_exportFields "exportfields" ///< Fields to export for staggered problems.
//@}

namespace oofem {
class DofDistributedPrimaryField;
class Function;

/**
 * Solves general nonlinear transient transport problems.
 * @author Mikael Öhman
 */
class AdditiveManufacturingProblem : public EngngModel
{
protected:
    SparseMtrxType sparseMtrxType = SMT_Skyline;
    std :: unique_ptr< DofDistributedPrimaryField > field;

    std :: unique_ptr< SparseMtrx > effectiveMatrix;

    FloatArray solution;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;

    double alpha = 0.5;
    int dtFunction = 0;
    FloatArray prescribedTimes;
    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    double deltaT = 1.;
    bool keepTangent = false, hasTangent = false;
    bool lumped = false;

    IntArray exportFields;

    // G-code filepath
    std::string gCodeFilePath = "";

public:
    AdditiveManufacturingProblem(int i, EngngModel *master=nullptr);

    void solveYourselfAt(TimeStep *tStep) override;
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    bool newDofHandling() override { return true; }
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    virtual void applyIC();

    int requiresUnknownsDictionaryUpdate() override;
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDomainLinks() override;

    Function *giveDtFunction();
    double giveDeltaT(int n);
    double giveDiscreteTime(int iStep);

    TimeStep *giveNextStep() override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;

    bool requiresEquationRenumbering(TimeStep *tStep) override;
    int forceEquationNumbering() override;

    void updateYourself(TimeStep *tStep) override;
    
    int checkConsistency() override;
    FieldPtr giveField (FieldType key, TimeStep *tStep) override;
    // identification
    const char *giveInputRecordName() const { return _IFT_AdditiveManufacturingProblem_Name; }
    const char *giveClassName() const override { return "AdditiveManufacturingProblem"; }
    fMode giveFormulation() override { return TL; }
};
} // end namespace oofem
#endif // AdditiveManufacturingProblem_h
