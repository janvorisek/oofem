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

#include "additivemanufacturingproblem.h"
#include "printer.h"
#include "unknownnumberingscheme.h"
#include "masterdof.h"
#include "initialcondition.h"
#include "load.h"
#include "boundaryload.h"
#include "freeconstantsurfaceload.h"
#include "dynamicinputrecord.h"
#include "tm/EngineeringModels/transienttransportproblem.h"
#include "connectivitytable.h"
#include <vector>
#include <chrono>
#include <thread>

#include "engngm.h"
#include "timestep.h"
#include "function.h"
#include "metastep.h"
#include "exportmodulemanager.h"
#include "mathfem.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "verbose.h"
#include "classfactory.h"
#include "domain.h"

#include <stdlib.h>

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(AdditiveManufacturingProblem);

bool add_node_if_not_exists(EngngModel *emodel, const CellNode &cn) {
    Domain *d = emodel->giveDomain(1);

    int nNodesBefore = d->giveNumberOfDofManagers();
    int nodeId = cn.id;

    // Resize node list to accomodate new nodes
    //if(cn.id > nNodesBefore)
    //    d->resizeDofManagers(nodeId);

    // was already added?
    if(d->dofManagerList[cn.id-1] != nullptr) return false;
    
    std::unique_ptr< DofManager >dman( classFactory.createDofManager("node", nodeId, d ) );
    if ( !dman ) {
        OOFEM_ERROR("Couldn't create node %d\n", nodeId);
    }

    FloatArray coords = {cn.coords[0]/1000, cn.coords[1]/1000, cn.coords[2]/1000};
    //std::cout << "coords=" << coords[0] << ", " << coords[1] << ", " << coords[2] << "\n";

    dman->setCoordinates(coords);
    dman->setGlobalNumber(nodeId);
    dman->appendDof( new MasterDof(dman.get(), T_f) );

    //OOFEM_LOG_INFO("Added node %d\n", nodeId);

    //dman->giveDofWithID(T_f)->giveUnknowns()->add(0, 100.);
    //dman->giveDofWithID(T_f)->giveUnknowns()->add(1, 100.);
    d->setDofManager(nodeId, std::move(dman));

    return true;
}

void add_element_if_not_exists(EngngModel *emodel, const Cell &cn) {
    Domain *d = emodel->giveDomain(1);

    // Add new element
    int nElBefore = d->giveNumberOfElements();
    int elId = cn.id;

    //if(cn.id > nElBefore)
    //    d->resizeElements(elId);

    // was already added?
    if(d->elementList[cn.id-1] != nullptr) return;
    
    std::unique_ptr< Element >elem(classFactory.createElement("brick1ht", elId, d) );
    if ( !elem ) {
        OOFEM_ERROR("Couldn't create element %d\n", elId);
    } 
    
    IntArray enodes = {
        cn.nodes[0], cn.nodes[1], cn.nodes[2], cn.nodes[3],
        cn.nodes[4], cn.nodes[5], cn.nodes[6], cn.nodes[7]
    };

    elem->setDofManagers(enodes);
    elem->setGlobalNumber(elId);
    elem->setCrossSection(1);
    elem->setParallelMode(elementParallelMode::Element_local); // TODO: WHY SOME ELEMENTS HAVE SET REMOTE????
    elem->postInitialize(); // this is needed to allocate the integration point array

    d->setElement(elId, std::move(elem));

    // Add new node set
    // Used to apply initial boundary condition to the new DOFs
    int nSetsBefore = d->giveNumberOfSets();
    d->resizeSets(nSetsBefore + 6);

    // Add 6 loads (1 per element face)
    int nBCBefore = d->giveNumberOfBoundaryConditions();
    d->resizeBoundaryConditions(nBCBefore + 6);
    
    for(int i = 0; i < 6; i++) {
        int setId = nSetsBefore + 1 + i;
        
        std::unique_ptr<Set> set = std::make_unique<Set>(setId, d);
        IntArray bounds = {elId, 1 + i};
        set->setNumber(setId);
        set->setBoundaryList(bounds);
        //set->setNodeList(nNodeIds);
        d->setSet(setId, std::move(set));

        int bcid = nBCBefore + 1 + i;
        //std::cout << "new bc " << bcid << "\n";
        std::unique_ptr< GeneralBoundaryCondition >bc(classFactory.createBoundaryCondition("freeconstantsurfaceload", bcid, d) );   
        //auto surfLoad = dynamic_cast<FreeConstantSurfaceLoad*> (bc.get());

        // freeconstantsurfaceload 1 loadtimefunction 1 loadtype 3 set 2 components 1 500 properties 1 97 10

        DynamicInputRecord myInput = DynamicInputRecord(_IFT_FreeConstantSurfaceLoad_Name, bcid);
        myInput.setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
        myInput.setField(3, _IFT_BoundaryLoad_loadtype);
        //Dictionary props = Dictionary();
        //props.add('a', 10.); // 97 == 'a'
        //props.printYourself();
        //myInput.setField(props, _IFT_BoundaryLoad_properties);

        FloatArray comps = {25.};
        myInput.setField(comps, _IFT_Load_components);
        myInput.setField(setId, _IFT_GeneralBoundaryCondition_set);
        
        bc->initializeFrom(myInput);
        
        bc->setNumber(bcid);
        bc->postInitialize();
        //std::cout << "initialized BC type " << bc->giveType() << "\n" ;

        d->setBoundaryCondition(bcid, std::move(bc));
        //std::cout << "moved?" << d->giveBc(bcid)->giveNumber();
        //std::cout << "bc set";
    }
    

    //std::cout << "bcs set";

    //OOFEM_LOG_INFO("Added element %d (nodes: %d %d %d %d) \n", elId, enodes[0], enodes[1], enodes[2], enodes[3], enodes[4], enodes[5], enodes[6], enodes[7]);   

    /*for ( auto &ic : d->giveIcs() ) {
        emodel->giveField(FT_Temperature, tStep);

        //this->applyInitialCondition(*ic);
    }*/

    //d->giveNumberOfRegions();
    //d->elementList [ elId - 1 ] = std::move(elem);

    //OOFEM_ERROR("elementList error %d", d->elementList.size());

    // MUST FORCE RENUMBER!
    //if ( this->requiresEquationRenumbering( this->giveCurrentStep() ) ) {
    
}

// Add single HT element based on existing dofs
void
addElement(EngngModel *emodel, IntArray &existingNodes, std::vector<FloatArray> &newNodes, double temperature) {
        Domain *d = emodel->giveDomain(1);
        OOFEM_LOG_INFO("Going to add Brick1HT element\n");

        int nNodesBefore = d->giveNumberOfDofManagers();
        int nNodesToAdd = newNodes.size();

        // Resize node list to accomodate new nodes
        d->resizeDofManagers(nNodesBefore+nNodesToAdd);

        // Add new nodes
        int nNodesAdded = 0;
        IntArray nNodeIds(nNodesToAdd);
        for(auto &n : newNodes) {
            int nodeId = nNodesBefore + nNodesAdded + 1;
            std::unique_ptr< DofManager >dman( classFactory.createDofManager("node", nodeId, d ) );
            if ( !dman ) {
                OOFEM_ERROR("Couldn't create node %d\n", nodeId);
            }

            dman->setCoordinates(n);
            dman->setGlobalNumber(nodeId);
            dman->appendDof( new MasterDof(dman.get(), T_f) );

            OOFEM_LOG_INFO("Added node %d\n", nodeId);

            //dman->giveDofWithID(T_f)->giveUnknowns()->add(0, 100.);
            //dman->giveDofWithID(T_f)->giveUnknowns()->add(1, 100.);
            d->setDofManager(nodeId, std::move(dman));

            //d->dofManagerList [ nodeId - 1 ] = std::move(dman);
            nNodeIds[nNodesAdded] = nodeId;
            nNodesAdded++;
        }

        OOFEM_LOG_INFO( "All new nodes added successfully\n");
    
        // Add new node set
        // Used to apply initial boundary condition to the new DOFs
        int nSetsBefore = d->giveNumberOfSets();
        int setId = nSetsBefore+1;
        d->resizeSets(setId);
        std::unique_ptr<Set> set = std::make_unique<Set>(setId, d);
        set->setNodeList(nNodeIds);
        d->setSet(setId, std::move(set));
        OOFEM_LOG_INFO( "Set %d created\n", setId);

        // Add initial condition for the new nodes
        /*int nICBefore = d->giveNumberOfInitialConditions();
        int icId = nICBefore+1;
        d->resizeInitialConditions(icId);
        std::unique_ptr< InitialCondition >ic( classFactory.createInitialCondition("initialcondition", icId, d) );

        if ( !ic ) {
            OOFEM_ERROR("Couldn't create IC\n");
        }

        DynamicInputRecord myInput(_IFT_InitialCondition_Name, icId);
        myInput.setField(IntArray{14}, _IFT_InitialCondition_dofs);

        Dictionary conditions = Dictionary();
        conditions.add('u', temperature);

        myInput.setField(conditions, _IFT_InitialCondition_conditions);
        myInput.setField(setId, _IFT_InitialCondition_set);
        ic->initializeFrom(myInput);
        d->setInitialCondition(icId, std::move(ic));

        OOFEM_LOG_INFO("Added IC %d\n", icId);    */

        // Add new element
        int nElBefore = d->giveNumberOfElements();
        int elId = nElBefore+1;
        
        std::unique_ptr< Element >elem(classFactory.createElement("brick1ht", elId, d) );
        if ( !elem ) {
            OOFEM_ERROR("Couldn't create element %d\n", elId);
        } 

        /*IntArray enodes = {
            3, 4, 6, 5
        };*/
        IntArray enodes(8);

        // Add existing node IDs
        for(int i = 0; i < 8-nNodesToAdd; i++) {
            enodes[i] = existingNodes[i];
        }

        // Add added node ids
        int j = 1;
        for(int i = 8-nNodesToAdd; i < 8; i++, j++) {
            enodes[i] = nNodesBefore+j;
        }
        
        elem->setDofManagers(enodes);
        elem->setGlobalNumber(elId);
        elem->setCrossSection(1);
        elem->postInitialize(); // this is needed to allocate the integration point array
        
        d->resizeElements(nElBefore + 1);
        d->setElement(elId, std::move(elem));

        OOFEM_LOG_INFO("Added element %d (nodes: %d %d %d %d) \n", elId, enodes[0], enodes[1], enodes[2], enodes[3], enodes[4], enodes[5], enodes[6], enodes[7]);   


        /*for ( auto &ic : d->giveIcs() ) {
            emodel->giveField(FT_Temperature, tStep);

            //this->applyInitialCondition(*ic);
        }*/

        //d->giveNumberOfRegions();
        //d->elementList [ elId - 1 ] = std::move(elem);

        //OOFEM_ERROR("elementList error %d", d->elementList.size());

        // MUST FORCE RENUMBER!
        //if ( this->requiresEquationRenumbering( this->giveCurrentStep() ) ) {
            emodel->forceEquationNumbering();
        //}

        //TransientTransportProblem *emodelns = reinterpret_cast <TransientTransportProblem *>(emodel);
        //emodelns->fiel;
        //OOFEM_LOG_INFO("I added everything I should\n");
}

AdditiveManufacturingProblem :: AdditiveManufacturingProblem(int i, EngngModel *_master) : EngngModel(i, _master),
    adaptiveStepLength(false),
    minStepLength(0.),
    maxStepLength(0.),
    reqIterations(0.),
    adaptiveStepSince(0.),
    endOfTimeOfInterest(0.),
    prevStepLength(0.),
    currentStepLength(0.)
{
    ndomains = 1; // domain is needed to store the time step function

    dtFunction = 0;
    stepMultiplier = 1.;
    timeDefinedByProb = 0;

    //fclose(stdout);
    this->totalTimer.startTimer();
}

AdditiveManufacturingProblem :: ~AdditiveManufacturingProblem()
{
    this->totalTimer.stopTimer();
    OOFEM_LOG_INFO( "TOTAL TIME SOLVING: %.2fs\n", this->totalTimer.getUtime() );
}

///////////
int
AdditiveManufacturingProblem :: instanciateYourself(DataReader &dr, InputRecord &ir, const char *dataOutputFileName, const char *desc)
{
    int result;
    result = EngngModel :: instanciateYourself(dr, ir, dataOutputFileName, desc);
    ir.finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();
    return result;
}

int
AdditiveManufacturingProblem :: instanciateDefaultMetaStep(InputRecord &ir)
{
    if ( timeDefinedByProb ) {
        /* just set a nonzero number of steps;
         * needed for instanciateDefaultMetaStep to pass; overall has no effect as time stepping is deteremined by slave
         */
        this->numberOfSteps = 1;
    }
    EngngModel :: instanciateDefaultMetaStep(ir);
    //there are no slave problems initiated so far, the overall metaStep will defined in a slave problem instantiation
    return 1;
}

int
AdditiveManufacturingProblem :: instanciateSlaveProblems()
{
    //first instantiate master problem if defined
    //EngngModel *timeDefProb = NULL;
    emodelList.resize( inputStreamNames.size() );
    if ( timeDefinedByProb ) {
        OOFEMTXTDataReader dr(inputStreamNames [ timeDefinedByProb - 1 ]);
        std :: unique_ptr< EngngModel >prob( InstanciateProblem(dr, this->pMode, this->contextOutputMode, this) );
        //timeDefProb = prob.get();
        emodelList [ timeDefinedByProb - 1 ] = std :: move(prob);
    }

    for ( int i = 1; i <= ( int ) inputStreamNames.size(); i++ ) {
        if ( i == timeDefinedByProb ) {
            continue;
        }

        OOFEMTXTDataReader dr(inputStreamNames [ i - 1 ]);
        //the slave problem dictating time needs to have attribute master=NULL, other problems point to the dictating slave
        std :: unique_ptr< EngngModel >prob( InstanciateProblem(dr, this->pMode, this->contextOutputMode, this) );
        emodelList [ i - 1 ] = std :: move(prob);
    }

    return 1;
}


void
AdditiveManufacturingProblem :: initializeFrom(InputRecord &ir)
{
    OOFEM_LOG_INFO("Starting Additive Manufacturing solver\n");

    IR_GIVE_FIELD(ir, this->gCodeFilePath, _IFT_AdditiveManufacturingProblem_gcode);
    IR_GIVE_FIELD(ir, this->stepX, _IFT_AdditiveManufacturingProblem_stepx);
    IR_GIVE_FIELD(ir, this->stepY, _IFT_AdditiveManufacturingProblem_stepy);
    IR_GIVE_FIELD(ir, this->stepZ, _IFT_AdditiveManufacturingProblem_stepz);

    this->printer = Printer(this->stepX, this->stepY, this->stepZ);
    
    OOFEM_LOG_INFO( "\n\nG-code file: %s\n", this->gCodeFilePath.c_str() );

    
    this->printer.parse(this->gCodeFilePath);
    this->printer.begin_print();

    OOFEM_LOG_INFO("\nFinished G-code parsing\n\n");

    IR_GIVE_FIELD(ir, numberOfSteps, _IFT_EngngModel_nsteps);
    if ( numberOfSteps <= 0 ) {
        throw ValueInputException(ir, _IFT_EngngModel_nsteps, "nsteps must be > 0");
    }
    if ( ir.hasField(_IFT_AdditiveManufacturingProblem_deltat) ) {
        EngngModel :: initializeFrom(ir);
        IR_GIVE_FIELD(ir, deltaT, _IFT_AdditiveManufacturingProblem_deltat);
        dtFunction = 0;
    } else if ( ir.hasField(_IFT_AdditiveManufacturingProblem_prescribedtimes) ) {
        EngngModel :: initializeFrom(ir);
        IR_GIVE_FIELD(ir, discreteTimes, _IFT_AdditiveManufacturingProblem_prescribedtimes);
        dtFunction = 0;
    } else if ( ir.hasField(_IFT_AdditiveManufacturingProblem_dtf) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, dtFunction, _IFT_AdditiveManufacturingProblem_dtf);
    } else {
        IR_GIVE_FIELD(ir, timeDefinedByProb, _IFT_AdditiveManufacturingProblem_timeDefinedByProb);
    }

    if ( ir.hasField(_IFT_AdditiveManufacturingProblem_adaptiveStepLength) ) {
        adaptiveStepLength = true;
        this->minStepLength = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_AdditiveManufacturingProblem_minsteplength);
        this->maxStepLength = 1.e32;
        IR_GIVE_OPTIONAL_FIELD(ir, maxStepLength, _IFT_AdditiveManufacturingProblem_maxsteplength);
        this->reqIterations = 1;
        IR_GIVE_OPTIONAL_FIELD(ir, reqIterations, _IFT_AdditiveManufacturingProblem_reqiterations);
        this->endOfTimeOfInterest = 1.e32;
        IR_GIVE_OPTIONAL_FIELD(ir, endOfTimeOfInterest, _IFT_AdditiveManufacturingProblem_endoftimeofinterest);
        this->adaptiveStepSince = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, adaptiveStepSince, _IFT_AdditiveManufacturingProblem_adaptivestepsince);
    }


    IR_GIVE_OPTIONAL_FIELD(ir, stepMultiplier, _IFT_AdditiveManufacturingProblem_stepmultiplier);
    if ( stepMultiplier < 0 ) {
        throw ValueInputException(ir, _IFT_AdditiveManufacturingProblem_stepmultiplier, "stepMultiplier must be > 0");
    }

    //    timeLag = 0.;
    //    IR_GIVE_OPTIONAL_FIELD(ir, timeLag, _IFT_AdditiveManufacturingProblem_timeLag);

    inputStreamNames.resize(2);
    if ( ir.hasField(_IFT_AdditiveManufacturingProblem_prob3) ){
        inputStreamNames.resize(3);
    }
    
    IR_GIVE_FIELD(ir, inputStreamNames [ 0 ], _IFT_AdditiveManufacturingProblem_prob1);
    IR_GIVE_FIELD(ir, inputStreamNames [ 1 ], _IFT_AdditiveManufacturingProblem_prob2);
    IR_GIVE_OPTIONAL_FIELD(ir, inputStreamNames [ 2 ], _IFT_AdditiveManufacturingProblem_prob3);
    
    
    renumberFlag = true; // The staggered problem itself should always try to check if the sub-problems needs renumbering.

    coupledModels.resize(3);
    IR_GIVE_OPTIONAL_FIELD(ir, this->coupledModels, _IFT_AdditiveManufacturingProblem_coupling);


    if ( dtFunction < 1 ) {
        ndomains = 0;
        domainNeqs.clear();
        domainPrescribedNeqs.clear();
        domainList.clear();
    }

    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if ( suppressOutput ) {
        printf("Suppressing output.\n");
    } else {

        if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
            throw ValueInputException(ir, "None", "can't open output file: " + this->dataOutputFileName);
        }

        fprintf(outputStream, "%s", PRG_HEADER);
        fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
        fprintf(outputStream, "%s\n", simulationDescription.c_str());

#ifdef __PARALLEL_MODE
        if ( this->isParallel() ) {
            fprintf(outputStream, "Problem rank is %d/%d on %s\n\n", this->rank, this->numProcs, this->processor_name);
        }
#endif
    }
}


void
AdditiveManufacturingProblem :: updateAttributes(MetaStep *mStep)
{
    auto &ir = mStep->giveAttributesRecord();

    EngngModel :: updateAttributes(mStep);

    // update attributes of slaves
    for ( auto &emodel: emodelList ) {
        emodel->updateAttributes(mStep);
    }

    if ( !timeDefinedByProb ) {
        if ( ir.hasField(_IFT_AdditiveManufacturingProblem_deltat) ) {
            IR_GIVE_FIELD(ir, deltaT, _IFT_AdditiveManufacturingProblem_deltat);
            IR_GIVE_OPTIONAL_FIELD(ir, dtFunction, _IFT_AdditiveManufacturingProblem_dtf);
            IR_GIVE_OPTIONAL_FIELD(ir, stepMultiplier, _IFT_AdditiveManufacturingProblem_stepmultiplier);
            if ( stepMultiplier < 0 ) {
                OOFEM_ERROR("stepMultiplier must be > 0")
            }
        } else if ( ir.hasField(_IFT_AdditiveManufacturingProblem_prescribedtimes) ) {
            IR_GIVE_FIELD(ir, discreteTimes, _IFT_AdditiveManufacturingProblem_prescribedtimes);
        }
    }
}

Function *AdditiveManufacturingProblem :: giveDtFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtFunction ) {
        return NULL;
    }

    return giveDomain(1)->giveFunction(dtFunction);
}

double
AdditiveManufacturingProblem :: giveDeltaT(int n)
{
    if ( giveDtFunction() ) {
        return giveDtFunction()->evaluateAtTime(n);
    }

    //in the first step the time increment is taken as the initial, user-specified value
    if ( stepMultiplier != 1 && currentStep ) {
        if ( currentStep->giveNumber() >= 2 ) {
            return ( currentStep->giveTargetTime() * ( stepMultiplier ) );
        }
    }

    if ( discreteTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    }

    if ( adaptiveStepLength ) {
        EngngModel *sp;
        int nite = 1;
        double adjustedDeltaT = deltaT;

        if ( currentStep != NULL ) {
            if ( currentStep->giveNumber() != 0 ) {
                // return prescribed deltaT for times until time = adaptiveStepSince
                // can be used for consecutive force loading applied in a specified number of steps
                if ( !( currentStep->giveTargetTime() > this->adaptiveStepSince ) ) {
                    return adjustedDeltaT;
                }

                for ( int i = 1; i <= this->giveNumberOfSlaveProblems(); i++ ) {
                    sp = this->giveSlaveProblem(i);
                    nite = max(sp->giveCurrentNumberOfIterations(), nite);
                }

                if ( nite > reqIterations ) {
                    adjustedDeltaT =  this->prevStepLength * reqIterations / nite;
                } else {
                    adjustedDeltaT  =  this->prevStepLength * sqrt( sqrt( ( double ) reqIterations / ( double ) nite ) );
                }

                if ( adjustedDeltaT > maxStepLength ) {
                    adjustedDeltaT = maxStepLength;
                }

                if ( adjustedDeltaT < minStepLength ) {
                    adjustedDeltaT = minStepLength;
                }
            }
        }

        this->currentStepLength = adjustedDeltaT;

        return adjustedDeltaT;
    }

    return deltaT;
}

double
AdditiveManufacturingProblem :: giveDiscreteTime(int iStep)
{
    if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( discreteTimes.at(iStep) );
    }

    if ( ( iStep == 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( 0.0 );
    }

    OOFEM_ERROR("invalid iStep");
    return 0.0;
}

TimeStep *
AdditiveManufacturingProblem :: giveCurrentStep(bool force)
{
    if ( timeDefinedByProb ) {
        return emodelList [ timeDefinedByProb - 1 ].get()->giveCurrentStep(true);
    } else {
        return EngngModel :: giveCurrentStep();
    }
}

TimeStep *
AdditiveManufacturingProblem :: givePreviousStep(bool force)
{
    if ( timeDefinedByProb ) {
        return emodelList [ timeDefinedByProb - 1 ].get()->givePreviousStep(true);
    } else {
        return EngngModel :: givePreviousStep();
    }
}

TimeStep *
AdditiveManufacturingProblem :: giveSolutionStepWhenIcApply(bool force)
{
    if ( timeDefinedByProb ) {
        return emodelList [ timeDefinedByProb - 1 ].get()->giveSolutionStepWhenIcApply(true);
    } else {
        if ( !stepWhenIcApply ) {
            int nFirst = giveNumberOfFirstStep();
            //stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, -giveDeltaT(nFirst), giveDeltaT(nFirst), 0); //previous version for [-dt, 0]
            stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, 0., giveDeltaT(nFirst), 0); //now go from [0, dt]
        }

        return stepWhenIcApply.get();
    }
}


EngngModel *
AdditiveManufacturingProblem :: giveTimeControl(){
    if ( !timeDefinedByProb ) {
        return this;
    } else { //time dictated by slave problem
        return this->giveSlaveProblem(timeDefinedByProb);
    }
}


int
AdditiveManufacturingProblem :: giveNumberOfFirstStep(bool force)
{
    if ( timeDefinedByProb && !force) {
        return emodelList [ timeDefinedByProb - 1 ].get()->giveNumberOfFirstStep(true);
    } else {
        return EngngModel :: giveNumberOfFirstStep(force);
    }
}


TimeStep *
AdditiveManufacturingProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    double totalTime = 0;
    StateCounterType counter = 1;

    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    double dt = this->giveDeltaT(currentStep->giveNumber()+1);
    istep =  currentStep->giveNumber() + 1;
    totalTime = currentStep->giveTargetTime() + this->giveDeltaT(istep);
    counter = currentStep->giveSolutionStateCounter() + 1;
    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(*previousStep, dt);

    if ( ( totalTime >= this->endOfTimeOfInterest ) && this->adaptiveStepLength ) {
        totalTime = this->endOfTimeOfInterest;
        OOFEM_LOG_INFO("\n==================================================================\n");
        OOFEM_LOG_INFO( "\nAdjusting time step length to: %lf \n\n", totalTime - previousStep->giveTargetTime() );
        currentStep = std::make_unique<TimeStep>(istep, this, 1, totalTime, totalTime - previousStep->giveTargetTime(), counter);
        this->numberOfSteps = istep;
    } else {
        if ( this->adaptiveStepLength ) {
            OOFEM_LOG_INFO("\n==================================================================\n");
            OOFEM_LOG_INFO( "\nAdjusting time step length to: %lf \n\n", totalTime - previousStep->giveTargetTime() );
        }
        currentStep = std::make_unique<TimeStep>(istep, this, 1, totalTime, totalTime - previousStep->giveTargetTime(), counter);
    }

    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep.get();
}


void
AdditiveManufacturingProblem :: solveYourself()
{
    EngngModel *sp;
    sp = giveTimeControl();

    int smstep = 1, sjstep = 1;
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    if ( sp->giveCurrentStep() ) {
        smstep = sp->giveCurrentStep()->giveMetaStepNumber();
        sjstep = sp->giveMetaStep(smstep)->giveStepRelativeNumber( sp->giveCurrentStep()->giveNumber() ) + 1;
    }

    for ( int imstep = smstep; imstep <= sp->giveNumberOfMetaSteps(); imstep++ ) { //loop over meta steps
        MetaStep *activeMStep = sp->giveMetaStep(imstep);
        // update state according to new meta step in all slaves
        this->initMetaStepAttributes(activeMStep);

        int nTimeSteps = activeMStep->giveNumberOfSteps();
        for ( int jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { //loop over time steps
            this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
            this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
            sp->preInitializeNextStep();
            sp->giveNextStep();

            // renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
            if ( this->requiresEquationRenumbering( sp->giveCurrentStep() ) ) {
                this->forceEquationNumbering();
            }

            Timer timer;
            timer.startTimer();
            this->initializeYourself( sp->giveCurrentStep() );
            timer.stopTimer();
            printf("this->initializeYourself( sp->giveCurrentStep() ); took %.3f s\n", timer.getUtime());
            
            timer.startTimer();
            this->solveYourselfAt( sp->giveCurrentStep() );
            timer.stopTimer();
            printf("this->solveYourselfAt( sp->giveCurrentStep() ); took %.3f s\n", timer.getUtime());
            
            timer.startTimer();
            this->updateYourself( sp->giveCurrentStep() );
            timer.stopTimer();
            printf("this->updateYourself( sp->giveCurrentStep() ); took %.3f s\n", timer.getUtime());
            
            timer.startTimer();
            this->terminate( sp->giveCurrentStep() );
            timer.stopTimer();
            printf("this->terminate( sp->giveCurrentStep() ); took %.3f s\n", timer.getUtime());

            this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
            double _steptime = this->timer.getUtime(EngngModelTimer :: EMTT_SolutionStepTimer);
            OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n",
                           sp->giveCurrentStep()->giveNumber(), _steptime);

            if(!suppressOutput) {
            	fprintf(this->giveOutputStream(), "\nUser time consumed by solution step %d: %.3f [s]\n\n",
                        sp->giveCurrentStep()->giveNumber(), _steptime);
            }

#ifdef __PARALLEL_MODE
            if ( loadBalancingFlag ) {
                this->balanceLoad( sp->giveCurrentStep() );
            }

#endif

            if ( ( sp->giveCurrentStep()->giveTargetTime() >= this->endOfTimeOfInterest ) && this->adaptiveStepLength ) {
                break;
            }
        }
    }
}

void
AdditiveManufacturingProblem :: solveYourselfAt(TimeStep *tStep)
{
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT("Solving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime());
#endif
    
    EngngModel *emodel = this->emodelList.at(0).get();

    if(!tStep->isTheFirstStep()) {
        Timer timer;
        timer.startTimer();
        this->printer.print_to_time(tStep->giveTargetTime());
        timer.stopTimer();
        OOFEM_LOG_INFO( "Printer::print_to_time: %.4fs\n", timer.getUtime() );

        bool adding = this->printer.getGrid().getCellNodes().size() > emodel->giveDomain(1)->giveNumberOfDofManagers();

        std::cout << "step total nodes: " << this->printer.getGrid().getCellNodes().size() << "\n";
        std::cout << "step total voxels: " << this->printer.getGrid().getCells().size() << "\n";

        if(adding) {
            OOFEM_LOG_INFO("gonna add nodes\n");

            timer.startTimer();
            
            const int nToAdd = this->printer.getGrid().getCellNodes().size() - emodel->giveDomain(1)->giveNumberOfDofManagers();
            emodel->giveDomain(1)->resizeDofManagers(emodel->giveDomain(1)->giveNumberOfDofManagers() + nToAdd);

            for(auto &n : this->printer.getGrid().getCellNodes()) {
                bool added = add_node_if_not_exists(emodel, n.second); 

                if(added) {
                    emodel->giveDomain(1)->giveDofManager(n.second.id)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep, VM_Total, 235.0);
                    emodel->giveDomain(1)->giveDofManager(n.second.id)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep, VM_TotalIntrinsic, 235.0);
                    emodel->giveDomain(1)->giveDofManager(n.second.id)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 235.0);
                    emodel->giveDomain(1)->giveDofManager(n.second.id)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep->givePreviousStep(), VM_TotalIntrinsic, 235.0);
                }           
            }
            timer.stopTimer();
            OOFEM_LOG_INFO( "Nodes added: %.4fs\n", timer.getUtime() );

            OOFEM_LOG_INFO("gonna add els\n");

            timer.startTimer();

            const int eToAdd = this->printer.getGrid().getCells().size() - emodel->giveDomain(1)->giveNumberOfElements();
            emodel->giveDomain(1)->resizeElements(emodel->giveDomain(1)->giveNumberOfElements() + eToAdd);

            for(auto &n : this->printer.getGrid().getCells()) {
                add_element_if_not_exists(emodel, n.second);
            }
            timer.stopTimer();
            OOFEM_LOG_INFO( "Elements added: %.4fs\n", timer.getUtime() );

            std::cout << "total dofmans " << emodel->giveDomain(1)->giveNumberOfDofManagers() << "\n";
            std::cout << "total elems " << emodel->giveDomain(1)->giveNumberOfElements() << "\n";

            OOFEM_LOG_INFO("added\n");
            if(adding) {
                timer.startTimer();
                emodel->forceEquationNumbering();
                timer.stopTimer();
                printf("emodel->forceEquationNumbering(); took %.3f s\n", timer.getUtime());
                
                // todo: is needed? yes it is!
                timer.startTimer();
                emodel->giveDomain(1)->giveConnectivityTable()->reset();
                timer.stopTimer();
                printf("emodel->giveDomain(1)->giveConnectivityTable()->reset(); took %.3f s\n", timer.getUtime());

                //emodelList.at(0)->giveExportModuleManager()->giveModule(1)->clear
                if(emodelList.at(0)->giveExportModuleManager()->giveNumberOfModules() > 0) {
                    emodelList.at(0)->giveExportModuleManager()->giveModule(1)->initialize();
                }
            }
        }
    }

    std::cout << "LETS SOLVE\n";

    /*if(emodel->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) == 0 && !tStep->isTheFirstStep()) {
        OOFEM_LOG_INFO("gonna add nodes\n");
        IntArray existingNodes = {};//{3,4};
        // std::vector<FloatArray> newNodes{ FloatArray{0.04, 0.24, 0.00}, FloatArray{0.00, 0.24, 0.00} };
        std::vector<FloatArray> newNodes{
            FloatArray{0.0/10, 0.0/10, 1.0/10},
            FloatArray{0.0/10, 1.0/10, 1.0/10},
            FloatArray{1.0/10, 1.0/10, 1.0/10},
            FloatArray{1.0/10, 0.0/10, 1.0/10},
            FloatArray{0.0/10, 0.0/10, 0.0/10},
            FloatArray{0.0/10, 1.0/10, 0.0/10},
            FloatArray{1.0/10, 1.0/10, 0.0/10},
            FloatArray{1.0/10, 0.0/10, 0.0/10},
        };

        addElement(emodel, existingNodes, newNodes, 1000.0);

        emodelList.at(0)->giveExportModuleManager()->giveModule(1)->initialize();

        //this->forceEquationNumbering();
        //OOFEM_LOG_INFO("added stuff\n");
    
        for(int i = 1; i <= 8; i++) {
            emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 235.0 + i);
            //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep, VM_Total, 235.0*i);
            //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep, VM_Total, 235.0*i);
        }
    }*/

    // TODO: THIS IS NOT WORKING! HOW TO SET VALUES DIRECTLY?
    /*for(int i = 1; i <= 8; i++) {
        //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->set`
        emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 235.0*i);
        //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->updateUnknownsDictionary(tStep, VM_Total, 235.0*i);   
        //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->unkn
        //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->giveUnknowns()->add(VM_Total, 235.);
        //emodel->giveDomain(1)->giveDofManager(i)->giveDofWithID(T_f)->giveUnknowns()->add(0, 235.);
    }*/

    /*if(tStep->giveTargetTime() >= 7200) {
        OOFEM_LOG_INFO("=====> init value %f\n", emodel->giveDomain(1)->giveDofManager(6)->giveDofWithID(T_f)->giveUnknowns()->at(0));
        OOFEM_LOG_INFO("=====> init value %f\n", emodel->giveDomain(1)->giveDofManager(6)->giveDofWithID(T_f)->giveUnknowns()->at(VM_Total));
    }*/

    for ( auto &emodel: emodelList ) {
        emodel->solveYourselfAt(tStep);
    }

    tStep->incrementStateCounter();

    //std::cout << "step total nodes: " << this->printer.getGrid().getCellNodes().size() << "\n";
    //std::cout << "step total voxels: " << this->printer.getGrid().getCells().size() << "\n";
}

int
AdditiveManufacturingProblem :: forceEquationNumbering()
{
    int neqs = 0;
    for ( auto &emodel: emodelList ) {
        // renumber equations if necessary
        if ( emodel->requiresEquationRenumbering( emodel->giveCurrentStep() ) ) {
            neqs += emodel->forceEquationNumbering();
        }
    }

    return neqs;
}

void
AdditiveManufacturingProblem :: updateYourself(TimeStep *tStep)
{
    if ( adaptiveStepLength ) {
        this->prevStepLength = this->currentStepLength;
    }

    for ( auto &emodel: emodelList ) {
        emodel->updateYourself(tStep);
    }

    EngngModel :: updateYourself(tStep);
}

void
AdditiveManufacturingProblem :: terminate(TimeStep *tStep)
{
    for ( auto &emodel: emodelList ) {
        emodel->terminate(tStep);
    }
}

void
AdditiveManufacturingProblem :: doStepOutput(TimeStep *tStep)
{
    for ( auto &emodel: emodelList ) {
        emodel->giveExportModuleManager()->doOutput(tStep);
    }
}


void
AdditiveManufacturingProblem :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // Subproblems handle the output themselves.
}


void
AdditiveManufacturingProblem :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);
    for ( auto &emodel: emodelList ) {
        emodel->saveContext(stream, mode);
    }
}


void
AdditiveManufacturingProblem :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);
    for ( auto &emodel: this->emodelList ) {
        emodel->restoreContext(stream, mode);
    }
}


EngngModel *
AdditiveManufacturingProblem :: giveSlaveProblem(int i)
{
    if ( ( i > 0 ) && ( i <= this->giveNumberOfSlaveProblems() ) ) {
        return this->emodelList [ i - 1 ].get();
    } else {
        OOFEM_ERROR("Undefined problem");
    }

    return NULL;
}


int
AdditiveManufacturingProblem :: checkProblemConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int result = 1;
    for ( auto &emodel: emodelList ) {
        result &= emodel->checkProblemConsistency();
    }

#  ifdef VERBOSE
    if ( result ) {
        OOFEM_LOG_DEBUG("Consistency check:  OK\n");
    } else {
        VERBOSE_PRINTS("Consistency check", "failed")
        exit(1);
    }

#  endif

    return result;
}

void
AdditiveManufacturingProblem :: updateDomainLinks()
{
    for ( auto &emodel: emodelList ) {
        emodel->updateDomainLinks();
    }
}

void
AdditiveManufacturingProblem :: setRenumberFlag()
{
    for ( auto &emodel: emodelList ) {
        emodel->setRenumberFlag();
    }
}

#ifdef __OOFEG
void AdditiveManufacturingProblem :: drawYourself(oofegGraphicContext &gc)
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem(ap)->drawYourself(gc);
    }
}

void AdditiveManufacturingProblem :: drawElements(oofegGraphicContext &gc)
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem(ap)->drawElements(gc);
    }
}

void AdditiveManufacturingProblem :: drawNodes(oofegGraphicContext &gc)
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem(ap)->drawNodes(gc);
    }
}
#endif
} // end namespace oofem
