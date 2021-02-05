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

#include "connectivitytable.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "intarray.h"
#include "timer.h"
#include <set>

namespace oofem {

void
ConnectivityTable :: reset()
{
    nodalConnectivityFlag = 0;
}

void
ConnectivityTable :: instanciateConnectivityTable()
//
// assembles table of nodal connectivities
//
{
    int ndofMan = domain->giveNumberOfDofManagers();
    int nelems = domain->giveNumberOfElements();
    IntArray dofManConnectivity(ndofMan);

    if ( nodalConnectivityFlag ) {
        return;                     // already initialized
    }

//    OOFEM_LOG_INFO("ConnectivityTable: initializing\n");

    for ( auto &elem : domain->giveElements() ) {
        int nnodes = elem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnodes; j++ ) {
            int jnode = elem->giveDofManager(j)->giveNumber();
            dofManConnectivity.at(jnode)++;
        }
    }

    // allocate Nodal connectivity table for domain
    nodalConnectivity.resize(ndofMan);
    for ( int i = 0; i < ndofMan; i++ ) {
        nodalConnectivity[i].resize( dofManConnectivity[i] );
    }

    // build Nodal connectivity table for domain
    dofManConnectivity.zero();

    for ( int i = 1; i <= nelems; i++ ) {
        Element *ielem = domain->giveElement(i);
        int nnodes = ielem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnodes; j++ ) {
            int jnode = ielem->giveDofManager(j)->giveNumber();
            nodalConnectivity[jnode-1].at( ++dofManConnectivity.at(jnode) ) = i;
        }
    }

    nodalConnectivityFlag = 1;
}

const IntArray *
ConnectivityTable :: giveDofManConnectivityArray(int dofman)
{
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    return &this->nodalConnectivity[dofman-1];
}


void
ConnectivityTable :: giveElementNeighbourList(IntArray &answer, IntArray &elemList)
{
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    answer.resize(0);

    for ( auto &el_num : elemList ) {
        Element *ielem = domain->giveElement( el_num );
        int nnode = ielem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnode; j++ ) {
            int jnode = ielem->giveDofManager(j)->giveNumber();
            for ( auto &val : this->nodalConnectivity[jnode-1] ) {
                answer.insertSortedOnce( val );
            }
        }
    }
}


void
ConnectivityTable :: giveNodeNeighbourList(IntArray &answer, IntArray &nodeList)
{
    int nnodes = nodeList.giveSize();
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    answer.resize(0);

    for ( int i = 1; i <= nnodes; i++ ) {
        int inode = nodeList.at(i);
        for ( auto &val : this->nodalConnectivity[inode-1] ) {
            answer.insertSortedOnce( val );
        }
    }
}


void
ConnectivityTable::giveElementColoring(std::vector<std::vector<int> > &elementColorSets)
{
  /* strategy == 1 the minimum colors alloacated, the lowest color assigned
     strategy == 2 the maxcolors created; element gets assigned color which has been assigned to the lowest number of elements
  */ 
  int strategy = 2;
  int strategy2_numberOfColors = 30;
  int nelem = this->domain->giveNumberOfElements();

  //coloring of elements
  IntArray cElemArray(1); //initial parameters of all ielem =0
  IntArray elementColors(nelem);
  IntArray colorCounter;
  IntArray neighboursArray;
  std::set<int> elementNeighborColors;
  Timer timer;

  if (strategy == 2) {
    colorCounter.resize(strategy2_numberOfColors);
  }
  
  timer.startTimer();
  elementColorSets.clear();
  if (nelem == 0) return;
  for (int ielem = 1; ielem <= nelem; ielem++)
    elementColors.at(ielem) = -1; //all members of elementColors assign initial value
  
  for (int ielem = 1; ielem <= nelem; ielem++) {
    // Element *element = domain->giveElement(ielem);
    
    cElemArray.at(1) = ielem; //first member [0] in cElemArray is set as ielem
    this->giveElementNeighbourList(neighboursArray, cElemArray);
    // to neigrhboursArray wrote neighbours ID of execute ielem
    elementNeighborColors.clear();
    
    for (int n = 1; n <= neighboursArray.giveSize(); n++) { // loop over element neighbors
      if (elementColors.at(neighboursArray.at(n)) > 0) {
	elementNeighborColors.insert(elementColors.at(neighboursArray.at(n)));
      }
    }


    // assign color code to element
    if (strategy==1) {
      // default strategy -> assign lowest unused color
      int _elementColor=1;
      while (elementNeighborColors.find(_elementColor) != elementNeighborColors.end()) {
	  _elementColor ++;
      }
      elementColors.at(ielem) = _elementColor;
      colorCounter.resizeWithValues(_elementColor);
      colorCounter.at(_elementColor)++;
    } else if (strategy ==2) {
      // determine the color with lowest occurence not in neighbour list
      bool init = true;
      int colorIndex=0, numberOfColors = colorCounter.giveSize();
      for (int c=1; c<=numberOfColors; c++) {
	if (elementNeighborColors.find(c) == elementNeighborColors.end()) {
	  // color is allowable
	  if (init) {
	    init=false;
	    colorIndex = c;
	  } else {
	    if (colorCounter.at(c) < colorCounter.at(colorIndex)) {
	      colorIndex = c;
	    }
	  }
	}
      }
      elementColors.at(ielem) = colorIndex;
      colorCounter.at(colorIndex)++;
    } // end strategy == 2
    
  }
  int numberOfColors = elementColors.maximum();
  elementColorSets.resize(numberOfColors);
  for (int e = 1; e <= nelem; e++) {
    elementColorSets[elementColors.at(e)-1].push_back(e);
  }
  
  timer.stopTimer();
  OOFEM_LOG_INFO("Element coloring done in %.2f[s] using %d colors \n", timer.getWtime(), numberOfColors);
  if (false) {
    elementColors.printYourself("Element colors");
    for (int c = 1; c <= numberOfColors; c++) {
      printf ("Color %d elements:", c);
      for (int v:elementColorSets[c-1]) {
	printf (" %d", v);
      }
      printf("\n");
    }
  }
}
  



} // end namespace oofem
