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

#ifndef oofemtxtdatareader_h
#define oofemtxtdatareader_h

#include "datareader.h"
#include "oofemtxtinputrecord.h"

#include <fstream>

namespace oofem {
/**
 * Class representing the implementation of plain text date reader.
 * It reads a sequence of input records from data file
 * and creates the corresponding input records.
 * There is no check for record type requested, it is assumed that records are
 * written in correct order, which determined by the coded sequence of
 * component initialization and described in input manual.
 */
class OOFEM_EXPORT OOFEMTXTDataReader : public DataReader
{
protected:
    OOFEMTXTInputRecord ir;
    std :: ifstream inputStream;
    std :: string dataSourceName;

    /// Keep track of read line from stream. Used for error repors.
    int lineNumber;

public:
    /// Constructor.
    OOFEMTXTDataReader(const char *inputfilename);
    OOFEMTXTDataReader(const OOFEMTXTDataReader & x);
    virtual ~OOFEMTXTDataReader();

    /// Return a line number, which is helpful for tracking errors.
    int giveLineNumber() { return lineNumber; }

    virtual InputRecord *giveInputRecord(InputRecordType, int recordId);
    virtual std :: string giveLine();
    virtual void finish();
    virtual const char *giveDataSourceName() const { return dataSourceName.c_str(); }

protected:
    void giveLineFromInput(std :: string &line);
    void giveRawLineFromInput(std :: string &line);
};
} // end namespace oofem
#endif // oofemtxtdatareader_h
