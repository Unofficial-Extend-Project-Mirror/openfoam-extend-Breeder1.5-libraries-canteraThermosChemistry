/*---------------------------------------------------------------------------*\
This file written by Institute of Energy Process Enineering and Chemical
	Engineering TU Freiberg  http://www.iec.tu-freiberg.de
and ICE Stroemungsfoschungs GmbH http://www.ice-sf.at
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    canteraLocalTimeChemistryModel

Description

SourceFiles
    canteraLocalTimeChemistryModel.C

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraLocalTimeChemistryModel.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef canteraLocalTimeChemistryModel_H
#define canteraLocalTimeChemistryModel_H

#include "canteraChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class canteraLocalTimeChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class canteraLocalTimeChemistryModel
:
    public canteraChemistryModel
{

        //- Disallow default bitwise assignment
        void operator=(const canteraLocalTimeChemistryModel&);

public:

    TypeName("canteraLocalTimeChemistryModel");

    // Constructors

        //- Construct from components
        canteraLocalTimeChemistryModel
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        );


    // Destructor

        virtual ~canteraLocalTimeChemistryModel();

        //- Solve the reaction system for the given start time and time-step
        // and return the characteristic time
        virtual scalar solve
        (
            const scalar t0,
            const scalar deltaT
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
