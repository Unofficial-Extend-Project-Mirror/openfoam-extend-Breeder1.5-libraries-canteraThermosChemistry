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
    canteraGasMixWrapper

Description

SourceFiles
    canteraGasMixWrapper.C

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraGasMixWrapper.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef canteraGasMixWrapper_H
#define canteraGasMixWrapper_H

// #include "dictionary.H"

//namespace cantera {
    // using namespace std;

#include <cstdlib>

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>
//}

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class canteraGasMixWrapper Declaration
\*---------------------------------------------------------------------------*/

class canteraGasMixWrapper
:
    public IdealGasMix
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
    //        canteraGasMixWrapper(const canteraGasMixWrapper&);

        //- Disallow default bitwise assignment
    //        void operator=(const canteraGasMixWrapper&);

    autoPtr<Cantera::Transport> trans_;

public:

    // Constructors

        //- Construct null
        canteraGasMixWrapper();

    canteraGasMixWrapper(string infile, string id);

    // Destructor

        ~canteraGasMixWrapper();


    // Member Functions

        // Access

    Cantera::Transport &transport()
        { return trans_(); }

        // Check

        // Edit

        // Write


    // Member Operators


    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
