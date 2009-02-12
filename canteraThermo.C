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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraThermo.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraThermo.H"

#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::canteraThermo::canteraThermo(const dictionary& dict)
{
    fileName fName(dict.lookup("gasFile"));
    fName.expand();
    word name(dict.lookup("gasId"));
    Info << " Reading from Cantera-File " << fName << " the mixture " << name << endl;
    try {
        gas_.set(
            new canteraGasMixWrapper(
                fName.c_str(),
                name.c_str()
            )
        );
    } 
    catch(Cantera::CanteraError) {
        Cantera::showErrors(std::cerr);
        
        FatalErrorIn("canteraThermo::canteraThermo(const dictionary& dict)")
            << "Cantera has a problem: " << endl
                << Cantera::lastErrorMessage() << endl
                << abort(FatalError);
    }
}

Foam::wordList Foam::canteraThermo::getNames() 
{
    wordList result(gas().nSpecies());

    forAll(result,i) {
        result[i]=gas().speciesName(i);
    }

    return result;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::canteraThermo::~canteraThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os,const canteraThermo &t) {
    const canteraGasMixWrapper &g=t.gas();
    label nsp=g.nSpecies();

    os << "T:" << g.temperature() << " p:" << g.pressure() 
        << " ";
    
    for(int i=0;i<nsp;i++) {
        Info << " " << g.speciesName(i) << ":" << g.massFraction(i); 
    }

    return os;
}

// ************************************************************************* //
