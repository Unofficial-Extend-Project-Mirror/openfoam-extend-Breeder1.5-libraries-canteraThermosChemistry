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
    canteraThermo

Description

SourceFiles
    canteraThermo.C

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraThermo.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef canteraThermo_H
#define canteraThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "canteraGasMixWrapper.H"

#include "scalarList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class canteraThermo Declaration
\*---------------------------------------------------------------------------*/

class canteraThermo
{
    // Private data

    mutable autoPtr<canteraGasMixWrapper> gas_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
    canteraThermo(const canteraThermo&);

        //- Disallow default bitwise assignment
    void operator=(const canteraThermo&);

    
public:

    // Static data members

    
    // Constructors

        //- Construct null
    //        canteraThermo();

        //- Construct from Istream
        canteraThermo(const dictionary& dict);


    // Selectors


    // Destructor

        ~canteraThermo();


    // Member Functions

    wordList getNames();

    inline void setState(const scalar &T,const scalar &p,const scalarList &Y) const;

        // Access

    canteraGasMixWrapper &gas() const { return gas_(); }

    // Transport stuff 
    inline scalar mu(const scalar T) const;

    inline scalar alpha(const scalar T) const;

    inline scalar kappa(const scalar T) const;

    // Gas relations
    inline scalar psi(const scalar p,const scalar T) const;
    
    inline scalar rho(const scalar p,const scalar T) const;
    

    // Thermo stuff
    inline scalar H(const scalar T) const;

    inline scalar TH(const scalar H,const scalar T0) const;

    inline scalar Cp(const scalar T) const;

    inline scalar Cv(const scalar T) const;

        // Check

        // Edit

        // Write


    // Member Operators

    inline scalar operator[](label i) const;

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&,const canteraThermo &);

    // IOstream Operators

};

Ostream& operator<<(Ostream&,const canteraThermo &);
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "canteraThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
