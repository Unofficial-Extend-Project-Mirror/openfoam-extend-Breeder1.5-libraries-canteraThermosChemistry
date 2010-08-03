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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraThermoI.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalar Foam::canteraThermo::H(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }
    return gas().enthalpy_mass();
}

Foam::scalar Foam::canteraThermo::TH(const scalar H,const scalar T0) const
{
//     if(T0!=gas().temperature()) {
//         gas().setTemperature(T0);
//     }

    try {
        gas().setState_HP(H,gas().pressure());
    } catch(Cantera::CanteraError) {
        scalarList comp(gas().nSpecies());
        gas().getMassFractions(comp.begin());

        Cantera::showErrors(std::cerr);

        FatalErrorIn("canteraThermo::TH")
            << " With the state: " << comp
                << " and the enthalpy " << H 
                << " and the T0 " << T0
                << " Cantera complained " << endl
                << abort(FatalError) ;
    }
    //    Info << T0 << " - " << H << " -> " << gas().temperature() << endl;

    return gas().temperature();
}

Foam::scalar Foam::canteraThermo::Cp(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }

    return gas().cp_mass();
}

Foam::scalar Foam::canteraThermo::Cv(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }

    return gas().cv_mass();
}

Foam::scalar Foam::canteraThermo::mu(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }

    return gas().transport().viscosity();
}

Foam::scalar Foam::canteraThermo::kappa(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }

    return gas().transport().thermalConductivity();
}

Foam::scalar Foam::canteraThermo::alpha(const scalar T) const
{
    if(T!=gas().temperature()) {
        gas().setTemperature(T);
    }

    return  gas().transport().thermalConductivity()/(gas().cp_mass());
}

Foam::scalar Foam::canteraThermo::rho(const scalar p,const scalar T) const
{
    if(T!=gas().temperature() || p!=gas().pressure()) {
        gas().setState_TP(T,p);
    }

    return gas().meanMolecularWeight()*gas().density();
}

Foam::scalar Foam::canteraThermo::psi(const scalar p,const scalar T) const
{
    if(T!=gas().temperature() || p!=gas().pressure()) {
        gas().setState_TP(T,p);
    }

    return gas().meanMolecularWeight()/gas()._RT();
}

void Foam::canteraThermo::setState(const scalar &T,const scalar &p,const scalarList &Y) const
{
    gas().setState_TPY(T,p,Y.begin());
}

Foam::scalar Foam::canteraThermo::operator[](label i) const
{
    return gas().massFraction(i);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
