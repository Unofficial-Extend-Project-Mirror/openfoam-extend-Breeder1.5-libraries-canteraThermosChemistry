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
    canteraChemistryModel

Description

SourceFiles
    canteraChemistryModel.C

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraChemistryModel.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef canteraChemistryModel_H
#define canteraChemistryModel_H

#include "alternateChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CanteraZeroD {
    class ReactorNet;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class canteraChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class canteraChemistryModel
:
    public alternateChemistryModel
{

        //- Disallow default bitwise assignment
        void operator=(const canteraChemistryModel&);

protected:

        //- the chemical source term
        PtrList<scalarField> RR_;

        //- Relative tolerance to control CVode
        scalar relTol_;
        //- Absolute tolerance to control CVode
        scalar absTol_;
        //- Relative tolerance sensitivity to control CVode
        scalar relTolSens_;
        //- Absolute tolerance sensitivity to control CVode
        scalar absTolSens_;
        //- Maximum number of steps to control CVode
        label maxSteps_;
        //- Minimum number stepsize to control CVode
        scalar minH_;

    void setNumerics(CanteraZeroD::ReactorNet &r);
public:

    TypeName("canteraChemistryModel");

    // Constructors

        //- Construct from components
        canteraChemistryModel
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        );


    // Destructor

        virtual ~canteraChemistryModel();


        //- The chemical time scale
        virtual tmp<volScalarField> tc() const;

        //- Return the chemical source terms
        virtual tmp<volScalarField> RR(const label i) const;

        //- Solve the reaction system for the given start time and time-step
        // and return the characteristic time
        virtual scalar solve
        (
            const scalar t0,
            const scalar deltaT
        );

        //- Calculate dQ
        virtual void calcDQ(volScalarField &);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
