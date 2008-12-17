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
    canteraMixture

Description

SourceFiles
    canteraMixture.C

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraMixture.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef canteraMixture_H
#define canteraMixture_H

#include "combustionMixture.H"
#include "HashPtrTable.H"

#include "canteraThermo.H"

#include "sutherlandTransport.H"
#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class canteraMixture Declaration
\*---------------------------------------------------------------------------*/

class canteraMixture
:
    public combustionMixture
{
    // Private data

        //- Species data
    canteraThermo speciesData_;

    mutable scalarList yTemp_;

    const volScalarField &T_;

    const volScalarField &p_;

        //- Temporary storage for the cell/face mixture thermo data
    //        mutable ThermoType mixture_;


    // Private member functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
    //        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
    //        void correctMassFractions();

        //- Construct as copy (not implemented)
        canteraMixture(const canteraMixture&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef canteraThermo thermoType;

    typedef sutherlandTransport<specieThermo<janafThermo<perfectGas> > > speciesType;

private:
    
    //- A mirror copy of the species data to make certain OF-codes happy
    PtrList<speciesType> speciesDataMirror_;

public:
    // Constructors

        //- Construct from dictionary, specie names, thermo database and mesh
//         canteraMixture
//         (
//             const dictionary&,
//             const wordList& specieNames,
//             const HashPtrTable<ThermoType>& specieThermoData,
//             const fvMesh&
//         );

        //- Construct from dictionary and mesh
        canteraMixture(const dictionary&, const fvMesh&);


    // Member functions

    inline const canteraThermo& cellMixture(const label celli) const;

    inline const canteraThermo& patchFaceMixture
    (
        const label patchi,
        const label facei
    ) const;

    //- Return the raw specie thermodynamic data
    const PtrList<speciesType>& speciesData() const
        {
            return speciesDataMirror_;
        }

        //- Read dictionary
    void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "canteraMixtureI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
