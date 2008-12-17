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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraMixture.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraMixture.H"
#include "fvMesh.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// const ThermoType& canteraMixture<ThermoType>::constructSpeciesData
// (
//     const dictionary& thermoDict
// )
// {
//     forAll(species_, i)
//     {
//         speciesData_.set
//         (
//             i,
//             new ThermoType(thermoDict.lookup(species_[i]))
//         );
//     }

//     return speciesData_[0];
// }


// void canteraMixture<ThermoType>::correctMassFractions()
// {
//     volScalarField Yt = Y_[0];

//     for(label n=1; n<Y_.size(); n++)
//     {
//         Yt += Y_[n];
//     }

//     forAll (Y_, n)
//     {
//         Y_[n] /= Yt;
//     }
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// canteraMixture<ThermoType>::canteraMixture
// (
//     const dictionary& thermoDict,
//     const wordList& specieNames,
//     const HashPtrTable<ThermoType>& specieThermoData,
//     const fvMesh& mesh
// )
// :
//     combustionMixture(thermoDict, specieNames, mesh),
//     speciesData_(species_.size()),
//     mixture_("mixture", *specieThermoData[specieNames[0]])
// {
//     forAll(species_, i)
//     {
//         speciesData_.set
//         (
//             i,
//             new ThermoType(*specieThermoData[species_[i]])
//         );
//     }

//     correctMassFractions();
// }


canteraMixture::canteraMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    combustionMixture(thermoDict, canteraThermo(thermoDict.subDict("cantera")).getNames(), mesh),
    speciesData_(thermoDict.subDict("cantera")),
    yTemp_(Y().size()),
    T_(mesh.db().lookupObject<volScalarField>("T")),
    p_(mesh.db().lookupObject<volScalarField>("p"))
    //    speciesData_(species_.size()),
    //    mixture_("mixture", constructSpeciesData(thermoDict))
{
    //    correctMassFractions();

    const Cantera::XML_Node& xml=*(speciesData_.gas().speciesData());
    //    Info << "Anzahl: " << xml.nChildren() << endl;
    //    xml.write(std::cout);

    speciesDataMirror_.setSize(xml.nChildren()/2);

    fileName thermoStandinFile(thermoDict.lookup("canteraStandinThermoFile"));
    thermoStandinFile.expand();

    WarningIn("canteraMixture::canteraMixture")
        << " The thermophysical properties of CANTERA are currently not"
            << " converted to OpenFOAM. Instead the properties from " 
            << thermoStandinFile << " are used" 
            << " This can lead to errors if the data is inconsistent with the Cantera-data" 
            << endl;

    IFstream tf(thermoStandinFile);

    HashPtrTable<speciesType> standinThermo(tf);

    //    Info << standinThermo.toc() << endl;

    for(label i=0;i<xml.nChildren()/2;i++) {
        Cantera::XML_Node& spec=xml.child(2*i+1);
        word name(spec["name"]);
        //        Info << i << " " << name << endl;
        //        const speciesType *tmp=standinThermo[name];
        speciesDataMirror_.set(
            i,
            new speciesType(name,*standinThermo[name])
        );
    }
    //    Info << "Fertig " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void canteraMixture::read(const dictionary& thermoDict)
{
    Info << " Read stuff " << endl;
    //  speciesData_ = canteraThermo(thermoDict.subDict("cantera"));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
