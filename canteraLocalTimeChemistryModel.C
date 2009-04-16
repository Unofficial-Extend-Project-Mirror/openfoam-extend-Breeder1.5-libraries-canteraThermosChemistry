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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraLocalTimeChemistryModel.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraLocalTimeChemistryModel.H"

#include "zeroGradientFvPatchFields.H"

#include "canteraMixture.H"

#include "hMixtureThermo.H"

#include <cantera/Cantera.h>
#include <cantera/kernel/CVodesIntegrator.h>
#include <cantera/zerodim.h>

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(canteraLocalTimeChemistryModel, 0);
addToRunTimeSelectionTable(alternateChemistryModel,canteraLocalTimeChemistryModel,steadyChemistry);

// Construct from components
canteraLocalTimeChemistryModel::canteraLocalTimeChemistryModel
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    canteraChemistryModel(thermo,rho)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

canteraLocalTimeChemistryModel::~canteraLocalTimeChemistryModel()
{}



// ************************************************************************* //

scalar canteraLocalTimeChemistryModel::solve(
    const scalar t0,
    const scalar defaultDeltaT
)
{
    scalar dtMin=defaultDeltaT;
    scalar Ns = RR_.size(); //number of species

    const canteraMixture &mix=static_cast<const canteraMixture&>(thermo().composition());

    Reactor react;
    ConstPressureReactor pReact;

    const scalarField localTime=characteristicTime()().internalField();

    forAll(rho_,cellI) {	
      try {

	scalarField y0(Ns), y1(Ns);
        const canteraThermo &gas=mix.cellMixture(cellI);  //reacting gas

	gas.gas().getMassFractions(y0.begin());
	scalar rhomean = gas.gas().density();
        //set up reactor
        Reactor react; //Isochoric reactor
        ConstPressureReactor pReact; //putting this outside the cellI loop speeds up calc but makes it unstable ?

    	ReactorNet sim;

        if(useOldImplementation_) {
            sim.addReactor(&react);
            react.insert(gas.gas());
        } else {
            sim.addReactor(&pReact);
            pReact.insert(gas.gas());
        }
 
        setNumerics(sim);

        scalar deltaT=min(defaultDeltaT,localTime[cellI]);


// 	scalar now=0;
//         while(now < deltaT) { //use this loop togeter with step
            try {
                sim.advance(deltaT);
//                  now = sim.step(deltaT);
            } catch(const Cantera::CVodesErr& e) {
                Cantera::showErrors(std::cerr);

                FatalErrorIn("canteraLocalTimeChemistryModel::solve")
                    << " With the state: " << gas 
                        << " and the mixture " << y0 
                        << " Cantera complained in cell " << cellI << endl
                    //                   << e << endl
                        << abort(FatalError) ;
            }
//         }

	gas.gas().getMassFractions(y1.begin());

        for(label i=0; i<Ns; i++)
        {
	    RR_[i][cellI] = (y1[i] - y0[i])*rhomean/deltaT; 
        }
      } catch(const Cantera::CanteraError& e) {
          Cantera::showErrors(std::cerr);

          FatalErrorIn("canteraLocalTimeChemistryModel::solve")
                  << " Cantera complained in cell " << cellI
                  << " with a Cantera::CanteraError"  << endl
                  << abort(FatalError) ;
      }
    }

    return dtMin; 

}

} // namespace Foam
