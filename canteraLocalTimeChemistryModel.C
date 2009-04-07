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
    //    Info << "deltaT=  " << deltaT << endl;
    const canteraMixture &mix=static_cast<const canteraMixture&>(thermo().composition());

    scalarList RRtmp(RR_.size());
    Reactor react;
    scalar Ns = RR_.size(); //number of species
    const canteraThermo &gas1=mix.cellMixture(0);	
    scalarList mw(RR_.size());
    //delete &gas1;
    gas1.gas().getMolecularWeights(mw.begin());

    const scalarField localTime=characteristicTime()().internalField();

    //to make reactor isobaric we need a reservoir for expansion
    //set up inert reservoir
    Reservoir inert;

    //install a wall to make the reactor isobaric
    //the wall area should be 1x-10x larger than the reactor volume
    //this is a tradeoff between stability and accuracy
    //wrong choice is a popular source for cantera crashes
    Wall w;  w.setArea(1); 
    w.setExpansionRateCoeff(10); 
    w.install(react,inert);

    forAll(rho_,cellI) {	
      try {
        ReactorNet sim;
        //void setTolerances(doublereal rtol, doublereal atol) 
// 	Info << "rtol=  " << sim.rtol() << "; atol=  " << sim.atol() << endl;
 	sim.setTolerances(1e-8, 1e-9);
// 	setNumerics(sim);
// 	Info << "rtol=  " << sim.rtol() << "; atol=  " << sim.atol() << endl;
        sim.addReactor(&react);
        //sim.setInitialTime(0);


        //        Info << cellI << endl;
        scalarField c0(Ns), c1(Ns);
        const canteraThermo &gas=mix.cellMixture(cellI);  //reacting gas
    	const canteraThermo &gas0=mix.cellMixture(cellI); //nonreacting reservoir

        gas.gas().getConcentrations(c0.begin());

 	//set up reactor
        react.insert(gas.gas());

	inert.insert(gas0.gas());

        scalar deltaT=min(defaultDeltaT,localTime[cellI]);
//     	Info << "defaultDeltaT=  " << defaultDeltaT << "localTime=  "<<localTime[cellI]< "deltaT= "<< deltaT <<endl;


//         while(now < deltaT) { //use this loop togeter with step
            try {
                sim.advance(deltaT);
//                  now = sim.step(deltaT);
            } catch(const Cantera::CVodesErr& e) {
                Cantera::showErrors(std::cerr);

                FatalErrorIn("canteraLocalTimeChemistryModel::solve")
                    << " With the state: " << gas 
                        << " and the mixture " << c0 
                        << " Cantera complained in cell " << cellI << endl
                    //                   << e << endl
                        << abort(FatalError) ;
            }
//         }


        gas.gas().getConcentrations(c1.begin());
        scalarField dc = c1 - c0;
        //        Info << cellI << " " << RRtmp << endl;
        for(label i=0; i<Ns; i++)
        {
            RR_[i][cellI] = dc[i]*mw[i]/deltaT;
        }
      } catch(const Cantera::CanteraError& e) {
          Cantera::showErrors(std::cerr);

          FatalErrorIn("canteraLocalTimeChemistryModel::solve")
                  << " Cantera complained in cell " << cellI
                  << " with a Cantera::CanteraError"  << endl
                  << abort(FatalError) ;
      }
    }

    return dtMin; //das ist nur für die Schrittweitenstrg da, denk ich

}

} // namespace Foam
