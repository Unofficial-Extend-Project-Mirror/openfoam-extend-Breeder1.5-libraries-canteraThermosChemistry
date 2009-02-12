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

    forAll(rho_,cellI) {	
        ReactorNet sim;
        sim.addReactor(&react);

        //        Info << cellI << endl;
        scalarField c0(Ns), c1(Ns);
        const canteraThermo &gas=mix.cellMixture(cellI);
        gas.gas().getConcentrations(c0.begin());

        react.insert(gas.gas());
        //        sim.initialize(0);
        sim.setInitialTime(0);
        setNumerics(sim);

        scalar deltaT=min(defaultDeltaT,localTime[cellI]);
        scalar timeLeft=deltaT;

        scalar old=0;

//         if(gas[0]>SMALL) {
//             Info << "Vorher:  " << gas << endl;
//         }
        while(timeLeft>SMALL) {
            try {
                // This doesn't seem to work yet - overshots the time
                sim.integrator().setMaxStepSize(timeLeft);
                //            sim.initialize(old);
                scalar now=sim.step(timeLeft);
                scalar dt=now-old;
                //            Info << timeLeft << " " << now << " " << dt << endl;
                if(old>0) {
                    // discard the first timestep
                    dtMin=min(dt,dtMin);
                }
                old=now; 
                timeLeft=deltaT-now;
            } catch(const Cantera::CVodesErr& e) {
                Cantera::showErrors(std::cerr);

                FatalErrorIn("canteraLocalTimeChemistryModel::solve")
                    << " With the state: " << gas 
                        << " and the mixture " << c0 
                        << " Cantera complained in cell " << cellI << endl
                    //                   << e << endl
                        << abort(FatalError) ;
            }
        }
//         if(gas[0]>SMALL) {
//             Info << "Nachher: " << gas << endl;
//         }
        
        gas.gas().getConcentrations(c1.begin());
        scalarField dc = c1 - c0;
        //        Info << cellI << " " << RRtmp << endl;
        for(label i=0; i<Ns; i++)
        {
            RR_[i][cellI] = dc[i]*mw[i]/deltaT;
        }
    }

    return dtMin; //das ist nur für die Schrittweitenstrg da, denk ich

}

} // namespace Foam
