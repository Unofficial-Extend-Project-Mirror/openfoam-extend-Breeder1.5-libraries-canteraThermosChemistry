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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraChemistryModel.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraChemistryModel.H"

#include "zeroGradientFvPatchFields.H"

#include "canteraMixture.H"

#include "hMixtureThermo.H"

#include <cantera/Cantera.h>
#include <cantera/kernel/CVodesIntegrator.h>
#include <cantera/zerodim.h>

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

defineTypeNameAndDebug(canteraChemistryModel, 0);
addToRunTimeSelectionTable(alternateChemistryModel,canteraChemistryModel,steadyChemistry);
addToRunTimeSelectionTable(alternateChemistryModel,canteraChemistryModel,transientChemistry);

// Construct from components
canteraChemistryModel::canteraChemistryModel
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    alternateChemistryModel(thermo,rho),
    relTol_(1e-9),
    absTol_(1e-15),
    relTolSens_(1e-5),
    absTolSens_(1e-4),
    maxSteps_(20000),
    minH_(SMALL),
    useOldImplementation_(false)
{
    if(!isA<hMixtureThermo<canteraMixture> >(thermo)){
        FatalErrorIn("canteraChemistryModel::canteraChemistryModel")
            << " thermo is required to be a hMixtureThermo<canteraMixture>"
                << endl << abort(FatalError);
    }

    RR_.setSize(thermo.composition().Y().size());
    forAll(RR_,i)
    {
        RR_.set
        (
            i,
            new scalarField(rho_.size(), 0.0)
        );
    }

    IOdictionary chemistryPropertiesDict
        (
            IOobject
            (
                "chemistryProperties",
                rho.time().constant(),
                rho.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    if(chemistryPropertiesDict.found("canteraNumerics")) {
        const dictionary &cn=chemistryPropertiesDict.subDict("canteraNumerics");

        relTol_=cn.lookupOrDefault<scalar>("relativeTolerance",1e-9);
        absTol_=cn.lookupOrDefault<scalar>("absoluteTolerance",1e-15);
        relTolSens_=cn.lookupOrDefault<scalar>("relativeSensitivityTolerance",1e-5);
        absTolSens_=cn.lookupOrDefault<scalar>("absoluteSensitivityTolerance",1e-4);
        maxSteps_=cn.lookupOrDefault<label>("maximumSteps",20000);
        minH_=cn.lookupOrDefault<scalar>("minimumStepSize",SMALL);
        useOldImplementation_=cn.lookupOrDefault<bool>("useOldImplementation",false);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

canteraChemistryModel::~canteraChemistryModel()
{}



// ************************************************************************* //

tmp<volScalarField> canteraChemistryModel::tc() const
{

    const canteraMixture &mix=static_cast<const canteraMixture &>(thermo().composition());
    scalar NSp = RR_.size();

    //Simple imlementation of tc() acording to OF1.4	
    /*	scalarList Yi(RR_.size());
	scalarField t(rho_.size(), GREAT);
	forAll(rho_, celli)
	{
        const canteraThermo &gas=mix.cellMixture(celli);
        gas.gas().getMassFractions(Yi.begin());
        forAll(RR_, i)
        {
        if(RR_[i][celli] < -SMALL)
        {
        t[celli] = min(t[celli],
        -Yi[i]/(rho_[celli]*RR_[i][celli]));
        }
        }
	}*/
	
    //Imlementation of tc() acording to OF1.5 - see also posting from Niklas Nordin
    const canteraThermo &gas=mix.cellMixture(0);
    scalar noReac = gas.gas().nReactions();
    scalarList Ci(NSp);
    scalarList Di(noReac);
    scalar CSum;

    scalarField t(rho_.size(), SMALL);

    forAll(rho_, celli)
    {	CSum = 0.0;
        const canteraThermo &gas=mix.cellMixture(celli);
        gas.gas().getConcentrations(Ci.begin());	//[kmol/kmol]	
        gas.gas().getFwdRatesOfProgress(Di.begin());

        forAll(RR_, i)
        {
            CSum += Ci[i];
        }

        for (scalar i=0;i<noReac;i++)
        {
            t[celli] += Di[i];
        }
        t[celli] = (0.5*noReac*CSum)/t[celli];
    }

    tmp<volScalarField> tsource
        (
            new volScalarField
            (
                IOobject
                (
                    "tc",
                    rho_.time().timeName(),
                    rho_.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_.mesh(),
                dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

    tsource().internalField() = t;
    tsource().correctBoundaryConditions();

    return tsource;
}

tmp<volScalarField> canteraChemistryModel::RR(const label i) const
{
    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + thermo().composition().Y()[i].name() + ')',
                rho_.time().timeName(),
                rho_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho_.mesh(),
            dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tRR().internalField() = RR_[i];
    tRR().correctBoundaryConditions();
    
    return tRR;
}

scalar canteraChemistryModel::solve(
    const scalar t0,
    const scalar deltaT
)
{
    scalar dtMin=deltaT;
    //    Info << "deltaT=  " << deltaT << endl;
    const canteraMixture &mix=static_cast<const canteraMixture&>(thermo().composition());

    scalarList RRtmp(RR_.size());
    Reactor react;
    scalar Ns = RR_.size(); //number of species
    const canteraThermo &gas1=mix.cellMixture(0);	
    scalarList mw(RR_.size());
    //delete &gas1;
    gas1.gas().getMolecularWeights(mw.begin());

    forAll(rho_,cellI) {
      try {
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

                FatalErrorIn("canteraChemistryModel::solve")
                    << " With the state: " << gas 
                        << " and the mixture " << c0 
                        << " Cantera complained in cell " << cellI
                        << " with a Cantera::CVodesErr"  << endl
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
      } catch(const Cantera::CanteraError& e) {
          Cantera::showErrors(std::cerr);
          
          FatalErrorIn("canteraChemistryModel::solve")
                  << " Cantera complained in cell " << cellI
                  << " with a Cantera::CanteraError"  << endl
                  << abort(FatalError) ;
      }
    }

    return dtMin; //das ist nur für die Schrittweitenstrg da, denk ich

}

void canteraChemistryModel::setNumerics(ReactorNet &sim) {
    sim.setTolerances(relTol_,absTol_);
    sim.setSensitivityTolerances(relTolSens_,absTolSens_);

    // these two calls don't seem to work
    sim.integrator().setMaxSteps(maxSteps_);
    sim.integrator().setMinStepSize(minH_);

    sim.initialize();
}

void canteraChemistryModel::calcDQ(volScalarField &dQ)
{
    WarningIn("canteraChemistryModel::calcDQ(volScalarField &dQ)")
        << "Strange error concerning the RR-Field. dQ is currently not computed"
            << endl;


    return;

    const canteraMixture &mix=static_cast<const canteraMixture &>(thermo().composition());

    scalarList cp(RR_.size());

    forAll(rho_, cellI)
    {
        const canteraThermo &gas=mix.cellMixture(cellI);
        
        gas.gas().getPartialMolarEnthalpies(cp.begin()); 
        
        Info << cp << endl;

        dQ[cellI]=0;

        forAll(cp,i) {
            dQ[cellI]+=RR_[i][cellI]*cp[i];
        }
    }
    
}

} // namespace Foam
