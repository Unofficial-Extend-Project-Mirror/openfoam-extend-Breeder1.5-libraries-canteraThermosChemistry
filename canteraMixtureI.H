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

 ICE Revision: $Id: /local/openfoam/Libraries/canteraThermosChemistry/canteraMixtureI.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "canteraMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const canteraThermo& canteraMixture::cellMixture
(
    const label celli
) const
{
    for (label n=0; n<Y_.size(); n++)
    {
        yTemp_[n] = Y_[n][celli];
        // make sure that none of the species becomes exactly zero: otherwise Cantera Chemistry solver might crash
        if(yTemp_[n]<SMALL) {
            yTemp_[n]=SMALL;
        }
    }
    
    scalar T=T_[celli];
    scalar p=p_[celli];

    speciesData_.setState(T,p,yTemp_);

    return speciesData_;
}


const canteraThermo& canteraMixture::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    for (label n=0; n<Y_.size(); n++)
    {
        yTemp_[n] = Y_[n].boundaryField()[patchi][facei];
    }
    
    scalar T=T_.boundaryField()[patchi][facei];
    scalar p=p_.boundaryField()[patchi][facei];

    speciesData_.setState(T,p,yTemp_);

    return speciesData_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
