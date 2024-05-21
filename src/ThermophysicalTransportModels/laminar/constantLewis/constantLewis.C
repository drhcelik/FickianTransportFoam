/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | FickianTransportFoam 
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  | 2024, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of FickianTransportFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/FickianTransportFoam

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "constantLewis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
constantLewis<laminarThermophysicalTransportModel>::
constantLewis
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    basicFickianTransportModel<unityLewisFourier<laminarThermophysicalTransportModel>>
    (
        typeName,
        momentumTransport,
        thermo,
        true
    )
{
    read();
    this->correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool
constantLewis<laminarThermophysicalTransportModel>::read()
{
    return basicFickianTransportModel
    <
        unityLewisFourier<laminarThermophysicalTransportModel>
    >::read();
}

template<class laminarThermophysicalTransportModel>
tmp<volScalarField> constantLewis<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();

    return volScalarField::New
    (
        "DEff",
        this->alphaEff()
       /this->Le_[composition.index(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField> constantLewis<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();

    return
        this->alphaEff(patchi)
       /this->Le_[composition.index(Yi)];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
