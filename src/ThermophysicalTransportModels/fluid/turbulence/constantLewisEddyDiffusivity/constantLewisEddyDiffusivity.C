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

#include "constantLewisEddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
constantLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
constantLewisEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    basicFickianTransportModel<unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>>
    (
        typeName,
        momentumTransport,
        thermo,
        true
    ),

    Sct_("Sct", dimless, this->coeffDict_)
{
    read();
    this->printCoeffs(typeName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
constantLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if
    (
        basicFickianTransportModel
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::read()
    )
    {
        Sct_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volScalarField>
constantLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->thermo().kappa()/this->thermo().Cp()
        /this->Le_[this->thermo().specieIndex(Yi)]
      + (this->Prt_/Sct_)*this->alphat()
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
constantLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->thermo().kappa().boundaryField()[patchi]
       /this->thermo().Cp().boundaryField()[patchi]
       /this->Le_[this->thermo().specieIndex(Yi)]
      + this->Prt_.value()/Sct_.value()*this->alphat(patchi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
