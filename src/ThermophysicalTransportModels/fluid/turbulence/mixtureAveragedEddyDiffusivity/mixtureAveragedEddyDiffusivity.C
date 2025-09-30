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

#include "mixtureAveragedEddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
mixtureAveragedEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
mixtureAveragedEddyDiffusivity
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
        false
    ),

    Sct_("Sct", dimless, this->coeffDict())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
mixtureAveragedEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
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
mixtureAveragedEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *this->Dm()[this->thermo().specieIndex(Yi)]
      + (this->Prt_/Sct_)*this->alphat()
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{

    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *this->Dm()[this->thermo().specieIndex(Yi)].boundaryField()[patchi]
      + this->Prt_.value()/Sct_.value()*this->alphat(patchi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
