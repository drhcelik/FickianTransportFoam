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
Class
    Foam::turbulenceThermophysicalTransportModels::constantLewis

Description
    Specie-specific constant Lewis number laminar transport model
    
SourceFiles
    constantLewis.C

\*---------------------------------------------------------------------------*/

#include "basicFickianTransportModel.H"
#include "unityLewisFourier.H"

#ifndef constantLewis_H
#define constantLewis_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class constantLewis Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class constantLewis
:
    public basicFickianTransportModel
    <
        unityLewisFourier<laminarThermophysicalTransportModel>
    >
{


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("constantLewis");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        constantLewis
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~constantLewis()
    {}


    // Member Functions
        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;
        //- Read thermophysicalTransport dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "constantLewis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
