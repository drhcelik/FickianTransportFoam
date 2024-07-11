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
    Foam::turbulenceThermophysicalTransportModels::basicFickianTransportModel

Description
    Base class for Fickian transport
    
SourceFiles
    basicFickianTransportModel.C

\*---------------------------------------------------------------------------*/


#ifndef basicFickianTransportModel_H
#define basicFickianTransportModel_H
#include "Function2.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class basicFickianTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class basicFickianTransportModel
:
    public BasicThermophysicalTransportModel
{
    // Private data
            bool constantLewis_;
            
            // Correction flux
            mutable surfaceScalarField Jc_;
	   	    
	    bool implicitFlux_; 

            //- Array of specie binary mass diffusion coefficient functions
            //  [m^2/s]
            List<PtrList<Function2<scalar>>> DFuncs_;
            
            scalar DmLimit_;
            
            void setDm(scalar Wi, scalarField& Dii, scalarField& Dmi, const scalarField& Yi, const scalarField& Wm) const;
            	    
 	    void correctJc() const;
 	    
            //- Update Dm_
            void updateDm() const;
  
            //- Mass diffusion coefficients in the mixture
            mutable PtrList<volScalarField> Dm_; 
            
protected:

            //- Update Dm_ if not yet constructed and return
            const PtrList<volScalarField>& Dm() const;            
            scalarList Le_;	    
 	    
public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;
    // Constructors

        //- Construct from a momentum transport model and a thermo model
        basicFickianTransportModel
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo,
            bool constantLewis
        );


    //- Destructor
    virtual ~basicFickianTransportModel()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

         //- Update the diffusion coefficients
         virtual void predict();
         
           // Mesh changes
  
             //- Update for mesh motion
             virtual bool movePoints();
  
             //- Update topology using the given map
             virtual void topoChange(const polyTopoChangeMap& map);
  
             //- Update from another mesh using the given map
             virtual void mapMesh(const polyMeshMap& map);
  
             //- Redistribute or update using the given distribution map
             virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicFickianTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
