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

#include "basicFickianTransportModel.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
basicFickianTransportModel<BasicThermophysicalTransportModel>::basicFickianTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    bool constantLewis
)
:
    BasicThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),
    constantLewis_(constantLewis),
    Jc_
     (
         IOobject
         (
             thermo.phasePropertyName("Jc"),
             this->thermo().T().mesh().time().name(),
             this->thermo().T().mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->thermo().T().mesh(),
         dimensionedScalar(dimMass/dimArea/dimTime, 0)
     ),

    implicitFlux_(true),
 
    DFuncs_(this->thermo().species().size()),

    DmLimit_(1.0 - 1e-6),
    
    Le_(this->thermo().species().size())
    
{
    if(this->thermo().he().name() == "e") {
        FatalErrorInFunction
            << "Internal energy is not supported as a solution variable. Please use enthalpy instead."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool basicFickianTransportModel<BasicThermophysicalTransportModel>::read()
{
    if
    (
        BasicThermophysicalTransportModel::read()
    )
    {
        const speciesTable& species = this->thermo().species();

	if(constantLewis_) {
	   Info << "Setting Lewis numbers " << endl;
	   dictionary LewisNumberDict(this->coeffDict_.subDict("Le"));

           const PtrList<volScalarField>& Y = this->thermo().Y();
           for (label i=0; i<species.size(); i++)
           {
              if (LewisNumberDict.found(Y[i].member()))
              {
                Le_[i] = readScalar(LewisNumberDict.lookup(Y[i].member()));
                Info<<"Lewis number of specie "<<Y[i].name()<<" is: "<<Le_[i]<<endl;
              }
              else {
           	Info<<"Setting Lewis number of specie "<<Y[i].name()<<" to default value 1.0 "<<endl;
           	Le_[i] = 1;
              }
           }
	}
	else 
	{
            const dictionary& Ddict = this->coeffDict_.subDict("D");

            // Read the array of specie binary mass diffusion coefficient
            // functions
            forAll(species, i)
            {
                DFuncs_[i].setSize(species.size());

                forAll(species, j)
                {
                    if (j >= i)
                    {
                        const word nameij(species[i] + '-' + species[j]);
                        const word nameji(species[j] + '-' + species[i]);

                        word Dname;

                        if (Ddict.found(nameij) && Ddict.found(nameji))
                        {
                            if (i != j)
                            {
                                WarningInFunction
                                    << "Binary mass diffusion coefficients "
                                       "for Both " << nameij
                                    << " and " << nameji << " provided, using "
                                    << nameij << endl;
                            }

                            Dname = nameij;
                        }
                        else if (Ddict.found(nameij))
                        {
                            Dname = nameij;
                        }
                        else if (Ddict.found(nameji))
                        {
                            Dname = nameji;
                        }
                        else
                        {
                            FatalIOErrorInFunction(Ddict)
                                << "Binary mass diffusion coefficient for pair "
                                << nameij << " or " << nameji << " not provided"
                                << exit(FatalIOError);
                        }

                        DFuncs_[i].set
                        (
                            j,
                            Function2<scalar>::New(Dname,                                
                                 dimPressure,
                                 dimTemperature,
                                 dimKinematicViscosity,
                                 Ddict).ptr()
                        );
                    }
                }
            }
            DmLimit_ = this->coeffDict_.lookupOrDefault("selfDiffusionLimit",1.0 - 1e-6);
        }
        
	implicitFlux_ = this->coeffDict_.lookupOrDefault("implicitHeatFlux",true);
	Info << "Selecting "<< (implicitFlux_ ? "implicit" : "explicit") << " formulation for the heat flux" << endl;
	
        return true;
    }
    else
    {
        return false;
    }
}



template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> basicFickianTransportModel<BasicThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    
    if (Y.size())
    {
        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {

                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(BasicThermophysicalTransportModel::j(Y[i]));
                sumJh += ji*fvc::interpolate(hi);
            
        }



        tmpq.ref() += sumJh;
    }
    tmpq.ref() -= Jc_*fvc::interpolate(this->thermo().he());
    return tmpq;
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> basicFickianTransportModel<BasicThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
     correctJc();
     
     tmp<fvScalarMatrix> tmpDivq
     (
        implicitFlux_ ? -fvm::laplacian(this->alpha()*this->alphaEff(), he) :
        fvm::Su 
         (
             -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
             he
         )
     );
     
     const PtrList<volScalarField>& Y = this->thermo().Y();
     const volScalarField& p = this->thermo().p();
     const volScalarField& T = this->thermo().T();

     if (!implicitFlux_) 
     {
     tmpDivq.ref() -=
         correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));
     }

 
     surfaceScalarField sumJh
     (
         surfaceScalarField::New
         (
             "sumJh",
             he.mesh(),
             dimensionedScalar(dimMass/dimArea/dimTime*he.dimensions(), 0)
         )
     );
 
     forAll(Y, i)
     {

             const volScalarField hi(this->thermo().hsi(i, p, T));
 
             const surfaceScalarField ji(BasicThermophysicalTransportModel::j(Y[i]));
 
             sumJh += ji*fvc::interpolate(hi);
             if (implicitFlux_) {
                 sumJh += fvc::interpolate(this->alpha()*this->alphaEff()*hi)*fvc::snGrad(Y[i]);
            
             }
         
     }
 
 
     tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());
     tmpDivq.ref() -= fvm::div(Jc_*he.mesh().magSf(),he,"div(phic,Yi_h)");
     return tmpDivq;
 }

template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> basicFickianTransportModel<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
       return BasicThermophysicalTransportModel::j(Yi) - Jc_*fvc::interpolate(Yi); 
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> basicFickianTransportModel<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const volScalarField& T = this->thermo().T();
    return BasicThermophysicalTransportModel::divj(Yi) - fvm::div(Jc_*T.mesh().magSf(),Yi,"div(phic,Yi_h)");  
}


template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::updateDm() const
{
    if (constantLewis_) return correctJc();
    
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());
    Dm_.setSize(Y.size());
    forAll(Dm_, i)
    {
     Dm_.set
	 (
	     i,
	     volScalarField::New
	     (
		 "Dm_" + Y[i].name(),
		 T.mesh(),
		 dimensionedScalar(dimKinematicViscosity,scalar(0))
	     ));
     }
     volScalarField Dji
	 (
	     volScalarField::New
	     (
		 "Dji",
		 T.mesh(),
		 dimless/dimKinematicViscosity
	     )
	  );  

     forAll(Y, i)
     {    
       for(label j = 0; j < i; j++)
       {
	   Dji = evaluate(DFuncs_[j][i], dimless/dimKinematicViscosity, p, T);
	   Dm_[i] += Y[j]/Dji*(1/this->thermo().Wi(j).value() + (1/this->thermo().Wi(i).value() - 1/this->thermo().Wi(j).value())*Y[i]);
	   Dm_[j] += Y[i]/Dji*(1/this->thermo().Wi(i).value() + (1/this->thermo().Wi(j).value() - 1/this->thermo().Wi(i).value())*Y[j]);
       }
     }         

     forAll(Dm_, i)
     { 

	  volScalarField & Dmi = Dm_[i];
	  Dji.dimensions().reset(dimKinematicViscosity);
	  Dji = evaluate(DFuncs_[i][i], dimKinematicViscosity, p, T);
	  setDm(this->thermo().Wi(i).value(),Dji,Dmi,Y[i](),Wm());
	  forAll(Dmi.boundaryFieldRef(), patchi)
	  { 
	     setDm(this->thermo().Wi(i).value(),Dji.boundaryFieldRef()[patchi],Dmi.boundaryFieldRef()[patchi],Y[i].boundaryField()[patchi],Wm.boundaryField()[patchi]);
	  }
      }     
       
      correctJc();
}

 template<class BasicThermophysicalTransportModel>
 void basicFickianTransportModel<BasicThermophysicalTransportModel>::predict()
 {
     BasicThermophysicalTransportModel::predict();
     updateDm();
 }
  
  
 template<class BasicThermophysicalTransportModel>
 bool basicFickianTransportModel<BasicThermophysicalTransportModel>::movePoints()
 {
     return true;
 }
  
  
 template<class BasicThermophysicalTransportModel>
 void basicFickianTransportModel<BasicThermophysicalTransportModel>::topoChange
 (
     const polyTopoChangeMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::mapMesh
 (
     const polyMeshMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::distribute
 (
     const polyDistributionMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }

 template<class BasicThermophysicalTransportModel>
 const PtrList<volScalarField>&
 basicFickianTransportModel<BasicThermophysicalTransportModel>::Dm() const
 {
     if (!Dm_.size())
     {
         updateDm();
     }
  
     return Dm_;
 }
   
template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::correctJc() const
{
   
   const PtrList<volScalarField>& Y = this->thermo().Y();
   Jc_ *= scalar(0);
   forAll(Y, i)
   {
     Jc_ -= fvc::interpolate(this->alpha()*this->DEff(Y[i]))*fvc::snGrad(Y[i]);
   }  

}

template<class BasicThermophysicalTransportModel>
void basicFickianTransportModel<BasicThermophysicalTransportModel>::setDm(scalar Wi, scalarField& Dii, scalarField& Dmi, const scalarField& Yi, const scalarField& Wm) const {
  forAll(Dmi,i) {
      if(Yi[i] > DmLimit_){
          Dmi[i] = Dii[i];
      } else {
          Dmi[i] = (1 - Yi[i])/Dmi[i]/Wm[i];
      } 
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
