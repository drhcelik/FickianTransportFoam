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
    Foam::Function2s::binaryDiffusionCoefficientsPolynomial

Description
    Function for polynomial fit of binary diffusion coefficients with respect
    to logarithm of temperature.

SourceFiles
    binaryDiffusionCoefficientsPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusionCoefficientsPolynomial_H
#define binaryDiffusionCoefficientsPolynomial_H

#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
          Class binaryDiffusionCoefficientsPolynomial Declaration
\*---------------------------------------------------------------------------*/

class binaryDiffusionCoefficientsPolynomial
:
    public FieldFunction2<scalar, binaryDiffusionCoefficientsPolynomial>
{
    // Private Data

        FixedList<scalar, 5> polynomialCoeffs_;


public:

    //- Runtime type information
    TypeName("binaryDiffusionCoefficientsPolynomial");


    // Constructors

        //- Construct from components
        binaryDiffusionCoefficientsPolynomial
        (
            const word& name,
            const FixedList<scalar, 5> polynomialCoeffs
        );

        //- Construct from name and dictionary
        binaryDiffusionCoefficientsPolynomial
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function2<scalar>> clone() const
        {
            return tmp<Function2<scalar>>(
                new binaryDiffusionCoefficientsPolynomial(*this)
            );
        }


    // Member Functions

        //- Inherit base class value methods to prevent clang warnings
        using FieldFunction2<scalar, binaryDiffusionCoefficientsPolynomial>::value;

        //- Binary diffusion coefficient function using pressure and temperature
        virtual scalar value(scalar p, scalar T) const
        {
            scalar one_over_p = 1/p;
            scalar logT = log(T);
            scalar sqrtT = sqrt(T);
            return one_over_p*T*sqrtT*(
                polynomialCoeffs_[0]
              + polynomialCoeffs_[1]*logT
              + polynomialCoeffs_[2]*logT*logT
              + polynomialCoeffs_[3]*logT*logT*logT
              + polynomialCoeffs_[4]*logT*logT*logT*logT
          );
        }

        //- Write the function coefficients
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
