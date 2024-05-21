# FickianTransportFoam 
**Stable** mixture-averaged and specie-specific constant Lewis transport models that are currently missing from the current OpenFOAM release. Supported OpenFoam versions are: 10.

Developed by Aleksi Rintanen & Ilya Morev, Aalto University, Finland
## Theory
### Mass diffusion fluxes
FickianTransportFoam implements a transport model based on Fick's law using the formulation, where diffusion fluxes are evaluated with respect to the mass fraction gradient as
```math
\tag{1}
\mathbf{j}_k = -\rho D_k \nabla Y_k,
```
where $D_k$ is a model specific diffusion coefficient.   
To ensure mass conservation, a correction velocity is introduced to the diffusion fluxes [1] as
```math
 \mathbf{j}_k = -\rho D_k \nabla Y_k + \rho Y_k \sum_{j=1}^{N} D_j \nabla Y_j $$
```


Diffusion coefficients are evaluted based on two assumptions:
1. **Mixture-averaged diffusion coefficients**  
   Diffusion coefficients are given by Eq. 12.178 in [1]
   $$D_k = \left(\sum_{j\neq k} \frac{X_j}{D_{jk}} + \frac{X_k}{1-Y_k}\sum_{j\neq k} \frac{Y_j}{D_{jk}}\right)^{-1}$$
3. **Constant Lewis number diffusion coefficients**
   $$D_k = \frac{\lambda}{c_p}\frac{1}{\mathrm{Le}}$$

In LES-models, additional subgrid scale diffusion is added using eddy diffusivity concept similarly as in the OpenFOAM's native transport models. Turbulent mass diffusion coefficients are calculated from the turbulent thermal diffusivity using constant turbulent Prantl and Schmidt numbers defined in the model dictionary. 

### Heat fluxes
In the enthalpy equation, heat flux $\mathbf{q}$ is given assuming Fourier's law as
$$\mathbf{q} = -\lambda \nabla T + \sum_{k=1}^N h_k \mathbf{j}_k   $$
Calculating $\nabla\cdot\bf q$ is problematic, since the term $-\nabla\cdot \lambda \nabla T$ cannot added directly explicitly for stability reasons.
This library implements two models for evaluating the term $\nabla\cdot\bf q$ in entalphy equation:
1. **Add the $-\nabla\cdot \lambda \nabla T$ explicitly**  
    A correction term is also included that stabilizes the equation.  
    This approach is used in the native implementation of FickianFourier in Openfoam 
2. **Reformulate in terms of entalphy (recommended)**
    $$-\lambda \nabla T = -\frac{\lambda}{c_p} \nabla h_s  + \frac{\lambda}{c_p}\sum_{k=1}^N h_k \nabla Y_k$$
    The implicit formulation is not compatible with the default implementation of a coupled temperature boundary condition in OpenFOAM, which      balances the heat fluxes using temperature gradients [2]. Thus, for conjugate heat transfer applications, the native explicit formulation should be used, when the default coupled temperature boundary condition is used.

## Why OpenFOAM's native implementation "FickianFourier" is unstable? 
#### 1. Correction velocity is missing
The disparities in diffusion velocities are dumped to the inert specie to ensure mass conservation. This works fine, when the mixture is diluted greatly (for example the inlet is premixed mixture of air and fuel and N2 is thus the dominant specie), but not when there is a separate fuel inlet.

#### 2. Mixture-averaged diffusion coefficients do not behave well
The formula for mixture-averaged diffusion coefficients is undefined when $Y_k\rightarrow 1$ (leads to $\frac{0}{0}$).   
(Side note: OpenFoam uses a wrong formulation for mixture-averaged diffusion coefficients)
$$D_k = \frac{1-X_k}{\sum_{j\neq k} X_j/D_{jk}}$$
In FickianFourier this is "fixed" by adding $\epsilon$ to the denominator, which is defined by default as "small" ($\approx10^{-16}$)
$$D_k = \frac{1-X_k}{\sum_{j\neq k} X_j/D_{jk} + \epsilon}$$
This causes the diffusion coefficient to go zero when $X_k\rightarrow 1$ and negative when $X_k\gt1$. This is both unphysical and numerically unstable leading to divergent simulations. 

In this implementation, the issue is handled by setting $D_k = D_{kk}$, when the denominator is smaller than a threshold (DmLimit), where $D_{kk}$ is the theoretical limit.

## Usage
Include in controlDict by libFickianTransport.so

Available models:

```
Laminar:

mixtureAveraged
constantLewis

LES/RAS:

mixtureAveragedEddyDiffusivity
constantLewisEddyDiffusivity
```

Diffusion coefficients are given pairwise for mixture-averaged formulation, either as constants or polynomials.
In constantLewis, Lewis numbers are included in the subdictionary as
```
Le {
  H2 1.0;
  O2 1.0;
  N2 1.0;
}
```

Example for dictionary entry in thermophysicalTransport: 

```
laminar
{
    model           mixtureAveraged;
    D
    {
      O-H2 {
        type binaryDiffusionCoefficientsPolynomial;
        polynomialCoeffs (-9.65147826e-03  5.71246016e-03 -1.09656624e-03  9.84367103e-05 -3.22064414e-06);
        }
      O-H {
        type binaryDiffusionCoefficientsPolynomial;
        polynomialCoeffs (-2.64601309e-02  1.37960208e-02 -2.44482561e-03  1.99897138e-04 -5.97357847e-06);
        }
    // ....     
    }
    implicitHeatFlux true; //Default true
    DmLimit 1e-10;
}
```

## References
[1] Kee, R. J., Coltrin, M. E. & Glarborg, P. Chemically Reacting Flow: Theory and Practice (John Wiley & Sons, 2003).

[2] R. Tuominen, Coupling Serpent and OpenFOAM for neutronics - CFD multi-physics calculations. Master's thesis, Aalto university, Espoo, Helsinki, Aug. 2015
