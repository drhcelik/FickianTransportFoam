# OpenFOAM's counter-flow diffusion flame tutorial example  

A simple demo of running OpenFOAM's standard tutorial with this library. Allrun script copies tutorial files and runs a case with mixtureAveraged transport model


## [optional] Running with [DLBFoam](https://github.com/Aalto-CFD/DLBFoam)

If you have DLBFoam compiled, uncomment corresponding section in the `Allrun` script and set `$DLBPATH` in the script to your DLBFoam installation location.  

Note that DLBFoam does not affect the simulation results, it only affects the simulation time (for this case, simulation will be approximately three times faster).
