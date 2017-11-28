# vumat

Working draft for material model implementations in Abaqus/Explicit.

Written in Fortran 90/95

Compiles with Abaqus 6.14-4, Visual Studio 2010 x64, Intel Fortran 2016 using Abaqus in windows command line using

abaqus double job=modelabq user=kinematic Int

For Abaqus to accept Fortran 90/95, it might be necessary to add "/free" compiler flag in abaqus environment file.
