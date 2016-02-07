# planetesimalDrift

A python tool for integrating the decaying orbit of a planetesimal.  Gas drag and "pebble" drag is included.

## Use

The python library can be imported using

> import planetesimalDrift as pd

which gives access to all of the tools. The basic disc configuration can be setup using

> d=pd.disc()

creating a disc model in object d.  The current conditions are for a distance of 1 AU about a solar-mass star.

The planetesimal object can be created using something like

> p=pd.pseed(Radius=1e5)

for a planetesimal with radius 1e5 cm (1 km).  

## Disc model

The disc model assumes a power law in midplane gas density (rho = rho0 R^-nrho
for disc distance R in AU) and a power law in midplane temperature (T = T0
R^-ntemp) The value T0 is the temperature at 1 AU and rho0 is the density at 1
AU.  Options for the disc model are as follows:

* ntemp=0.75 // index
* nrho=2.5   // index
* T0=300  // K
* rho0=1e-9 // g/cc
* R=1.0  // AU
* g2d=0.01 // fraction
* mu=2.3 // g/mol
* mstar=1. // solar masses

where g2d is the gas-to-dust ratio, mu is the mean weight of the gas (g/mol), and mstar is the mass of the star in solar masses.

## Planetesimal model

The planetesimal is spherical, characterized by radius Radius and average internal density rhom.  Defaults are

* Radius=100e5 // cm
* rhom=3. // g/cc

## Integration

Integration uses a simple predictor-corrector for now.  Integration requires the disc object dObj, planetesimal object pObj.
Additional options and defaults are

* tstart=0 // yr
* tend=1e6 // yr
* tout=1e5 // yr
* ifac // fraction

where tstart, tend, and tout are the starting, ending, and output time interval, respectively.  The fraction ifac is a 
factor that controls the size of the time step determined from the derivatives (smaller value takes smaller time steps). 

Example:

>pd.integrate(d,p,tend=500000,tout=10000)

## Example script and output

A simple example script is included in the repo, called "example.py".  It can be invoked by 

> python3 example.py > output.dat

where output.dat is just some file.  Output is in the format time (yr), disc distance (AU), planetesimal size (cm), and planetesimal mass (g).



