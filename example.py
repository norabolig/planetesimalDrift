# A simple example of 

import  planetesimalDrift as pd

d=pd.disc(R=0.4)
p=pd.pseed(Radius=1e5)

pd.integrate(d,p,tend=500000,tout=10000)

