import math
import copy

class disc:
    ''' class to setup basic disc model '''
  
    def __init__(self, ntemp=0.75, nrho=2.5,T0=300,rho0=1e-9,R=1.0,g2d=0.01,mu=2.3,mstar=1.):
        self.ntemp = ntemp
        self.nrho  = nrho
        self.T0    = T0
        self.rho0  = rho0
        self.g2d = g2d
        self.mu = mu
        self.R  = R


        self.AU = 1.496e13
        self.G  = 6.67e-8
        self.MSUN = 1.989e33
        self.RGAS = 8.314e7

        self.mstar = mstar*self.MSUN
        self.Rcm   = R*self.AU

        self.set_rho(R)
        self.set_T(R)
        self.set_vcirc(R)
        self.set_vrel()

    def set_rho(self,R): 
        self.rho  = self.rho0 * R**(-self.nrho)
        self.rhod = self.rho*self.g2d

    def set_T(self,R):   self.T  = self.T0   * R**(-self.ntemp)

    def set_vcirc(self,R): self.vcirc = math.sqrt(self.G*self.mstar/(R*self.AU))

    def set_vrel(self):
        vorbit = math.sqrt( self.vcirc**2 - (self.nrho+self.ntemp)*self.RGAS*self.T/(self.mu) )
        self.vrel = vorbit - self.vcirc 


    def set_disc_position(self,R):
        self.R=R
        self.Rcm = R*self.AU
        self.set_rho(R)
        self.set_T(R)
        self.set_vcirc(R)
        self.set_vrel()


class pseed:
    ''' planet seed (planetesimal) parameters, such as mass, density, and radius '''

    def __init__(self, Radius=100e5, rhom=3.):
        self.Radius=Radius
        self.rhom=rhom
        self.set_mass(Radius)

    def set_mass(self,Radius): self.mass=4*math.pi*self.Radius**3*self.rhom/3.

    def set_radius(self,mass): self.Radius = (3.*mass/(4.*math.pi*self.rhom))**(1./3.)

    def set_mass_radius(self,mass):
        self.mass=mass
        self.set_radius(mass)

def drag_torque(dObj,pObj,Cd=0.5):
    ''' set rdot from head wind with gas '''

    return 0.75 * Cd * (dObj.rho * dObj.Rcm * abs(dObj.vrel)*dObj.vrel)/(pObj.rhom * pObj.Radius * dObj.vcirc)

def grav_focus(dObj,pObj):
    ''' return gravitational focusing factor for relative wind '''

    return 1. + (8.*math.pi*dObj.G*pObj.Radius**2*pObj.rhom)/(3.*dObj.vrel**2)

def accr_torque(dObj,pObj):
     ''' set rdot from accretion of slightly lower-angular momentum pebbles '''
     focus = grav_focus(dObj,pObj)
     return 1.5 * (dObj.rhod * dObj.Rcm * dObj.vrel * abs(dObj.vrel) * focus ) / (pObj.rhom * pObj.Radius * dObj.vcirc )

def accr_mass(dObj,pObj):
    ''' planetesimal accretes mass, return mdot '''

    focus = grav_focus(dObj,pObj)
    return abs(dObj.vrel)*math.pi*pObj.Radius**2*focus*dObj.rhod

def integrate(dObj,pObj,tstart=0,tend=1e6,tout=1e5,ifac=0.1):
    ''' integrate particle in radius and mass.
        tstart, tend, and tout should be in years.
        Just use simple predictor-corrector for now.'''

    yr2s=365.25*24*3600.
    tstart*=yr2s
    tend  *=yr2s
    tout  *=yr2s

    t=tstart*1.
    tdump=tout+tstart
    print("time distance radius mass")
    print(t,tstart,tend,tout)
    while ( t<tend):

        rdot_d = drag_torque(dObj,pObj)
        rdot_a = accr_torque(dObj,pObj)
        mdot   = accr_mass(dObj,pObj)

        dt_rdot = dObj.Rcm/(abs(rdot_d)+abs(rdot_a))*ifac
        dt_mdot = pObj.mass/abs(mdot)*ifac

        dt = min(dt_rdot,dt_mdot)
        
        print_values = 0
        if( t+dt >= tdump ): 
           dt = tdump-t
           print_values = 1
           tdump+=tout

        Rprime = dObj.Rcm + (rdot_d + rdot_a)*dt
        Mprime = pObj.mass + mdot * dt

        # debug
        #print(Rprime,Mprime,dt,dObj.Rcm,pObj.mass,rdot_d,rdot_a,mdot,t/yr2s)

        # this is inefficient.  Moving on for now.  Should fix
        pObj_trial = copy.deepcopy(pObj)
        dObj_trial = copy.deepcopy(dObj)

        dObj_trial.set_disc_position(Rprime/dObj.AU)
        pObj_trial.set_mass_radius(Mprime)

        # now get rates again with new objects
        rdot_d = drag_torque(dObj_trial,pObj_trial)
        rdot_a = accr_torque(dObj_trial,pObj_trial)
        mdot   = accr_mass(dObj_trial,pObj_trial)

        R = 0.5*(Rprime + dObj.Rcm  +(rdot_d + rdot_a)*dt )
        M = 0.5*(Mprime + pObj.mass + mdot *dt )


        dObj.set_disc_position(R/dObj.AU)
        pObj.set_mass_radius(M)
      
        t+=dt

        if(print_values):
          print(t/yr2s,dObj.R,pObj.Radius,pObj.mass)

    #print(t/yr,dObj.R,pObj.Radius,pObj.mass)
