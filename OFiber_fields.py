'''
Copyright (C) 2022 Sergio G Rodrigo sergut@unizar.es

OFiber is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as published by 
the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

OFiber is distributed in the hope that it will be useful for 
research or/and academic purpouses, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License 
for more details. You should have received a copy of the GNU Affero General 
Public License along with OFiber. If not,see http://www.gnu.org/licenses/.

Commercial applications may also acquire a commercial license. 
Please contact Sergio G Rodrigo sergut@unizar.es
'''

from OFiber_class import OFiber
import numpy as np

class EMfield(OFiber):
   '''
   EM fields class:
   Snyder&Love: Section 12-8 Table 12-3 Field components of all the modes.
   Class EMfield inherits class OFiber as a means to have simple access to the 
   electromagnetic fields. Therefore, the electromagnetic fields and the quantities 
   necesary to define them are implented in OFiber class.      
   '''
   def __init__(self,U=None,V=None,nu=None,rho=None,nco=None,ncl=None,delta=None,mode='hybrid'):
    super().__init__(rho,nco,ncl,delta) 
    self.V=V
    self.U=U
    self.nu=nu   
    self.mode=mode    
    pass
   
   def __str__(self):
    super().__str__()
    msg= f'\nEMfields for mode= {self.mode} V = {self.V}; U= {self.U}; order(nu)= {self.nu}'
    return msg

   def E_intensity(self):
    U=self.U; V=self.V; nu=self.nu
    A1=self.ex()
    A2=self.ey()
    A3=self.ez()
    def Intensity(R,phi):
        E2=np.abs(A1(R,phi))**2+np.abs(A2(R,phi))**2+np.abs(A3(R,phi))**2 
        return E2
    return Intensity

   def ex(self):
    U=self.U; V=self.V; nu=self.nu
    #Warning: no arguments are passed because is self for this class
    er=self.er() 
    ephi=self.ephi()
    def Ex(R,phi):
      return er(R,phi)*np.cos(phi)-ephi(R,phi)*np.sin(phi)      
    return Ex

   def ey(self):
    U=self.U; V=self.V; nu=self.nu
    #Warning: no arguments are passed because is self for this class
    er=self.er()    
    ephi=self.ephi()
    def Ey(R,phi):
      return er(R,phi)*np.sin(phi)+ephi(R,phi)*np.cos(phi)      
    return Ey

   def ez(self):
     '''
     Ez-field calculations for given {V,U,nu,mode}
     '''
     func=super().ez(self.U,self.V,self.nu,self.mode)
     return func

   def er(self):
     '''
     E_r-field calculations for given {V,U,nu,mode}
     '''
     func=super().er(self.U,self.V,self.nu,self.mode)
     return func

   def ephi(self):
     '''
     Ephi-field calculations for given {V,U,nu,mode}
     '''
     func=super().ephi(self.U,self.V,self.nu,self.mode)
     return func

   def hx(self):
    U=self.U; V=self.V; nu=self.nu
    #Warning: no arguments are passed because is self for this class
    hr=self.hr() 
    hphi=self.hphi()
    def Hx(R,phi):
      return hr(R,phi)*np.cos(phi)-hphi(R,phi)*np.sin(phi)      
    return Hx

   def hy(self):
    U=self.U; V=self.V; nu=self.nu
    #Warning: no arguments are passed because is self for this class
    hr=self.hr()    
    hphi=self.hphi()
    def Hy(R,phi):
      return hr(R,phi)*np.sin(phi)+hphi(R,phi)*np.cos(phi)      
    return Hy

   def hz(self):
     '''
     Hz-field calculations for given {V,U,nu,mode}
     '''
     func=super().hz(self.U,self.V,self.nu,self.mode)
     return func

   def hr(self):
     '''
     Hr-field calculations for given {V,U,nu,mode}
     '''
     func=super().hr(self.U,self.V,self.nu,self.mode)
     return func

   def hphi(self):
     '''
     Hphi-field calculations for given {V,U,nu,mode}
     '''
     func=super().hphi(self.U,self.V,self.nu,self.mode)
     return func
 
   def Sz_analytic(self):
     '''
     sz-field calculations for given {V,U,nu,mode}
     '''
     func=super().Sz(self.U,self.V,self.nu,self.mode)
     return func
 
   def Sz_numeric(self):
        ex=self.ex()
        ey=self.ey()
        hx=self.hx()
        hy=self.hy()
        def curr(R,phi):
          return 0.5*(ex(R,phi)*hy(R,phi)-ey(R,phi)*hx(R,phi))
        return curr
    
   def sum_fields(self,fields):
        def superpositon(R,phi):
            fs=0.0
            for f in fields:                
                fs+=f(R,phi)
            return fs
        return superpositon
    
   def integrate_Sz(self,Rmax,U,V,nu,mode,nrcore=100,nrcladding=50,nphi=100,analytic=True):        
        from math import pi
        #Field(R,phi)    
        if(not(analytic)):
            sz= self.Sz_numeric()  # To much slow (due to loop)              
        else:
            sz=self.Sz_analytic()
        Nco=0.0      
        radial_distance=1.0 # 1.0 for normalized radius R=r/rho rho=radius
        drdphi=radial_distance*2.0*pi/(nrcore*nphi)
        for r in np.linspace(0, 1, nrcore):
            for phi in np.radians(np.linspace(0, 360, nphi)):
              Nco+=r*drdphi*sz(r,phi)        
        print('Nco=',Nco)              
        Ncl=0.0 
        radial_distance=Rmax-radial_distance
        drdphi=radial_distance*2.0*pi/(nrcladding*nphi)
        for r in np.linspace(1, Rmax, nrcladding):
            for phi in np.radians(np.linspace(0, 360, nphi)):
              Ncl+=r*drdphi*sz(r,phi)        
        print('Ncl=',Ncl)    
        frac_power=Nco/(Nco+Ncl)
        return frac_power
           
           
