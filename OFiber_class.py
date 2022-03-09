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

from math import pi
from scipy.special import jv,kv,jvp,kvp,jn_zeros
import numpy as np
from numpy import cos,sin
import pandas as pd

'''
Physical constants
'''
c0= 299792458.0 #m/s
eps0= 8.854187817e-12    #F·m-1 
mu0= 1.0/(eps0*(c0**2)) #N·A-2
'''
Conversion factors
'''
con_wltoeV=1239.828330

class OFiber(object):  
  '''
  ---------------------------------------------------------------------------
  Step-profile Optical Fiber class. 
  The object catch the main characteristics of an optical fiber. The notation 
  from Snyder&Love book is followed ("Optical Waveguide Theory", A.W. Snyder 
  and J. Love Springer, Boston, MA (1983) 1st Ed 
  https://doi.org/10.1007/978-1-4613-2813-1).)
  
  It includes:
  - Eigenvalue equations and cuttoff expressions (Table 12-4)
  - Field components of the exact modes (HE,EH,TE and TM - Table 12-3)
  - Modal properties of the step-profile fiber (Table 12-5)
  ---------------------------------------------------------------------------
  '''
  def __init__(self,rho_microns=1.0,nco=None,ncl=None,delta=None):
    self.rho=rho_microns # in microns
    self.nco=nco
    self.ncl=ncl
    self.delta=delta    
    self.dr=pd.DataFrame()
    if(self.delta!=None):
      if(nco==None):
        self.nco=self.nco_from_delta()
      if(ncl==None):
        self.ncl=self.ncl_from_delta()  
    else:
      self.delta=self.Delta()
    pass

  def __str__(self):    
    if(not(self.dr.empty)):
      msg= f'rho (microns)= {self.rho}; nco= {self.nco}; ncl= {self.ncl}; $\Delta$= {self.delta};\n V-U=\n {self.dr.head()}'        
    else:
      msg= f'rho (microns)= {self.rho}; nco= {self.nco}; ncl= {self.ncl}; $\Delta$= {self.delta}'        
    return msg
  
  def import_dispersion_relation(self,file_name):
    import pandas as pd
    self.dr=pd.read_csv(file_name,sep=';',index_col=False)    
    pass    

  '''
  Definitions of adimensional parameters
  '''
  def Delta(self):
    return 0.5*(1-self.ncl**2/self.nco**2)

  def V(self,lmbda):    
    k=2.0*pi/lmbda
    return k*self.rho*np.sqrt(self.nco**2-self.ncl**2)
  
  def U(self,lmbda,beta):   
    k=2.0*pi/lmbda
    return self.rho*np.sqrt(k**2*self.nco**2-beta**2)

  def W(self,lmbda,beta):    
    k=2.0*pi/lmbda
    return self.rho*np.sqrt(beta**2-k**2*self.ncl**2)
  
  def b(U,V): # From Ghatak book "Introduction to fiber optics"
    return 1.0-(U/V)**2
  
  '''
  Radius,beta,k...from adimensinal paramters
  '''  
  def rho_from_V(self,V,lmbda):    
    k=2.0*pi/lmbda
    return V/(k*np.sqrt(self.nco**2-self.ncl**2))

  def beta_from_U(self,U,lmbda):    
    k=2.0*pi/lmbda
    return np.sqrt(k**2*self.nco**2-(U/self.rho)**2)

  def k_from_V(self,V):        
    return V/(self.rho*np.sqrt(self.nco**2-self.ncl**2))

  def omega_from_V(self,V):# eV
    lmbda=2.0*np.pi*(self.rho*np.sqrt(self.nco**2-self.ncl**2))/V
    lmbda=lmbda/1e-3 #nm
    return con_wltoeV/lmbda

  def lambda_from_V(self,V):        
    return 2.0*np.pi*(self.rho*np.sqrt(self.nco**2-self.ncl**2))/V

  def ncl_from_delta(self):
    return self.nco*np.sqrt(1.0-2*self.delta)

  def nco_from_delta(self):
    return self.ncl/np.sqrt(1.0-2*self.delta)

  def beta_from_UV(self,k,U,V):
     return k*self.nco*np.sqrt(1.0-2.0*self.delta*(U/V)**2)
  '''
  Auxiliary definitions for EM modes
  '''
  def fnu(self,nu):    
      def cosnu(phi): return cos(nu*phi)
      def sinnu(phi): return sin(nu*phi)            
      if(nu%2 == 0):
        return cosnu #even
      else:
        return sinnu #odd

  def gnu(self,nu):
      def cosnu(phi): return cos(nu*phi)
      def sinnu(phi): return -sin(nu*phi)      
      if(nu%2 == 0):
        return sinnu
      else:
        return cosnu

  def b1(self,nu,U): 
    return  (jv(nu-1,U)/jv(nu,U)-jv(nu+1,U)/jv(nu,U))/(2.0*U)

  def b2(self,nu,W): 
    return -(kv(nu-1,W)/kv(nu,W)+kv(nu+1,W)/kv(nu,W))/(2.0*W)

  def F1(self,nu,U,V,W): 
    b1=self.b1(nu,U)
    b2=self.b2(nu,W)
    return ((U*W/V)**2)*(b1+(1-2.0*self.delta)*b2)/nu

  def F2(self,nu,U,V,W): 
    b1=self.b1(nu,U)
    b2=self.b2(nu,W)
    return ((V/(U*W))**2)*(nu/(b1+b2))

  def a1(self,nu,U,V,W):
     return (self.F2(nu,U,V,W)-1.0)/2.0
  
  def a2(self,nu,U,V,W):
     return (self.F2(nu,U,V,W)+1.0)/2.0

  def a3(self,nu,U,V,W):
     return (self.F1(nu,U,V,W)-1.0)/2.0

  def a4(self,nu,U,V,W):
     return (self.F1(nu,U,V,W)+1.0)/2.0
  
  def a5(self,nu,U,V,W):
    return self.a3(nu,U,V,W)+self.delta

  def a6(self,nu,U,V,W):
    return self.a4(nu,U,V,W)-self.delta

  '''
  Dispersion relation transcendental functions
  '''  
  def em_modes(self,case,verbose=False):
    '''
    Definitions from Snyder and Love "Optical Waveguide Theory" Section 12
    for step-profile fibers
    '''
    def dispersion_relation_exact_hybrid(U,V,nu):  
        '''
        Snyder&Love: Section 12-9 Table 12-4 (a) Eigenvalue equations for the 
        step-profile fiber (HEnum and EHnum modes)        
        '''
        nco  =self.nco
        ncl  =self.ncl
        k    =self.k_from_V(V)    
        beta =self.beta_from_UV(k,U,V)
        W=np.sqrt(V**2-U**2)         
        fac_left1= (jvp(nu,U)/(U*jv(nu,U))+kvp(nu,W)/(W*kv(nu,W))) 
        fac_left2= (jvp(nu,U)/(U*jv(nu,U))+(ncl**2/nco**2)*kvp(nu,W)/(W*kv(nu,W)))
        fac_rigth1=(nu*beta/(k*nco))**2
        fac_rigth2=(V/(U*W))**4
        return fac_left1*fac_left2 -fac_rigth1*fac_rigth2

    def dispersion_relation_exact_TE(U,V,nu):    
        '''
        Snyder&Love: Section 12-9 Table 12-4 (a) Eigenvalue equations for the 
        step-profile fiber (TE0m modes)        
        '''              
        W=np.sqrt(V**2-U**2)   
        return jv(1,U)/(U*jv(0,U))+kv(1,W)/(W*kv(0,W))

    def dispersion_relation_exact_TM(U,V,nu):
        '''
        Snyder&Love: Section 12-9 Table 12-4 (a) Eigenvalue equations for the 
        step-profile fiber (TM0m modes)        
        '''
        nco  =self.nco
        ncl  =self.ncl        
        W=np.sqrt(V**2-U**2)   
        return (self.nco**2)*jv(1,U)/(U*jv(0,U))+(self.ncl**2)*kv(1,W)/(W*kv(0,W))                

    def dispersion_relation_exact_hybrid_alt(U,V,nu):  
        '''
        Snyder&Love: Section 12-9 Table 12-4 (a) Eigenvalue equations for the 
        step-profile fiber (HEnum and EHnum modes alternative form)        
        '''
        nco  =self.nco
        ncl  =self.ncl
        k    =self.k_from_V(V)    
        beta =self.beta_from_UV(k,U,V)
        W=np.sqrt(V**2-U**2) 

        if(nu==0): nu=1e-12 # To avoid 0/0 errors        
        fac_left1= (k*nco)**2            
        fac_left2= self.F1(nu,U,V,W)
        fac_rigth1=beta**2        
        fac_rigth2=self.F2(nu,U,V,W)
        return  fac_left1*fac_left2-fac_rigth1*fac_rigth2

    def dispersion_relation_weakly_guiding(U,V,nu):
        '''
        Snyder&Love: Section 12-11 Table 12-6 Weak-guidance expansion for the even
        HE11 mode eigenvalue equation in a step-profile fiber
        '''
        # Order number nu is included for consistency with the rest of dispersion 
        # relations but unusued in this    
        try:        
          W=np.sqrt(V**2-U**2)        
        except:
          print("W=np.sqrt(V**2-U**2) cannot be imaginary") 
        try:   
          fac_left1= (U*jv(1,U)/(jv(0,U))-W*kv(1,W)/(kv(0,W))) #weakly guiding  
        except:
          fac_left1=1e12
        return fac_left1              

    func={'hybrid':dispersion_relation_exact_hybrid,
          'TE':dispersion_relation_exact_TE,
          'TM':dispersion_relation_exact_TM,
          'weak':dispersion_relation_weakly_guiding,
          'hybrid-alt':dispersion_relation_exact_hybrid_alt}
    return func[case]

  '''
  Cut-off of modes
  '''
  def cutoff_em_modes_exact(self,case='TE0m',verbose=False):
      '''
      Snyder&Love: Section 12-9 Table 12-4 (b) Limiting values of the modal 
      parameter U (cutoff W=0, U=V)
      '''
      def TE0m(nt):
        return jn_zeros(0, nt) # First zero order nu; nt=number of zeros to return
      
      def TM0m(nt):
        return jn_zeros(0, nt) # First zero order nu; nt=number of zeros to return

      def HE1m(nt):
        return jn_zeros(1, nt) # First zero order nu; nt=number of zeros to return

      def EHnum(nu,nt):
        return jn_zeros(nu, nt) # First zero order nu; nt=number of zeros to return

      def HEnum(U,nu):
        delta=self.delta                
        left=U*jv(nu-2,U)/((nu-1)*jv(nu-1,U))
        right=-2.0*delta/(1-2.0*delta)
        return left-right

      func={'TE0m':TE0m,
            'TM0m':TM0m,
            'HE1m':HE1m,
            'EHnum':EHnum,
            'HEnum':HEnum}
      return func[case]

  '''
  EM fields
  '''
  def ez(self,U,V,nu,mode):
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      A=-1j*(U/(rho*beta))
      f=self.fnu(nu) 
      def Ez(R,phi):            
        if(np.max(R)<=1.0):
            return np.imag(A*jv(nu,U*R)*f(phi)/jv(nu,U))
        else:
            return np.imag(A*kv(nu,W*R)*f(phi)/kv(nu,W))      
      def Ez_TE(R,phi):         
        return R-R
      def Ez_TM(R,phi):                      
        B=-1j*(W/(rho*beta)) 
        if(np.max(R)<=1.0):
            return np.imag(-A*jv(0,U*R)/jv(1,U))
        else:
            return np.imag(B*(self.nco**2/self.ncl**2)*kv(0,W*R)/kv(1,W))
      func={'hybrid':Ez,
            'TE':Ez_TE,
            'TM':Ez_TM}
      return func[mode]

  def er(self,U,V,nu,mode ):  
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      a1=self.a1(nu,U,V,W)
      a2=self.a2(nu,U,V,W)
      f=self.fnu(nu)    
      def Er(R,phi):         
        if(np.max(R)<=1.0):
            return       -(a1*jv(nu-1,U*R)+a2*jv(nu+1,U*R))*f(phi)/jv(nu,U)
        else:
            return -(U/W)*(a1*kv(nu-1,W*R)-a2*kv(nu+1,W*R))*f(phi)/kv(nu,W)      
      def Er_TE(R,phi):         
        return R-R
      def Er_TM(R,phi):         
        if(np.max(R)<=1.0):
            return jv(1,U*R)/jv(1,U)
        else:
            return (self.nco**2/self.ncl**2)*kv(1,W*R)/kv(1,W)
      func={'hybrid':Er,
            'TE':Er_TE,
            'TM':Er_TM}
      return func[mode]

  def ephi(self,U,V,nu,mode): 
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      ''' 
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      a1=self.a1(nu,U,V,W)
      a2=self.a2(nu,U,V,W)
      g=self.gnu(nu)
      def Ephi(R,phi):                     
        if(np.max(R)<=1.0):
            return       -(a1*jv(nu-1,U*R)-a2*jv(nu+1,U*R))*g(phi)/jv(nu,U)
        else:
            return -(U/W)*(a1*kv(nu-1,W*R)+a2*kv(nu+1,W*R))*g(phi)/kv(nu,W)            
      def Ephi_TM(R,phi):         
        return R-R
      def Ephi_TE(R,phi):         
        if(np.max(R)<=1.0):
            return -jv(1,U*R)/jv(1,U)
        else:
            return -kv(1,W*R)/kv(1,W)
      func={'hybrid':Ephi,
            'TE':Ephi_TE,
            'TM':Ephi_TM}
      return func[mode]

  def ex(self,U,V,nu,mode):
    er=self.er(U,V,nu,mode)
    ephi=self.ephi(U,V,nu,mode)
    def Ex(R,phi):
      return er(R,phi)*cos(phi)-ephi(R,phi)*sin(phi)      
    return Ex

  def ey(self,U,V,nu,mode):
    er=self.er(U,V,nu,mode)
    ephi=self.ephi(U,V,nu,mode)
    def Ey(R,phi):
      return er(R,phi)*sin(phi)+ephi(R,phi)*cos(phi)      
    return Ey

  def hz(self,U,V,nu,mode):  
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      '''
      k  = self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W  = np.sqrt(V**2-U**2)  
      F2 = self.F2(nu,U,V,W)      
      rho= self.rho
      B  =-1j*np.sqrt(eps0/mu0)*(U*F2/(rho*k))
      g=self.gnu(nu) 
      def Hz(R,phi):            
        if(np.max(R)<=1.0):
            return np.imag(B*jv(nu,U*R)*g(phi)/jv(nu,U))
        else:
            return np.imag(B*kv(nu,W*R)*g(phi)/kv(nu,W))
      def Hz_TM(R,phi):         
        return R-R
      def Hz_TE(R,phi):                      
        C=-1j*np.sqrt(eps0/mu0)/(k*beta)
        if(np.max(R)<=1.0):
            return np.imag(C*U*jv(0,U*R)/jv(1,U))
        else:
            return np.imag(-C*W*kv(0,W*R)/kv(1,W))

      func={'hybrid':Hz,
            'TE':Hz_TE,
            'TM':Hz_TM}
      return func[mode]

  def hr(self,U,V,nu,mode): 
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      ''' 
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      a3=self.a3(nu,U,V,W)
      a4=self.a4(nu,U,V,W)
      a5=self.a5(nu,U,V,W)
      a6=self.a6(nu,U,V,W)
      Eps=np.sqrt(eps0/mu0)*k*self.nco**2/beta
      g=self.gnu(nu)          
      def Hr(R,phi):         
        if(np.max(R)<=1.0):
            return Eps*(a3*jv(nu-1,U*R)-a4*jv(nu+1,U*R))*g(phi)/jv(nu,U)
        else:
            return Eps*(U/W)*(a5*kv(nu-1,W*R)+a6*kv(nu+1,W*R))*g(phi)/kv(nu,W)
      def Hr_TM(R,phi):         
        return R-R
      def Hr_TE(R,phi):                      
        C=np.sqrt(eps0/mu0)*beta/k
        if(np.max(R)<=1.0):
            return C*jv(1,U*R)/jv(1,U)
        else:
            return C*kv(1,W*R)/kv(1,W)

      func={'hybrid':Hr,
            'TE':Hr_TE,
            'TM':Hr_TM}
      return func[mode]

  def hphi(self,U,V,nu,mode):  
      '''
      Snyder&Love: Section 12-8 Table 12-3  
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      a3=self.a3(nu,U,V,W)
      a4=self.a4(nu,U,V,W)
      a5=self.a5(nu,U,V,W)
      a6=self.a6(nu,U,V,W)
      Eps=np.sqrt(eps0/mu0)*k*self.nco**2/beta
      f=self.fnu(nu)          
      def Hphi(R,phi):         
        if(np.max(R)<=1.0):
            return -Eps*(a3*jv(nu-1,U*R)+a4*jv(nu+1,U*R))*f(phi)/jv(nu,U)
        else:
            return -Eps*(U/W)*(a5*kv(nu-1,W*R)-a6*kv(nu+1,W*R))*f(phi)/kv(nu,W)
      def Hphi_TE(R,phi):         
        return R-R
      def Hphi_TM(R,phi):                      
        C=np.sqrt(eps0/mu0)*(self.nco**2)*k/beta
        if(np.max(R)<=1.0):
            return C*jv(1,U*R)/jv(1,U)
        else:
            return C*kv(1,W*R)/kv(1,W)

      func={'hybrid':Hphi,
            'TE':Hphi_TE,
            'TM':Hphi_TM}
      return func[mode]

  def hx(self,U,V,nu,mode):
    hr=self.hr(U,V,nu,mode)
    hphi=self.hphi(U,V,nu,mode)
    def Hx(R,phi):
      return hr(R,phi)*cos(phi)-hphi(R,phi)*sin(phi)      
    return Hx

  def hy(self,U,V,nu,mode):
    hr=self.hr(U,V,nu,mode)
    hphi=self.hphi(U,V,nu,mode)
    def Hy(R,phi):
      return hr(R,phi)*sin(phi)+hphi(R,phi)*cos(phi)      
    return Hy

  '''
  Modal properties of the step-profile fiber: Nco,Ncl,eta and group velocity  
  '''
  def v_group(self,U,V,nu,mode):
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)  
      delta=self.delta          
      fac1=c0*beta/(self.nco**2*k)
      fac2=1.0-2.0*delta*(1-self.frac_power_core(U,V,nu,mode))
      return fac1/fac2

  def frac_power_core(self,U,V,nu,mode):
    Nco_=self.Nco(U,V,nu,mode)
    Ncl_=self.Ncl(U,V,nu,mode)
    return Nco_/(Nco_+Ncl_)

  def Nco(self,U,V,nu,mode):  
      '''
      Snyder&Love: Section 12-10 Table 12-5 
      (see also Table 11-1 Properties of bound modes)
      Nco=(P/|a|**2) -> Normalized modal power in the core
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      a1=self.a1(nu,U,V,W)
      a2=self.a2(nu,U,V,W)      
      a3=self.a3(nu,U,V,W)
      a4=self.a4(nu,U,V,W)
      if(mode=='hybrid'):      
         a=(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*k*self.nco**2/(beta*jv(nu,U)**2)
         s1=a1*a3*(jv(nu-1,U)**2-jv(nu,U)*jv(nu-2,U))
         s2=a2*a4*(jv(nu+1,U)**2-jv(nu,U)*jv(nu+2,U))
         return a*(s1+s2)                 
      if(mode=='TE'):
         a=(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*beta/k
         s1=1.0
         s2=-jv(0,U)*jv(2,U)/jv(1,U)**2
         return a*(s1+s2)                          
      if(mode=='TM'):
         a=(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*k*self.nco**2/beta
         s1=1.0
         s2=-jv(0,U)*jv(2,U)/jv(1,U)**2         
         return a*(s1+s2)      
      pass

  def Ncl(self,U,V,nu,mode):  
      '''
      Snyder&Love: Section 12-10 Table 12-5 
      (see also Table 11-1 Properties of bound modes)
      Ncl=(P/|a|**2) -> Normalized modal power in the core
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      delta=self.delta  
      a1=self.a1(nu,U,V,W)
      a2=self.a2(nu,U,V,W)      
      a5=self.a5(nu,U,V,W)
      a6=self.a6(nu,U,V,W)
      if(mode=='hybrid'):
         a=-(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*k*self.nco**2/(beta*kv(nu,W)**2)
         a=a*(U/W)**2
         s1=a1*a5*(kv(nu-1,W)**2-kv(nu,W)*kv(nu-2,W))
         s2=a2*a6*(kv(nu+1,W)**2-kv(nu,W)*kv(nu+2,W))
         return a*(s1+s2)                 
      if(mode=='TE'):
         a=-(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*beta/k
         s1=1.0
         s2=-kv(0,W)*kv(2,W)/kv(1,W)**2
         return a*(s1+s2)                    
      if(mode=='TM'):   
         a=-(np.pi*rho**2/2.0)*np.sqrt(eps0/mu0)*k*self.nco**2/(beta*(1.0-2.0*delta))
         s1=1.0
         s2=-kv(0,W)*kv(2,W)/kv(1,W)**2
         return a*(s1+s2)      
      pass

  def Sz(self,U,V,nu,mode):
      '''
      Snyder&Love: Section 12-10 Table 12-5 
      (see also Table 11-1 Properties of bound modes)
      Sz/|a|^2
      '''
      k    =self.k_from_V(V)    
      beta =self.beta_from_UV(k,U,V)
      W=np.sqrt(V**2-U**2)  
      rho=self.rho
      delta=self.delta  
      a1=self.a1(nu,U,V,W)
      a2=self.a2(nu,U,V,W)      
      a3=self.a3(nu,U,V,W)
      a4=self.a4(nu,U,V,W)
      a5=self.a5(nu,U,V,W)
      a6=self.a6(nu,U,V,W)      
      F1=self.F1(nu,U,V,W)
      F2=self.F2(nu,U,V,W)
      def sz(R,phi):            
        if(np.max(R)<=1.0):
         a=0.5*np.sqrt(eps0/mu0)*k*self.nco**2/(beta*jv(nu,U)**2)
         s1=a1*a3*jv(nu-1,U*R)**2
         s2=a2*a4*jv(nu+1,U*R)**2         
         s3=0.5*(1.0-F1*F2)*jv(nu-1,U*R)*jv(nu+1,U*R)*np.cos(2.0*nu*phi) #even
         if(nu%2 != 0): s3=-s3 #odd
         return a*(s1+s2+s3)   
        else:
         a=0.5*np.sqrt(eps0/mu0)*(k*self.nco**2)/(beta*kv(nu,W)**2)
         a=a*(U/W)**2
         s1=a1*a5*kv(nu-1,W*R)**2
         s2=a2*a6*kv(nu+1,W*R)**2         
         s3=0.5*(1.0-2.0*delta-F1*F2)*kv(nu-1,W*R)*kv(nu+1,W*R)*np.cos(2.0*nu*phi) #even
         if(nu%2 != 0): s3=-s3 #odd
         return a*(s1+s2+s3)      
      def sz_TE(R,phi):   
        a=0.5*np.sqrt(eps0/mu0)*beta/k
        if(np.max(R)<=1.0):         
         s1=jv(1,U*R)**2/jv(1,U)**2         
        else:         
         s1=kv(1,W*R)**2/kv(1,W)**2
        return a*s1
      def sz_TM(R,phi):    
        a=0.5*np.sqrt(eps0/mu0)*k*self.nco**2/beta                          
        if(np.max(R)<=1.0):          
          s1=jv(1,U*R)**2/jv(1,U)**2        
        else:          
          a=a/(1.0-2.0*delta)
          s1=kv(1,W*R)**2/kv(1,W)**2
        return a*s1

      func={'hybrid':sz,
            'TE':sz_TE,
            'TM':sz_TM}
      return func[mode]