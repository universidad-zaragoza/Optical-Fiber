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

import numpy as np
import scipy
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

class RootFinder:
    '''
    Simple zero/root finder for functions with one variable. 
    Finds all the zeroes within a range of values, using user-defined searchstep and tolerance.
    Utilizes the scipy.optimize.fsolve: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html
    '''
    def __init__(self, start, stop, step=0.01, root_dtype="float64", xtol=1e-9):

        self.start = start
        self.stop = stop
        self.step = step
        self.xtol = xtol
        self.roots = np.array([], dtype=root_dtype)

    def add_to_roots(self, x):
        if (x < self.start) or (x > self.stop):
            return  # outside range
        if any(abs(self.roots - x) < self.xtol):
            return  # root already found.
        self.roots = np.append(self.roots, x)
    
    def find_root(self, f, x0, *args):
        '''
        Uses scipy.optimize.fsolve
        '''
        x, _, ier, _ = fsolve(f, x0=x0, args=args, full_output=True, xtol=self.xtol)
        if ier == 1:
            return x[0]
        return None

    def find(self, f, *args):
        current = self.start
        for x0 in np.arange(self.start, self.stop + self.step, self.step):
            if x0 < current:
                continue
            x = self.find_root(f, x0, *args)
            if x is None:  # no root found.
                continue
            current = x
            self.add_to_roots(x)
        return self.roots
         
def find_EM_modes(func,args,Uini,Ufin,resolution,verbose=False):              
    r = RootFinder(Uini,Ufin,resolution)        
    roots = r.find(func,*args)    
    # plot results
    if(verbose==True):           
      print("U: ", roots)    
      u = np.linspace(Uini,Ufin,num=500)
      fu=func(u,*args)      
      fu_roots=func(np.array(roots),*args)   
      fig, ax = plt.subplots(figsize=(20,6))
      ax.plot(u, np.log(np.abs(fu)))
      ax.scatter(roots, np.log(np.abs(fu_roots)), color="r", s=10)
      ax.grid(color="grey", ls="--", lw=0.5)
      plt.show()
    return roots

def find_EM_modes_loop(func,args,Uini,steps,verbose=True):             
  U=[]
  Vloop=args['Vloop']
  tolerance=1e-3
  Ufin=max(Vloop)-tolerance  
  order=args['order']
  of=args['ofiber']
  resolution=(Ufin-Uini)/steps
  print(of)
  print("Order=",order)
  print(Ufin)
  for v in Vloop:
      args_local=(v,order)    
      u=find_EM_modes(func,args_local,Uini,Ufin,resolution)
      U.append(u)   
      if(verbose):         
        msg= f'V= {v}; U= {u}'
        print(msg)
  return U