---
Material for the course "Photonics systems and devices"

The following notes written by SERGIO G RODRIGO (2021-2022)</font> 

---
---

  **Step-profile Optical Fiber class (OFiber_class.py):**

  The object catch the main characteristics of a step-profile optical fiber. The notation is
  from Snyder&Love book is followed ("Optical Waveguide Theory", A.W. Snyder 
  and J. Love Springer, Boston, MA (1983) 1st Ed 
  https://doi.org/10.1007/978-1-4613-2813-1).)
  
  It includes:
  - Eigenvalue equations and cuttoff expressions (Table 12-4)
  - Field components of the exact modes (Table 12-3)
  - Modal properties of the step-profile fiber (Table 12-5)

---
---
   **EM fields class (OFiber_fields.py):**

   Snyder&Love: Section 12-8 Table 12-3 Field components of all the modes.
   Class EMfield inherits class OFiber as a means to have simple access to the 
   electromagnetic fields. The electromagnetic fields and the quantities necesary to define them are implented in OFiber class.      

---
---

   **OFiber_plot.py and OFiber_find.py:**

   Methods to find the solution of the Optical Fiber dispersion relation.
   Specific methods for plotting the results.  


---------------------------------------------------------------------------

License:
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