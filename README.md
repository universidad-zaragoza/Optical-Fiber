---
Material for the course "Photonics systems and devices"

Universidad de Zaragoza, Zaragoza (Spain)

The following notes written by SERGIO G RODRIGO (2021-2022)

sergut@unizar.es

---
---

# User notes
- Go to Github option Code-> Download ZIP.
- Save in your Google Drive*.
- Open Jupyter notebooks with Jupyter or Colab.
- See below the description of python and jupyter notebooks of this package.

*By creating a folder in path "/content/drive/MyDrive/optical_fiber" the jupyter notebooks will run in Colab.*

---
---

# Python files
  **Step-profile Optical Fiber class (OFiber_class.py):**

  The class contain the main modal properties of a step-profile optical fiber. We followed the notation from Snyder&Love book is followed ("Optical Waveguide Theory", A.W. Snyder and J. Love Springer, Boston, MA (1983) 1st Ed (https://doi.org/10.1007/978-1-4613-2813-1).
  
  It includes:
  - Eigenvalue equations and cuttoff expressions (Table 12-4).
  - Field components of the exact modes (Table 12-3).
  - Modal properties of the step-profile fiber (Table 12-5).


   **EM fields class (OFiber_fields.py):**

   Field components of all the modes (see Snyder&Love: Section 12-8 Table 12-3).
   Class EMfield inherits class OFiber as a means to have simple access to the 
   electromagnetic fields. The electromagnetic fields and the quantities necesary to define them are implented in OFiber class.      

   **Plotting and finding (OFiber_plot.py and OFiber_find.py):**

   Methods to find the solution of the Optical Fiber dispersion relation.
   Specific methods for plotting the results.  


---
---

# Jupyter notebooks
+ Two examples are provided: 
 - High index contrast optical fiber ($n_{co}$ >> $n_{cl}$) to compare with Fig. 12-4 Snyder & Love (see fig12-4_Snyder_Love.csv)
 - Weak guidance example ($n_{co}$ ≈ $n_{cl}$) to compare with Fig. 14-4 Snyder & Love (see fig14-4_Snyder_Love.csv)

**ofiber_calc_dispersion_relation.ipynb**
+ Solver of the eigenvalue equations for the step-profile fiber from Snyder&Love, Section 12-9 Table 12-4. 

**ofiber_plot_dispersion_relation.ipynb**

+ Plot dispersion relation either as: $U$ vs $V$, $\omega (eV)$ vs $\beta(\mu m^{-1})$ and $\lambda (nm)$ vs $\beta(\mu m^{-1})$, following the analytical expressions in *Snyder & Love, Section 12-8 Table 12-4*

**ofiber_em_fields.ipynb**

+ Plot the electromagnetic modes of a step-profile optical fiber following the analytical expressions in *Snyder & Love, Section 12-8 Table 12-3*. 
+ Plot the fraction of power in the core of a step-profile optical fiber. It takes analytical expressions applied to several modes. The power for the electromagnetic mode chosen for plotting is also calculated from the its electromagnetic fields.

**ofiber_modal_properties.ipynb**
+ Several modal properties of the step-profile fiber from analytical expressions in Snyder&Love, Section 12-10 Table 12-5. Namely, fraction of power in the core and group velocity can be obtained.

---
---

*License:*

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