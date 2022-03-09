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

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from math import pi

def plot_1d(xlabel,ylabel,*args): 
    no_axis=len(args)         
    fig, ax = plt.subplots(1,no_axis,figsize=(12,6))
    for i,graph in enumerate(args):            
        colors=graph[1]
        ax[i].scatter(graph[0],graph[1],c=colors, cmap='seismic', alpha=0.75)    
        ax[i].plot(graph[0],graph[1])
        ax[i].set_xlabel(xlabel)
        ax[i].set_ylabel(ylabel[i])
        #ax[i].set_ylim(0.0,1.0)   

def plot_field_polar(field_func,Rmax=10.0):   
    # Using linspace so that the endpoint of 360 is included
    phi_core_l = np.radians(np.linspace(0, 360, 100))
    R_core_l = np.linspace(0, 1, 100)
    phi_cladding_l = np.radians(np.linspace(0, 360, 100))
    R_cladding_l = np.linspace(1, Rmax, 100)
    # Creates a mesh (R,phi)
    R_core, phi_core = np.meshgrid(R_core_l, phi_core_l)
    R_cladding, phi_cladding = np.meshgrid(R_cladding_l, phi_cladding_l)    
    
    #Field(R,phi)    
    core = field_func(R_core,phi_core)
    cladding= field_func(R_cladding,phi_cladding)
    vmaxcore=np.max(np.abs(core))
    vmaxcladding=np.max(np.abs(cladding))
    vmax=max(vmaxcore,vmaxcladding)
    
    R=np.hstack((R_core,R_cladding))    
    phi=np.hstack((phi_core,phi_cladding))    
    field=np.hstack((core,cladding))
    
    #Plot
    print("max=",vmax)
    fig, ax = plt.subplots(figsize=(3,3),dpi=100,
                           subplot_kw=dict(projection='polar'))
    img=ax.contourf(phi, R,field,cmap='bwr',vmin=-vmax,vmax=vmax)        
    ax.set_theta_offset(pi/2.0)
    ax.set_title(field_func.__name__,loc= 'left')          
    plt.colorbar(img)
    plt.show() 
    pass

def plot_field_cartesian(field_func,Rmax=10.0):   
    # Using linspace so that the endpoint of 360 is included
    phi_core_l = np.radians(np.linspace(0, 360, 100))
    R_core_l = np.linspace(0, 1, 100)
    phi_cladding_l = np.radians(np.linspace(0, 360, 100))
    R_cladding_l = np.linspace(1, Rmax, 100)
    # Creates a mesh (R,phi)
    R_core, phi_core = np.meshgrid(R_core_l, phi_core_l)
    R_cladding, phi_cladding = np.meshgrid(R_cladding_l, phi_cladding_l)    
    
    #Field(R,phi)    
    core = field_func(R_core,phi_core)
    cladding= field_func(R_cladding,phi_cladding)
    vmaxcore=np.max(core)
    vmaxcladding=np.max(cladding)
    vmax=max(vmaxcore,vmaxcladding)
    
    R=np.hstack((R_core,R_cladding))    
    phi=np.hstack((phi_core,phi_cladding))    
    field=np.hstack((core,cladding))
    
    X=R*np.cos(phi)
    Y=R*np.sin(phi)
    #Plot
    fig, ax = plt.subplots(figsize=(5,5))
    img=ax.contourf(X, Y,field,cmap='bwr',vmin=-np.max(field),vmax=np.max(field))        
    ax.set_title('$|F|^2$ (a.u.)')          
    plt.colorbar(img)
    plt.show() 
    pass
'''
https://matplotlib.org/stable/gallery/color/colormap_reference.html
'''

def step_array_to_regular(x):
    x=np.array(x)        
    no=len(x)  
    print('Number of V-values=',no)  
    len_x=[len(xi) for xi in x]      
    mx_len_x=max(len_x)
    print("Max len of U-values=",mx_len_x)        
    y=np.zeros(shape=(no,mx_len_x))           
    for i in range(1,no):      
      x[i].resize(mx_len_x,refcheck=False)      
      y[i]=y[i]+x[i]      
    print(y.shape)
    #print(y)          
    return np.transpose(y)

def plot_dispersion_relation(x,y,col_names,light_cone_x,light_cone_y,
                             column_name=None,labelsxy=['x','y'],cutoff=[],
                             xmin=None,xmax=None,ymin=0.01,ymax=None,
                             invert_yaxis=False,linewidth=0.0,markersize=4):    
                             
    plt.rcParams["figure.figsize"] = (5,5)        
    plt.figure(dpi=120)
    # Plot optical fiber dispersion relation
    if(x.ndim<=y.ndim):      
      if(column_name==None):
        for name in col_names:                                        
            plt.plot(x,y[name],label=name,
                     linewidth=linewidth,marker='o',markersize=markersize) 
      else:          
          plt.plot(x,y[column_name].values,label=column_name,
                     linewidth=linewidth,marker='o',markersize=markersize) 
    else:
      if(column_name==None):
        for name in col_names:                 
            plt.plot(x[name],y,label=name,
                     linewidth=linewidth,marker='o',markersize=markersize) 
      else:
          plt.plot(x[column_name].values,y,label=column_name,
                     linewidth=linewidth,marker='o',markersize=markersize) 
    # Plot the light cone  
    for _ in light_cone_y:
        plt.plot(light_cone_x,_,linewidth=4)
    # Plot the cutoff (if provided)
    if(cutoff!=[]):
        plt.scatter(cutoff,cutoff,s=80,alpha=0.8)

    plt.xlabel(labelsxy[0],fontsize=20)
    plt.ylabel(labelsxy[1],fontsize=20)
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if(invert_yaxis): plt.gca().invert_yaxis()
    plt.show()    
    pass


def EM_field_boundaries(field,Rmax):
    R_core_l = np.linspace(-1, 1, 100)
    R_cladding_l1=np.linspace(1, Rmax, 100)
    R_cladding_l2=np.linspace(-1,-Rmax, 100)
    from numpy import pi
    xcore=R_core_l
    xcladding1=R_cladding_l1
    xcladding2=R_cladding_l2
    ycore=np.abs(field(xcore,pi/4.0))
    ycladding=np.abs(field(xcladding1,pi/4.0))

    # Along a given direction
    fig, ax = plt.subplots(figsize=(4,4))
    ax.scatter(xcore,ycore)
    ax.scatter(xcladding1,ycladding)
    ax.scatter(xcladding2,ycladding)
    ax.set_xlabel(r'$R (=r/\rho)$')    
    plt.show()

    # At rho=R+-delta_R
    phi_l = np.radians(np.linspace(0, 360, 100))
    ycore=np.abs(field(0.999,phi_l))
    ycladding=np.abs(field(1.001,phi_l))

    if(ycore!=[]):
      fig, ax = plt.subplots(figsize=(4,4))
      ax.scatter(phi_l,ycore,alpha=0.5,marker='o')
      ax.scatter(phi_l,ycladding,alpha=0.5,marker='^')
      plt.show()
    pass

def plot_frac_energy(df,frac_power,mode_labels,
                     frac_power_fields=None,Vmax=5.0,
                     markersize=1,linewidth=1.0):
    fig, ax = plt.subplots(figsize=(8,8)) 
    for i,power in enumerate(frac_power):
      ax.plot(df['V'],power,label=mode_labels[i],
              linewidth=linewidth,marker='o',markersize=markersize)
    
    try:
        for value in frac_power_fields:
          ax.scatter(value[0],value[1])
    except:
        print("No individual values of power fraction given.")            
    
    ax.set_xlabel('V')
    ax.set_ylabel('$\eta$')
    ax.set_xlim(0.0,Vmax)
    ax.set_ylim(0.0,1.0)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show() 
    pass

def plot_vg(of,df,vg,mode_labels,
            vg_fields=None,Vmax=5.0,
            markersize=1,linewidth=1.0):
    c0= 299792458.0 #m/s
    fig, ax = plt.subplots(figsize=(8,8))
    for i,v in enumerate(vg):
      ax.plot(df['V'],v/c0,label=mode_labels[i],
              linewidth=linewidth,marker='o',markersize=markersize)   
    try:
        for value in vg_fields:
          ax.scatter(value[0],value[1]/c0)
    except:
        print("No individual values of power fraction given.")       
        
    ax.hlines(1.0/of.nco,0.0,Vmax,colors='r',linestyles='--', label='$v_{co}/c$')
    ax.hlines(1.0/of.ncl,0.0,Vmax,linestyles='--', label='$v_{cl}/c$')
    ax.set_xlabel('V')
    ax.set_ylabel('$v_g/c$')
    ax.set_xlim(0.0,Vmax)
    ax.set_ylim(0.0,1.0)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show() 
    pass

def plotpolar_mode_fields(field,Rmax,E=True,H=False,Sz=True):
    if(E):
        plot_field_polar(field.ex(),Rmax)
        plot_field_polar(field.ey(),Rmax)
        plot_field_polar(field.ez(),Rmax)
        plot_field_polar(field.er(),Rmax)
        plot_field_polar(field.ephi(),Rmax)
        
    if(H):
        plot_field_polar(field.hx(),Rmax)
        plot_field_polar(field.hy(),Rmax)              
        plot_field_polar(field.hz(),Rmax)  
        plot_field_polar(field.hr(),Rmax)
        plot_field_polar(field.hphi(),Rmax)
    
    if(Sz):
        plot_field_polar(field.Sz_numeric(),Rmax)    
        plot_field_polar(field.Sz_analytic(),Rmax)    
    
    pass