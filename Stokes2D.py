#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:22:13 2020

@author: abauville
"""

import numpy as np
from numpy import array as arr
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import time
plt.figure(1)
plt.clf()
## ========================================================================
#                           Numerical parameters
# =========================================================================
AX = arr([-1, 1, -1, 1]);
xmin = AX[0];
xmax = AX[1];
ymin = AX[2];
ymax = AX[3];
nx = 300;
ny = 300;
nNodes = nx*ny;

Kappa = 1.0;

# Reference grid
# =========================================================================
XX,YY = np.meshgrid(np.linspace(xmin,xmax,nx) , np.linspace(ymin,ymax,ny))
XX = XX.T
YY = YY.T
dx = (XX[2,1]-XX[1,1])/2.0;
dy = (YY[1,2]-YY[1,1])/2.0;

# staggered grid x
# =========================================================================
[XXsx,YYsx] = np.meshgrid(np.linspace(xmin,xmax,nx) , np.linspace(ymin-dy,ymax+dy,ny+1));
XXsx = XXsx.T
YYsx = YYsx.T
# staggered grid y
# =========================================================================
[XXsy,YYsy] = np.meshgrid(np.linspace(xmin-dx,xmax+dx,nx+1) , np.linspace(ymin,ymax,ny));
XXsy = XXsy.T
YYsy = YYsy.T
# staggered grid P
# =========================================================================
[XXsP,YYsP] = np.meshgrid(np.linspace(xmin+dx,xmax-dx,nx-1) , np.linspace(ymin+dy,ymax-dy,ny-1));
XXsP = XXsP.T
YYsP = YYsP.T

nVx = XXsx.shape[0]*XXsx.shape[1]
nVy = XXsy.shape[0]*XXsy.shape[1]
nP = XXsP.shape[0]*XXsP.shape[1]

Num_Vx = np.arange(0,nVx)#1:nVx;
Num_Vy = Num_Vx[-1]+1+np.arange(0,nVy)# Num_Vx(end)+[1:nVy];
Num_P = Num_Vy[-1]+1+np.arange(0,nP)#Num_Vy(end) +[1:nP];
no_eq_ini = Num_P[-1]+1;

## ========================================================================
#                           Physical parameters
# =========================================================================
Eta_n = np.ones(nP)             # viscosity for normal stress nodes, defined on the P grid
Eta_s = np.ones(nNodes)         # visc. for shear stress nodes, defined on the reference grid
rho = 1.0*np.ones(nP)
g = [0.0, 0.0]

## ========================================================================
#                   Debugging plt.plot - numbering etc...
# =========================================================================
if False:
    plt.plot(XX,YY,'-k',c=[.5,.5,.5])
    plt.plot(XX.T,YY.T,'-k',c=[.5,.5,.5])
    plt.plot(XX,YY,'sr',markersize=11,markeredgewidth=0,markerfacecolor=[.8,.8,.0])
    plt.plot(XXsx,YYsx,'>g',markersize=11,markeredgewidth=0.0,markerfacecolor=[0,.8,0])
    plt.plot(XXsy,YYsy,'^r',markersize=11,markeredgewidth=0.0,markerfacecolor=[.8,0,0])
    plt.plot(XXsP,YYsP,'ob',markersize=11,markeredgewidth=0.0,markerfacecolor=[0,0,.8])

    # # text(XXsx(:),YYsx(:),num2str([1:nVx]'),'Color',[0 0 0],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # # text(XXsy(:),YYsy(:),num2str([1:nVy]'),'Color',[1 1 1],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # # text(XXsP(:),YYsP(:),num2str([1:nP ]'),'Color',[1 1 1],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # # text(XX(:),YY(:),num2str([1:nNodes ]'),'Color',[0 0 0],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)

    # text(XXsx(:),YYsx(:),num2str(Num_Vx([1:nVx])'),'Color',[0 0 0],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # text(XXsy(:),YYsy(:),num2str(Num_Vy([1:nVy])'),'Color',[1 1 1],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # text(XXsP(:),YYsP(:),num2str(Num_P( [1:nP ])'),'Color',[1 1 1],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)
    # # text(XX(:),YY(:),num2str([1:nNodes ]'),'Color',[0 0 0],'Fontweight','bold','HorizontalAlignment','center','Fontsize',9)



    # axis(AX+2*[-dx dx -dy dy])
    # axis equal




# ## ========================================================================
# #                       Define boundary conditions
# # =========================================================================

# # Get boundary numbering
# # =========================================================================

# Vx boundaries
bound_Ind_VxBot  = (Num_Vx[1:nx-1]);                                        # inner bottom
bound_Ind_VxTop  = (Num_Vx[nVx-1+np.arange(-(nx-2),-1+1)])                              # inner top
bound_Ind_VxLeft = (Num_Vx[0+nx*np.arange(0,ny+1)]);                                   # left
bound_Ind_VxRight= (Num_Vx[nx-1+nx*np.arange(0,ny+1)]);

# Vy boundaries
bound_Ind_VyBot = (Num_Vy[0:nx+1]);                                         # bottom
bound_Ind_VyTop = (Num_Vy[nVy-1+np.arange(-nx,0+1)]);                                    # top
bound_Ind_VyLeft = (Num_Vy[0+(nx+1)*np.arange(1,ny-2+1)]);                             # inner left
bound_Ind_VyRight = (Num_Vy[nx+(nx+1)*np.arange(1,ny-2+1)]);                         # inner right

# P boundaries
# bound_Ind_P = (Num_P[arr([0,nx-2,nP-1,nP-1-(nx-1)+1])]);
bound_Ind_P = (Num_P[arr([nP-1-(nx-1)+1])]);
Xf = XXsP.T.flatten()
Yf = YYsP.T.flatten()
plt.plot(Xf[bound_Ind_P-Num_Vy[-1]-1],Yf[bound_Ind_P-Num_Vy[-1]-1],'.w')

# Set boundary values and type
# =========================================================================
# Define BC types
BC_free = -1
BC_dirichlet = 0
BC_dirichlet_ghost = 1
BC_numann_ghost = 2
BC_neumann_ghost = 3

# Vx values
bound_Val_VxBot      =  0.0*np.ones(bound_Ind_VxBot.shape);
bound_Val_VxTop      =  0.0*np.ones(bound_Ind_VxTop.shape);
bound_Val_VxLeft     =  1.0*np.ones(bound_Ind_VxLeft.shape);
bound_Val_VxRight    = -1.0*np.ones(bound_Ind_VxRight.shape);
# Vx types
bound_Type_VxBot     = BC_dirichlet * np.ones(bound_Ind_VxBot.shape,dtype=int);
bound_Type_VxTop     = BC_dirichlet * np.ones(bound_Ind_VxTop.shape,dtype=int);
bound_Type_VxLeft    = BC_dirichlet * np.ones(bound_Ind_VxLeft.shape,dtype=int);
bound_Type_VxRight   = BC_dirichlet * np.ones(bound_Ind_VxRight.shape,dtype=int);


# Vy values
bound_Val_VyBot      = -1.0*np.ones(bound_Ind_VyBot.shape);
bound_Val_VyTop      =  1.0*np.ones(bound_Ind_VyTop.shape);
bound_Val_VyLeft     =  0.0*np.ones(bound_Ind_VyLeft.shape);
bound_Val_VyRight    =  0.0*np.ones(bound_Ind_VyRight.shape);
# Vy types
bound_Type_VyBot     = BC_dirichlet * np.ones(bound_Ind_VyBot.shape,dtype=int);
bound_Type_VyTop     = BC_dirichlet * np.ones(bound_Ind_VyTop.shape,dtype=int);
bound_Type_VyLeft    = BC_dirichlet * np.ones(bound_Ind_VyLeft.shape,dtype=int);
bound_Type_VyRight   = BC_dirichlet * np.ones(bound_Ind_VyRight.shape,dtype=int);

# P values
bound_Val_P          = np.zeros(bound_Ind_P.shape);
# P types
bound_Type_P         = BC_dirichlet * np.ones(bound_Ind_P.shape,dtype=int);

# Fill the IndAll and ValAll boundary vectors
# =========================================================================
bound_Ind_All = np.zeros(no_eq_ini,dtype=bool);
bound_Val_All = np.zeros(no_eq_ini);
bound_Type_All = BC_free * np.ones(no_eq_ini,dtype=int);

bound_Ind_list = [bound_Ind_VxBot, bound_Ind_VxTop, bound_Ind_VxLeft, bound_Ind_VxRight,
                  bound_Ind_VyBot, bound_Ind_VyTop, bound_Ind_VyLeft, bound_Ind_VyRight,
                  bound_Ind_P]
bound_Val_list = [bound_Val_VxBot, bound_Val_VxTop, bound_Val_VxLeft, bound_Val_VxRight,
                  bound_Val_VyBot, bound_Val_VyTop, bound_Val_VyLeft, bound_Val_VyRight,
                  bound_Val_P]
bound_Type_list = [bound_Type_VxBot, bound_Type_VxTop, bound_Type_VxLeft, bound_Type_VxRight,
                   bound_Type_VyBot, bound_Type_VyTop, bound_Type_VyLeft, bound_Type_VyRight,
                   bound_Type_P]
no_eq = no_eq_ini

#%%
Num_nodir_All = np.ones(no_eq_ini,dtype=int);
for iF in range(len(bound_Ind_list)):
    Ind = bound_Ind_list[iF]; #Ind = bound_Ind_(Fields{iF});
    bound_Ind_All[Ind]  = True;
    bound_Val_All[Ind]  = bound_Val_list[iF];
    bound_Type_All[Ind] = bound_Type_list[iF];
    Num_nodir_All[Ind] = np.logical_not(np.logical_or(bound_Type_list[iF] == 0, bound_Type_list[iF] == 1));
    no_eq -= len(Ind)

Num_nodir_All = np.cumsum(Num_nodir_All)-1
Num_nodir_All[np.logical_or(bound_Type_All == 0, bound_Type_All == 1)] = -1
# print(Num_nodir_All)
#%%

# Get number of non-zeros
# =========================================================================
bound_nVx = len(bound_Ind_VxBot) + len(bound_Ind_VxTop) + len(bound_Ind_VxLeft ) + len(bound_Ind_VxRight); # number of Vx Dirichlets equations
bound_nVy = len(bound_Ind_VyBot) + len(bound_Ind_VyTop) + len(bound_Ind_VyRight) + len(bound_Ind_VyRight);
bound_nP  = len(bound_Ind_P);
no_nz = (nVx-bound_nVx)*11 + (nVy-bound_nVy)*11 + (nP-bound_nP)*4;
# no_nz += -bound_nVx - bound_nVy - bound_nP
# no_nz += 2*(-bound_nVx-4 - bound_nVy-4) - 2*bound_nP + 2

## Create numbering without dirichlet
Num_Vx_nodir = np.ones(Num_Vx.shape,dtype=int)
Num_Vx_nodir[np.logical_or(bound_Type_All[:nVx] == 0,
                           bound_Type_All[:nVx] == 1)] = 0
Num_Vx_nodir = np.cumsum(Num_Vx_nodir) - 1
Num_Vx_nodir[np.logical_or(bound_Type_All[:nVx] == 0,
                           bound_Type_All[:nVx] == 1)] = -1

Num_Vy_nodir = np.ones(Num_Vy.shape,dtype=int)
Num_Vy_nodir[np.logical_or(bound_Type_All[nVx:nVx+nVy] == 0,
                           bound_Type_All[nVx:nVx+nVy] == 1)] = 0
Num_Vy_nodir = np.cumsum(Num_Vy_nodir) - 1
Num_Vy_nodir[np.logical_or(bound_Type_All[nVx:nVx+nVy] == 0,
                           bound_Type_All[nVx:nVx+nVy] == 1)] = -1

Num_P_nodir = np.ones(Num_P.shape,dtype=int)
Num_P_nodir[np.logical_or(bound_Type_All[nVx+nVy:nVx+nVy+nP] == 0,
                          bound_Type_All[nVx+nVy:nVx+nVy+nP] == 1)] = 0
Num_P_nodir = np.cumsum(Num_P_nodir) - 1
Num_P_nodir[np.logical_or(bound_Type_All[nVx+nVy:nVx+nVy+nP] == 0,
                           bound_Type_All[nVx+nVy:nVx+nVy+nP] == 1)] = -1

# 1736
# no_nz = no_nz + bound_nVx  + bound_nVy + bound_nP; # add diagonal terms

## ========================================================================
#                           Allocate memory
# =========================================================================

# Allocate sparse triplet
# =========================================================================
# Isp = np.zeros(no_nz);
Isp = np.zeros(no_eq+1,dtype=int); Isp[-1] = no_nz
Jsp = np.zeros(no_nz,dtype=int);
Vsp = np.zeros(no_nz);
b   = np.zeros(no_eq);


# Initialize structures for local numbering and coeffs
# =========================================================================
nzC = 0; #counter of non-zeros
eqC = 0; # counter of equations
class Num:
    def __init__(self):
        self.C = 0
        self.N = 0
        self.S = 0
        self.E = 0
        self.W = 0
        self.NE = 0
        self.NW = 0
        self.SE = 0
        self.SW = 0

Vxl     = Num()
Vyl     = Num()
Pl      = Num()
Eta_l   = Num()

Vyl_Coeffs  = Num()
Vxl_Coeffs  = Num()
Pl_Coeffs   = Num()
Eta_l.Val   = Num()



## ========================================================================
#                           Fill triplets
# =========================================================================
tic = time.time()
# =========================================================================
#                           Fill Vx equations
# =========================================================================
no_eq_l   = 11; # number of local equations
loc_J = np.zeros(no_eq_l,dtype=int);
loc_V = np.zeros(no_eq_l);
for i in range(nVx):
    # if  not bound_Ind_All[Num_Vx[i]]: # Free nodes
    if  bound_Type_All[Num_Vx[i]] == BC_free: # Free nodes
        ## ================================================================
        #                           Get local infos
        # =================================================================
        shift = int(np.floor((i)/nx)); # used to index Vy taking into account that a row vy is +1 long than a vx line


        # Get numbering
        # =================================================================
        Vxl.C= i;             Vxl.N= i+nx;      Vxl.S= i-nx;         Vxl.E= i+1;            Vxl.W= i-1; # numbering of the local Vx equations (C,N,S,E,W)
        Vyl.NE= shift+i+1;    Vyl.NW= shift+i;  Vyl.SE= shift+i-nx;  Vyl.SW= shift+i-nx-1;
        Pl.E= -nx-shift+i+1;  Pl.W= -nx-shift+i;
        Eta_l.N= i;           Eta_l.S= i-nx;     Eta_l.E= Pl.E;      Eta_l.W= Pl.W;


        # Get values of viscosity
        # =================================================================
        etaN = Eta_s[Eta_l.N];        etaS = Eta_s[Eta_l.S];
        etaE = Eta_n[Eta_l.E];        etaW = Eta_n[Eta_l.W];


        # Get Coeffs
        # =================================================================
        Vxl_Coeffs.C= (-2.0*etaE/(dx**2)  -2.0*etaW/(dx**2)  -etaN/(dy**2)  -etaS/(dy**2));
        Vxl_Coeffs.N= etaN/(dy**2);       Vxl_Coeffs.S= etaS/(dy**2);
        Vxl_Coeffs.E= 2.0*etaE/(dx**2);     Vxl_Coeffs.W= 2.0*etaW/(dx**2);

        Vyl_Coeffs.NE=  etaN/(dy*dx);     Vyl_Coeffs.NW= -etaN/(dy*dx);
        Vyl_Coeffs.SE= -etaS/(dy*dx);     Vyl_Coeffs.SW=  etaS/(dy*dx);

        Pl_Coeffs.E= -Kappa/dx;           Pl_Coeffs.W= Kappa/dx;


        # Fill local sparse triplets
        # =================================================================
        loc_J[0:5]      = Num_Vx[[Vxl.C,   Vxl.N,   Vxl.S,   Vxl.E,   Vxl.W]];
        loc_J[5:9]      = Num_Vy[[Vyl.NE,  Vyl.NW,  Vyl.SE,  Vyl.SW]];
        loc_J[9:11]     = Num_P[ [Pl.E,    Pl.W]];

        loc_V[0:5]      = [Vxl_Coeffs.C,   Vxl_Coeffs.N,   Vxl_Coeffs.S,   Vxl_Coeffs.E,   Vxl_Coeffs.W];
        loc_V[5:9]      = [Vyl_Coeffs.NE,  Vyl_Coeffs.NW,  Vyl_Coeffs.SE,  Vyl_Coeffs.SW];
        loc_V[9:11]     = [Pl_Coeffs.E,    Pl_Coeffs.W];


        # Fill global sparse triplets
        # =================================================================
        I = bound_Type_All[loc_J] == BC_free
        l = len(loc_J[I])

        Jsp[nzC:nzC+l] = Num_nodir_All[loc_J[I]];
        Isp[eqC]  = nzC#Num_Vx[i]#*np.ones(no_eq_l);
        Vsp[nzC:nzC+l]  = loc_V[I];

        nzC = nzC + l;


        # Add gravity term to right-hand side
        # =================================================================
        b[eqC] -= 0.5*( rho[Pl.E] + rho[Pl.W] )*g[0]

        # Check and apply BC
        # =================================================================
        for loc_i in range(no_eq_l):
            j = loc_J[loc_i]
            if bound_Type_All[j] == BC_dirichlet:
                b -= loc_V[loc_i] * bound_Val_All[j]
        eqC += 1
#   end if
# end for



# =========================================================================
#                           Fill Vy equations
# =========================================================================
loc_J = np.zeros(no_eq_l,dtype=int);
loc_V = np.zeros(no_eq_l);
for i in range(nVy):
    if bound_Type_All[Num_Vy[i]] == BC_free: # Free nodes
        ## ================================================================
        #                           Get local infos
        # =================================================================
        shift = int(np.floor((i)/(nx+1))); # used to index Vy taking into account that a row vy is +1 long than a vx line

        # Get numbering
        # =================================================================
        Vyl.C= i;                    Vyl.N= i+(nx+1);         Vyl.S= i-(nx+1);    Vyl.E= i+1;          Vyl.W= i-1; # numbering of the local Vx equations (C,N,S,E,W)
        Vxl.NE= -shift+i+nx;         Vxl.NW= -shift+i+nx-1;   Vxl.SE= -shift+i;   Vxl.SW= -shift+i-1;
        Pl.N= -nx-2*shift+i+(nx-1);  Pl.S= -nx-2*shift+i;
        Eta_l.N= Pl.N;               Eta_l.S= Pl.S;           Eta_l.E= i-shift;   Eta_l.W= i-shift-1;

        # Get values of viscosity
        # =================================================================
        etaN = Eta_n[Eta_l.N];        etaS = Eta_n[Eta_l.S];
        etaE = Eta_s[Eta_l.E];        etaW = Eta_s[Eta_l.W];

        # Get Coeffs
        # =================================================================
        Vyl_Coeffs.C= (-2.0*etaN/(dy**2)  -2.0*etaS/(dy**2)  -etaE/(dx**2)  -etaW/(dx**2));
        Vyl_Coeffs.N= 2.0*etaN/(dy**2);   Vyl_Coeffs.S= 2.0*etaS/(dy**2);
        Vyl_Coeffs.E=   etaE/(dx**2);     Vyl_Coeffs.W=   etaW/(dx**2);

        Vxl_Coeffs.NE=  etaE/(dy*dx);     Vxl_Coeffs.NW= -etaW/(dy*dx);
        Vxl_Coeffs.SE= -etaE/(dy*dx);     Vxl_Coeffs.SW=  etaW/(dy*dx);

        Pl_Coeffs.N= -Kappa/dy;           Pl_Coeffs.S= Kappa/dy;


        # Fill local sparse triplets
        # =================================================================
        loc_J[0:5]      = Num_Vy[[Vyl.C,   Vyl.N,    Vyl.S,    Vyl.E,   Vyl.W]];
        loc_J[5:9]      = Num_Vx[[Vxl.NE,  Vxl.NW,   Vxl.SE,   Vxl.SW]];
        loc_J[9:11]     = Num_P[[Pl.N,    Pl.S]];

        loc_V[0:5]      = [Vyl_Coeffs.C,   Vyl_Coeffs.N,   Vyl_Coeffs.S,   Vyl_Coeffs.E,   Vyl_Coeffs.W];
        loc_V[5:9]      = [Vxl_Coeffs.NE,  Vxl_Coeffs.NW,  Vxl_Coeffs.SE,  Vxl_Coeffs.SW];
        loc_V[9:11]     = [Pl_Coeffs.N,    Pl_Coeffs.S];


        # Fill global sparse triplets
        # =================================================================
        I = bound_Type_All[loc_J] == BC_free #
        l = len(loc_J[I])
        Jsp[nzC:nzC+l] = Num_nodir_All[loc_J[I]];
        Isp[eqC]  = nzC#Num_Vx[i]#*np.ones(no_eq_l);
        Vsp[nzC:nzC+l]  = loc_V[I];

        nzC = nzC + l;


        # Add gravity term to right-hand side
        # =================================================================
        b[eqC] -= 0.5*( rho[Pl.N] + rho[Pl.S] )*g[1]


        # Check and apply BC
        # =================================================================
        for loc_i in range(no_eq_l):
            j = loc_J[loc_i]
            if bound_Type_All[j] == BC_dirichlet:
                b -= loc_V[loc_i] * bound_Val_All[j]
    # elif bound_Type_All[Num_Vy[i]] == BC_dirichlet: # Boundary nodes
    #     Jsp[nzC] = Num_Vy[i];
    #     Isp[eqC] = nzC#Num_Vy[i];
    #     Vsp[nzC] =     1;
    #     nzC = nzC+1;
        eqC += 1
#   end if
# end for



# =========================================================================
#                           Fill P equations
# =========================================================================
no_eq_l = 4
loc_J = np.zeros(no_eq_l,dtype=int);
loc_V = np.zeros(no_eq_l);
for i in range(nP):
    if bound_Type_All[Num_P[i]] == BC_free: # Free nodes
        ## ================================================================
        #                           Get local infos
        # =================================================================
        shift = int(np.floor((i)/(nx-1))); # used to index Vy taking into account that a row vy is +1 long than a vx line
        # Get numbering
        # =================================================================
        Vxl.E= nx+shift+i+1;      Vxl.W= nx+shift+i;
        Vyl.N= nx+1+2*shift+i+1;  Vyl.S= nx+1+2*shift+i-nx;

        # Get Coeffs
        # =================================================================
        Vxl_Coeffs.E = 1.0/dx;    Vxl_Coeffs.W = -1.0/dx;
        Vyl_Coeffs.N = 1.0/dy;    Vyl_Coeffs.S = -1.0/dy;


        # Fill local sparse triplets
        # =================================================================
        loc_J[:2]  = Num_Vx[[Vxl.E,  Vxl.W]];
        loc_J[2:4] = Num_Vy[[Vyl.N,  Vyl.S]];

        loc_V[:] = [Vxl_Coeffs.E,   Vxl_Coeffs.W,   Vyl_Coeffs.N,   Vyl_Coeffs.S];


        # Fill global sparse triplets
        # =================================================================
        I = bound_Type_All[loc_J] == BC_free #
        l = len(loc_J[I])
        Jsp[nzC:nzC+l] = Num_nodir_All[loc_J[I]];
        Isp[eqC]  = nzC#Num_Vx[i]#*np.ones(no_eq_l);
        Vsp[nzC:nzC+l]  = loc_V[I];

        nzC = nzC + l;

        # Check and apply BC
        # =================================================================
        for loc_i in range(no_eq_l):
            j = loc_J[loc_i]
            if bound_Type_All[j] == BC_dirichlet:
                b -= loc_V[loc_i] * bound_Val_All[j]
        eqC += 1
  # end if
# end for

# Update no_nz and trim the end of Vsp etc... (initially overallocated)
no_nz = nzC
Isp[-1] = no_nz
Jsp = Jsp[:no_nz]
Vsp = Vsp[:no_nz]
print('Triplet filling: %.1fs' % (time.time()-tic))

## ========================================================================
#                           Assemble stiffness
# =========================================================================
tic = time.time()

A = sp.csr_matrix((Vsp, Jsp, Isp), shape=(no_eq, no_eq))
print('Assembly: %.1fs\n' % (time.time()-tic))
# Isp = [];
# Jsp = [];
# V = [];
plt.clf()
plt.spy(A,markersize=3)

# ## ========================================================================
# #                       Apply boundary conditions
# # =========================================================================
# b[bound_Ind_All] += bound_Val_All[bound_Ind_All];


# ========================================================================
                      ######     SOLVE     ######
# =========================================================================
tic = time.time()
solve = linalg.factorized(A)
print('Factorization: %.1fs' % (time.time()-tic))
tic = time.time()
x = solve(b)
print('Solve: %.1fs' % (time.time()-tic))






# #%%  ========================================================================
# #                               Visualization
# # =========================================================================

# Vx = x[Num_Vx];
# Vy = x[Num_Vy];
# P  = x[Num_P];

# ## plt.plot
# plt.clf()
# Vxplot = np.reshape(Vx,XXsx.T.shape).T;
# Vyplot = np.reshape(Vy,XXsy.T.shape).T;

# Vxplot = 0.5*(Vxplot[:,:-1] + Vxplot[:,1:])
# Vyplot = 0.5*(Vyplot[:-1,:] + Vyplot[1:,:])



# # X = np.reshape(XXsx,nVx,1);
# # Y = np.reshape(YYsx,nVx,1);
# # Pplot = np.reshape(P,size(XXsP));
# # X = np.reshape(XXsP,nP,1);
# # Y = np.reshape(YYsP,nP,1);

# # plt.imshow(X,Y,Vxplt.plot')
# plt.imshow(Vxplot.T,origin='lower',extent=(xmin,xmax,ymin,ymax))
# r = 20
# plt.quiver(XX[::r,::r],YY[::r,::r],Vxplot[::r,::r],Vyplot[::r,::r])
# # surface(XXsx,YYsx,VxP)
# # shading interp
# plt.axis('equal')



# # pcolor(XXsx(1:end,:),YYsx(1:end,:),VxP(1:end,:))
# # shading interp
# # colorbar

# # caxis([-1 1])
# # end




