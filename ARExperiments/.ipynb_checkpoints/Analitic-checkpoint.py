#!/usr/bin/env python
# coding: utf-8

# # Analisi analitic

# In[31]:


import numpy
from   thermo             import Mixture
from   scipy.constants    import atmosphere, zero_Celsius, mmHg, g
import scipy           
import math
from   skaero.gasdynamics import isentropic, shocks
from   sympy.solvers      import solve
from   sympy              import Symbol
import gasTables
isen = Isentropic() 


# In[34]:


def pressure_from_mach_p0 ( M , p0 , gamma ):
    x = Symbol('x')
    
    a = gamma - 1
    b = gamma / a
    
    pressure_from_mach_p0 = solve( - p0 + x * (1 + a * M*M /2 ) ** b , x )
    just = (1 + a * M*M /2 ) ** - b
    return pressure_from_mach_p0, just


# In[35]:





# In[22]:


airN           = Mixture( 'air' , T = 273.15 , P = atmosphere )
rhoN           = airN.rho
# Variables geometriques

d_chicle                  = 2.75e-3           #m
d_sortidaChicle           = 7e-3              #m

# Parametres de fluids
p_0                       = 6e5               #pa
temperatura_entrada       = 25                #ºC
gamma                     = 1.40

# Parametres calculables

area_coll          = area(d_chicle)
area_sortidaChicle = area (d_sortidaChicle)
t_entradaConsum    = temperatura_entrada + 273.15
p_alimentacio      = atmosphere + p_0


cabal = chocked_mass_flow(p_alimentacio,t_entradaConsum,area_coll,gamma)
print ("cabal màssic en Nl/min:", cabal*60*1000/rhoN)


# In[23]:


p , just = pressure_from_mach_p0 ( 1 , p_alimentacio , gamma )

print ( p , just )


# In[24]:


# psortida = p_alimentacio* just         # aquest en teoria hauria de clavar el mach 1

psortida = p_alimentacio * 0.1
mach     = mach_from_p0_p ( p_alimentacio , psortida , gamma )

print ( psortida , mach)


# In[25]:


rel_areafinal = area_from_astar_mach ( mach , gamma)
areafinal = area_coll * rel_areafinal

diam = diametre (areafinal)

print ( diam*1e3)


# In[30]:


machh = isen.get_M_from_A_by_Astar (mach)
mstaaaar = isen.get_Mstar_from_M(2)
print (machh,mstaaaar)


# In[ ]:




