#!/usr/bin/env python
# coding: utf-8

# # Functions definition
# This functions calculates flow rate from a long radious nozzle. Which in our case it means: the secondary flow or the aspiration flow 

# In[1]:


from   fluids               import *
from   thermo               import Mixture
from   scipy.constants      import atmosphere, zero_Celsius, mmHg, g
import scipy                as     sci
import uncertainties 
from   uncertainties        import *
from   uncertainties.unumpy import *
import uncertainties        as     unumpy
from   numpy                import *
from   pandas               import *
import pandas               as     pd
import matplotlib.pyplot    as     plt
import pickle


# In[2]:


# For the calculation of flow rate with uncertainties:
wrappedFlowMeter = uncertainties.wrap( differential_pressure_meter_solver)
# Firsts definitions
airN               = Mixture('air',T=293,P=101325)
rhoN               = airN.rho


# In[3]:


#Here it is defined all our ambi 
def set_of_const  (P_amb , T_amb):
            
    P_amb     = P_amb * mmHg           # Pressure in Pa
    T_amb     = T_amb + zero_Celsius   # T in K 

# Here we calculate all what it is important with the air at the ambient pressure
    
    air       = Mixture( 'air', T = T_amb , P = P_amb)
    k         = air.isentropic_exponent
    mu        = air.mu
    rho       = air.rho
    
    return P_amb, T_amb, air, k , mu , rho


# In[4]:


# here we define the function that will calculate the secondary flow
def flow_meter_calculator (flow_meter_calculator, mesures1 , mesures2 , manRelRho, DeltaH , Do , D1 , rho , k , mu , P_amb ):  

# Experimental setup
    taps                            = 'flange'
    meter_type                      = 'long radius nozzle'
    
# Inci tractament de les dades
#     flow_meter_calculator           = DataFrame()
   
    flow_meter_calculator['h1']     = uarray(mesures1 , DeltaH *ones(len (mesures1)))
    flow_meter_calculator['h2']     = uarray(mesures2 , DeltaH *ones(len (mesures2)))
    
# ManRelRho is the relative density of the manometric fluid
    flow_meter_calculator['p1']     = flow_meter_calculator.h2 * g * manRelRho # p1 is the bigger pressure, and 
    flow_meter_calculator['p2']     = flow_meter_calculator.h1 * g * manRelRho # p2 is the low pressure in the nozzle
    flow_meter_calculator['DeltaP'] = flow_meter_calculator.p1 - flow_meter_calculator.p2

# First we have to initialize the vector
    mdot_s_1                        = uarray( zeros(len (mesures1)) , zeros(len (mesures1)) )

    for i,row in flow_meter_calculator.iterrows():
        
        if nominal_value(row.DeltaP) != 0.0:             # Just to avoid division by 0 when DeltaP=0. .
            P2                      = P_amb - row.DeltaP
            mdot_s_1[i]             = wrappedFlowMeter(D =D1, D2 =Do, P1 =P_amb, P2 =P_amb-row.DeltaP, rho=rho, mu=mu, k=k, meter_type=meter_type, taps=taps)
        else:
            mdot_s_1[i]             = ufloat(0.0,0.0)

    flow_meter_calculator['mdot_s'] = mdot_s_1
    flow_meter_calculator['QN_s']   = flow_meter_calculator.mdot_s * 60000/rhoN #flow_meter_calculator.mdot_s*60000/rhoN 
    
    return flow_meter_calculator


# In[5]:


# This function calculates the pressure in the vessel. It uses mercury.
def pressure_v  (pressure_vessel, mesures1, mesures2, DeltaH_Hg , P_amb):
    
    pressure_vessel['h3']           = uarray( mesures1 , DeltaH_Hg * ones(len (mesures1)) )
    pressure_vessel['h4']           = uarray( mesures2 , DeltaH_Hg * ones(len (mesures2)) )
    pressure_vessel['p3']           = pressure_vessel.h3 * mmHg  #mmHg that it is been taken from scipy.constants
    pressure_vessel['p4']           = pressure_vessel.h4 * mmHg
    pressure_vessel['P_vessel']     = pressure_vessel.p3 - pressure_vessel.p4
    pressure_vessel['P_vessel_dim'] = (P_amb - pressure_vessel.P_vessel) /P_amb
    pressure_vessel['P_vess_AR']    = 100 - pressure_vessel.P_vessel_dim*100           # Aquí ho transportem al percentatge normal, per AR


# In[6]:


#Here it is calculated the flow consumtion, or the primary flow
def mdot_f ( mdot_p , mesures1 , mesures2 , DeltaP , DeltaT , mesuresT , D1 , D2 ):

    # Experimental setup
    taps              = 'flange'
    meter_type        = 'ISO 5167 orifice'  # meter_type = 'ISA 1932 nozzle'

# First we have to initialize the vector   
    mdot              = uarray( zeros(len (mesures1)) , zeros(len (mesures1)) )  

# Here we inicialize the vectors with their uncertainties
    p1                = uarray( mesures1 ,   DeltaP * ones(len (mesures1)) )
    p2                = uarray( mesures2 ,   DeltaP * ones(len (mesures2)) )
    T1                = uarray( mesuresT ,   DeltaT * ones(len (mesures1)) )
# Here we start to carry the units as we need them     
    T1                = T1 + zero_Celsius    # Temperature in KELVIN
    p1                = ( p1 + 1 ) * 1e5     # Pressió en Pa, absoluts
    p2                = ( p2 + 1 ) * 1e5

    for i in range( len(mesures1) ):
# l'aire es calcula amb cada valor de P i de T. IMPORTANT: el nominal_value aconsegueix fer-hi passar només el valor correcte.
        air           = Mixture( 'air' , T = float(nominal_value (T1[i]) ) , P = float( nominal_value( p1[i]) ) )
        r             = sci.constants.R / air.MW * 1000
        k             = air.isentropic_exponent
        mu            = air.mu
        rho           = p1[i] / ( r * T1[i] )    
#  Here is where we write the vector mdot
        mdot[i]       = wrappedFlowMeter(D=D1, D2=D2, P1=p1[i], P2=p2[i], rho=rho, mu=mu, k=k, meter_type=meter_type, taps=taps)

    mdot_p['mdot_p']  = mdot
    mdot_p['QN_p']    = mdot_p.mdot_p * 60000 / rhoN
    mdot_p['mu']      = mdot_p.QN_s / mdot_p.QN_p


# In[7]:


# Here I define the function that we will use for the numerical results
def numDat    (mdot_s , mdot_p, p0_s , errNumSec , errNumPri ):

#Here we first have to inicialize the variables. In this case has to be multiplied by 360/5=72 
    mdot_s               = mdot_s*72
    mdot_p               = mdot_p*72
    mdot_f_s             = uarray( mdot_s , mdot_s * errNumSec)
    mdot_f_p             = uarray( mdot_p , mdot_p * errNumPri)

    # First we have to create the panda dataFrame
    numDat               = DataFrame()
    numDat['p0_s_dim']   = p0_s                     # These are the total pressure in the vessel
    numDat['mdot_s']     = mdot_f_s
    numDat['mdot_p']     = mdot_f_p  
    numDat['p0_s_AR']    = 100 - numDat.p0_s_dim*100
    
# We transform the calculations into Nl/min
    numDat['QN_s']       = numDat.mdot_s*60000/rhoN
    numDat['QN_p']       = numDat.mdot_p*60000/rhoN
# Finally we calculate the mu
    numDat['mu']         = numDat.QN_s / numDat.QN_p

    return numDat     


# In[8]:


def firstExperimentLabson (N, P_amb, T_amb ,manRelRho, DeltaH, DeltaH_Hg, DeltaP, DeltaT , D_s, D_p, Do_s, Do_p ):
    nom = 'experiment'
    Data = DataFrame()
    Data = nom+'{}'.format(N)
    Data = pd.read_csv ( '{}{}.csv'.format(nom,N) , sep = ',' )
    
    P_amb, T_amb, air , k , mu , rho = set_of_const (P_amb , T_amb)
    
    flow_meter_calculator (Data, Data.mesures1_s , Data.mesures2_s , manRelRho, DeltaH , Do_s , D_s , rho , k , mu, P_amb  )
    pressure_v (Data , Data.mesures1_pv , Data.mesures2_pv , DeltaH_Hg, P_amb )
    mdot_f (Data, Data.mesures1_p, Data.mesures2_p , Data.mesuresT , DeltaP , DeltaT  , D_p , Do_p )
    
    with open('{}.pickle' .format(N), 'wb')  as f:
        pickle.dump( Data , f , pickle.HIGHEST_PROTOCOL )       # Pickle the 'data' dictionary using the highest protocol available.
    


# In[9]:


def read_dades (dades):
    with open('smallN1.pickle', 'rb') as f:
        smallNozzleData1     = pickle.load(f)
    

