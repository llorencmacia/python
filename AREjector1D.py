# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy

import math
import scipy.optimize as opt

from gastables.Isentropic import Isentropic
isen = Isentropic()

from gastables.ObliqueShock import ObliqueShock
os = ObliqueShock()

from gastables.NormalShock import NormalShock
ns = NormalShock()

def f(variables):
    # This is the function that we want to make 0
    eq=np.zeros(16)
    # We take the unknowns
    (Mc1,Tc1,Mc2,Tc2,pc,Ac1,mdot2,Mm,Tm,pm,theta,Mo,po,To,Mw,pw)=variables
    print "mdot2 = ",mdot2
    print "Mm = ",Mm  
    vc1 = Mc1*np.sqrt(g*r*Tc1)
    vc2 = Mc2*np.sqrt(g*r*Tc2)    
    vm = Mm*np.sqrt(g*r*Tm)
    print "vm = ", vm    
    eq[0] = Tc1 - isen.get_T_by_To_from_M(Mc1)*T01
    eq[1] = pc - isen.get_P_by_Po_from_M(Mc1)*pt01
    eq[2] = Mc2 - isen.get_M_from_P_by_Po(pc/p02)
    eq[3] = Tc2 - isen.get_T_by_To_from_M(Mc2)*T02
    eq[4] = Ac1 - phi2*At*Mt11/Mc1*pow((2+(g-1)*Mc1**2)/(2+(g-1)*Mt11**2), \
            0.5*(g+1)/(g-1))
    eq[5] = mdot2 - pc*(Ac-Ac1)*vc2/(r*Tc2) #In the paper, there is a "pi"
                                                    # probably a mistake
    # Here we estimate the conditions after the mixing zone, applying 
    # momentum and energy conservation and continuity
    eq[6] = pc*Ac - pm*Am + vm*(mdot1+mdot2) - phi3*(mdot1*vc1+mdot2*vc2)
    eq[7] = (mdot1 + mdot2)*(cp*Tm+vm**2/2) - \
            mdot1*(cp*Tc1+vc1**2/2) + mdot2*(cp*Tc2+vc2**2/2)
    eq[8] = mdot1 + mdot2 - vm*pm*Am/(r*Tm)
    # Here we simulate expansion effects with an oblique shock wave 
    if (Mm > 1):     
        eq[9] = theta - os.get_Wave_Angle_from_Mx_and_Turn_Angle(Mm,delta)[0]
        eq[10] = Mo - os.get_My_from_Mx_and_Turn_Angle(Mm,delta)
        eq[11] = po - os.get_Py_by_Px_from_Mx_and_Turn_Angle(Mm,delta)*pm
        eq[12] = To - os.get_Ty_by_Tx_from_Mx_and_Turn_Angle(Mm,delta)*Tm
    else:
        eq[9] = theta - 0
        eq[10] = Mo - Mm
        eq[11] = po - pm
        eq[12] = To - Tm        

    # There is no friction, so far
    Ms = Mo
    ps = po
    # Now, the final normal shock wave
    if (Ms > 1):
        eq[13] = Mw - ns.get_My_from_Mx(Ms)
        eq[14] = pw - ns.get_Py_by_Px_from_Mx(Ms)*ps
    else:
        eq[13] = Mw - Ms
        eq[14] = pw -ps
    #And the isentropic flow to outlet
    eq[15] = pe - pw/isen.get_P_by_Po_from_M(Mw)
    
    return eq

if __name__=="__main__":
    
    #Some parameters
    g = 1.4
    cp = 1005
    r = 287
    delta = np.deg2rad(10)
    phi1 = 0.95
    phi2 = 1
    phi3 = 1    

    #input data
    p01 = 7.0e5
    T01 = 300.0
    p02 = .3e5
    T02 = 300.0
    dg = 2.75e-3 # This is the diameter in the throat
    dt = 6.5e-3 # diameter in the exit of first nozzle
    dc = 20.0e-3 # diameter of entrance of secondary flow. Beginning of mixing
    dm = 8.5e-3 # diameter of the second nozzle. End of mixing zone
    Ag = np.pi*dg**2/4
    At = np.pi*dt**2/4
    Ac = np.pi*dc**2/4
    Am = np.pi*dm**2/4
    pe = 1.0e5
    Mt1 = isen.get_M_from_A_by_Astar(At/Ag)[1] #The component 0 is subsonic
    pt1 = isen.get_P_by_Po_from_M(Mt1)*p01
    Tt1 = isen.get_T_by_To_from_M(Mt1)*T01
    Mt11 = ns.get_My_from_Mx(Mt1)
    pt11 = ns.get_Py_by_Px_from_Mx(Mt1)*pt1
    Tt11 = ns.get_Ty_by_Tx_from_Mx(Mt1)*Tt1
    mdot1 = p01*Ag*np.sqrt(phi1*g/(r*T01)*pow(2/(g+1),(g+1)/(g-1)))
    pt01 = ns.get_Poy_by_Pox_from_Mx(Mt1)*p01
    
    #Initial guesses
    Mc10 = 4*Mt11
    Tc10 = isen.get_T_by_To_from_M(Mc10)*T01
    Ac10 = phi2*At*Mt11/Mc10*pow((2+(g-1)*Mc10**2)/(2+(g-1)*Mt11**2), \
            0.5*(g+1)/(g-1))
    pc0 = isen.get_P_by_Po_from_M(Mc10)*pt01
    Mc20 = isen.get_M_from_P_by_Po(pc0/p02)
    Tc20 = isen.get_T_by_To_from_M(Mc20)*T02
    vc10 = Mc10*np.sqrt(g*r*Tc10)
    vc20 = Mc20*np.sqrt(g*r*Tc20)
    mdot20 = pc0*(Ac-Ac10)*vc20/(r*Tc20)
    vm0 = phi3*(mdot1*vc10+mdot20*vc20)/(mdot1+mdot20)
    Tm0 = ((mdot1*(cp*Tc10+vc10**2/2) + mdot20*(cp*Tc20+vc20**2/2))/ \
            (mdot1 + mdot20) - vm0**2/2)/cp
    pm0 = (mdot1+mdot20)*r*Tm0/(vm0*Am)
    Mm0 = vm0/np.sqrt(g*r*Tm0)
    #theta0 = os.get_Wave_Angle_from_Mx_and_Turn_Angle(Mm0,delta)[0]
    #Mo0 = os.get_My_from_Mx_and_Turn_Angle(Mm0,delta)
    #po0 = os.get_Py_by_Px_from_Mx_and_Turn_Angle(Mm0,delta)*pm0
    theta0 = delta
    Mo0 = Mm0
    po0 = pm0    
    To0 = Tm0
    #Mw0 = ns.get_My_from_Mx(Mm0)
    #pw0 = ns.get_Py_by_Px_from_Mx(Mm0)*po0
    Mw0 = Mm0
    pw0 = po0    
    
    # We just call opt passing the function and an initial value
    solution = opt.root(f,(Mc10,Tc10,Mc20,Tc20,pc0,Ac10,mdot20, \
            Mm0,Tm0,pm0,theta0,Mo0,po0,To0,Mw0,pw0))
    Mc1 = solution.x[0]
    Tc1 = solution.x[1]
    Mc2 = solution.x[2]
    Tc2 = solution.x[3]
    pc = solution.x[4]
    Ac1 = solution.x[5]
    mdot2 = solution.x[6]
    Mm = solution.x[7]
    Tm = solution.x[8]
    pm = solution.x[9]
    theta = solution.x[10]
    Mo = solution.x[11]
    po = solution.x[12]
    To = solution.x[13]
    Mw = solution.x[14]
    pw = solution.x[15]

    print "The solution is :"
    print "mdot1 = ",mdot1
    print "Mc1 = ",Mc1
    print "Mc2 = ",Mc2
    print "Tc1 = ",Tc1
    print "Tc2 = ",Tc2
    print "pc = ",pc
    print "Ac1 = ",Ac1
    print "mdot2 = ",mdot2
    print "Mm = ",Mm
    print "Tm = ",Tm
    print "pm = ",pm
    print "theta = ",theta
    print "Mo = ",Mo
    print "po = ",po
    print "To = ",To
    print "Mw = ",Mw
    print "pw = ",pw
    print "Success : ",solution.success
    