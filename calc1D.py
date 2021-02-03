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

    # We take the unknowns
    (mdot1,Mc1,Mc2,Tc1,Tc2,pc,Ac1,mdot2,Mm,Tm,theta,Mo,po,To,Mw,pw)=variables
    vc1 = Mc1*np.sqrt(g*r*Tc1)
    vc2 = Mc2*np.sqrt(g*r*Tc2)
    Mm2 = Mm**2
    Mw2 = Mw**2
    vm = Mm*np.sqrt(g*r*Tm)
    eq=np.zeros(16)
    eq[0] = mdot1 - p01*Ag*np.sqrt(phi1*g/(r*T01)*pow(2/(g+1),(g+1)/(g-1)))
    #eq[1] = pc - p01/pow(1+(g-1)/2*Mc12,g/(g-1))
    eq[1] = pc - isen.get_P_by_Po_from_M(Mc1)*p01
    #eq[2] = Tc1 - T01/(1+(g-1)/2*Mc12)
    eq[2] = Tc1 - isen.get_T_by_To_from_M(Mc1)*T01  
    #eq[3] = Ac1 - phi2*At*Mc1/Mt1*pow((2+(g-1)*Mc12)/(2+(g-1)*Mt1**2), \
    #        0.5*(g+1)/(g-1))
    eq[3] = Ac1 - phi2*isen.get_A_by_Astar_from_M(Mc1)*Ag
    #eq[4] = pc - p02/pow(1+(g-1)/2*Mc22,g/(g-1))
    eq[4] = pc - isen.get_P_by_Po_from_M(Mc2)*p02
    #eq[5] = Tc2 - T02/(1+(g-1)/2*Mc22)
    eq[5] = Tc2 - isen.get_T_by_To_from_M(Mc2)*T02
    vc2 = Mc2*np.sqrt(g*r*Tc2)
    eq[6] = mdot2 - pc*(Ac-Ac1)*vc2/(r*Tc2) #In the paper, there is a "pi"
                                                    # probably a mistake
    # Here we estimate the conditions after the mixing zone, applying 
    # momentum and energy conservation    
    eq[7] = vm*(mdot1+mdot2) - phi3*(mdot1*vc1+mdot2*vc2)
    eq[8] = (mdot1 + mdot2)*(cp*Tm+vm**2/2) - \
            mdot1*(cp*Tc1+vc1**2/2) + mdot2*(cp*Tc2+vc2**2/2)
    # Here we simulate expansion effects with an oblique shock wave        
#    eq[9] = np.tan(delta) - 2/np.tan(theta)*(Mmst2-1) \
#            /(2+Mm2*(g+np.cos(2*theta)))
    eq[9] = theta - os.get_Wave_Angle_from_Mx_and_Turn_Angle(Mm,delta)[0]
#    eq[10] = Mo - 1/np.sin(theta-delta)*np.sqrt((2+(g-1)*Mmst2) \
#            /(1-g+2*g*Mmst2))
    eq[10] = Mo - os.get_My_from_Mx_and_Turn_Angle(Mm,delta)
#    eq[11] = po - pc*(1+2*g/(g+1)*(Mmst2-1))
    eq[11] = po - os.get_Py_by_Px_from_Mx_and_Turn_Angle(Mm,delta)*pc
#    eq[12] = To - Tm*(1+2*g/(g+1)*(Mmst2-1)) * ((g-1)*Mmst2+2)/((g+1)*Mmst2)
    eq[12] = To - os.get_Ty_by_Tx_from_Mx_and_Turn_Angle(Mm,delta)*Tm
    # There is no friction, so far
    Ms = Mm
    ps = pc
    Ms2 = Ms**2
    # Now, the final normal shock wave
    eq[13] = Mw - ns.get_My_from_Mx(Mm)
    eq[14] = pw - ns.get_Py_by_Px_from_Mx(Mm)*po
    #And the isentropic flow to outlet
    eq[15] = pe - pw/isen.get_P_by_Po_from_M(Mw)

    
    
    return eq

if __name__=="__main__":
    
    #Some parameters
    g = 1.4
    cp = 1005
    r = 287
    delta = np.deg2rad(15)
    phi1 = 1
    phi2 = 1
    phi3 = 1    

    #input data
    p01 = 7.0e5
    T01 = 300.0
    p02 = 1.0e5
    T02 = 300.0
    dg = 2.75e-3
    dt = 5.5e-3
    dc = 8.5e-3
    Ag = np.pi*dg**2/4
    At = np.pi*dt**2/4
    Ac = np.pi*dc**2/4
    pe = 1.0e5
    Mt1 = isen.get_M_from_A_by_Astar(At/Ag)[1] #The component 0 is subsonic
    
    #Initial guesses
    mdot10 = 0.01
    Mc10 = 3.5
    #Tc10 = T01/(1+(g-1)/2*Mc10**2)
    Tc10 = isen.get_T_by_To_from_M(Mc10)*T01
    pc0 = p01/pow(1+(g-1)/2*Mc10**2,g/(g-1))
    Mc20 = np.sqrt(2/(g-1)*(pow((p02/pc0),(g-1)/g)-1))
    Tc20 = T02/(1+(g-1)/2*Mc20**2)
    #Ac10 = phi2*At*Mc10/Mt1*pow((2+(g-1)*Mc10**2)/(2+(g-1)*Mt1**2), \
    #        0.5*(g+1)/(g-1))
    Ac10 = phi2*isen.get_A_by_Astar_from_M(Mc10)*Ag
    vc10 = Mc10*np.sqrt(g*r*Tc10)
    vc20 = Mc20*np.sqrt(g*r*Tc20)
    mdot20 = pc0*(Ac-Ac10)*vc20/(r*Tc20)
    vm0 = phi3*(mdot10*vc10+mdot20*vc20)/(mdot10+mdot20)
    Tm0 = ((mdot10*(cp*Tc10+vc10**2/2) + mdot20*(cp*Tc20+vc20**2/2))/ \
            (mdot10 + mdot20) - vm0**2/2)/cp
    Mm0 = vm0/np.sqrt(g*r*Tm0)
    theta0 = os.get_Wave_Angle_from_Mx_and_Turn_Angle(Mm0,delta)[0]
    Mo0 = os.get_My_from_Mx_and_Turn_Angle(Mm0,delta)
    po0 = os.get_Py_by_Px_from_Mx_and_Turn_Angle(Mm0,delta)*pc0
    To0 = Tm0
    Mw0 = ns.get_My_from_Mx(Mm0)
    pw0 = ns.get_Py_by_Px_from_Mx(Mm0)*po0
    
    
    # We just call opt passing the function and an initial value
    solution = opt.root(f,(mdot10,Mc10,Mc20,Tc10,Tc20,pc0,Ac10,mdot20, \
            Mm0,Tm0,theta0,Mo0,po0,To0,Mw0,pw0))
    mdot1 = solution.x[0]
    Mc1 = solution.x[1]
    Mc2 = solution.x[2]
    Tc1 = solution.x[3]
    Tc2 = solution.x[4]
    pc = solution.x[5]
    Ac1 = solution.x[6]
    mdot2 = solution.x[7]
    Mm = solution.x[8]
    Tm = solution.x[9]
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
    print "theta = ",theta
    print "Mo = ",Mo
    print "po = ",po
    print "To = ",To
    print "Mw = ",Mw
    print "pw = ",pw
    print "Success : ",solution.success
    