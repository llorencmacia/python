# -*- coding: utf-8 -*-
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

def SIGMA(M):
    return isen.get_A_by_Astar_from_M(M)



if __name__=="__main__":

    
