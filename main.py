#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#import CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI

print(">Main<")

#8<------------------------------------------------------------------------
## Para pasar de grados Celsius a grados Kelvin
def c2K(T):
  return T + 273.15


#8<------------------------------------------------------------------------
## Para encontrar propiedades termodinámicas, por ejemplo entalpía del agua a 3MPa y 350°C.

P = 3*10**6
T = c2K(350)
h = PropsSI("H","P",P,"T",T,"Water")
print(h)