#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#import CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI


#8<------------------------------------------------------------------------
## Para pasar de grados Celsius a grados Kelvin
def c2K(T):
  return T + 273.15


#8<------------------------------------------------------------------------

print(">Main<")
