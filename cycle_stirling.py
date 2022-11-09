#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

#import CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import StateContainer
from CoolProp.Plots import SimpleRankineCycle

#from CoolProp.Plots import SimpleCompressionCycle

print('CoolProp version: ', CoolProp.__version__)
print('CoolProp gitrevision: ', CoolProp.__gitrevision__)
print('CoolProp fluids: ', CoolProp.__fluids__)

def c2K(T):
    return T + 273.15;

## YouTube
#DataIn = {"Fluid" : "Air", "TLow" : c2K(27), "THigh" : c2K(527), "PLow" : 0.1*10**3, "PHigh" : 2*10**6, "RendIseTurb" : 1, "RendIsePump" : 1}

# fuente :http://laplace.us.es/wiki/index.php/Caso_pr%C3%A1ctico_de_ciclo_de_Stirling
## Diatómico
DataIn = {"Fluid" : "Air", "TLow" : 300, "THigh" : 2000, "PLow" : 150*10**3, "PHigh" : 3*10**6, "RendIseTurb" : 1, "RendIsePump" : 1, "moles":100}

stat = StateContainer()
##M / R / Cp / Cv

Mair = PropsSI("MOLARMASS",DataIn["Fluid"]);print("Mair: ",Mair)
R = PropsSI("GAS_CONSTANT",DataIn["Fluid"]);print("R: ",R)
Rair = R/(Mair*1000); print("Rair: ",Rair) ## En kJ/(kgK)
## Para diatómico
cv = 5/2*Rair
#Cp = PropsSI("CP0MASS","T",300,DataIn["Fluid"]);print(Cp)
#Cv = PropsSI("CVMASS","",0,"",0,DataIn["Fluid"]);print(Cv)

n = DataIn["moles"]
MassAir = Mair*n
print("Masa de Aire: ",MassAir)
## Pv = RairT
## v = RairT/P
## ESTADO 0 :
stat[0,"P"] = DataIn["PLow"]
stat[0,"T"] = DataIn["TLow"]
v0 = Rair*stat.get_point(0).T/(stat.get_point(0).P/1000)
stat[0,"D"] = 1/v0
#print(1/stat.get_point(0).D*MassAir)

## ESTADO 2 :
stat[2,"T"] = DataIn["THigh"]
stat[2,"P"] = DataIn["PHigh"]
v2 = Rair*stat.get_point(2).T/(stat.get_point(2).P/1000)
stat[2,"D"] = 1/v2
#print(1/stat.get_point(2).D*MassAir)

## ESTADO 1 :
stat[1,"T"] = stat.get_point(0).T
stat[1,"D"] = stat.get_point(2).D
P1 = Rair*stat.get_point(1).T*stat.get_point(1).D*1000
#print(P1/10**6)
stat[1,"P"] = P1

## ESTADO 3 :
stat[3,"T"] = stat.get_point(2).T
stat[3,"D"] = stat.get_point(0).D
P3 = Rair*stat.get_point(3).T*stat.get_point(3).D*1000
#print(P3/10**6)

#8<------------------------------------------------------------------------

print("Estado 1: ",stat.get_point(0))
print("Estado 2: ",stat.get_point(1))
print("Estado 3: ",stat.get_point(2))
print("Estado 4: ",stat.get_point(3))

#8<------------------------------------------------------------------------
## Balance, calor y tabajo
## Compresion isotermica 0-1
w01 = 1000*Rair *stat.get_point(0).T*np.log(stat.get_point(1).D/stat.get_point(0).D)
#print(w01*MassAir/10**6)
q01 = -w01

## Calentamiento isocórico 1-2
w12 = 0
q12 = cv*(stat.get_point(2).T-stat.get_point(1).T)*1000
#print(q12*MassAir/10**6)

## Expansion Isotermica 2-3
w23=1000*Rair*stat.get_point(2).T*np.log(stat.get_point(3).D/stat.get_point(2).D)
#print(w23*MassAir/10**6)
q23=-w23

## Enfriamiento Isocórico 3-0
w30 = 0
q30 = cv*(stat.get_point(0).T-stat.get_point(3).T)*1000
#print(q30*MassAir/10**6)

win = w01
wout = w23
wnet = np.abs(wout + win)

#print(win*MassAir/10**6)
#print(wout*MassAir/10**6)
#print(wnet*MassAir/10**6)

print("Trabajo In: ",round(win/1000,1),"[kJ/kg]")
print("Trabajo Out: ",round(wout/1000,1),"[kJ/kg]")
print("Trabajo Net: ",round(wnet/1000,1),"[kJ/kg]")

qin = q12 + q23
#print(qin*MassAir/10**6)
#qout = stat.get_point(1).H - stat.get_point(2).H
print("Calor In: ",round(qin/1000,1),"[kJ/kg]")
#print("Calor Out: ",round(qout/1000,1),"[kJ/kg]")

rend_term = (wnet/qin)*100
print("Rend Termico: ", round(rend_term,1))

#8<------------------------------------------------------------------------
## Obtención de entropias

## Compresion isotermica 0-1
s01 = q01/stat.get_point(0).T
print("ds01: ",round(s01/10**3,3),"[kJ/kgK]")
print("dS01: ",round(s01*MassAir/10**3,3),"[kJ/K]")

## Calentamiento isocórico 1-2
s12 =cv*np.log(stat.get_point(2).T/stat.get_point(1).T)*1000
print("ds12: ",round(s12/10**3,3),"[kJ/kgK]")
print("dS12: ",round(s12*MassAir/10**3,3),"[kJ/K]")

## Expansion Isotermica 2-3
s23 = q23/stat.get_point(2).T
print("ds23: ",round(s23/10**3,3),"[kJ/kgK]")
print("dS23: ",round(s23*MassAir/10**3,3),"[kJ/K]")

## Enfriamiento Isocórico 3-0
s30 =cv*np.log(stat.get_point(0).T/stat.get_point(3).T)*1000
print("ds30: ",round(s30/10**3,3),"[kJ/kgK]")
print("dS30: ",round(s30*MassAir/10**3,3),"[kJ/K]")

#8<------------------------------------------------------------------------

## Realizar rutinas para graficar en matplotlib