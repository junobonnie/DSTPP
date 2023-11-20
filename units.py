# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:16:35 2023

@author: replica
"""

#time
Gyr = 1

s = 3.1710e-17*Gyr
ms = 1e-3*s
ns = 1e-6*s

minute = 60*s
hour = 3600*s
day = 86400*s
year = 3.154e7*s

#length
kpc = 1
pc = 1e-3*kpc
Gpc = 1e3*kpc

m = 3.24078e-20*kpc
km = 1e3*m
cm = 1e-2*m
mm = 1e-3*m
nm = 1e-6*m

#mass
M_sun = 1

kg = 5.0276521e-31*M_sun
g = 1e-3*kg

if __name__=="__main__":
    G = 6.67384e-11*(m**3)*(kg**-1)*(s**-2)
    print(G)
    print(5e8*M_sun*2.17*(cm**2)/g/25/3.14)



