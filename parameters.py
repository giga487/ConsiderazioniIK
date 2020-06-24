# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 15:44:05 2020

@author: ariannagasparri
"""

# Finite Elements
N = 250

### Solver Parameters
myMaxIter = 1000

### Numerical Constants
# lunghezza dei link
a1 = 1.3
a2 = 1.6
# Masse dei link approssimate
m1 = 1500
m2 = 1350
# Accelerazione di gravità
g = 9.81
# Coefficienti della funzione costo
kh = 0.
kj = 1.
kt = 0.

# Momenti di inerzia principali rispetto ai frame locali
Ixx1 = (m1*((a1/2)**2))/3
Izz1 = (m1*((a1/2)**2))/3
Ixx2 = (m2*((a2/2)**2))/3
Izz2 = (m2*((a2/2)**2))/3