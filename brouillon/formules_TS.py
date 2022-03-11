# -*- coding: utf-8 -*-
# The MIT License (MIT)

# Copyright (c) 2018 Damien KUSS

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import math
from numba import jit


######################################################
#### Lois de transport solide#########################
######################################################
@jit('f8(f8,f8,f8,f8,f8,f8)')
def Qsrick1990(Q,d90,d30,d50,I,b):
    Q=max(Q,0.001)
    I=max(I,0.001)
    q=Q/b
    qcr=0.065*1.65**1.67*9.81**0.5*d50**1.5*I**(-1.12)

    
    qs=12.6*(d90/d30)**0.2*(q-qcr)*I**2.*1.65**(-1.5)
#     qs=max(qs,0.)
    Qs=qs*b/0.75 #on divise par 0.75 pour passer en débit apparent
    Qs=max(Qs,0.)
    return Qs

#@jit('f8(f8,f8,f8,f8)')
@jit(nopython=True)
def Qsrick1991(Q,d50,I,b):
    Q=max(Q,0.001)
    I=max(I,0.001)
    q=Q/b
    qcr=0.065*1.65**1.67*9.81**0.5*d50**1.5*I**(-1.12)
    qs=1.5*(q-qcr)*I**1.5
#     qs=max(qs,0.)
    Qs=qs*b/0.75 #on divise par 0.75 pour passer en débit apparent
    Qs=max(Qs,0.)
    return Qs

#Loi de transport solide de Lefort
@jit(nopython=True)
def Qslef2015(Q,dm,Gr,L,J):
    Q=max(Q,0.001)
    J=max(J,0.001)
    
    dms=dm*(9.81*1.65/0.000001**2)**(1./3.)
    Cdms=0.0444*(1+15/(1+dms)-1.5*math.exp(-dms/75))
    qs=(Q/L)/(9.81*J*dm**3)**0.5
    if qs<200 :
        rkskr=0.75*(qs/200)**0.23
    else :
        rkskr=0.75
    n=1.6+0.06*math.log10(J)    
    m=1.8+0.08*math.log10(J)
    q0=(9.81*(1.65*dm)**3)**0.5*Cdms*(dm/L)**(1./3.)*rkskr**-0.5*J**-n
    if (dms<14 and rkskr<0.63) :
        cor=1-1.4*math.exp(-0.9*rkskr**2*(Q/L/q0)**0.5)
    else :
        cor=1.
    M=(qs+2.5)/200.
    Z=1+0.38/dms**0.45*(Q/L/(9.81*dm**3)**0.5)**0.192

    if (Q/L)<q0:
        F=0.06*M*(Q/L)/q0
    else:
        F=(6.1*(1-0.938*(q0*L/Q)**0.284)**1.66)**Z
    
    Cp=1700000*J**m*2.65/1.65**1.65*Gr**0.2*cor*F
    Qs=Cp/1000.*Q/2650*(2650./2000.)
    return Qs

@jit(nopython=True)
def Qslef1991(Q,d90,d30,dm,I):
    Q=max(Q,0.001)
    I=min(0.83,max(I,0.001))
    Qcrit=0.0776*(9.81*dm**5)**0.5*1.65**(8./3.)/I**(13./6.)*(1-1.2*I)**(8./3.)
    Qs=4.45*Q*math.pow(I,1.5)/1.65*math.pow(d90/d30,0.2)*(1-math.pow((Qcrit/Q),0.375))
    return Qs

@jit(nopython=True)
def Qsmeunier(Q,I):
    Qs=8.2*Q*I**2
    return Qs

@jit(nopython=True)
def Qspitonrecking2017a(Q,d84tb,d84bs,L,J):  
    Q=max(Q,0.001)
    J=max(J,0.001)
    q=Q/L
    qe=q/(9.81*J*d84bs**3)**0.5
    if qe<100:
        p=0.24
    else:
        p=0.31
    taue=0.015*q**(2*p)*d84bs**(1-3*p)*J**(1-p)/(p**2.5*9.81**p*1.65*d84tb)
    taume=1.5*J**0.75
    phi=14*taue**2.5/(1+(taume/taue)**4)
    qsv=phi*(9.81*1.65*d84tb**3)**0.5
    Qsapp=qsv*L*2650/2000
    return Qsapp

@jit(nopython=True)
def Qspitonrecking2017b(Q,d84tb,d84bs,d50,L,J):  
    Q=max(Q,0.001)
    J=max(J,0.001)
    q=Q/L
    qe=q/(9.81*J*d84bs**3)**0.5
    if qe<100:
        p=0.24
    else:
        p=0.31
    taue=0.015*q**(2*p)*d84bs**(1-3*p)*J**(1-p)/(p**2.5*9.81**p*1.65*d84tb)
    taume=(5*J+0.06)*(d84bs/d50)**(4.4*J**0.5*-1.5)
    phi=14*taue**2.5/(1+(taume/taue)**4)
    qsv=phi*(9.81*1.65*d84tb**3)**0.5
    Qsapp=qsv*L*2650/2000
    return Qsapp

@jit(nopython=True)
def Qspitonrecking2017c(Q,d84tb,d84bs,d50,L,J):  
    Q=max(Q,0.001)
    J=max(J,0.001)
    q=Q/L
    qe=q/(9.81*J*d84bs**3)**0.5
    if qe<100:
        p=0.24
    else:
        p=0.31
    taue=0.015*q**(2*p)*d84bs**(1-3*p)*J**(1-p)/(p**2.5*9.81**p*1.65*d84tb)
    taume=(5*J+0.06)*(d84bs/d50)
    phi=14*taue**2.5/(1+(taume/taue)**4)
    qsv=phi*(9.81*1.65*d84tb**3)**0.5
    Qsapp=qsv*L*2650/2000
    return Qsapp

@jit(nopython=True)
def Qspitonrecking2016(Q,d84bs,L,J):
    Q=max(Q,0.001)
    J=max(J,0.001) 
    q=Q/L
    qe=q/(9.81*J*d84bs**3)**0.5
    if qe<100:
        p=0.24
    else:
        p=0.31
    Qsv=0.00058*L*d84bs**(1.5-7.5*p)*q**(5*p)*J**(2.5*(1-p))/p**6.25/9.81**(2.5*p)
    Qsapp=Qsv*2650/2000
    return Qsapp


