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
"""
Created on Mon Sep 17 08:53:02 2018

@author: dk33085
"""


import numpy as np
import math
import xlwings as xw
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lib import formules_TS as TS
from scipy import integrate



def main():
    ##### Définition des feuilles
    sht = xw.Book.caller().sheets[u'Sédimentogramme']
    sht1 = xw.Book.caller().sheets['Hydrogramme']
    
    ##### User Inputs
    data=sht.range('B2:B17').value
    [loi,vide,vide,vide,b,J,vide,vide,vide,dm,d50,d90,d30,d84tb,d84bs,Gr]=data
    dt=sht1.range('G4').value-sht1.range('G3').value
    [t,Q]=np.transpose(sht1.range('G3:H3').expand('down').options(np.array,ndim=2).value)

    ##### Préallocation   
    Qsin=np.zeros(len(Q))
    Qscumul=np.zeros(len(Q))
    Ql=10. #Valeurbidon pour dictionnaire d'aguments de fonctions TS

    ##### DICTIONNAIRe Fonction TS   
#    fndict = {"LEF2015": TS.Qslef2015, "LEF1990":  TS.Qslef1990, "RICK":  TS.Qsrick, "MEUNIER": TS.Qsmeunier}
    fndict = {"LEF2015": TS.Qslef2015, "LEF1991":  TS.Qslef1991, "RICK1990":  TS.Qsrick1990, "RICK1991": TS.Qsrick1991,"MEUN1989": TS.Qsmeunier,
              "PIRECK2017a":TS.Qspitonrecking2017a,  "PIRECK2017b":TS.Qspitonrecking2017b,"PIRECK2016":TS.Qspitonrecking2016}
#    argsdict = {"LEF2015": [Ql,dm,Gr,b,J], "LEF1990":[Ql,d90,d30,dm,J], "RICK": [Ql,d90,d30,d50,J,b], "MEUNIER": [Ql,J]}
    argsdict = {"LEF2015": [Ql,dm,Gr,b,J], "LEF1991":[Ql,d90,d30,dm,J], "RICK1990": [Ql,d90,d30,d50,J,b],"RICK1991": [Ql,d50,J,b],"MEUN1989":[Ql,J],
                "PIRECK2017a":[Ql,d84tb,d84bs,b,J],"PIRECK2017b":[Ql,d84tb,d84bs,d50,b,J],"PIRECK2016":[Ql,d84bs,b,J]}
    
    
    ##### Effacer valeurs
    sht.range('d3').expand().clear_contents()
    
    ##### Calcul Qs et Qs_cumul
    for i in range(len(Q)):
        argsdict[loi][0]=Q[i]
        Qsin[i]=fndict[loi](*argsdict[loi])        
        Qscumul=integrate.cumtrapz(Qsin, t,initial=0)*3600
        
        Vs=Qscumul[-1]
        Vl=np.trapz(Q, t)*3600
        Cv=Vs/(Vs+Vl)
    
    ##### Ecriture résultats
    sht.range('D3').options(transpose=True).value=np.array([t,Qsin,Qscumul])
    sht.range('P2').options(transpose=True).value=np.array([Vl,Vs,Cv])
    
    
#    
    ##### Figure SEDIMENTOGRAMMES 1
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1,1) #height_ratios=[4,2,2])
    ax1 = fig.add_subplot(gs[0])
    
    s1=ax1.plot(t,Q,'k',label=u'débit liquide (m$^{3}$/s)')
    
    ax2=ax1.twinx()
    s2=ax2.plot(t,Qscumul,'r',label=u'Volume solide cumulé (m$^{3}$)')

        
    ax1.grid()
    ax1.set_xlabel(u'Temps (h)',fontsize=22)
    ax1.set_ylabel('Débit (m$^{3}$/s)',fontsize=22)
    ax2.set_ylabel('Volume solide (m$^{3}$)',fontsize=22,color='r')
    ax2.tick_params('y', colors='r')

    # added these three lines
    courbes = s1+s2
    labs = [l.get_label() for l in courbes]
    ax1.legend(courbes, labs, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize=16)
#    ax1.legend(loc=1,fontsize=16)

    matplotlib.rcParams.update({'font.size':18})
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 18)
    plt.tight_layout()

    plot=sht.pictures.add(fig, name='S2',left=sht.range('H20').left, top=sht.range('H20').top,update=True)
    plot.height = 8*30
    plot.width = 12*30
#    
    #### Figure SEDIMENTOGRAMMES 2
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1,1) #height_ratios=[4,2,2])
    ax1 = fig.add_subplot(gs[0])
    
    s1=ax1.plot(t,Q,'k',label=u'débit liquide (m$^{3}$/s)')
    
    ax2=ax1.twinx()
    s2=ax2.plot(t,Qsin,'r',label=u'débit solide (m$^{3}$/s)')

        
    ax1.grid()
    ax1.set_xlabel(u'Temps (h)',fontsize=22)
    ax1.set_ylabel('Débit (m$^{3}$/s)',fontsize=22)
    ax2.set_ylabel('Débit solide (m$^{3}$/s)',fontsize=22,color='r')
    ax2.tick_params('y', colors='r')

    # added these three lines
    courbes = s1+s2
    labs = [l.get_label() for l in courbes]
    ax1.legend(courbes, labs, loc=1,fontsize=16)
#    ax1.legend(loc=1,fontsize=16)

    matplotlib.rcParams.update({'font.size':18})
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 18)
    plt.tight_layout()

    plot=sht.pictures.add(fig, name='S1',left=sht.range('H2').left, top=sht.range('H2').top,update=True)
    plot.height = 8*30
    plot.width = 12*30
    
def effacer():
    sht = xw.Book.caller().sheets['Sédimentogramme']
    sht.range('D3').expand().clear_contents()
#    sht.range('D3').expand().clear_contents()
    try:
        sht.pictures('S1').delete()
    except :
        pass
    try:
        sht.pictures('S2').delete()
    except :
        pass

