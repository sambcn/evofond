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
from lib import visvalingamwyatt as vw
import win32api
import pandas as pd
from tkinter import filedialog
from tkinter import Tk
import os
from scipy import integrate
 


def hydcrue(Qmax,tm,alpha,Qbase,d,dt):
    tlav=np.arange(0,d+dt,dt)
    Q=np.array([(Qmax-Qbase)*2*math.pow(x/tm,alpha)/(1+math.pow(x/tm,2*alpha))+Qbase for x in tlav])
    Q=[0.0001 if x==0 else x for x in Q]
    return [tlav,Q]

def hydrogramme():
    ##### Definition feuilles
    sht = xw.Book.caller().sheets['Hydrogramme']
    
    ##### User Inputs
    [Qmax,Qbase,tm,alpha,d,dt] = sht.range('B2:B7').value
 
    ##### Excel: clear output, write out initial values of percentiles/sample path and set chart source
    ##### and x-axis values
    sht.range('G3').expand().clear_contents()

    ###### Preallocation
    t=hydcrue(Qmax,tm,alpha,Qbase,d,dt)[0]
    Q=hydcrue(Qmax,tm,alpha,Qbase,d,dt)[1]
    V=integrate.cumtrapz(Q, t,initial=0)*3600
    sht.range('G3').options(transpose=True).value=np.array([t,Q,V])
    sht.range('M19').options(transpose=True).value=V[-1]
    ##### Figure HYDROGRAMMES 
    plot_hydrogramme(t,Q,V,sht)   
    

def interpolation():
    ##### Definition feuilles
    sht = xw.Book.caller().sheets['Hydrogramme']
    
    try:
        ##### User inputs
        [t,Q]=sht.range('D3').expand('table').options(np.array,ndim=2,transpose=True).value
        dt=sht.range('B10').value
        
        ##### Effacer donnees existantes
        sht.range('G3').expand().clear_contents()
        
        ##### Interpolation et calcul V cumulé
        [tnew,Qnew]=interpol(t,Q,dt)
        V=integrate.cumtrapz(Qnew, tnew,initial=0)*3600
        
        ##### Ecriture résultats
        sht.range('G3').options(transpose=True).value=np.array([tnew,Qnew,V])
        sht.range('M19').options(transpose=True).value=V[-1]
        
        ##### Graphique
        plot_hydrogramme(tnew,Qnew,V,sht)           
    except:
        wb = xw.Book.caller()
        win32api.MessageBox(wb.app.hwnd, "Erreur : t et Q probablement non renseignés")
        

def interpol(t,Q,dt):
    #Interpolation de t,Q tous les dt
    tnew=[]
    for i in range(len(t)-1):  #len(x0)-1
        a=np.arange(t[i],t[i+1],dt)
        if len(a)==1 and a==[t[i]]:
            a=[t[i],t[i+1]]
        tnew=np.append(tnew,a)
    tnew=np.append(tnew,t[-1])
    tnew=np.unique(tnew)
    Qnew=np.interp(tnew,t,Q)
    return [tnew,Qnew]

def plot_hydrogramme(t,Q,V,sht):
    
    ### Figure HYDROGRAMMES 
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1,1) #height_ratios=[4,2,2])
    ax1 = fig.add_subplot(gs[0])
    
    s1=ax1.plot(t,Q,'k',label=u'débit (m$^{3}$/s)')
    ax2=ax1.twinx()
    s2=ax2.plot(t,V,'r',label=u'Volume liquide cumulé (m$^{3}$)')

        
    ax1.grid()
    ax1.set_xlabel(u'Temps (h)',fontsize=22)
    ax1.set_ylabel('Débit (m$^{3}$/s)',fontsize=22)
    ax2.set_ylabel('Volume liquide (m$^{3}$)',fontsize=22,color='r')
    ax2.tick_params('y', colors='r')

    # added these three lines
    courbes = s1+s2
    labs = [l.get_label() for l in courbes]
    ax1.legend(courbes, labs, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize=16)

    matplotlib.rcParams.update({'font.size':18})
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 18)
    plt.tight_layout()
            
    
    plot=sht.pictures.add(fig, name='Hydrogramme',left=sht.range('K2').left, top=sht.range('K2').top,update=True)
    plot.height = 8*30
    plot.width = 12*30      

 
def simplification():
    sht = xw.Book.caller().sheets['Hydrogramme']
#    t=sht.range('G3').expand('down').value
#    Q=sht.range('H3').expand('down').value
    
    try:
        assert np.array(sht.range('G3').expand().value).size>2
        if np.array(sht.range('G3').expand().value).size>2:
            t,Q=np.array(sht.range('G3').expand().value).T
            pts=[(x,y) for x,y in zip(t,Q)]
            r = sht.range('B13').value
            
            sht.range('J3').expand().clear_contents()
        
            
            #Methode simplification
            simplifier = vw.Simplifier(pts)
            
            # Simplify by percentage of points to keep
            pts1=simplifier.simplify(ratio=r)
            # pts1=simplifier.simplify(threshold=30)
            ts=pts1[:,0]
            Qs=pts1[:,1]
        
            sht.range('J3').options(transpose=True).value=np.array([ts,Qs])
#            sht.range('J3').options(transpose=True).value=ts
#            sht.range('K3').options(transpose=True).value=Qs
            
            #### Figure HYDROGRAMMES 
            fig = plt.figure(figsize=(12,8))
            gs = gridspec.GridSpec(1,1) #height_ratios=[4,2,2])
            ax1 = fig.add_subplot(gs[0])
            
            ax1.plot(t,Q,'k',label=u'Hydrogramme (m$^{3}$/s)')
            ax1.plot(ts,Qs,'r--',label=u'Hydrogramme simplifiié (m$^{3}$/s)')
        
                
            ax1.grid()
            ax1.set_xlabel(u'Temps (h)',fontsize=22)
            ax1.set_ylabel('Débit (m$^{3}$/s)',fontsize=22)
        
        
            ax1.legend(loc=1,fontsize=16)
        
            matplotlib.rcParams.update({'font.size':18})
            ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
            plt.tight_layout()
        
        
            
            plot1=sht.pictures.add(fig, name='H_simplifie',left=sht.range('M20').left, top=sht.range('M20').top,update=True)
            plot1.height = 8*30
            plot1.width = 12*30
        
    except AssertionError:
        wb = xw.Book.caller()
        win32api.MessageBox(wb.app.hwnd, "Vous devez calculer l'hydrogramme au préalable")
 
        
    
    
def effacer():
    sht = xw.Book.caller().sheets['Hydrogramme']
    sht.range('G3').expand().clear_contents()
    sht.range('D3').expand().clear_contents()
    sht.range('J3').expand().clear_contents()
    try:
        sht.pictures('H_simplifie').delete()
    except :
        pass
    try:
        sht.pictures('Hydrogramme').delete()
    except :
        pass
        
    #sht.pictures('H_simplifie').delete()
    #sht.pictures('Hydrogramme').delete()
    
    #☻plot1.clear_contents()
    

def importer():
    """
    Importe fichier .txt d'hydrogramme
    fichier comportant entête t(h),Q(m3/s)
    puis sur chaque ligne les valeurs ti Qi
    """
    sht = xw.Book.caller().sheets['Hydrogramme']
    ##### Adresse du répertoire parent
    pardir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    
    ##### Récupère le symbole de séparation des décimales (, ou .)
    dec=sht.range("B13").value
    
    try:
        Tk().withdraw()  
        txt_filename =  filedialog.askopenfilename(title = "Selectionner fichier",filetypes = [("txt file","*.txt")],initialdir=pardir)
        data = pd.read_csv(txt_filename,sep="\t",skiprows=1,header=None,decimal=dec)
        sht.range('D3:E3').expand('down').clear_contents()
        sht.range('D3').value=np.array(data)
    except:
        pass    

    
#    #### Figure HYDROGRAMMES 
#    fig = plt.figure(figsize=(12,8))
#    gs = gridspec.GridSpec(1,1) #height_ratios=[4,2,2])
#    ax1 = fig.add_subplot(gs[0])
#    
#    ax1.plot(tl,Ql,'k',label=u'débit (m$^{3}$/s)')
#
#        
#    ax1.grid()
#    ax1.set_xlabel(u'Temps (h)',fontsize=22)
#    ax1.set_ylabel('Débit (m$^{3}$/s)',fontsize=22)
#
#
#    ax1.legend(loc=1,fontsize=16)
#
#    matplotlib.rcParams.update({'font.size':18})
#    ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
#    plt.tight_layout()
#
#
#    
#    plot=sht.pictures.add(fig, name='Hydrogramme',left=sht.range('M2').left, top=sht.range('M2').top,update=True)
#    plot.height = 8*30
#    plot.width = 12*30

    
    
    
    

