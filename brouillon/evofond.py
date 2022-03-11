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
Created on Tue Sep 18 10:07:01 2018

@author: dk33085
"""
import matplotlib.pyplot as plt
import numpy as np
import timeit
import matplotlib
import matplotlib.gridspec as gridspec
import xlwings as xw
from lib import formules_TS as TS
from lib import hydrogramme as hyd
import os
import copy

def main():

    #############################################
    ############## User inputs ##################
    #############################################
    sht = xw.Book.caller().sheets[u'Profil EVOFOND']
    sht1 = xw.Book.caller().sheets[u'Calculs']
    sht2 = xw.Book.caller().sheets[u'Hydrogramme']
    sht3 = xw.Book.caller().sheets[u'Sédimentogramme']
    sht4 = xw.Book.caller().sheets[u'Granulometrie']
    sht5 = xw.Book.caller().sheets[u'zfond']
    sht6 = xw.Book.caller().sheets[u'B']
    sht9 = xw.Book.caller().sheets[u'zfond + hc']
    sht10 = xw.Book.caller().sheets[u'zfond + Hsc']
    sht11 = xw.Book.caller().sheets[u'Vol']

    
    ##### User Inputs
    [x0,z0,b0] = (np.transpose(sht.range('A3:C3').expand('down').options(np.array,ndim=2).value))
    [xmin0,zmin0] = (np.transpose(sht.range('F3:G3').expand('down').options(np.array,ndim=2).value))
    [xseuil0,zseuil0] = (np.transpose(sht.range('I3:J3').expand('down').options(np.array,ndim=2).value))
    [dx0,dt,kborda,loi,tenr]=sht1.range('B2:B6').value
    [t0,Q0] = (np.transpose(sht2.range('G3:H3').expand('down').options(np.array,ndim=2).value))
    [tsin0,Qsin0] = (np.transpose(sht3.range('D3:E3').expand('down').options(np.array,ndim=2).value))
    [[zmingr1,zmingr2,zmingr3,zmingr4],[zmaxgr1,zmaxgr2,zmaxgr3,zmaxgr4]]=sht.range('M3:P4').value
    [dm_Gr1,d50_Gr1,d90_Gr1,d30_Gr1,d84tb_Gr1,d84bs_Gr1,Gr_Gr1]=sht4.range('D4:D10').value
    
    if sht4.range('E4').expand('down').value is not None:
        [dm_Gr2,d50_Gr2,d90_Gr2,d30_Gr2,d84tb_Gr2,d84bs_Gr2,Gr_Gr2]=sht4.range('E4').expand('down').value
    if sht4.range('F4').expand('down').value is not None:
        [dm_Gr3,d50_Gr3,d90_Gr3,d30_Gr3,d84tb_Gr3,d84bs_Gr3,Gr_Gr3]=sht4.range('F4').expand('down').value
    if sht4.range('G4').expand('down').value is not None:
        [dm_Gr4,d50_Gr4,d90_Gr4,d30_Gr4,d84tb_Gr4,d84bs_Gr4,Gr_Gr4]=sht4.range('G4').expand('down').value
    
    
    ##### DICTIONNAIRe Fonction TS   
    fndict = {"LEF2015": TS.Qslef2015, "LEF1991":  TS.Qslef1991, "RICK1990":  TS.Qsrick1990, "RICK1991":  TS.Qsrick1991, "MEUN1989": TS.Qsmeunier,
              "PIRECK2017a":TS.Qspitonrecking2017a,"PIRECK2017b":TS.Qspitonrecking2017b,"PIRECK2016":TS.Qspitonrecking2016}
    [Ql,dm,d50,d90,d30,d84tb,d84bs,Gr,J,b]=np.zeros(10)

    
    ############################################
    ####### Profil en long #####################
    ############################################
    #Interpolation des abscisses tous les dx0
    [x,z]=hyd.interpol(x0,z0,dx0)
    zinit=np.interp(x,x0,z0)             #Pour garder en mémoire le profil initial
    zmin=np.interp(x,xmin0,zmin0)
    b=np.interp(x,x0,b0)
    
    if xseuil0 is not None:
        zseuil=np.interp(x,xseuil0,zseuil0)
    
    
    ############################################
    ####### Granulométrie  #####################
    ############################################
    dm=np.ones(len(b))*dm_Gr1
    d50=np.ones(len(b))*d50_Gr1
    d90=np.ones(len(b))*d90_Gr1
    d30=np.ones(len(b))*d30_Gr1
    d84tb=np.ones(len(b))*d84tb_Gr1
    d84bs=np.ones(len(b))*d84bs_Gr1
    Gr=np.ones(len(b))*Gr_Gr1
    if zmingr2 is not None :
        dm[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=dm_Gr2
        d50[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=d50_Gr2
        d90[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=d90_Gr2
        d30[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=d30_Gr2
        d84tb[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=d84tb_Gr2
        d84bs[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=d84bs_Gr2
        Gr[np.argwhere((zinit>zmingr2) & (zinit<zmaxgr2))]=Gr_Gr2
        
    if zmingr3 is not None :
        dm[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=dm_Gr3
        d50[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=d50_Gr3
        d90[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=d90_Gr3
        d30[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=d30_Gr3
        d84tb[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=d84tb_Gr3
        d84bs[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=d84bs_Gr3
        Gr[np.argwhere((zinit>zmingr3) & (zinit<zmaxgr3))]=Gr_Gr3
        
    if zmingr4 is not None :
        dm[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=dm_Gr4
        d50[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=d50_Gr4
        d90[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=d90_Gr4
        d30[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=d30_Gr4
        d84tb[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=d84tb_Gr4
        d84bs[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=d84bs_Gr4
        Gr[np.argwhere((zinit>zmingr4) & (zinit<zmaxgr4))]=Gr_Gr4
    
    
    ############################################
    ####### Hydrogramme ########################
    ############################################
    tli=np.arange(0,max(t0)*3600.+dt,dt)
    Q=np.interp(tli,np.array(t0)*3600.,Q0)
       
    ############################################
    ####### Injection Transport solide ######## 
    ############################################
    Qsin=np.interp(tli,np.array(tsin0)*3600.,Qsin0)
       
    ##############################################
    ### Initialisation de vecteuret variables ####
    ##############################################
    dx = (np.diff(x)) # b
    Qsout=np.zeros(len(Q))
    tsaved=[]
    Qsinsaved=0
    elapsed=0.
    Qs,Qsmin,Qsmin1,Qssaved,J,jborda=(np.zeros(len(dx)) for i in range(6))
    Vs,h,u,H,zsaved,hsaved,Hsaved=(np.array(zinit) for i in range(7))
    Bsaved=copy.deepcopy(b)
    Vsii=np.array(zinit)
    Vsi=np.array(zinit)
    Vsi=0
    Vs=0
    Vss=[]
    
    ###########################################
    ###########       CALCULS    ##############
    ###########################################
    start_time = timeit.default_timer()
    for j in range(len(Q)):
           
        ##### Calcul de la hauteur, vitesse, charge à chaque noeud
        h[0:]=(Q[j]/b[0:]/9.81**0.5)**(2./3.)                   # Hauteur critique
        u[0:]=(9.81*h[0:])**0.5                                 # Vitesse critique
        H=h+u**2/2/9.81+z                                       # Charge
        jborda=kborda*abs(u[1:]**2/2/9.81-u[0:-1]**2/2/9.81)    # Perte de charge singulière à chaque noeud
        H[1:]=H[1:]+jborda                                      # Charge + perte de charges sing.

        for i in range(1,len(H)):                               # Modification de la ligne de cjarge de telle façon que
            if H[i]<H[i-1]:                                     # la charge soit strictement décroissante avec pente min =0
                H[i]=H[i-1]
        J=(H[1:]-H[0:-1])/dx[0:]                                # Recalcul de la pente de la ligne de charge


        ##### Calcul du débit solide potentiel
        for k in range(len(J)):
            argsdict = {"LEF2015": [Q[j],dm[k],Gr[k],(b[k]+b[k+1])/2,J[k]], "LEF1991":[Q[j],d90[k],d30[k],dm[k],J[k]], "RICK1990": [Q[j],d90[k],d30[k],d50[k],J[k],(b[k]+b[k+1])/2],
                                    "RICK1991": [Q[j],d50[k],J[k],(b[k]+b[k+1])/2],"MEUN1989":[Q[j],J[k]], "PIRECK2017a":[Q[j],d84tb[k],d84bs[k],(b[k]+b[k+1])/2,J[k]],
                                    "PIRECK2017b":[Q[j],d84tb[k],d84bs[k],d50[k],(b[k]+b[k+1])/2,J[k]],"PIRECK2016":[Q[j],d84bs[k],(b[k]+b[k+1])/2,J[k]]}     
            Qs[k]=fndict[loi](*argsdict[loi])

        ##### Limitation du transport soilde avec profil minimum (radier, substratum) défini par zmin
        ##### NE PAS TOUCHER !!!!
        Qsmin1[-1]=Qsin[j]+(z[-1]-zmin[-1])/2/dt*(1.75*b[-1]+0.25*b[-2])
        Qsmin1[0:-1]=Qsmin[1:]+(z[1:-1]-zmin[1:-1])/2/dt*((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1])  
        Qs=np.minimum(Qs,Qsmin1)
           
        ##### Condition limite aval transport solide   
        Qsout[j]=Qs[0]
    
        
        ##### Gestion d'éventuels seuils. Niveaux de crête de seuil définis par zseuil
        try:
            a=(zinit-zseuil)
            indseuil=np.argwhere(abs(a)<0.01)
            
            for k in range(len(indseuil)):
                if z[int(indseuil[k])-1]<zinit[int(indseuil[k])]-(z[int(indseuil[k])+1]-zinit[int(indseuil[k])])/dx[int(indseuil[k])]*dx[int(indseuil[k])-1]: #+(Qs[26]-Qs[25])/b/dx[25]*dt:
                    Qs[int(indseuil[k])-1]=Qs[int(indseuil[k])]
        except:
            pass
        
       ##### Calcul évolution altitude au niveau de chaque noeud (Exner)
       ##### NE PAS TOUCHER !!!!    
        zold=copy.deepcopy(z)
        z[0]=z[0]+(Qs[0]-Qsout[j])*dt*2/dx[0]/(1.75*b[0]+0.25*b[1])
        z[1:-1]=z[1:-1]+(Qs[1:]-Qs[0:-1])*dt*2/((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1]) 
        z[-1]=z[-1]+(Qsin[j]-Qs[-1])*dt*2/dx[-1]/(1.75*b[-1]+0.25*b[-2])
    

       #Enregistrement de certaines variables tous les 10 pas de temps pour animation (voir code plus bas à adapter)
        if j%round(tenr/dt,-1)==0:
            tsaved=np.append(tsaved,tli[j]/3600.)
            if j>0:
                hsaved=np.vstack([hsaved, h])
                zsaved=np.vstack([zsaved, z])
                Bsaved=np.vstack([Bsaved, b])
                Qsinsaved=np.vstack([Qsinsaved,Qsin[j]])
                Hsaved=np.vstack([Hsaved, H])
                Qssaved=np.vstack([Qssaved, Qs])

        Vsii[0]=(z[0]-zold[0])*dx[0]/2*(1.75*b[0]+0.25*b[1])
        Vsii[-1]=(z[-1]-zold[-1])*dx[-1]/2*(1.75*b[-1]+0.25*b[-2])
        Vsii[1:-1]=(z[1:-1]-zold[1:-1])*((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1])/2
        Vsi=Vsi+Vsii
        Vs=sum(Vsi)

    #Vérification bilan de masse
    #NE PAS TOUCHER !!!!  
    Vin=sum(Qsin)*dt                        # Volume sédiments entrant
    Vout=sum(Qsout)*dt                      # Volume sédiment sortant
    
    Vsim=sum([i for i in Vsi if i<0])       # Volume érodé au sein du modèle
    Vsip=sum([i for i in Vsi if i>0])       # Volume déposé au sein du modèe

    
    ##### Ecriture résultats
    ##### Temps de calcul
    elapsed = timeit.default_timer() - start_time             # Résulat du temps de calcul
    sht1.range('B16').options(transpose=True).value=elapsed   # Ecriture temps de calcul dans Excel
    
    ##### Zfond
    sht5.clear_contents()
    sht5.range('b4').options(transpose=True).value=zsaved
    sht5.range('a4').options(transpose=True).value=x
    sht5.range('B2').value=tsaved
    sht5.range('b3').value=['z(t='+str(round(t,2))+' h)' for t in tsaved]
    
    ##### Zfond + hc
    sht9.clear_contents()
    sht9.range('b4').options(transpose=True).value=zsaved+hsaved
    sht9.range('a4').options(transpose=True).value=x
    sht9.range('B2').value=tsaved
    sht9.range('b3').value=['z+hc(t='+str(round(t,2))+' h)' for t in tsaved]
    
    ##### Zfond + Hsc
    sht10.clear_contents()
    sht10.range('b4').options(transpose=True).value=Hsaved
    sht10.range('a4').options(transpose=True).value=x
    sht10.range('B2').value=tsaved
    sht10.range('b3').value=['z+Hsc(t='+str(round(t,2))+' h)' for t in tsaved] 

    ##### B
    sht6.clear_contents()
    sht6.range('b4').options(transpose=True).value=Bsaved
    sht6.range('a4').options(transpose=True).value=x
    sht6.range('B2').value=tsaved
    sht6.range('b3').value=['B(t='+str(round(t,2))+' h)' for t in tsaved]
    
    ##### Vol
    sht11.clear_contents()
    sht11.range('a5').options(transpose=True).value=x
    sht11.range('b5').options(transpose=True).value=Vsi
    sht11.range('a1').options(transpose=True).value=['Somme Er','Somme Dep']
    sht11.range('b1').options(transpose=True).value=[Vsim,Vsip]
    sht11.range('a4').value=['X','Delta V']

    sht1.range('B11').options(transpose=True).value=[Vin,Vs,Vout,Vin-Vout-Vs]  # Bilan de vol. dans feuille calcul


    ##### Figure de résulats (Feuille calcul)
    ##### Annotation sur le graphique
    annotation_string = r'Vol. sédiments entrant : %.1f m$^{3}$' %(Vin)
    annotation_string += "\n"
    annotation_string += r'Vol. sédiments stocké  : %.1f m$^{3}$' %(Vs)
    annotation_string += "\n"
    annotation_string += r'Vol. sédiments sortant : %.1f m$^{3}$' %(Vout)
    annotation_string += "\n"
    annotation_string += r'Bilan de volume            : %.1f m$^{3}$' %(Vin-Vout-Vs)
    
    ##### Définition figure
    fig=plt.figure(figsize=(16,10))
    gs = gridspec.GridSpec(1, 1)
    ax1=plt.subplot(gs[0])
    ax1.grid()    
    ax1.plot(x,zinit,'-k',label=u'Fond t=0s')
    ax1.plot(x,z,'r--',label=u'Fond fin de crue',linewidth=2)
    ax1.plot(x,zmin,'--m',linewidth=0.5,label=u'Profil minimum')
    ax1.set_xlim(min(x),max(x))
    ax1.legend(loc=2,fontsize=14)
    ax1.set_xlabel('Abscisse (m)',fontsize=18)
    ax1.set_ylabel('Altitude (m)',fontsize=18)

    ax2=ax1.twinx()
    ax2.set_ylabel(r'z$_{final}$-z$_{initial}$ (m)',fontsize=18)
    ax2.plot(x,z-zinit,'o',color='grey',alpha=0.5,markersize=6,label=u'Z final - Z initial')

    ###### Légende
    ax1.legend(loc=2,fontsize=14)

    ###### Tracé de l'annotation
    fig.text(0.88, 0.2,annotation_string,
             horizontalalignment='right',style='normal',bbox={'facecolor':'b', 'alpha':0.2, 'pad':10},
             multialignment='left',transform=ax1.transAxes,fontsize=14)
    

    ###### Paramètres de figure additionnels
    plt.rcParams.update({'font.size': 18})
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 18)
    plt.tight_layout()
    
   
    ###### Enregistrement de la figure
    directory=os.getcwd()
    fig.savefig(os.path.dirname(directory)+'\Profil.jpg',dpi=150)
    sht1.pictures.add(os.path.dirname(directory)+'\Profil.jpg',name='Profil',left=sht.range('E1').left,
                            top=sht.range('E1').top,height=16*30,width=24*30,update=True)


    
    

    
    