#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import formules_TS as TS

def hydrogrammeLavabre(Qmax,tm,alpha,Qbase,t):
    #t=np.arange(0,d+dt,dt)
    Q=np.array([(Qmax-Qbase)*2*np.power(ti/tm,alpha)/(1+np.power(ti/tm,2*alpha))+Qbase for ti in t])
    Q=[0.0001 if x==0 else x for x in Q]
    return np.array(Q)


def main():
    
    # CONSTANTE

    g = 9.81
    kborda = 0.5

    # TIME GRID
    d = 24*3600
    dt = 10
    
    t=np.arange(0,d+dt,dt)
    lenT = len(t)

    # PROFIL BIEF
    xParam = np.array([0, 295, 300, 350, 355, 700])
    zInitParam = 0.02*xParam
    zMinParam = zInitParam-2
    bParam = np.array([10, 10, 5, 5, 10, 10])
    stepX = 50

    x = []
    for i in range(len(xParam)-1):
        xInter = np.arange(xParam[i], xParam[i+1], stepX)
        x += list(xInter)
    x.append(xParam[-1])
    zInit = np.interp(x, xParam, zInitParam)
    zMin = np.interp(x, xParam, zMinParam)
    b = np.interp(x, xParam, bParam)
    dx = np.diff(x)
    lenX = len(x)

    # HYDROGRAMME 
    Qmax = 50
    tm = 8*3600
    alpha = 3
    Qbase = 1
    
    Q = hydrogrammeLavabre(Qmax, tm, alpha, Qbase, t)

    # GRANULOMÉTRIE
     
    dm=np.ones(lenX)*0.065
    d50=np.ones(lenX)*0.03
    d90=np.ones(lenX)*0.09
    d30=np.ones(lenX)*0.01
    d84tb=np.ones(lenX)*0.04
    d84bs=np.ones(lenX)*2.2
    Gr=np.ones(lenX)*2.5

    # SEDIMENTOGRAMME
    loi = "MEUN1989"
    b0 = 10
    I0 = 0.02
    fndict = {"LEF2015": TS.Qslef2015, "LEF1991":  TS.Qslef1991, "RICK1990":  TS.Qsrick1990, "RICK1991":  TS.Qsrick1991, "MEUN1989": TS.Qsmeunier,
              "PIRECK2017a":TS.Qspitonrecking2017a,"PIRECK2017b":TS.Qspitonrecking2017b,"PIRECK2016":TS.Qspitonrecking2016}
    argsdict = {"LEF2015": [Q,dm,Gr,b0,I0], "LEF1991":[Q,d90,d30,dm,I0], "RICK1990": [Q,d90,d30,d50,I0,b0],"RICK1991": [Q,d50,I0,b0],"MEUN1989":[Q,I0],
                "PIRECK2017a":[Q,d84tb,d84bs,b0,I0],"PIRECK2017b":[Q,d84tb,d84bs,d50,b0,I0],"PIRECK2016":[Q,d84bs,b0,I0]}
    Qsin = np.zeros(len(Q))
    for i in range(len(Q)):
        args = argsdict[loi]
        args[0] = Q[i]
        Qsin[i] = fndict[loi](*args)

    # Nommage des tailles de la grille 
    lenT = len(t)
    lenX = len(x)

    Qs = np.zeros(lenX-1)
    QsMin = np.zeros(lenX-1)
    Qsout = np.zeros(lenT)
    Vd = np.zeros(lenX)
    VdTot = 0
    z = np.copy(zInit)

    for ti in range(lenT):
        ##### Calcul de la hauteur, vitesse, charge à chaque noeud
        h=(Q[ti]/(b*np.sqrt(g)))**(2./3.)                       # Hauteur critique
        u=np.sqrt(g*h)                                          # Vitesse critique
        H=h+z+u**2/(2*g)                                        # Charge
        jborda=kborda*abs(u[1:]**2-u[0:-1]**2)/(2*g)    # Perte de charge singulière à chaque noeud
        # H[1:]=H[1:]+jborda                                      # Charge + perte de charges sing.

        for i in range(1,lenX):                               # Modification de la ligne de cjarge de telle façon que
            if H[i]<H[i-1]:                                     # la charge soit strictement décroissante avec pente min =0
                H[i]=H[i-1]
        J=(np.diff(H))/dx                                # Recalcul de la pente de la ligne de charge


        ##### Calcul du débit solide potentiel
        for k in range(len(J)):
            argsdict = {"LEF2015": [Q[ti],dm[k],Gr[k],(b[k]+b[k+1])/2,J[k]], "LEF1991":[Q[ti],d90[k],d30[k],dm[k],J[k]], "RICK1990": [Q[ti],d90[k],d30[k],d50[k],J[k],(b[k]+b[k+1])/2],
                                    "RICK1991": [Q[ti],d50[k],J[k],(b[k]+b[k+1])/2],"MEUN1989":[Q[ti],J[k]], "PIRECK2017a":[Q[ti],d84tb[k],d84bs[k],(b[k]+b[k+1])/2,J[k]],
                                    "PIRECK2017b":[Q[ti],d84tb[k],d84bs[k],d50[k],(b[k]+b[k+1])/2,J[k]],"PIRECK2016":[Q[ti],d84bs[k],(b[k]+b[k+1])/2,J[k]]}     
            Qs[k]=fndict[loi](*argsdict[loi])

        ##### Limitation du transport soide avec profil minimum (radier, substratum) défini par zmin
        QsMin[-1]=Qsin[ti]+(z[-1]-zMin[-1])/2/dt*(1.75*b[-1]+0.25*b[-2])
        QsMin[0:-1]=+(z[1:-1]-zMin[1:-1])/2/dt*((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1])  
        Qs=np.minimum(Qs,QsMin)
           
        ##### Condition limite aval transport solide   
        Qsout[ti]=Qs[0]
        
        ##### Calcul évolution altitude au niveau de chaque noeud (Exner)
        ##### NE PAS TOUCHER !!!!    
        zold=np.copy(z)
        z[0]=z[0]+(Qs[0]-Qsout[ti])*dt*2/dx[0]/(1.75*b[0]+0.25*b[1])
        z[1:-1]=z[1:-1]+(Qs[1:]-Qs[0:-1])*dt*2/((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1]) 
        z[-1]=z[-1]+(Qsin[ti]-Qs[-1])*dt*2/dx[-1]/(1.75*b[-1]+0.25*b[-2])

        Vd[0]=(z[0]-zold[0])*dx[0]/2*(1.75*b[0]+0.25*b[1])
        Vd[-1]=(z[-1]-zold[-1])*dx[-1]/2*(1.75*b[-1]+0.25*b[-2])
        Vd[1:-1]=(z[1:-1]-zold[1:-1])*((0.75*b[1:-1]+0.25*b[2:])*dx[1:]+(0.75*b[1:-1]+0.25*b[0:-2])*dx[0:-1])/2
        VdTot += sum(Vd)

    
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(t/3600, Q, 'b-')
    ax2.plot(t/3600, Qsin, 'r-')

    ax1.set_xlabel('temps')
    ax1.set_ylabel('débit liquide', color='b')
    ax2.set_ylabel('débit solide', color='r')

    print("Volume entrant = " + str(sum(Qsin*dt)))
    print("Volume sortant = " + str(sum(Qsout*dt)))
    print("Volume déposé = " + str(VdTot))
    print("Bilan de volume = " + str(sum(Qsin*dt)-sum(Qsout*dt)-VdTot))
    plt.figure()
    plt.plot(x, zInit)
    plt.plot(x, z)
    plt.show()

    return


if __name__=="__main__":
    main()