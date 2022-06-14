from sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
import math

class Lefort1991(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs(self, section, Q, y, y_down):
        Q=max(Q,0.001)
        if section.is_downstream():
            I = section.get_Sf(Q, y)
        else:
            I = (section.get_H(Q, y) - section.get_down_section().get_H(Q, y_down)) / (section.get_down_section().get_x() - section.get_x())
        I=min(0.83,max(I,0.001))
        # I=max(section.get_S0(),0.001)
        b = 0.5*(section.get_b(y) + section.get_down_section().get_b(y_down))
        dm = section.get_granulometry().dm
        d90 = section.get_granulometry().d90
        d30 = section.get_granulometry().d30

        Qcrit=0.0776*(9.81*dm**5)**0.5*1.65**(8./3.)/I**(13./6.)*(1-1.2*I)**(8./3.)
        Qs=4.45*Q*math.pow(I,1.5)/1.65*math.pow(d90/d30,0.2)*(1-math.pow((Qcrit/Q),0.375))
        return max(Qs, 0.01) # return Qs a la base


    def __str__(self):
        return "Lefort 1991"