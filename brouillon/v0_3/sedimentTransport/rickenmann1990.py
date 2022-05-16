from sedimentTransport.sedimentTransportLaw import SedimentTransportLaw

class Rickenmann1990(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs(self, section, Q, y, y_down):
        Q=max(Q,0.001)
        if section.is_downstream():
            I = section.get_Sf(Q, y)
        else:
            I = (section.get_H(Q, y) - section.get_down_section().get_H(Q, y_down)) / (section.get_down_section().get_x() - section.get_x())
        I=max(I,0.001)
        # I=max(section.get_S0(),0.001)
        b = 0.5*(section.get_b(y) + section.get_down_section().get_b(y_down))
        d30 = section.get_granulometry().d30
        d50 = section.get_granulometry().d50
        d90 = section.get_granulometry().d90
        q=Q/b
        qcr=0.065*1.65**1.67*9.81**0.5*d50**1.5*I**(-1.12)        
        qs=12.6*(d90/d30)**0.2*(q-qcr)*I**2.*1.65**(-1.5)
        Qs=qs*b/0.75 
        Qs=max(Qs,0.)
        return Qs

    def __str__(self):
        return "Rickenmann 1990"