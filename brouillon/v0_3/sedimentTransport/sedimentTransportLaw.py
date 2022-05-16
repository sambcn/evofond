from abc import ABC, abstractmethod

class SedimentTransportLaw(ABC):
    """
    This class is like an interface of Object-Oriented programming.
    It is only used to hide all the transport law under an unique class.
    It allows to avoid the 
    It could be done better with the right modules, but in the first place it will work well.
    """
     
    # computation functions

    @abstractmethod
    def compute_Qs(self, section, Q, y, y_next): 
        pass

    # operators overloading

    def __str__(self):
        return "unknown SedimentTransportLaw"
