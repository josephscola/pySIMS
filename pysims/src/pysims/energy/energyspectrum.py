#====================================================================#
#                                                                    #
#      This file contains the code for energy spectrum analysis      #
#                                                                    #
#====================================================================#

from pysims.datamodel import Crater
from pysims.utils import * 

class EnergySpectrum(Crater):
    """
    Main class used for energy analysis, contains the raw data and the
    results of the processings applied on data.

    :param path: the path of the depthprofile ascii file
    :type path: str
    """
    
    @property
    def energy(self, elem: str):
        """
        Getter method to acces the energy values for an element

        :param elem: the element
        :type elem: str

        :rtype: list
        """
        return self._get_elem_attr(elem, ENERGY)
    
    @property
    def intensity(self, elem: str):
        """
        Getter method to acces the intensity values for an element

        :param elem: the element
        :type elem: str

        :rtype: list
        """
        return self._get_elem_attr(elem, INTENSITY)
