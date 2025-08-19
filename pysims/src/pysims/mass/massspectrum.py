#====================================================================#
#                                                                    #
#       This file contains the code for mass spectrum analysis       #
#                                                                    #
#====================================================================#

from pprint import pprint
import numpy as np

from pysims.datamodel import Crater
from pysims.utils import * 

from .isotopes import *


class MassSpectrum(Crater):
    def __init__(self, path: str):
        Crater.__init__(self,path)
        mass = self._raw_data["1"][MASS] + self._raw_data["101"][MASS]
        intens = self._raw_data["1"][INTENSITY] + self._raw_data["101"][INTENSITY]
        self._raw_data[MASS] = mass 
        self._raw_data[INTENSITY] = intens

    @property
    def mass(self) :
        return self._get_attr(MASS)
    
    @property
    def intensity(self) :
        return self._get_attr(INTENSITY)

    def local_max(self, m_ref, n=16, DEBUG=False):
        """
        Calculate the local max intensity around a reference mass
        (+/- 0.5 at. unit)

        :param m_ref: the reference mass
        :type m_ref: int

        :param n: local intervall size
        :type n: int

        :param DEBUG: if set to true, prints additionnal debug infos
        :type DEBUG: bool

        :return: local max intensity and associated mass
        :rtype: tuple
        """
    
        idx = int(self.mass.index(m_ref) - n / 2)
        local_spectrum = self.intensity[idx : idx + n]
        localmax_i = max(local_spectrum)
        localmax_m = self.mass[local_spectrum.index(localmax_i) + idx]
    
        if DEBUG:
            print(f"Integer mass {m_ref} located at index {idx}")
            print(f"exp. intensity local max. intensity {localmax_i}\n")
        return(localmax_m, localmax_i)


    def deviation_to_natural_abundance(self, ref: str, relevance_threshold=100, n=16, DEBUG=False):
        """
        Calculate the standard deviation of experimental mass spectrum
        to natural abundance for a specific element

        :param ref: the reference element
        :type ref: str

        :param relevance_threshold: threshold below which the
            experimental intensity is not relevant
        :type relevance_threshold: float | int

        :param n: local intervall size for local max evaluation
        :type n: int

        :param DEBUG: if set, prints additionnal debug infos
        :type DEBUG: bool
        """
        int_mass, elem = read_isotope_ref(ref)
    
        # reference abundance
        ref_intens_theo = get_isotope_abundance(ref)
        ref_int_mass_expe, ref_intens_expe = self.local_max(int_mass,
                                                            n=n,
                                                            DEBUG=DEBUG)
        if DEBUG:
            print(ref, ref_intens_theo, ref_intens_expe)

        # test intensity relevance
        if ref_intens_expe <= relevance_threshold:
            return -1 

        minors = get_minors_isotopes(ref)
        if DEBUG:
            pprint(minors)
            
        var_ionic_interference = 0
        for iso in minors:
            minor_int_mass_theo = np.round(iso.mass)
            minor_intens_theo = iso.abundance
            minor_int_mass_expe, minor_intens_expe = self.local_max(
                minor_int_mass_theo,
                n=n,
                DEBUG=DEBUG
            )
            r_theo = minor_intens_theo / ref_intens_theo
            r_expe = minor_intens_expe / ref_intens_expe
            var_ionic_interference += ((r_expe - r_theo) / r_theo) ** 2
            
            if DEBUG:
                print(str(iso.mass_number)+elem, minor_intens_theo, minor_intens_expe)
                print(f" : r_theo = {r_theo:.2e}" + f" : r_expe = {r_expe:.2e}")
                print("relative err : ", abs(r_expe - r_theo) / r_theo, "\n")
                
        sigma_ionic_interference = np.sqrt(var_ionic_interference / len(minors))
        if DEBUG:
            print(f"std dev = {sigma_ionic_interference * 100:.1f} %")
            
        return sigma_ionic_interference
