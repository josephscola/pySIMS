#====================================================================#
#                                                                    #
# This file contains the code for multilayer depth profiles analysis #
#                                                                    #
#====================================================================#

import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d

from pysims.utils import *
from pysims.datamodel import Crater


class DepthProfiles(Crater) :
    """
    Main class used for multilayer analysis, contains the raw data and
    the results of the processings applied on data.

    :param path: the path of the depthprofile ascii file
    :type path: str
    """
    def __init__(self, path):
        super().__init__(path)
        self._properties = {
            "plateaux indices": {},
            "plateaux" : {},
            "standard deviation": {},
            "interfaces": {},
            "ideal time": {},
            "ideal depth": {},
            "ideal intens in time": {},
            "ideal intens in depth": {},
            "calculated depth": {}
        }

    def get_list_elem(self) -> list:
        """
        Returns the list of the elements that are in the analysis

        :rtype: list
        """
        return list(self.data.keys())
    
    def intensity(self, elem: str) -> list:
        """
        Getter method to acces the intensities of an element's
        profile

        :param elem: the element
        :type elem: str

        :rtype: list
        """
        return self._get_elem_attr(elem, INTENSITY)

    def time(self, elem: str) -> list:
        """
        Getter method to acces the times of an element's profile

        :param elem: the element
        :type elem: str

        :rtype: list
        """
        return self._get_elem_attr(elem, TIME)

    def depth(self, elem: str) -> list:
        """
        Getter method to acces the depths of an element's profile

        :param elem: the element
        :type elem: str

        :rtype: list
        """
        if elem in self.properties["calculated depth"]:
            return self.properties["calculated depth"][elem]
        else:
            return self._get_elem_attr(elem, DEPTH)

    #====================== Plateaux detection ======================#

    def locate_interfaces(
                self,
                elem: str,
                prominence: float = 0.2
    ) -> list:
        """
        Detects interfaces between plateaux.  The method used is to
        filter the noise, calculate the gradient of the filtered
        intensity and detect it's peaks.

        :param elem: the element which interfaces we want to locate
        :type elem: str

        :param prominence: sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float

        :rtype: list
        """
        
        y = self.intensity(elem)
        # smoothed_y = uniform_filter1d(y, size=noise_filter_size)
        gradient = normalize(np.abs(np.gradient(y)))
        interfaces = list(find_peaks(gradient, prominence=prominence, distance=10)[0])
        
        self.properties["interfaces"][elem] = interfaces
        return interfaces
    
    def get_plateaux_indices(
                self,
                elem: str,
                interfaces_margin: int,
                prominence: float = 0.2
    ) -> list:
        """
        Returns the indices of the limits of all plateaux for the
        given element, +/- the interfaces_margin.  The results is a
        list of tuples containing the start and end indices of each
        plateau.

        .. note::

            ``prominence`` is only used for interfaces detection

        :param elem: the element
        :type elem: str

        :param interfaces_margin: the margin the apply to the detected
            interfaces
        :type interfaces_margin: int

        :param prominence: sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float

        :rtype: list
        """
        
        interfaces = self.locate_interfaces(elem, prominence)
        interfaces.append(len(self.intensity(elem)) - 1)

        list_indices = []
        for i in range(len(interfaces) -1):
            start = interfaces[i] + interfaces_margin
            end = interfaces[i+1] - interfaces_margin
            list_indices.append((start, end))

        
        self.properties["plateaux indices"][elem] = list_indices
        return list_indices

    def get_plateaux(
                self,
                elem: str,
                interfaces_margin: int,
                prominence: float = 0.2
    ) -> list:
        """
        Calculate the value of each plateau for the given element

        .. note::

            ``prominence`` is only used for interfaces detection

        :param elem: the element
        :type elem: str

        :param interfaces_margin: the margin to apply to the detected
            plateaux
        :type interfaces_margin: int

        :param prominence: sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float

        :rtype: list
        """
        
        plateaux_indices = self.get_plateaux_indices(elem,
                                                     interfaces_margin,
                                                     prominence)
        
        plateaux = []
        for indices in plateaux_indices:
            plateaux.append(
                calculate_plateau_value(self.intensity(elem), indices)
            )
        self.properties["plateaux"][elem] = plateaux
        return plateaux

    def get_plateaux_std(
                self,
                elem: str,
                interfaces_margin: int,
                prominence: float = 0.2
    ) -> list:
        """
        Calculate the standard deviation of each plateau for the given
        element

        .. note::

            ``prominence`` is only used for interfaces/plateaux detection

        :param elem: the element
        :type elem: str

        :param interfaces_margin: the margin to apply to the detected
            plateaux
        :type interfaces_margin: int

        :param prominence:  sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float

        :rtype: list
        """
        
        std = []
        plateaux_indices = self.get_plateaux_indices(elem,
                                                     interfaces_margin,
                                                     prominence)
        for indices in plateaux_indices:
            std.append(
                calculate_std(self.intensity(elem), indices)
            )
        self.properties["standard deviation"][elem] = std
        return std

    #=================== Ideal profile generation ===================#

    def calculate_profile_depth(
                self,
                elem: str,
                mining_speeds: list,
                min_interface_idx: int = 20,
                prominence: float = 0.2
    ) -> list:
        """
        Given a list containing the mining speed for each layers,
        calculates the depths for the given element.

        :param elem: the element
        :type elem: str

        :param mining_speeds: the list of the mining speeds of each
            layers
        :type mining_speeds: list

        :param min_interface_idx: minimum interface index below which
            it is considered to be at 0
        :type min_interface_idx: int

        :param prominence: sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float
        """

        time = self.time(elem)
        interfaces = self.locate_interfaces(elem, prominence)
        
        if interfaces[0] < min_interface_idx:
            interfaces[0] = 1
        else:
            interfaces.insert(0,1)
        interfaces.append(len(time)-1)

        depth = [0]
        for inter in range(len(interfaces) - 1):
            for i in range(interfaces[inter], interfaces[inter+1]+1):
                depth.append(depth[-1] + mining_speeds[inter] * (time[i] - time[i-1]))

        self.properties["calculated depth"][elem] = depth
        return depth
    
    def generate_ideal_profile(
                self,
                profile_type: str,
                elem: str,
                n_ideal: int, 
                shift: float = 0,
                interfaces_margin: int = 5,
                cancellation_thresh: float = 1e-3,
                prominence: float = 0.2
    ) -> tuple:
        """
        Generate the ideal profile intensity, and time or depth for the
        given element.  Returns the tuple (Tideal|Dideal, Iideal)

        .. note::

            ``prominence`` is only used for interfaces/plateaux detection

        :param profile_type: the type of profile to generate.  Must be
            either ``'time'``or ``'depth'``
        :type profile_type: str

        :param elem: the element
        :type elem: str

        :param n_ideal: the ideal profile resolution
        :type n_ideal: int

        :param shift: factor by which the experimental
            values will be shifted.

        :type shift_Time: float

        :param interfaces_margin: the margin the apply to the detected
            interfaces
        :type interfaces_margin: int

        :param cancellation_thresh: threshold bellow which a plateau
            is considered to be null
        :type cancellation_thresh: float

        :param min_plateau_width: minimum width of a plateau, used to
            clamp the plateaux to the limits of the time/depth region
        :type min_plateau_width: int

        :param prominence: sensitivity of interface detection. Must be
            between 0 and 1, lower value means more sensitivity.
        :type prominence: float

        :rtype: tuple
        """

        interfaces = self.locate_interfaces(elem, prominence)
        plateaux = self.get_plateaux(elem, interfaces_margin, prominence)
        
        exp = self.__getattribute__(profile_type)(elem)
        ideal = np.linspace(exp[0], exp[-1], n_ideal)

        # get indices of interfaces in Tideal
        idx_ideal = get_ideal_interface_indices(interfaces, exp, ideal, shift)
        idx_ideal.append(n_ideal) # add end indices to fully delimit the plateaux
        
        # calculate ideal intensity
        Iideal = np.ones(n_ideal)
        for i, plateau in enumerate(plateaux): 
            if plateau > cancellation_thresh:
                start = idx_ideal[i]
                end = idx_ideal[i+1]
                Iideal[start : end] *= plateau

        self.properties["ideal " + profile_type][elem] = ideal
        self.properties["ideal intens in " + profile_type][elem] = Iideal
        return ideal, Iideal


    
#===================== calculus/utils functions =====================#


def normalize(array: np.ndarray) -> np.ndarray:
    """
    normalize a numpy array

    :param array: array to be normalized
    :type array: np.ndarray

    :rtype: np.ndarray
    """
    m = np.min(array)
    M = np.max(array)
    return (array - m)/(M - m)


def calculate_plateau_value(array: list, indices: tuple) -> float:
    """
    Calculate the mean of the array values between the given indices

    :param array: a list of float or int values
    :type array: list

    :param indices: the indices of the plateaux
    :type indices: tuple

    :rtype: float
    """
    return np.mean(array[indices[0] : indices[1]])


def calculate_std(array: list, indices: tuple) -> float:
    """
    Calculate the standard deviation of the intensity between the given indices

    :param array: a list of float or int values
    :type array: list

    :param indices: the indices of the plateaux
    :type indices: tuple

    :rtype: float
    """
    return np.std(array[indices[0] : indices[1]])



#============= Ideal Profiles Generation utils functions ============#

def index_of_closest_element(value: float, array: list) -> int:
    """
    Return the index of the closest element to value in the given
    array

    :param value: the value which index is to be found
    :type value: float

    :param array: a list of float
    :type array: list
    """
    return min(range(len(array)), key=lambda i: abs(array[i] - value))


def get_ideal_interface_indices(
    interfaces: list,
    exp: list,
    ideal: list,
    shift: float = 0,
) -> list:
    """
    Find the index of the interfaces in the generated ideal profile
    for the given element, with the possibility to shift the
    experimental values if necessary.

    :param exp: the eperimental values
    :type exp: list

    :param ideal: the generated ideal values
    :type ideal: list

    :param interfaces: the list of interfaces indices in the
        experimental data

    :param shift: factor by which the experimental
        values will be shifted
    :type shift: float
    """
    idx_ideal = []
    for i in interfaces:
        idx = index_of_closest_element(exp[i] * (1 + shift), ideal)
        idx_ideal.append(idx)
    return idx_ideal
