from dataclasses import dataclass
from functools import wraps
from typing import Any, Dict
from json import dumps
import tatsu
import importlib.resources

from .semantic import Semantic

# Type Aliases
Data_t = Dict[str, list | Dict[str, list]]


class FileFormatError(Exception):
    pass

    
class Crater:
    """
    This is the generic data model class, which is inherited by
    DepthProfiles, MassSpectrum and EnergySpectrum.  It contains the
    raw data and metadata along with the results of processings
    applied on data.  Those results are stored in the _properties
    dictionnary.

    .. note::

        This class is not meant to be used directly.

    :param path: path to data file (must be in the CAMECA .ms, .dp or
        .nrj ascii file format)
    :type path: str

    :raises FileFormatError: if the given file does not follow the
        right format
    """
    def __init__(self, path: str):
        with open(path, "r", encoding="iso-8859-1") as f:
            raw = f.read()

        grammar = importlib.resources.files("pysims.datamodel").joinpath("grammar.ebnf").read_text(encoding='iso-8859-1')
        parser = tatsu.compile(grammar, semantics = Semantic())

        try:
            extract = parser.parse(raw)
        except tatsu.exceptions.FailedParse as exc:

            err_msg = "Bad file format"            
            err  = FileFormatError(err_msg)
            raise err from exc

        self._raw_data: Data_t = extract["data"]
        self._raw_metadata: Data_t = extract["metadata"]
        self._properties: Dict[str, Any] = {}

    def __str__(self):
        return dumps(self._data)

    @property
    def data(self) -> Data_t:
        """
        Returns a copy of the internal ``_data`` dictionnary.

        :rtype: Data_t
        """
        return self._raw_data.copy()

    @property
    def metadata(self) -> Data_t:
        """
        Returns a copy of the internal ``_metadata`` dictionnary.

        :rtype: Data_t
        """
        return self._raw_metadata.copy()
    
    @property
    def properties(self) -> Dict[str, Any]:
        """
        Returns the internal _properties dictionnary.

        :rtype: Dict[str, Record]
        """
        return self._properties

    def _get_attr(self, attr) -> list | dict:
        """
        Internal helper method to acces data attributes.  The return
        object is a copy of the one stored in the data dictionnary.

        :param attr: the attribute to be accessed
        :type attr: str

        :rtype: list | dict
        """
        return self._raw_data[attr].copy()    
    
    def _get_elem_attr(self, elem, attr) -> list:
        """
        Internal helper method to acces and element's attribute.
        The return object is a copy the one stored in the data
        dictionnary.

        :param elem: the element of data to be accessed
        :type elem: str

        :param attr: the element's attribute to be accessed
        :type attr: str

        :rtype: list
        """
        return self._raw_data[elem][attr].copy()
