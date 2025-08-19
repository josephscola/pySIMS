import mendeleev

_ALPHA = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
_NUM = "0123456789"

def read_isotope_ref(ref):
    int_mass = int(ref.rstrip(_ALPHA))
    name = ref.lstrip(_NUM)
    return int_mass, name

def get_isotope_abundance(ref):
    int_mass, name = read_isotope_ref(ref)
    isotope = mendeleev.isotope(name, int_mass)
    return isotope.abundance    

def get_minors_isotopes(ref):
    int_mass, name = read_isotope_ref(ref)
    elem = mendeleev.element(name)
    
    minors = []
    for iso in elem.isotopes:
        if iso.abundance and iso.mass_number != int_mass:
            minors.append(iso)
    return minors


    
    
