from mendeleev import element
import re
import itertools
import os
import csv
from math import prod
import matplotlib.ticker as ticker

class Molecularion:
    def __init__(self, name, mass, abundance):
        self.name = name
        self.mass = mass
        self.abundance = abundance

    def __str__(self):
        return f"{self.name},{self.mass},{self.abundance}" 

    def dissociate_constituent_atoms(self):
        return re.findall(r'[A-Za-z]+', self.name)

    def dissociate_constituent_isotopes(self):
        return re.findall(r'([0-9]+[A-Za-z]*)', self.name)
    
def generate_single_isotopes(atoms_selection, UNSTABLE_ISOTOPES = False):
    single_molecularion = []
    for elem in atoms_selection:
        for iso in element(elem).isotopes:
            if iso.abundance is not None:
                single_molecularion.append(Molecularion(f'{iso.mass_number}{elem}', iso.mass, iso.abundance/100))
            if iso.abundance == None and UNSTABLE_ISOTOPES:
                single_molecularion.append(Molecularion(f'{iso.mass_number}{elem}', iso.mass, 1e10))
    return single_molecularion

def generate_molecularion_combination(single_molecularion, multiplicity):
    molecularion_combination = []
    all_combinations = list(itertools.chain.from_iterable(itertools.combinations_with_replacement(single_molecularion, n) for n in range(1, multiplicity+ 1)))
    for comb in all_combinations:
          combined_name = "".join(molion.name for molion in comb)
          combined_mass = sum(molion.mass for molion in comb)
          combined_abundance = prod(molion.abundance for molion in comb)
          molecularion_combination.append(Molecularion(combined_name, combined_mass, combined_abundance))
    return molecularion_combination

def multiplicity_increment(molecularion_list, atoms_selection, multiplicity) :  #multiplicity est l'actuel, objectif +1
    single_isotopes = generate_single_isotopes(atoms_selection)
    molecularion_list_incremented = []
    for molion in molecularion_list:
        for iso in single_isotopes:
            if iso.name not in molion.name:
                molecularion_list_incremented.append(Molecularion(molion.name + iso.name, molion.mass + iso.mass, molion.abundance + iso.abundance))
    return molecularion_list_incremented

def get_weight_abundance_combination(iso_combination):
    weight = 0
    consituent_isoptopes = iso_combination.dissociate_constituent_isotopes()
    for iso in consituent_isoptopes:
        mass_number, symbol = re.match(r"(\d+)([A-Za-z]+)", iso).groups()
        weight += next((iso.mass for iso in element(symbol).isotopes if iso.mass_number == int(mass_number)))
        abundance += next((iso.abundance for iso in element(symbol).isotopes if iso.mass_number == int(mass_number)))
    return weight, abundance

def sort_molecularion_list(molecularion_list, criteria):
    molecularion_list.sort(key=lambda Molecularion: getattr(Molecularion, criteria))
    return molecularion_list


####################################################
###                FILES GESTION                 ###
####################################################

"""
def str2frame(estr, sep = ',', lineterm = '\n', set_header = True):
    dat = [x.split(sep) for x in estr.strip(lineterm).split(lineterm)]
    df = pd.DataFrame(dat)
    if set_header:
        df = df.T.set_index(0, drop = True).T # flip, set ix, flip back
    return df
"""

def molecularion_list_to_table (molecularion_list):
    # table is a list of list ready for DataFrame conversion
    table = [['name', 'mass', 'abundance']]
    for molecularion in molecularion_list:
        table.append ([molecularion.name, molecularion.mass, molecularion.abundance])
    return table

def check_file_atoms_selection(path_dossier,file_name): 
    #permet de vérifier si notre sélection d'atome à déjà un fichier de combinaison existant
    file_find=False
    fichiers = [f for f in os.listdir(path_dossier) if file_name in f]  

    if fichiers:
        print("Fichier(s) trouvé(s) :", fichiers)
        file_find=True
    else:
        print("Aucun fichier trouvé.") 
    return fichiers,file_find  

#permet de vérifier l'existance d'un fichier et d'en vérifier la multiplicité 

def check_file_with_multiplicity (atoms_selection, multiplicity,path_fichier_combinaison): #selection_atomes,multipliciter,path_fichier_combinaison, 
    
     atoms_selection_sorted=file_name_in_correct_order(atoms_selection)#donne le nom du fichier dans le bon format
     liste_fichier,file_find=check_file_atoms_selection(path_fichier_combinaison,atoms_selection_sorted) # liste des noms des fichiers correspondant à notre sélection d atomes
     if file_find==True: # condition d'existance d'un fichier de combinaison avec notre sélection d'atome
         
          file_name_useable=''
          file_usable=False
          
          for file_name in liste_fichier:
               
               file_multiplicity=''
               
               for i in range(len(file_name)):
        
                    try : # vérifie si le caractère dans le nom du fichier est un chiffre
                         if int(file_name[i]) in [i for i in range(10)]:
                              file_multiplicity+=file_name[i]
                         else:
                          break
                    except : None

               if int(file_multiplicity)<=multiplicity: # compare la multiplicité du fichier et celle attendu
                   
                    print(f'fichier adapté à une multiplicité de {multiplicity}: {file_name}')
                    file_usable=True
                    file_name_useable=file_name
                    break
          
          if file_usable==False:
               print('multiplicité demandé trop grande, création du nouveau fichier')
     
     else:
          print('pas de fichiers existant')
          
          
     return  file_name_useable,file_usable
          

# (selection d'atomes, chemin(s) d'accès des fichiers de
# combinaison, multiplicité)
#      - recherche de la présence d'un fichier résultat (nommage, contenu)
#      -> si oui -> affichage


def file_name_in_correct_order(atoms_selection): #donne le nom du fichier dans l'ordre des numéros atomiques
    
    liste_atom_sort_by_atomic_number=list()
    liste_atomic_number=list()
    dico_atom=dict()
    
    for atom in atoms_selection:
        dico_atom[element(atom).atomic_number]=atom
        liste_atomic_number.append(element(atom).atomic_number)
    liste_atomic_number.sort()
    
    for atomic_number in liste_atomic_number:
        
        liste_atom_sort_by_atomic_number.append(dico_atom[atomic_number])
        
    file_name=''
    for atom in liste_atom_sort_by_atomic_number:
        file_name+=str(atom)+'_'
    
    return file_name[:-1]

def file_name_with_multiplicity(atoms_selection, multiplicity): # donne le nom complet du fichier sous le format: multiplicité_atome1_atome2_..._atomeN.csv
    
    return f"{str(multiplicity)}_{file_name_in_correct_order(atoms_selection)}.csv"

########################################################
#AFFICHAGE
########################################################

def add_points(liste_xaxe_point,liste_yaxe_point,liste_label,ax):
    
    for i in range(len(liste_xaxe_point)):
        ax.scatter(liste_xaxe_point[i],liste_yaxe_point[i],label=liste_label[i])
    

def synthax_graph(ax):
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        
    return ax
    
####################################################
###                   TEST                       ###
####################################################


def test_classe():
    toto = Molecularion("O2H8Pg26", "20", "99")
    print(f'la masse de {toto.name} est {toto.mass} et son abondance est {toto.abundance}') 
    print(f'toto est composer des atomes : {toto.dissociate_constituent_atoms()}')
    print(f'toto est composer des isotopes : {toto.dissociate_constituent_isotopes()}')

def test_generate_single_isotopes():
    for elem in  generate_single_isotopes(['Fe','O','Au']):
        print(elem)
        
def test_generate_molecularion_combination():
    single_molecularion = generate_single_isotopes(['Fe','O','Au'])
    for elem in generate_molecularion_combination(single_molecularion, 3):
        print(elem)
        
def test_multiplicity_increment():
    single_molecularion = generate_single_isotopes(['Fe','O','Au'])
    molecularion_list = generate_molecularion_combination(single_molecularion, 2)
    atoms_list = ['O','H','Fe','Si']
    multiplicity = 4
    for elem in multiplicity_increment(molecularion_list, atoms_list, multiplicity):
       print(elem.__str__()) 

def test_sort_molecularion_list():
    single_molecularion = generate_single_isotopes(['Fe','O','Au'])
    molecularion_list = generate_molecularion_combination(single_molecularion, 2)
    for elem in molecularion_list:
        elem.mass = get_weight_abundance_combination(elem)
    criteria = 'mass'
    molecularion_list.sort(key=lambda Molecularion: getattr(Molecularion, criteria))
    for elem in molecularion_list:
        print(elem.__str__())

def test_check_file_atoms_selection(): 
    
    path_dossier=r'./molecularions'
    selection_atome=r'Fe'
    print(check_file_atoms_selection(path_dossier,selection_atome))

def test_check_file_with_multiplicity():
    path_dossier=r'./molecularions'
    atoms_selection=['H','O']
    multiplicity=5
    check_file_with_multiplicity(atoms_selection,multiplicity)
     
def test_read_file_csv(): # répertoire à changer
    
    path_file=r'./molecularions\10_O_Fe_Au.csv'    
    for element in read_csv(path_file):
        print(element.__str__())

def test_write_file_csv():  # execution assez longue, répertoire à modifier 
    multiplicity=4
    path_dossier=r'./molecularions' 
    atoms_selection=['H','Ti']
    single_molecularion = generate_single_isotopes(atoms_selection)
    molecularion_list = generate_molecularion_combination(single_molecularion, multiplicity)
    write_csv(molecularion_list ,path_dossier+file_name_with_multiplicity(atoms_selection,multiplicity),file_name_with_multiplicity(atoms_selection,multiplicity))

def test_read_write():
    try :
        liste_molecularion=read_csv(r'./molecularions\2_O_Fe.csv')
        print('lecture ok')
    except:
        'lecture_echec'
    try :
        write_csv(liste_molecularion,r'./molecularions\test_2_O_Fe.csv','test_2_O_Fe.csv')
        print('ecriture ok')
    except:
        print('echec ecriture')
        
    try :
        liste_test=read_csv(r'./molecularions\2_O_Fe.csv')
        for elem in liste_test:
            print(elem.__str__())
    except:
        print('echec lecture')
        
##########################################
#   Execution des fonctions de test
###########################################

#test_classe()
#test_generate_single_isotopes()
#test_generate_molecularion_combination() 
#test_multiplicity_increment()
#test_sort_molecularion_list()
# est_check_file_atoms_selection()
#test_check_file_with_multiplicity()
#test_read_file_csv()
#test_write_file_csv()
