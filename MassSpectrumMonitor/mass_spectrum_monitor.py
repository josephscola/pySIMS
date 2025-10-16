import tkinter as tk
import ast
from pysims.mass import *
from tkinter import ttk, Scale, HORIZONTAL, Label, StringVar, Toplevel, Frame, Button, messagebox
from tkinter import filedialog
from tkinter import END
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from mendeleev import element
import re
import itertools
import os
import csv
from math import prod
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import matplotlib as mpl
from operator import itemgetter, attrgetter
from src.MassInterference import Molecularion, generate_single_isotopes, generate_molecularion_combination, file_name_with_multiplicity, molecularion_list_to_table

debug = None
debug_bis = False

FOLDER_SEPARATION_CHAR = '/'

fig_parameters = {              #set les paramètres par défaut des settings du plot
     'axes.titlesize' : 16,
     'axes.labelsize' : 16,
     'lines.linewidth' : 1.5,
     'lines.markersize' : 5,
     'xtick.labelsize' : 16,
     'ytick.labelsize' : 16,
     'legend.fontsize' : 13,
      'font.size' : 16
}

for p in fig_parameters :
    mpl.rcParams [p] = fig_parameters [p]

# --- Periodic Table GUI Logic ---
selected_atoms = set()
buttons = {}
is_data = 0         # is_data = 0 si pas de données du Sims plotées /// is_data = 1 si des données du Sims sont plotées

periodic_table_layout = [
    ["H", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "He"],
    ["Li", "Be", "", "", "", "", "", "", "", "", "", "", "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", "", "", "", "", "", "", "", "", "", "", "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "", "", ""],
    ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],
    ["", "", "", "", "La", "Ce", "Pr", "Nd", "", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"],
    ["", "", "", "", "", "Th", "Pa", "U", "", "", "", "", "", "", "", "", "", ""]
]

def handle_click(symbol):
    # Permet de gérer le click des boutons : Si un bouton est cliqué il devient enfoncé et si on reclique cela redevient normal. Le bouton reste également enfoncé après la sélection
    if symbol not in selected_atoms:
        selected_atoms.add(symbol)
        buttons[symbol].config(relief="sunken")
    else:
        selected_atoms.remove(symbol)
        buttons[symbol].config(relief="raised")
    update_selected_atoms_label()

def create_periodic_table_window(root):
    # Créer un tableau périodique interactive qui revoie les atomes cliqués
    try:
        if table_window:
            table_window.destroy()
    except:
        pass
    
    table_window = tk.Toplevel(root)
    frame = tk.Frame(table_window)
    frame.pack(expand=True, fill="both")

    global buttons
    buttons = {}
    for row_index, row in enumerate(periodic_table_layout):
        for col_index, elem in enumerate(row):
            if elem:
                b = tk.Button(frame, text=elem, width=4, height=2,
                              command=lambda e=elem: handle_click(e))
                b.grid(row=row_index, column=col_index, padx=1, pady=1)
                buttons[elem] = b
                if elem in selected_atoms:
                    b.config(relief="sunken")
    
    
    
    #ALLOW_UNSTABLE_ISOTOPES = tk.IntVar ()
    checkUnstable = tk.Checkbutton (table_window, text = 'Allow unstable isotopes',
                                    variable = ALLOW_UNSTABLE_ISOTOPES)#, onvalue=1, offvalue=0)
    checkUnstable.pack (padx=10, pady=1)

    
    
    validate_button = tk.Button(table_window, text="Valider", command=table_window.destroy)
    validate_button.pack(pady=10, padx=10)
    validate_button.pack(side="bottom", fill="x")
    
    
    return table_window


def get_selected_atoms_pt():
    return selected_atoms

# --- Main GUI Logic ---
current_canvas = None
selected_atoms_label = None
csv_frame = None
root = None
slider_multiplicite = None
right_frame = None
graph_app = None
center_frame = None


selected_atoms_list = []
selected_names_in_csv_list = []
selected_masses_in_csv_list = []
selected_abundances_in_csv_list = []


def afficher_csv():
    # Génère ou lit un fichier CSV contenant des combinaisons d'ions moléculaires et affiche son contenu dans un widget Treeview pour pouvoir intéragire sur l'interface
    global selected_atoms, entry_multiplicite, csv_frame, tree

    if not selected_atoms:
        messagebox.showwarning("Warning", "Please select at least one element")
        return

    multiplicity = slider_multiplicite.get()
    path_folder = f".{FOLDER_SEPARATION_CHAR}molecularions{FOLDER_SEPARATION_CHAR}"
    os.makedirs(path_folder, exist_ok=True)

    file_path = path_folder + file_name_with_multiplicity(sorted(list(selected_atoms)), multiplicity)
    if debug: print(file_path)
    
    """
    if not os.path.exists(file_path):
    """
    try:
        single_isotopes = generate_single_isotopes(selected_atoms, UNSTABLE_ISOTOPES = ALLOW_UNSTABLE_ISOTOPES.get ())
        combination_list = generate_molecularion_combination(single_isotopes, multiplicity)
        table = molecularion_list_to_table (combination_list)
        csv_content = pd.DataFrame (table)
        csv_content = csv_content.T.set_index(0, drop = True).T # makes first row the index
        csv_content = csv_content.sort_values(by = "mass", ascending = True)
        csv_content.to_csv (file_path, index = False)
    except Exception as e:
        messagebox.showerror("Error", f"Molecular ion file generation error :\n{e}")
        return


    try:
        csv_content = pd.read_csv(file_path, delimiter=',')
    except Exception as e:
        messagebox.showerror("Error", f"Molecularion file reading error :\n{e}")
        return

    
    # Nettoyer le contenu du frame CSV
    for widget in csv_frame.winfo_children():
        widget.destroy()

    # Créer le tableau dans le csv_frame
    tree = ttk.Treeview(csv_frame, height=10)
    tree.pack(fill="both", expand=True)

    tree.bind("<<TreeviewSelect>>", on_tree_select)

    tree.tag_configure("selected", background="light blue")
    tree.tag_configure("not_selected", background = "white")

    vsb = ttk.Scrollbar(csv_frame, orient="vertical", command=tree.yview)
    vsb.place(relx=0.978, rely=0.175, relheight=0.713, relwidth=0.020)

    tree["columns"] = list(csv_content.columns)
    tree["show"] = "headings"

    for col in csv_content.columns:
        tree.heading(col, text=col)
        tree.column(col, width=80, anchor="center")
        if debug: print(col)

    for _, row in csv_content.iterrows():
        formatted_row = [row["name"], f'{row["mass"]:.4f}', f'{row["abundance"]:.4f}']
        tree.insert("", "end", values=formatted_row)


def selecting_atoms_in_csv():

    global tree, selected_item_id

    selected_item_id = tree.focus()
    selected_item = tree.item(selected_item_id)

    values = selected_item['values']
    
    selected_name_in_csv = values[0]   #les valeurs prises ici sont les valeurs approchées à 5 chiffres significatifs, ne sont donc pas exactement les mêmes que dans le fichier csv
    selected_mass_in_csv = float(values[1])
    selected_abundance_in_csv = float(values[2]) 

    molion = Molecularion(selected_name_in_csv, selected_mass_in_csv, selected_abundance_in_csv)

    app.update_data_to_plot(molion)
    app.draw_figure()
    
    return selected_atoms_list, molion


def on_tree_select(event):
    # Juste appeler la fonction, ne pas essayer de récupérer les retours ici
    selecting_atoms_in_csv()


def open_selection():
    global root
    create_periodic_table_window(root)

        
def update_selected_atoms_label():
    global selected_atoms_label, selected_atoms
    if selected_atoms_label:
        selected_atoms_label.config(text=f"Atomes sélectionnés : {', '.join(sorted(list(selected_atoms)))}")


######################################################################
#                Fonctions lancées par les boutons                   #
######################################################################


def open_graph_button_command() :
    app.SIMS_data()


def erase_graph():
    app.clear_graph()
    app.forget_data()


def close_interface():
    root.quit()


def fct_for_export():
    app.fct_valider_export_donner()


def open_export_menu():
    app.export()


def select_path_file():
    file = filedialog.askopenfilename()
    return file


def fct_linear_graph() :
    app.graph_style_from_linear_to_log()


def fct_log_graph():
    app.graph_syle_from_log_to_linear()


def fct_ajuster_echelle():
    app.update_i_norm()


def zoom_fct():
    try : 
        valeur_entry_min_mass = entry_min_mass.get()
        valeur_entry_max_mass = entry_max_mass.get()

    except Exception as e:

        print(e)
        messagebox.showerror("Erreur",f"{e}")
        return


    app.zoom_management(min_mass_entry = valeur_entry_min_mass, max_mass_entry = valeur_entry_max_mass)


######################################################################
#                       Class GraphApp                               #
######################################################################  

class GraphApp:

    def __init__(self, parent_frame):

        """
        Version adaptée de GraphApp pour les spectres de masse
        
        """
        self.parent_frame = parent_frame

        # Création de la figure matplotlib
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.fig.patch.set_facecolor('white')
        self.ax.set_xlabel("Masse (a.m.u) ")
        self.ax.set_ylabel("Intensité (cps)")
        
        # Intégration de matplotlib dans tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.path_file_graph_entry = None
        self.title_entry = None
        self.graph_file_window = None
        self.path_entry = None
        self.nom_fichier_entry = None

        self.sims_data_x = None
        self.sims_data_y = None
        self.selected_point = None
        self.x_min = None
        self.x_max = None
        self.i_norm = 1
        self.delta_mass = ''

        self.displayed_molion_list = []
        self.data_to_plot = []
        
        self.current_graph_type = 0  # 0 = linear, 1 = log

        self.Mass_profile = None

        self.path_folder = r"./"

        self.user_colors_path = self.path_folder + 'user_colors.txt'
        self.default_colors_path = self.path_folder + 'default_colors.txt'
  

    def clear_graph(self):
        # efface le graphique et réinitialise les données enregistrées
        self.ax.clear()
        self.canvas.draw()


    def forget_data(self):

        self.current_graph_type = 0
        button_linear.config(relief = 'sunken')
        button_log.config(relief = "raised")

        if csv_frame.winfo_children():
            for item_id in tree.get_children():
                tree.item(item_id, tags=())

        self.path_file_graph_entry = None
        self.title_entry = None
        self.graph_file_window = None
        self.path_entry = None
        self.nom_fichier_entry = None

        self.sims_data_x = []
        self.sims_data_y = []
        self.selected_point = None
        self.x_min = None
        self.x_max = None
        self.i_norm = 1
        self.aimed_mass = ''
        self.aimed_intensity = ''
        self.delta_mass = ''
        self.Mass_profile = None

        label_min_mass.config(text = f'Masse minimale : ')
        label_max_mass.config(text = f'Masse minimale : ')
        mass_label.config(text = f'Masse visée : {self.aimed_mass}')
        delta_label.config(text = f' \u0394 Masse visée : {self.delta_mass}')

        self.displayed_molion_list = []
        self.data_to_plot = []  #[[mass, abundance, kwargs], [] ...]
        app.draw_figure()


    def draw_figure(self):
        self.clear_graph()

        for x in self.data_to_plot :
            for i in range(0,len(x[2])):               
                if type(x[1]) == list:
                    self.ax.plot(x[0], x[1], **x[2][i])
                if type(x[1]) == float:
                    self.ax.plot(x[0], x[1] * self.i_norm, **x[2][i])

        if self.current_graph_type == 0 :
            pass
        else : 
            self.ax.set_yscale('log')

        if self.x_min is not None and self.x_max is not None:
            self.ax.set_xlim(self.x_min, self.x_max)
            self.y_min = -100
            y_max = 0          
            converted_sims_data_x = np.array(self.sims_data_x)

            for x in converted_sims_data_x[(converted_sims_data_x >= self.x_min) & (converted_sims_data_x <= self.x_max)]:
                local_y_max = self.Mass_profile.local_max(x)[1]
                if local_y_max > y_max :
                    y_max = local_y_max
            self.y_max = y_max * 1.08
            self.ax.set_ylim(self.y_min, self.y_max)
        self.ax.set_xlabel("Masse (a.m.u) ")
        self.ax.set_ylabel("Intensité (cps)")
        # display credits
        credits_label = 'GEMaC, CNRS UMR 8635, U. Paris-Saclay/UVSQ, 2025'
        self.ax.text (.82, -.07, credits_label, fontsize = 10, fontstyle = 'italic',
                      horizontalalignment='center', verticalalignment='center',
                      transform=self.ax.transAxes)
        plt.tight_layout()
        self.canvas.draw()
        

    def get_current_data(self):
           #Retourne les données actuellement affichées
           return {
               'masses': self.current_masses,
               'abundances': self.current_abundances,
               'names': self.current_names
           }
    

    def read_isotope_name(self, name):
        elem = ''
        int_mass_str = ''
        inmass = True
        for c in name:
             if not c.isdigit () and not inmass :
                 elem += c
             if not c.isdigit () and inmass:
                 inmass = not inmass
                 elem += c
             if c.isdigit () and inmass :
                 int_mass_str += c
             if c.isdigit () and not inmass :
                 elem += c
        int_mass = int (int_mass_str)
        return int_mass, elem


    def set_marker_format(self, molion) :        
       # information = [[name_0, mass_0, abundance_0], [name_1, mass_1,abundance_1], ... ,  ]   
       #### fonction colors qui choisi les couleur pour chaque atome (entrée utilisateur ou dico par défaut)
        #max_size = slider_marker_size.get()
        max_size = 16
        dissociated_elements = molion.dissociate_constituent_isotopes()
        dissociated_molion_into_elements = sorted(list(molion.dissociate_constituent_atoms()))

        if not os.path.exists(self.path_folder + 'user_colors.txt') :
            # cas le fichier user_colors.txt n'existe pas dans le dossier du programme (cas défault)
            M = len(dissociated_elements)
            if debug_bis:
                print('len(dissociated_element) : ', M)
            kwargs = []

            for k in range (0,M):
                i = sorted(list(selected_atoms)).index((dissociated_molion_into_elements[k]))                
                size = max_size * (M - k) / M
                kwargs.append({'marker' : 'o',
                	  'markersize' : size,
                      'color' : f'C{i}'})    
            if debug_bis: print ("liste kwargs : ", kwargs)
                
        else :                      
            #cas user_colors.txt existe
            if debug_bis :
                print("Boucle user_colors")
            d = open('user_colors.txt', 'r')
            read_dico = ast.literal_eval(d.read ())
            d.close ()
    
            M = len(dissociated_elements)
            int_mass, elem = self.read_isotope_name(dissociated_elements[0])
            i = sorted(list(selected_atoms)).index((elem))
            if debug: print("i : ", i)
            kwargs = []

            for i in range (0,M):
                try :   
                    size = max_size * (M - i) / M
                    alpha = 1.0 - (i * 0.2)
                      # Transparence croissante vers l'intérieur   
                    kwargs.append({'marker' : 'o',
                    	  'markersize' : size,
                          'alpha' : alpha,
                          'color' : read_dico[dissociated_elements[i]]})  
                except KeyError as e:
                    messagebox.showerror("Erreur", f"L'isotope '{molion.name}' n'a pas été trouvé dans les couleurs définies.")
                    tree.item(selected_item_id, tag = ("not_selected"))                 
        return kwargs        
    
    
    def custom_marker(self):        #markers camemberts
        pass
        

    def update_data_to_plot(self, molion):
        tree.tag_configure("selected", background="light blue") 
        tree.tag_configure("not_selected", background = "white")
        name = molion.name
        mass = molion.mass
        abundance = molion.abundance
        kwargs = self.set_marker_format(molion)
        selected_info_to_plot = [mass, abundance, kwargs]

        if selected_info_to_plot not in self.data_to_plot :
            self.data_to_plot.append(selected_info_to_plot)
            tree.item(selected_item_id, tag = ("selected"))
            self.displayed_molion_list.append([name, mass, abundance])
        else:
            self.data_to_plot.remove(selected_info_to_plot)
            tree.item(selected_item_id, tag = ("not_selected"))
            self.displayed_molion_list.remove([name, mass, abundance])

    def SIMS_data(self):
    
        path_file = select_path_file()

        self.kwargs_sims = {
            #'linewidth': 1.5,
            'linestyle': '-',
            'color': 'k'
            }

        if debug_bis :
            print(path_file)

        #title = self.title_entry.get()
        # Récupère et trace un spectre de masse à partir d'un fichier utilisateur
        #if not self.path_file_graph_entry or not self.title_entry:
        #    print("Champs de saisie non initialisés.")
        #    return


        self.Mass_profile = MassSpectrum(path_file)
        
        #data_mass = Mass_profile.mass
        #data_intensity = Mass_profile.intensity

        self.sims_data_x = self.Mass_profile.mass
        self.sims_data_y = self.Mass_profile.intensity

        self.sims_data = [self.Mass_profile.mass, self.Mass_profile.intensity, [self.kwargs_sims]]

        self.data_to_plot.append(self.sims_data)

        if self.graph_file_window:
            self.graph_file_window.destroy()

        self.draw_figure()


    def zoom_management(self, min_mass_entry = None, max_mass_entry = None):

        if min_mass_entry == 'min':
            self.x_min = min(self.sims_data_x)
            self.y_min = min(self.sims_data_y)
        else:
            self.x_min = int(min_mass_entry)

        if max_mass_entry == 'max':
            self.x_max = max(self.sims_data_x)
            self.y_max = max(self.sims_data_y)
        else:
            self.x_max = int(max_mass_entry)

        label_min_mass.config(text = f'Masse minimale : {self.x_min}')
        label_max_mass.config(text = f'Masse minimale : {self.x_max}')

        entry_min_mass.delete(0, END)
        entry_max_mass.delete(0, END)

        self.draw_figure()


    def export(self):

        if self.sims_data_x is None:
            messagebox.showerror("Warning", 'Please load an experimental mass spectrum first')
            return
        
        self.nom_fichier_entry = 'ms-monitor_figure_data' #Ou demander le nom du fichier ?

        self.export_path_folder = filedialog.askdirectory()

        if not self.export_path_folder:
            return

        export_file = os.path.join(self.export_path_folder, self.nom_fichier_entry + '.csv')
        if os.path.exists(export_file):
            messagebox.showerror("Warning", 'An exported figure file already exists')
            return
        
        try:
            with open(export_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)

                #Ligne 1 : titres généraux
                header_row_1 = ['SIMS DATA', '', '', '']  
                for molion in self.displayed_molion_list:
                    header_row_1.append(molion[0])        
                    header_row_1.extend(['', '', ''])    
           
                writer.writerow(header_row_1)

                #Ligne 2 : noms des colonnes 
                header_row_2 = ['Mass', 'Intensity', '']  
                for _ in self.displayed_molion_list:
                    header_row_2.extend(['Mass', 'Abundance', 'Normalized', '']) 
                writer.writerow(header_row_2)

                # Lignes n : données
                for i in range(len(self.sims_data_x)):
                    
                    row = [self.sims_data_x[i], self.sims_data_y[i], '']  

                    for molion in self.displayed_molion_list:
                        if i == 0:
                            row.extend([molion[1], molion[2], molion[2] * self.i_norm])
                        else:
                            row.extend(['', '', ''])
                        row.append('')  

                    writer.writerow(row)

            messagebox.showinfo("Succès", "Le fichier a été enregistré avec succès.")

        except Exception as e:
            messagebox.showerror("Erreur", f"Une erreur est survenue lors de l'export : {e}")


    def update_i_norm(self):
        self.aimed_mass = entry_aimed_mass.get()
        self.aimed_intensity = entry_aimed_intensity.get()

       #if entry_min_mass.get() == '':
       #    self.x_min = -10
       #if entry_max_mass.get() == '':
       #    self.x_max = 210


        if self.aimed_mass:
            mass_label.config(text = f'Masse visée : {self.aimed_mass}')
            self.aimed_intensity == ''
            intensity_label.config(text = f'Intensité visée : {self.aimed_intensity}')
            try:
                masse_cible = float(self.aimed_mass)
                self.delta_mass = entry_delta.get()
                if debug: print('self.delta_mass : ', self.delta_mass)
                if self.delta_mass != '' :
                    dm = self.Mass_profile.mass [1] - self.Mass_profile.mass [0]
                    n = int (float (self.delta_mass) / dm)
                    local_max = self.Mass_profile.local_max(masse_cible, n = n)
                    delta_label.config(text = f' Recherche du max sur : {self.delta_mass} (a.m.u.)')
                else : 
                    local_max = self.Mass_profile.local_max(masse_cible)
                converted_x_to_plot = np.array(self.sims_data_x)
                converted_y_to_plot = np.array(self.sims_data_y)
                mask_plot = (converted_x_to_plot >= self.x_min) & (converted_x_to_plot <= self.x_max)
                
                max_abundance = 0

                for molion in self.displayed_molion_list:

                    if molion[2] > max_abundance:
                        max_abundance = molion[2]

                self.i_norm = local_max[1] / max_abundance
    
            except Exception as e:
                print(e)
                messagebox.showerror('Erreur', str(e))

        if self.aimed_intensity:

            max_abundance = 0
            for molion in self.displayed_molion_list:
                if molion[2] > max_abundance:
                    max_abundance = molion[2]

            self.i_norm = float(self.aimed_intensity) / max_abundance
            intensity_label.config(text = f'Intensité visée : {self.aimed_intensity}')
            self.aimed_mass = ''
            self.delta_mass = ''
            mass_label.config(text = f'Masse visée : {self.aimed_mass}')
            delta_label.config(text = f' \u0394 Masse visée : {self.delta_mass}')


        entry_aimed_mass.delete(0, END)
        entry_aimed_intensity.delete(0, END)
        entry_delta.delete(0, END)

        self.draw_figure()


    def graph_style_from_linear_to_log(self):
        if self.current_graph_type == 0 :
            pass
            
        else : 
            self.current_graph_type = 0
            button_linear.config(relief = 'sunken')
            button_log.config(relief = 'raised')
            self.draw_figure()

    def graph_syle_from_log_to_linear(self) :
        if self.current_graph_type == 1:
            pass

        else : 
            self.current_graph_type = 1
            button_linear.config(relief = 'raised')
            button_log.config(relief = 'sunken')
            self.draw_figure()
    

##################################################################
#                         Main Tkinter                           #
##################################################################



def main():  
    global app, root, left_frame, center_frame, right_frame, slider_multiplicite, selected_atoms_label, csv_frame, entry_min_mass, entry_max_mass, button_linear, button_log
    #global slider_marker_size   
    global entry_aimed_intensity, entry_aimed_mass, mass_label, intensity_label, label_min_mass, label_max_mass, entry_delta, delta_label
    global ALLOW_UNSTABLE_ISOTOPES

    root = tk.Tk()
    root.title("Mass Spectrum Monitor")
    root.geometry("1800x1000")

    # Zone gauche — boutons
    left_frame = tk.Frame(root, bg="lightgrey", width=250)
    left_frame.pack(side="left", fill="y")
    left_frame.pack_propagate(False)

    # Zone centrale — graphique
    center_frame = tk.Frame(root, bg="white", width=880)
    center_frame.pack(side="left", fill="both", expand=True)
    center_frame.pack_propagate(False)

    # Zone droite — CSV
    right_frame = tk.Frame(root, bg="white", width=270)
    right_frame.pack(side="left", fill="y")
    right_frame.pack_propagate(False)


    tk.Button(left_frame, text="Periodic table", command=open_selection).pack(pady=5, padx=10, fill="x")
    tk.Button(left_frame, text="Molecular ions", command=afficher_csv).pack(pady=5, padx=10, fill="x")
    tk.Button(left_frame, text="Open datafile", command=open_graph_button_command).pack(pady=5, padx=10, fill="x")
    tk.Button(left_frame, text="Export figure", command=open_export_menu).pack(pady=5, padx=10, fill="x")
    tk.Button(left_frame, text="Clear figure", command=erase_graph).pack(pady=5, padx=10, fill="x")

    scale_frame = tk.Frame(left_frame, bg="#f0f0f0")
    scale_frame.pack(pady=5, padx=10, fill="x")
    button_linear = tk.Button(scale_frame, text="Linear", command=fct_linear_graph)
    button_linear.pack(side="left", expand=True, fill="x", padx=2)
    button_log = tk.Button(scale_frame, text="Log", command=fct_log_graph)
    button_log.pack(side="left", expand=True, fill="x", padx=2)
    button_linear.config(relief="sunken")

    info_frame = tk.Frame(left_frame, bg="#f0f0f0")
    info_frame.pack(pady=10, padx=10, fill="x")

    ALLOW_UNSTABLE_ISOTOPES = tk.IntVar ()

    label_min_mass = tk.Label(info_frame, text="Min mass : ")
    label_min_mass.pack(fill="x")
    entry_min_mass = tk.Entry(info_frame)
    entry_min_mass.pack(fill="x")

    label_max_mass = tk.Label(info_frame, text="Max mass : ")
    label_max_mass.pack(fill="x")
    entry_max_mass = tk.Entry(info_frame)
    entry_max_mass.pack(fill="x")

    tk.Button(info_frame, text="Validate", command=zoom_fct).pack(pady=10, padx=10, fill="x")

    slider_multiplicite = tk.Scale(left_frame, from_=1, to=6, orient="horizontal", label="Multiplicité")
    slider_multiplicite.pack(pady=10, padx=10, fill="x")
    #slider_multiplicite.set(2)

    #slider_marker_size = tk.Scale(left_frame, from_=1, to=20, orient="horizontal", label="Taille marker")
    #slider_marker_size.pack(pady=10, padx=10, fill="x")

    mass_label=tk.Label(info_frame, text='Target mass : ')
    mass_label.pack(pady=10, padx=10)
    entry_aimed_mass = tk.Entry(info_frame)
    entry_aimed_mass.pack(fill='x')

    intensity_label = tk.Label(info_frame, text="Target intensity : ")
    intensity_label.pack(fill="x")
    entry_aimed_intensity = tk.Entry(info_frame)
    entry_aimed_intensity.pack(fill="x")

    delta_label = tk.Label(info_frame, text=' \u0394 Target Mass : ')
    entry_delta = tk.Entry(info_frame)
    delta_label.pack(fill="x")
    entry_delta.pack(fill="x")    

    tk.Button(info_frame, text = 'Validate', command = fct_ajuster_echelle).pack(pady=10, padx=10, fill="x")

    #tk.Label(info_frame, text="Taille du marqueur").pack(fill="x")
    #entry_marker_size = tk.Entry(info_frame)
    #entry_marker_size.pack(fill="x")

    tk.Button(left_frame, text="Quit", command=close_interface).pack(pady=5, padx=10, fill="x")
    

    app = GraphApp(center_frame) 

    selected_atoms_label = tk.Label(right_frame, text="Atomes sélectionnés :", font=("Arial", 12), wraplength=280, bg="#f0f0f0")
    selected_atoms_label.pack(pady=10, padx=10, fill="x")

    csv_frame = tk.Frame(right_frame, bd=1, relief="solid", bg="white")
    csv_frame.pack(pady=10, padx=10, fill="both", expand=True)

    root.mainloop()

if __name__ == "__main__":
    main()
