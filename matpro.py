from pymatgen.matproj.rest import MPRester
import matplotlib.pyplot as plt
import pandas as pd

a = MPRester("m390ovHWT4KfcLu2")
 
def make_database():
    atom_list = ["K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br"]
    data_list = []
    id_list = []
    for atom in atom_list:
        chemical_system = atom + "-O"
        entries = a.get_entries(chemical_system,property_data=["formation_energy_per_atom","volume","nsites","density"])
        for entry in entries:
            name = entry.name
            mpid = entry.entry_id
            fepa = entry.data["formation_energy_per_atom"]
            volume = entry.data["volume"]
            nsites = entry.data["nsites"]
            density = entry.data["density"]
            vpa = volume/nsites
            data = [name,fepa,vpa,density]
            data_list.append(data)
            id_list.append(mpid)
    frame = pd.DataFrame(data_list)
    frame.index = id_list
    frame.columns = ["chemical_formula","formation_energy","volume","density"]
    frame.to_csv("./database.csv")

def check_database():
    number_of_data = 0
    miss_id_list = []
    for i in range(100000):
        mpid = "mp-" + str(i)
        try:
            a.get_entry_by_material_id(mpid)
            number_of_data += 1    
        except IndexError:
            print i
            miss_id_list.append(i)
            continue
    print number_of_data

def gather_data():
    for j in [5000,10000]:
        make_database_by_id(j,j+5000)
def make_database_by_id(min_id,max_id):
    number_of_data = 0
    miss_id_list = []
    property_list = ["formation_energy_per_atom","volume","nsites","density","pretty_formula",\
"unit_cell_formula","icsd_id","energy","energy_per_atom","elements","nelements",\
"initial_structure","final_structure","structure","e_above_hull","is_hubbard",\
"hubbards","is_compatible","band_gap"]
    data_list = []
    for i in range(min_id,max_id):
        mpid = "mp-" + str(i)
        try:
            entry = a.get_entry_by_material_id(mpid,property_data=property_list)
        except IndexError:
            print i
            continue
        except MPRestError:
            print i,"MPRestError"
            continue
        data = entry.data
        data["mpid"] = entry.entry_id
        data["name"] = entry.name
        data_list.append(data)
    frame = pd.DataFrame(data_list)
    frame.to_csv("./database_by_id"+str(min_id)+".csv")
