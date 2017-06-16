from pymatgen.matproj.rest import MPRester
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

a = MPRester("m390ovHWT4KfcLu2")

# make oxides database---------------------------------------------------------------------------
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
    df = pd.DataFrame(data_list)
    df.index = id_list
    df.columns = ["chemical_formula","formation_energy","volume","density"]
    df.to_csv("matpro_oxides_database.csv")


# check materials project-------------------------------------------------------------------------
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


# make database by id------------------------------------------------------------------------------
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
    df = pd.DataFrame(data_list)
    df.to_csv("database_by_id"+str(min_id)+".csv")


# gather data by function 'make_database_by_id()'----------------------------------------------
def gather_data():
    for j in [5000,10000]:
        make_database_by_id(j,j+5000)


# get matpro data for regression---------------------------------------------------------------
def get_data_for_regression():
    df = pd.read_csv("compound_table_ave.csv").sort('formula').drop('Nval', axis=1).set_index('formula')
    data_list = []
    for i in np.array(df.index):
        entries = a.get_entries(i,property_data=["energy","energy_per_atom","volume"])
        for entry in entries:
            name = entry.name
            e = entry.data["energy"]
            epa = entry.data["energy_per_atom"]
            volume = entry.data["volume"]
            data = [name,e,epa,volume]
            data_list.append(data)
    df_regression = pd.DataFrame(data_list, columns = ["name","energy","energy_per_atom","volume"]).set_index("name")
    df_regression.to_csv("matpro_data_regression.csv")

# get_data_for_regression()
