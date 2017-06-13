from pymatgen.matproj.rest import MPRester
import matplotlib.pyplot as plt
import pandas as pd

a = MPRester("m390ovHWT4KfcLu2")
 
def make_database():
    data_list = []
    id_list = []
    entries = a.get_entries_in_chemsys(elements=["H","He","Li","Be"],\
property_data=["energy","energy_per_atom","e_above_hull","icsd_id",\
"band_gap","formation_energy_per_atom","volume","nsites","density",\
"nelements","spacegroup","unit_cell_formula"])
    for entry in entries:
        name = entry.name
        mpid = entry.entry_id
        fepa = entry.data["formation_energy_per_atom"]
        volume = entry.data["volume"]
        nsites = entry.data["nsites"]
        density = entry.data["density"]
        vpa = volume/nsites
        e = entry.data["energy"]
        epa = entry.data["energy_per_atom"]
        eah = entry.data["e_above_hull"]
        icsd = entry.data["icsd_id"]
        bg =  entry.data["band_gap"]
        nelements = entry.data["nelements"]
        spacegroup = entry.data["spacegroup"]
        unit_formula = entry.data["unit_cell_formula"]
        data = [name,e,epa,fepa,eah,vpa,density,nsites,nelements,icsd,bg,spacegroup,\
unit_formula]
        data_list.append(data)
        id_list.append(mpid)
    frame = pd.DataFrame(data_list)
    frame.index = id_list
    frame.columns = ["chemical_formula","energy","energy_per_atom","formation_energy_per_atom",\
"e_above_hull","volume","density","nsites","nelements","icsd_id","band_gap",\
"spacegroup","unit_cell_formula"]
    frame.to_csv("./all_database.csv")
