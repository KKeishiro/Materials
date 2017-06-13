from pymatgen.matproj.rest import MPRester
import matplotlib.pyplot as plt
import pandas as pd

a = MPRester("m390ovHWT4KfcLu2")
 
def make_database():
    data_list = []
    id_list = []
    entries = a.get_entries_in_chemsys(elements=["H","He","Li","Be","B","C","N",\
"O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",\
"Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",\
"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",\
"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",\
"Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Ac","Th","Pa","U","Np","Pu"],\
property_data=["energy","energy_per_atom","e_above_hull","icsd_id",\
"band_gap","formation_energy_per_atom","volume","nsites","density","nelements"])
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
        data = [name,e,epa,fepa,eah,volume,vpa,density,nsites,nelements,icsd,bg]
        data_list.append(data)
        id_list.append(mpid)
    frame = pd.DataFrame(data_list)
    frame.index = id_list
    frame.columns = ["chemical_formula","energy","energy_per_atom","formation_energy_per_atom",\
"e_above_hull","volume","volume/nsites","density","nsites",\
"nelements","icsd_id","band_gap"]
    frame.to_csv("./all_database.csv")
