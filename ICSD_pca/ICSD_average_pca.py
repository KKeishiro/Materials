# import module-------------------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
plt.style.use("ggplot")


# tidy up and make new atomic dataframe-------------------------------------------------------
def tidyup_atomicdata():
    atomicdata = pd.read_csv('/Users/user/Documents/kyoyu2/data_part1/atomic_alldata.csv').set_index('Symbol')
    atomicdata = atomicdata.drop(atomicdata.ix['Po':'Rg'].index).drop(['RareEarth','Nval','IE2','IE3','Cp-g','Cp-mol',\
        'E-vapor','ElectricalConductivity','MolerMagneticSusceptibility','ThermalExpansion'],axis=1).fillna(0)

    dataframe = pd.DataFrame(atomicdata)
    dataframe.to_csv('/Users/user/Documents/kyoyu2/data_part1/atomicdata_processed.csv')


# tidy up and make new compounds list----------------------------------------------------------
def tidyup_compoundslist(filename):
    # make unnecessary compounds list
    ld = open(filename)
    lines = ld.readlines()
    ld.close

    delete_list = []
    for i in ['He','Ne','Ar','Kr','Xe','Tc','Pm','Tb','Ne','Po','At','Rn','Fr','Ra','Ac',\
    'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']:
        for line in lines:
            if line.find(i) >= 0:
                delete_list.append(line[:-1])

    unnecessary_compounds = []
    for j in delete_list:
        unnecessary_compounds.append(j.split()[0])

    # delete unnecessary compounds from compounds list
    new_complist=pd.read_csv(filename,header=None,delim_whitespace=True)
    new_complist.columns=["Compounds","H","He","Li","Be","B","C","N","O","F",\
                                     "Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",\
                                     "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As",\
                                     "Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh",\
                                     "Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La",\
                                     "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",\
                                    "Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi"]
    new_complist = new_complist.set_index("Compounds")

    for i in unnecessary_compounds:
        new_complist= new_complist.drop(i)

    new_complist = new_complist.reset_index().sort("Compounds")
    new_complist = new_complist[~new_complist["Compounds"].duplicated(take_last=False)].set_index("Compounds")
    df = pd.DataFrame(new_complist)
    name, ext = os.path.splitext(filename)
    df.to_csv(name+"_processed.csv")

# filename = "complist_integer.csv"
# tidyup_compoundslist(filename)


# make average database of compounds------------------------------------------------------------
def make_ave_data(inputfile):
    A = pd.read_csv(inputfile).set_index("Compounds")
    B = pd.read_csv("atomicdata_processed.csv").set_index("Symbol")
    C = A.dot(B)
    for i in range(len(A.index)):
        C.ix[i , :] = C.ix[i , :]/(A.sum(axis=1)[i])
    name, ext = os.path.splitext(inputfile)
    C.to_csv(name+"_ave.csv")

# inputfile = "complist_integer_processed.csv"
# make_ave_data(inputfile)


# eliminate simple substances and unnecessary properties---------------------------------------
def elim_simp(filename):
    df = pd.read_csv(filename).set_index("Compounds")\
    .drop(['Z','Nf','CuKa','MoKa'], axis=1)\
    .drop(['H1','Li1','Be1',"B1","C1","N1","O1","Na1","Mg1","Al1","Si1","P1","S1","K1","Ca1","Sc1",\
            "Ti1","V1","Cr1","Mn1","Fe1","Co1","Ni1","Cu1","Zn1","Ga1","Ge1","As1",\
           "Se1","Rb1","Sr1","Y1","Zr1","Nb1","Mo1","Ru1","Rh1","Pd1","Ag1","Cd1","In1",\
           "Sn1","Sb1","Te1","I1","Cs1","Ba1","La1","Ce1","Pr1","Nd1","Sm1","Eu1","Gd1",\
           "Dy1","Ho1","Er1","Tm1","Yb1","Lu1","Hf1","Ta1","W1","Re1","Os1","Ir1","Pt1",\
           "Au1","Hg1","Tl1","Pb1","Bi1"])
    name, ext = os.path.splitext(filename)
    df.to_csv(name+"_nosimp.csv")

# filename = "complist_integer_processed_ave.csv"
# elim_simp(filename)


# nomalization---------------------------------------------------------------------------------
def nomalization(filename):
    df = pd.read_csv(filename).set_index("Compounds")
    scaled_list = []
    for i in np.array(df.columns):
        sc = StandardScaler()
        sc.fit(df.ix[:,i])
        array_scaled = sc.transform(df.ix[:,i])
        scaled_list.append(array_scaled)
    df_nom = pd.DataFrame(scaled_list).T
    df_nom.index = np.array(df.index)
    df_nom.columns = np.array(df.columns)
    name, ext = os.path.splitext(filename)
    df_nom.to_csv(name+"_nom.csv")

# filename = "complist_integer_processed_ave_nosimp.csv"
# nomalization(filename)


#  PCA(Principal Component Analysis)-------------------------------------------------------------
def pca(filename):
    df = pd.read_csv(filename)
    df.rename(columns={df.columns[0]:"Compounds"},inplace=True)
    df = df.set_index("Compounds")
    pca = PCA(n_components = 3, whiten = False)
    pca.fit(df)
    df_pca = pca.transform(df)
    E = pca.explained_variance_ratio_
    pca_components = pd.DataFrame(pca.components_ ,index = ["first component","second component","third component"], columns = np.array(df.columns))
    df_pca = pd.DataFrame(df_pca, columns = ["first component","second component","third component"], index = df.index)
    print pca.components_
    print np.cumsum(E)
    name, ext = os.path.splitext(filename)
    df_pca.to_csv(name+"_pca.csv")

# filename = "complist_integer_processed_ave_nosimp_nom.csv"
# pca(filename)


# visualize PCA data-------------------------------------------------------------------------------
def plot_2D(filename):
    df = pd.read_csv(filename).set_index("Compounds")
    plt.scatter(df.ix[:,0], df.ix[:,1])
    plt.xlabel("first component")
    plt.ylabel("second component")
    plt.show()


def plot_3D(filename):
    df = pd.read_csv(filename).set_index("Compounds")
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(df.ix[:,0], df.ix[:,1], df.ix[:,2])
    ax.set_xlabel("first component")
    ax.set_ylabel("second component")
    ax.set_zlabel("third component")
    plt.show()

# filename = "complist_integer_processed_ave_nosimp_nom_pca.csv"
# plot_2D(filename)
# plot_3D(filename)

def plot_3D_categorized_together():
    fig = plt.figure()
    ax = Axes3D(fig)
    file_list = ["nitrides_pca_only.csv","oxides_pca_only.csv","fluorides_pca_only.csv","phosphides_pca_only.csv","sulfides_pca_only.csv","clorides_pca_only.csv"]
    color_list = ["r","b","pink","g","olive","orange"]
    for i,j in zip(file_list, color_list):
        df = pd.read_csv(i,header=None,names=["Compounds","first component","second component","third component"],index_col="Compounds")
        name, ext = os.path.splitext(i)
        ax.plot(df.ix[:,0], df.ix[:,1], df.ix[:,2], 'o', color=str(j), ms=3, label=name)
        ax.set_xlabel("first component")
        ax.set_ylabel("second component")
        ax.set_zlabel("third component")
        ax.legend(numpoints=1)
    plt.show()

# plot_3D_categorized_together()

def plot_3D_categorized_separated():
    fig = plt.figure()
    file_list = ["nitrides_pca_only.csv","oxides_pca_only.csv","fluorides_pca_only.csv","phosphides_pca_only.csv","sulfides_pca_only.csv","clorides_pca_only.csv"]
    color_list = ["r","b","pink","g","olive","orange"]
    subplot_list = ["231","232","233","234","235","236"]
    for i,j,k in zip(file_list, color_list, subplot_list):
        df = pd.read_csv(i,header=None,names=["Compounds","first component","second component","third component"],index_col="Compounds")
        name, ext = os.path.splitext(i)
        ax = fig.add_subplot(k, projection='3d')
        ax.plot(df.ix[:,0], df.ix[:,1], df.ix[:,2], 'o', color=str(j), ms=3, label=name)
        ax.set_xlabel("first component")
        ax.set_ylabel("second component")
        ax.set_zlabel("third component")
        ax.legend(numpoints=1)
    plt.show()

# plot_3D_categorized_separated()
