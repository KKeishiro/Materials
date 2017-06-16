# import module----------------------------------------------------------------------------------
import pandas as pd
import numpy as np
from pymatgen.matproj.rest import MPRester
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
from sklearn import svm
from sklearn import cross_validation
from sklearn.grid_search import GridSearchCV
from sklearn.learning_curve import learning_curve
plt.style.use("ggplot")


# make average database-------------------------------------------------------------------------
def make_average_data():
    atomicdata = pd.read_csv("atomicdata_2015.5edition.csv").set_index("Symbol")
    atomicdata.ix['None',:] = 0

    compoundtable = pd.read_csv("compound_table.csv").dropna()\
    .drop(["Ecoh","BM","bandgap","Compound","type","st"],axis=1).sort("formula").drop('Nval', axis=1)
    compoundtable = compoundtable[~compoundtable["formula"].duplicated(take_last=False)].set_index("formula")

    average_data = []
    for i in np.array(compoundtable.index):
        a1 = compoundtable.ix[i,'a1']
        a2 = compoundtable.ix[i,'a2']
        a3 = compoundtable.ix[i,'a3']
        n1 = compoundtable.ix[i,'n1']
        n2 = compoundtable.ix[i,'n2']
        n3 = compoundtable.ix[i,'n3']
        row_data = [((atomicdata.ix[a1,j]*n1)+(atomicdata.ix[a2,j]*n2)+(atomicdata.ix[a3,j]*n3))/(n1+n2+n3)\
                    for j in np.array(atomicdata.columns)]
        average_data.append(row_data)
    df = pd.DataFrame(average_data, index = compoundtable.index, columns = atomicdata.columns)
    df.to_csv("compound_table_ave.csv")

# make_average_data()


# eliminate double counted compounds of gained matpro data-----------------------------------------
def tidyup_matpro_data():
    df = pd.read_csv("matpro_data_regression.csv")
    df = df[~df["name"].duplicated(take_last=False)].set_index("name")
    df.to_csv("matpro_data_regression_processed.csv")


# eliminate unnecessary compounds from compound_table_ave------------------------------------------
def elim_unne_comp():
    df = pd.read_csv("compound_table_ave.csv").set_index('formula')
    df_matpro = pd.read_csv("matpro_data_regression_processed.csv").set_index("name")
    df = df.ix[df_matpro.index,:].dropna()
    df.to_csv("compound_table_ave_regression.csv")


# standardization----------------------------------------------------------------------------------
def standardization():
    df_table = pd.read_csv("compound_table_ave_regression.csv").set_index('name')
    df_matpro = pd.read_csv("matpro_data_regression_processed.csv").set_index('name')
    df_matpro = df_matpro.ix[df_table.index,:]
    df = pd.concat([df_table, df_matpro], axis=1)
    scaled_list = []
    for i in np.array(df.columns):
        sc = StandardScaler()
        sc.fit(df.ix[:,i])
        array_scaled = sc.transform(df.ix[:,i])
        scaled_list.append(array_scaled)
    df_std = pd.DataFrame(scaled_list).T
    df_std.index = np.array(df.index)
    df_std.columns = np.array(df.columns)
    df_std.to_csv("compound_table_ave_regression_std.csv")


# support vector regression-----------------------------------------------------------------------
def SVR():
    df = pd.read_csv("compound_table_ave_regression_std.csv")
    df.rename(columns={df.columns[0]:"Compounds"},inplace=True)
    df = df.set_index("Compounds")
    df_x = df.drop(['energy','energy_per_atom','volume'], axis=1)
    x = np.array(df_x.values.tolist())
    y = np.array(df["energy_per_atom"])

    def main():
        tuned_parameters = [{'kernel': ['rbf'], 'gamma': [10**i for i in range(-4,1)], 'C': [10**i for i in range(0,4)]}]
        gscv = GridSearchCV(svm.SVR(), tuned_parameters, cv=4, scoring="mean_squared_error")
        gscv.fit(x, y)
        reg_best = gscv.best_estimator_
        x_line = range(-4,3)
        y_line = x_line
        plt.plot(reg_best.predict(x),y,'ro',label='SVR')
        plt.plot(x_line,y_line,'k-',label='ideal')
        plt.xlabel('Predicted value [a.u.]')
        plt.ylabel('True value [a.u.]')
        plt.title('SVR')
        plt.legend(loc='upper left',numpoints=1)
        plt.show()
        print ("RMSE: %f" % np.sqrt(-gscv.score(x,y)))
        print gscv.best_params_

    main()
SVR()
