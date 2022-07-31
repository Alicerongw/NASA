#%%
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pandas as pd
import json
import warnings
warnings.simplefilter('ignore', np.RankWarning)
from scipy.optimize import curve_fit



# Create the dictionary storing the split points for each species
def create_split_dict():
    # To improve the accuracy of polynomial fit, the temperature has been divided into several intervls
    # The split points are used for determining the interval
    split_list=[['CH4'],['H2O2'],['CLO'],['HCN','CO','HF','HI','N2','NH3','OCS','OH','PH3','C2H2','C2N2','C2H4','CH3CL','CH3F','OH','O2'],['SF6'],['CS','HO2','NO','NO+','SO'],['H2','HBr','HCL'],['SO3'],['C2H6'],['C4H2']]
    split=[[200,500,1300,1500],[200,1500],[200,4000],[200,1000],[200,1000,2000,3000,4000],[200,1000,4000],[200,1000,5000],[200,500,650],[200,1000,2000,3000],[200,1000,2000]]    
    split_dict=dict()
    for i in range(len(split_list)):
        for j in range(len(split_list[i])):
            if split_list[i][j] not in split_dict:
                split_dict[split_list[i][j]]=split[i]
    return split_dict

# Calculate specific heat, get and store the value of Tmax from HITRAN database
# This project is interested in the temperature from 200K to 6000K
# If Tmax > 6000K, then Tmax=6000K
def get_Cp(R,split_dict):
    Cp_dict=dict()
    Tmax_dict=dict()
    path = "data"
    files= os.listdir(path) 
    for file in tqdm(files): 
        if ".txt" in file:
            if not os.path.isdir(file): 
                filename=path+"/"+file
                species=file.split(".")[0]
                print(species)
                storepath="Cp_results"
                storename=storepath+"/"+species+".csv"
                # Read partition functions from database
                T,Q=np.loadtxt(filename,usecols=(0,1),unpack=True)

                # Get the value of split points
                # The default split point is 200 where the specific heat starts
                if species in split_dict:
                    split=split_dict[species]
                else:
                    split=[200]

                # Divide the original data into intervals
                T_interval=[]
                Q_interval=[]
                for i in range(len(split)):
                    if i < len(split)-1:
                        if split[i]<split[i+1]:
                            T_interval.append(T[(split[i]-1):split[i+1]])
                            Q_interval.append(Q[(split[i]-1):split[i+1]])
                    else:
                        T_interval.append(T[(split[i]-1):])
                        Q_interval.append(Q[(split[i]-1):])

                # Calculate Cp interval by interval
                for k in range(len(T_interval)):
                    x=T_interval[k]
                    y=Q_interval[k]

                    # Determine the degree of the fitting polynomial
                    n=0
                    error_min=10e20
                    for i in range(25):
                        t=i+1
                        z1 = np.polyfit(x,y,t) 
                        p1= np.poly1d(z1)
                        Qvals = p1(x)
                        error=np.mean(np.abs(Qvals-y))
                        if error < error_min:
                            error_min=error
                            n=t   

                    # Polyfit the value to get derivative
                    z1 = np.polyfit(x,y,n) 
                    p1= np.poly1d(z1)
                    Qvals = p1(x)
                    dfx = p1.deriv() 
                    ddfx = dfx.deriv()
                    dQ=dfx(x)
                    ddQ=ddfx(x)

                    # Calculate specific heat
                    Q1=x*dQ
                    Q2=(x**2)*ddQ+2*Q1
                    Cp_temp=R*((Q2/y)-(Q1/y)**2)+(5/2)*R
                    for i in range(len(x)):
                        t=float(x[i])    
                        cp=float(Cp_temp[i])  
                        if species in Cp_dict:
                            Cp_dict[species][t]=cp
                        else:
                            Cp_dict[species]=dict()
                            Cp_dict[species][t]=cp 

                # Store Cp results to csv files
                Cpdict_temp=Cp_dict[species]            
                list_sorted=sorted(Cpdict_temp.items(),key=lambda item:item[0])
                T_this=[]
                Cp_this=[]
                for i in range(len(list_sorted)):
                    T_this.append(list_sorted[i][0])
                    Cp_this.append(list_sorted[i][1])
                # Cp_df=pd.DataFrame({'T(K)': T_this, 'Cp_Thiswork': Cp_this})
                # Cp_df.to_csv(storename,index=False)

                # Store the maximum temperature in HITRAN database
                Tmax_this=T_this[-1]
                Tmax=min(Tmax_this,6000)
                Tmax_dict[species]=Tmax
            
    return Cp_dict, Tmax_dict


# Generate a dictionary storing the specific heat results from JANAF
def get_JANAF(): 
    path = "Other_Cp/JANAF"
    files= os.listdir(path) 
    Cp_JANAF=dict()
    for file in tqdm(files): 
        if ".txt" in file:
            if not os.path.isdir(file): 
                filename=path+"/"+file
                species=file.split(".")[0]
                print(species)
                # Read partition functions from database
                T,Cp=np.loadtxt(filename,usecols=(0,1),unpack=True)

                for i in range(len(T)):
                    t=T[i]    
                    cp=float(Cp[i])  
                    if species in Cp_JANAF:
                        Cp_JANAF[species][t]=cp
                    else:
                        Cp_JANAF[species]=dict()
                        Cp_JANAF[species][t]=cp 
    return Cp_JANAF

# Generate a dictionary storing the specific heat fit coefficients from Capitelli
def get_Capitelli(): 
    Capitelli=np.loadtxt("Other_Cp/Capitelli_fit.csv",delimiter=",",usecols=(3,4,5,6,7,8,9))
    Capitelli_specials=['CO','CO2','N2','NO','NO+','NO2','N2O','O2','O3']
    coe_Capitelli=dict()
    for i in range(len(Capitelli_specials)):
        n=i*3
        species=Capitelli_specials[i]
        coefficients_Ca = np.zeros([3, 7])
        for ii in range(3):
            for j in range(7):
                Capitelli[n][j]=Capitelli[n][j]*10**((-5)*(j-2))
            coefficients_Ca[ii,]=Capitelli[n]
            n+=1
        if species not in coe_Capitelli:
            coe_Capitelli[species]=coefficients_Ca
    return coe_Capitelli

# Read and store the coefficients from original nasa glenn polynomials
def get_nasa(Cp_dict):

    # The name of some species in NASA is different from this work
    # The list is used for unifying the names
    name_nasa=['C2H2,acetylene','CH3CL','CLO','COCL2','HCL','HOCL']
    name_this=['C2H2','CH3Cl','ClO','COCl2','HCl','HOCl']

    coe_nasa=dict()
    with open('Other_Cp/nasa9.dat', 'r') as file:
        stop_index = True

        while stop_index:
            line = file.readline()
            
            if line=='':
                stop_index = False
                
            elif (line.split()[0] in Cp_dict) or (line.split()[0] in name_nasa):
                print(line.split()[0])
                if line.split()[0] in name_nasa:
                    name_index = name_nasa.index(line.split()[0])
                    species=name_this[name_index]
                else:
                    species=line.split()[0]

                # Get reference for the species data:
                nasa_data = file.readline()
                num_intervals = int(nasa_data[0:2])

                coefficients = np.zeros([num_intervals, 7])
                
                for j in range(num_intervals):
                    order = file.readline()
                    coe_line1 = file.readline()
                    coe_line2 = file.readline()

                    coefficients[j,:] = [float(coe_line1[0:16].replace('D','E')), float(coe_line1[16:32].replace('D','E')),
                                    float(coe_line1[32:48].replace('D','E')), float(coe_line1[48:64].replace('D','E')),
                                    float(coe_line1[64:80].replace('D','E')), float(coe_line2[0:16].replace('D','E')),
                                    float(coe_line2[16:32].replace('D','E'))]
                
                coe_nasa[species]=coefficients
                        
            else:
                pass

    return coe_nasa

# Define the equation in the format of nasa polynomials
def get_polynomial_results(coefficients,x):
    Cp=8.314*(coefficients[0]*x**(-2)+coefficients[1]*x**(-1)+coefficients[2]+coefficients[3]*x+coefficients[4]*x**2+coefficients[5]*x**3+coefficients[6]*x**4)
    return Cp

def get_Tmax_this():
    Tmax_this=dict()
    species_list=['C2H2','CH3F','ClO','CO2','CS','H2O','H2S','HCN','HI','HO2','N2','N2O','NH3','O2','OCS','PH3','SO','SO2','SO3','NO','HCOOH']
    Tmax_list=[1200.,1900.,800.,2600.,3400.,3000.,2000.,1200.,3100.,1500.,2300.,1500.,1400.,450.,1500.,2000.,1200.,1300.,1100.,2500.,400.]
    for i in range(len(species_list)):
        if species_list[i] not in Tmax_this:
            Tmax_this[species_list[i]]=Tmax_list[i]
    return Tmax_this

def get_difference(Cp_dict,T,Cp,species):
    cp_Temp=Cp_dict[species]
    T_di=[]
    Cp_that_di=[]
    Cp_this_di=[]
    for i in range(len(T)):
        T[i]=round(T[i],2)
        if T[i]%1.0==0.:
            T_di.append(float(T[i]))
            Cp_that_di.append(float(Cp[i]))
            Cp_this_di.append(float(cp_Temp[T[i]]))  

    # Cp_di=np.abs(np.array(Cp_that_di)-np.array(Cp_this_di))
    Cp_di=(np.abs(np.array(Cp_that_di)-np.array(Cp_this_di))/np.array(Cp_this_di))*100
    return T_di, Cp_di

def get_mean_max_difference(T_di,Tmax):
    diff_temp=[]
    diff_mean_max=[]
    T=T_di[0]
    differ=T_di[1]
    for i in range(len(T)):
        if T[i]<=Tmax:
            diff_temp.append(float(differ[i]))
    diff_mean_max.append(np.mean(diff_temp))   
    diff_mean_max.append(np.max(diff_temp)) 
    diff_mean_max.append(np.var(diff_temp))    
    return diff_mean_max

def create_dict(diff_dict,species,sourcename,diff):
    if species in diff_dict:
        diff_dict[species][sourcename]=diff
    else:
        diff_dict[species]=dict()
        diff_dict[species][sourcename]=diff
    return diff_dict    

def compare_difference(di_J_dict, di_nasa_dict, di_Ca_dict, T_di_Furtenbacher, T_di_SS_N, T_di_SS_P,Cp_dict,Tmax_this,Tmax_HITRAN):
    diff_dict=dict()
    for species in tqdm(Cp_dict):
        Tmax=Tmax_HITRAN[species]
        if species in Tmax_this:
            Tmax=Tmax_this[species]
        if species in di_J_dict:
            sourcename='J'
            T_J=di_J_dict[species]
            diff_J=get_mean_max_difference(T_J,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_J)   
        if species in di_nasa_dict:
            sourcename='N'
            T_N=di_nasa_dict[species]
            diff_N=get_mean_max_difference(T_N,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_N) 
        if species in di_Ca_dict:
            sourcename='Ca'
            T_Ca=di_Ca_dict[species]
            diff_Ca=get_mean_max_difference(T_Ca,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_Ca) 
        if species=='H2O':
            sourcename='Fur'
            diff_Fur=get_mean_max_difference(T_di_Furtenbacher,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_Fur)
        if species=='NH3':
            sourcename='SS'
            diff_SS_N=get_mean_max_difference(T_di_SS_N,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_SS_N)
        if species=='PH3':
            sourcename='SS'
            diff_SS_P=get_mean_max_difference(T_di_SS_P,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_SS_P) 
    return diff_dict        
# Using plots to compare the results with existing data
def compare_and_plot_Cp(Cp_dict,Tmax_dict,Cp_JANAF,coe_Capitelli,coe_nasa,Tmax_this):
    di_J_dict=dict()
    di_nasa_dict=dict()
    di_Ca_dict=dict()
    for species in tqdm(Cp_dict):
        other_data=False
        cp_Temp=Cp_dict[species]
        plt.figure(figsize=(8,6))
        Tmax=Tmax_dict[species]
        # The interval starts at 200K for gases and 298.15K for ions
        Tmin=200.
        if species=='NO+':
            Tmin=298.15

        # Plot JANAF specific heat
        if species in Cp_JANAF:

            other_data=True
            J_temp=Cp_JANAF[species]
            T_J=[]
            Cp_J=[]

            for t in J_temp:
                if (float(t) <=Tmax) & (float(t) >=Tmin):
                    T_J.append(float(t))
                    Cp_J.append(float(J_temp[t])) 
            T_J_di,di_J=get_difference(Cp_dict,T_J,Cp_J,species)
            T_di_J=np.vstack((T_J_di,di_J))
            if species not in di_J_dict:
                di_J_dict[species]=T_di_J
            plt.plot(T_J,Cp_J,color="blue",label="JANAF")   

        # plot specific heat using original nasa glenn polynomials
        if species in coe_nasa:
            other_data=True
            if Tmax<=1000:
                x_nasa=np.arange(Tmin, Tmax+1., 1.0)
                coefficient_nasa=coe_nasa[species][0]
                Cp_fit_nasa=get_polynomial_results(coefficient_nasa,x_nasa)
            else:
                x_nasa1=np.arange(Tmin, 1000, 1.0)
                coefficient_nasa1=coe_nasa[species][0]
                Cp_fit_nasa1=get_polynomial_results(coefficient_nasa1,x_nasa1)
                x_nasa2=np.arange(1000, Tmax+1., 1.0)
                coefficient_nasa2=coe_nasa[species][1]
                Cp_fit_nasa2=get_polynomial_results(coefficient_nasa2,x_nasa2)
                x_nasa=np.hstack((x_nasa1,x_nasa2))
                Cp_fit_nasa=np.hstack(( Cp_fit_nasa1, Cp_fit_nasa2))

            T_nasa_di,di_nasa=get_difference(Cp_dict,x_nasa,Cp_fit_nasa,species)
            T_di_nasa=np.vstack((T_nasa_di,di_nasa))
            if species not in di_nasa_dict:
                di_nasa_dict[species]=T_di_nasa
            plt.plot(x_nasa,Cp_fit_nasa,color="black",label="NASA Glenn") 



        if species in coe_Capitelli:
            other_data=True
            if Tmax<=1000:
                x_Ca=np.arange(Tmin, Tmax+1., 1.0)
                coefficient_Ca=coe_Capitelli[species][0]
                Cp_fit_Ca=get_polynomial_results(coefficient_Ca,x_Ca)
            else:
                x_Ca1=np.arange(Tmin, 1000, 1.0)
                coefficient_Ca1=coe_Capitelli[species][0]
                Cp_fit_Ca1=get_polynomial_results(coefficient_Ca1,x_Ca1)
                x_Ca2=np.arange(1000, 3000, 1.0)
                coefficient_Ca2=coe_Capitelli[species][1]
                Cp_fit_Ca2=get_polynomial_results(coefficient_Ca2,x_Ca2)
                x_Ca3=np.arange(3000, Tmax+1., 1.0)
                coefficient_Ca3=coe_Capitelli[species][2]
                Cp_fit_Ca3=get_polynomial_results(coefficient_Ca3,x_Ca3)
                x_Ca=np.hstack((x_Ca1,x_Ca2,x_Ca3))
                Cp_fit_Ca=np.hstack(( Cp_fit_Ca1, Cp_fit_Ca2,Cp_fit_Ca3))

            T_Ca_di,di_Ca=get_difference(Cp_dict,x_Ca,Cp_fit_Ca,species)
            T_di_Ca=np.vstack((T_Ca_di,di_Ca))
            if species not in di_Ca_dict:
                di_Ca_dict[species]=T_di_Ca
            plt.plot(x_Ca,Cp_fit_Ca,color="green",label="Capitelli et.al")  

        if species=='H2O':
            T_VT,Cp_VT=np.loadtxt("Other_Cp/H2O_Vidler_Tennyson.txt",usecols=(0,1),unpack=True) 
            T_Harris,Cp_Harris=np.loadtxt("Other_Cp/H2O_Harris.txt",usecols=(0,1),unpack=True) 
            T_Furtenbacher,Cp_Furtenbacher=np.loadtxt("Other_Cp/H2O_Furtenbacher.txt",usecols=(0,1),unpack=True) 
            T_VT_di,di_VT=get_difference(Cp_dict,T_VT,Cp_VT,species)
            T_Harris_di,di_Harris=get_difference(Cp_dict,T_Harris,Cp_Harris,species)
            T_Furtenbacher_di,di_Furtenbacher=get_difference(Cp_dict,T_Furtenbacher,Cp_Furtenbacher,species)
            # T_di_VT=np.vstack((T_VT_di,di_VT))
            # T_di_Harris=np.vstack((T_Harris_di,di_Harris))
            T_di_Furtenbacher=np.vstack((T_Furtenbacher_di,di_Furtenbacher))
            plt.plot(T_VT,Cp_VT,color='orange',label="Vidler & Tennyson")
            plt.plot(T_Harris,Cp_Harris,color='aqua',label="Harris et.al")
            plt.plot(T_Furtenbacher,Cp_Furtenbacher,color='purple',label="Furtenbacher et.al")

        if species=='NH3':
            other_data=True
            T_SS,Cp_SS=np.loadtxt("Other_Cp/NH3_Sousa_Silva.txt",usecols=(0,1),unpack=True) 
            T_SS_di,di_SS_N=get_difference(Cp_dict,T_SS,Cp_SS,species)
            T_di_SS_N=np.vstack((T_SS_di,di_SS_N))
            plt.plot(T_SS,Cp_SS,color='deepskyblue',label="C. Sousa-Silva et al.")

        if species=='PH3':
            other_data=True
            T_SS,Cp_SS=np.loadtxt("Other_Cp/PH3_Sousa_Silva.txt",usecols=(0,1),unpack=True) 
            T_SS_di,di_SS_P=get_difference(Cp_dict,T_SS,Cp_SS,species)
            T_di_SS_P=np.vstack((T_SS_di,di_SS_P))
            plt.plot(T_SS,Cp_SS,color='deepskyblue',label="C. Sousa-Silva et al.") 

        # Plot specific heat results in this work 

        T=[]
        Cp=[]
        for t in cp_Temp:
            if (float(t) <=Tmax) & (float(t) >=Tmin):
                T.append(float(t))
                Cp.append(float(cp_Temp[t]))
        plt.plot(T,Cp,color="red",label="This work")
     

        if species in Tmax_this:
            tmax=float(Tmax_this[species])
            plt.axvline(x=tmax,linestyle='--')  

        storepath="pictures"
        storename=storepath+"/"+species


        if species=='N2O':
            plt.legend(loc='lower left')
        else:
            plt.legend()


        plt.ylabel('$C_{p}$')
        plt.xlabel('T(K)')
        plt.title('Specific Heat Fit For '+species)
        plt.savefig(storename)
        plt.show() 

        if other_data:
            plt.figure(figsize=(8,6))
            if species in Cp_JANAF:
                plt.plot(T_J_di,di_J,color="blue",label="JANAF")  
            if species in coe_nasa:
                plt.plot(T_nasa_di,di_nasa,color="black",label="NASA Glenn") 
            if species in coe_Capitelli:
                plt.plot(T_Ca_di,di_Ca,color="green",label="Capitelli et.al") 
            if species=='H2O':
                plt.plot(T_VT_di,di_VT,color='orange',label="Vidler & Tennyson")
                plt.plot(T_Harris_di,di_Harris,color='aqua',label="Harris et.al")
                plt.plot(T_Furtenbacher_di,di_Furtenbacher,color='purple',label="Furtenbacher et.al")
            if species=='NH3':
                plt.plot(T_SS_di,di_SS_N,color='deepskyblue',label="C. Sousa-Silva et al.")
            if species=='PH3':
                plt.plot(T_SS_di,di_SS_P,color='deepskyblue',label="C. Sousa-Silva et al.")


            # plt.axhline(y=1.0,c='gray')
            # plt.axhline(y=2.0,c='gray')
            # plt.axhline(y=3.0,c='gray')
            # plt.axhline(y=5.0,c='gray')
            # storepath="Cp_Difference"
            # storepath="Difference"
            # storename=storepath+"/"+species        
            plt.legend()
            plt.ylabel('Difference(%)')
            plt.xlabel('T(K)')
            plt.title('Specific Heat Difference(%) For '+species)
            # plt.savefig(storename)
            plt.show() 

    return di_J_dict, di_nasa_dict, di_Ca_dict, T_di_Furtenbacher, T_di_SS_N, T_di_SS_P
# Define the nasa polynomial format using for fitting
def objective(x,a1,a2,a3,a4,a5,a6,a7):
    return R*(a1*x**(-2)+a2*x**(-1)+a3+a4*x+a5*x**2+a6*x**3+a7*x**4)
#%%
if __name__ == '__main__':

    # Calculate the specific heat based on the partition functions from HITRAN database
    R=8.314 # Universal gas constant
    split_dict=create_split_dict()
    Cp_dict,Tmax_dict=get_Cp(R,split_dict)

    # Read the specific heat results from other sources
    Cp_JANAF=get_JANAF()
    coe_Capitelli=get_Capitelli()
    coe_nasa=get_nasa(Cp_dict)
    Tmax_this=get_Tmax_this()
#%%
di_J_dict, di_nasa_dict, di_Ca_dict, T_di_Furtenbacher, T_di_SS_N, T_di_SS_P=compare_and_plot_Cp(Cp_dict,Tmax_dict,Cp_JANAF,coe_Capitelli,coe_nasa,Tmax_this)
#%%
diff_dict=compare_difference(di_J_dict, di_nasa_dict, di_Ca_dict, T_di_Furtenbacher, T_di_SS_N, T_di_SS_P,Cp_dict,Tmax_this,Tmax_dict)
#%%
diff_dict['SO3']
# %%
diff_dict
# %%
Tmax_this
# %%
