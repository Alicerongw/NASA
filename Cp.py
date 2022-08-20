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

def create_split_dict_theorets():
    # To improve the accuracy of polynomial fit, the temperature has been divided into several intervls
    # The split points are used for determining the interval
    split_list=[['GeH4']]
    split=[[1]]    
    split_dict=dict()
    for i in range(len(split_list)):
        for j in range(len(split_list[i])):
            if split_list[i][j] not in split_dict:
                split_dict[split_list[i][j]]=split[i]
    return split_dict


def get_Cp_theorets(R,split_dict):
    Cp_dict=dict()
    Tmax_dict=dict()
    path = "Other_Cp/Theorets"
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
            
    return Cp_dict


# Create the dictionary storing the split points for each species
def create_split_dict():
    # To improve the accuracy of polynomial fit, the temperature has been divided into several intervls
    # The split points are used for determining the interval
    split_list=[['CH4'],['H2O2'],['ClO'],['HCN','CO','HF','HI','N2','NH3','OCS','OH','PH3','C2H2','C2N2','C2H4','CH3Cl','CH3F','OH','O2'],['SF6'],['CS','HO2','NO','NO+','SO'],['H2','HBr','HCl'],['SO3'],['C2H6'],['C4H2']]
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
    path = "Partition_functions"
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
    Capitelli_specials=['CO','CO2','N2','NO','NO+','NO2','N2O','O2','O3','H2']
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

# Generate a dictionary storing the specific heat fit coefficients from Burcat
def get_Burcat(): 
    Burcat=np.loadtxt("Other_Cp/Burcat_fit.csv",delimiter=",",usecols=(3,4,5,6,7,8,9))
    Burcat_specials=['CH3Br','CH3I','CH3OH','ClO','ClONO2','COF2','GeH4','H2O2','HBr','HNO3','HO2','HOBr','HOCl','NH3','PH3','OH']
    coe_Burcat=dict()
    for i in range(len(Burcat_specials)):
        n=i*2
        species=Burcat_specials[i]
        coefficients_Burcat = np.zeros([2, 7])
        for ii in range(2):
            coefficients_Burcat[ii,]=Burcat[n]
            n+=1
        if species not in coe_Burcat:
            coe_Burcat[species]=coefficients_Burcat
    return coe_Burcat

# Generate a dictionary storing the specific heat fit coefficients from Barklem
def get_Barklem(): 
    Barklem=np.loadtxt("Other_Cp/Barklem_fit.csv",delimiter=",",usecols=(1,2,3,4,5,6,7))
    Barklem_specials=['ClO','CO','CS','H2','HBr','HCl','HF','HI','N2','NO','O2','OH','SO']
    coe_Barklem=dict()
    for i in range(len(Barklem_specials)):
        n=i*2
        species=Barklem_specials[i]
        coefficients_Barklem = np.zeros([2, 7])
        for ii in range(2):
            coefficients_Barklem[ii,]=Barklem[n]
            n+=1
        if species not in coe_Barklem:
            coe_Barklem[species]=coefficients_Barklem
    return coe_Barklem

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

def get_Tmax_this(Tmax_hitran):
    Tmax_this=dict()
    species_list=['C2H2','ClO','CO2','H2O','H2S','HCN','HI','HO2','N2','N2O','NH3','O2','PH3','SO','SO3']
    Tmax_list=[1400.,1000.,3000.,3500.,2250.,1700.,3100.,1900.,2500.,1500.,1700.,500.,2200.,1500.,1400.]
    for species in Tmax_hitran:
        Tmax=Tmax_hitran[species]
        if species in species_list:
            index= species_list.index(species)
            Tmax=Tmax_list[index]
        if species not in Tmax_this:
            Tmax_this[species]=Tmax
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

def get_mean_max_var_difference(T_di,Tmax):
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

def compare_difference(di_J_dict, di_nasa_dict, di_Burcat_dict,di_Ca_dict, T_di_Furtenbacher, T_di_Furtenbacher_O2,T_di_SS_N, T_di_SS_P,Cp_dict,Tmax_this):
    diff_dict=dict()
    for species in tqdm(Cp_dict):
        Tmax=Tmax_this[species]
        if species in di_J_dict:
            sourcename='J'
            T_J=di_J_dict[species]
            diff_J=get_mean_max_var_difference(T_J,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_J)   
        if species in di_nasa_dict:
            sourcename='N'
            T_N=di_nasa_dict[species]
            diff_N=get_mean_max_var_difference(T_N,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_N) 
        if species in di_Burcat_dict:
            sourcename='B'
            T_B=di_Burcat_dict[species]
            diff_B=get_mean_max_var_difference(T_B,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_B) 
        if species in di_Ca_dict:
            sourcename='Ca'
            T_Ca=di_Ca_dict[species]
            diff_Ca=get_mean_max_var_difference(T_Ca,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_Ca) 
        if species=='H2O':
            sourcename='Fur'
            diff_Fur=get_mean_max_var_difference(T_di_Furtenbacher,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_Fur)
        if species=='O2':
            sourcename='Fur'
            diff_Fur_O2=get_mean_max_var_difference(T_di_Furtenbacher_O2,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_Fur_O2)
        if species=='NH3':
            sourcename='SS'
            diff_SS_N=get_mean_max_var_difference(T_di_SS_N,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_SS_N)
        if species=='PH3':
            sourcename='SS'
            diff_SS_P=get_mean_max_var_difference(T_di_SS_P,Tmax)
            diff_dict=create_dict(diff_dict,species,sourcename,diff_SS_P) 
    return diff_dict  

# Using plots to compare the results with existing data
def compare_and_plot_Cp(Cp_dict,Tmax_dict,Cp_JANAF,coe_Capitelli,coe_nasa,coe_Burcat,coe_Barklem,Tmax_this,Cp_theorets):
    di_J_dict=dict()
    di_nasa_dict=dict()
    di_Ca_dict=dict()
    di_Burcat_dict=dict()
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

        # plot specific heat using original nasa glenn polynomials
        if species in coe_Burcat:
            other_data=True
            if Tmax<=1000:
                x_Burcat=np.arange(Tmin, Tmax+1., 1.0)
                coefficient_Burcat=coe_Burcat[species][0]
                Cp_fit_Burcat=get_polynomial_results(coefficient_Burcat,x_Burcat)
            else:
                x_Burcat1=np.arange(Tmin, 1000, 1.0)
                coefficient_Burcat1=coe_Burcat[species][0]
                Cp_fit_Burcat1=get_polynomial_results(coefficient_Burcat1,x_Burcat1)
                x_Burcat2=np.arange(1000, Tmax+1., 1.0)
                coefficient_Burcat2=coe_Burcat[species][1]
                Cp_fit_Burcat2=get_polynomial_results(coefficient_Burcat2,x_Burcat2)
                x_Burcat=np.hstack((x_Burcat1,x_Burcat2))
                Cp_fit_Burcat=np.hstack(( Cp_fit_Burcat1, Cp_fit_Burcat2))

            T_Burcat_di,di_Burcat=get_difference(Cp_dict,x_Burcat,Cp_fit_Burcat,species)
            T_di_Burcat=np.vstack((T_Burcat_di,di_Burcat))
            if species not in di_Burcat_dict:
                di_Burcat_dict[species]=T_di_Burcat
            plt.plot(x_Burcat,Cp_fit_Burcat,color='green',label="Burcat") 

        # plot specific heat using original nasa glenn polynomials
        if species in coe_Barklem:
            other_data=True
            if Tmax<=1000:
                x_Barklem=np.arange(298.15, Tmax+1., 1.0)
                coefficient_Barklem=coe_Barklem[species][0]
                Cp_fit_Barklem=get_polynomial_results(coefficient_Barklem,x_Barklem)
            else:
                x_Barklem1=np.arange(298.15, 1000, 1.0)
                coefficient_Barklem1=coe_Barklem[species][0]
                Cp_fit_Barklem1=get_polynomial_results(coefficient_Barklem1,x_Barklem1)
                x_Barklem2=np.arange(1000, Tmax+1., 1.0)
                coefficient_Barklem2=coe_Barklem[species][1]
                Cp_fit_Barklem2=get_polynomial_results(coefficient_Barklem2,x_Barklem2)
                x_Barklem=np.hstack((x_Barklem1,x_Barklem2))
                Cp_fit_Barklem=np.hstack(( Cp_fit_Barklem1, Cp_fit_Barklem2))

            # T_Burcat_di,di_Burcat=get_difference(Cp_dict,x_Burcat,Cp_fit_Burcat,species)
            # T_di_Burcat=np.vstack((T_Burcat_di,di_Burcat))
            # if species not in di_Burcat_dict:
            #     di_Burcat_dict[species]=T_di_Burcat
            plt.plot(x_Barklem,Cp_fit_Barklem,color='orange',label="Barklem") 

        if species in coe_Capitelli:
            other_data=True
            T_middle=1000
            if species =='H2':
                T_middle=500
            if Tmax<=1000:
                x_Ca=np.arange(Tmin, Tmax+1., 1.0)
                coefficient_Ca=coe_Capitelli[species][0]
                Cp_fit_Ca=get_polynomial_results(coefficient_Ca,x_Ca)
            else:
                x_Ca1=np.arange(Tmin, T_middle, 1.0)
                coefficient_Ca1=coe_Capitelli[species][0]
                Cp_fit_Ca1=get_polynomial_results(coefficient_Ca1,x_Ca1)
                x_Ca2=np.arange(T_middle, 3000, 1.0)
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
            T_di_Furtenbacher=np.vstack((T_Furtenbacher_di,di_Furtenbacher))
            T_Ruscic,Cp_Ruscic=np.loadtxt("Other_Cp/H2O_Ruscic.txt",usecols=(0,1),unpack=True) 
            plt.plot(T_VT,Cp_VT,color='orange',label="Vidler & Tennyson")
            plt.plot(T_Harris,Cp_Harris,color='aqua',label="Harris et.al")
            plt.plot(T_Furtenbacher,Cp_Furtenbacher,color='purple',label="Furtenbacher et.al")
            plt.plot(T_Ruscic,Cp_Ruscic,color='green',label="Ruscic")

        if species=='O2':
            T_Furtenbacher,Cp_Furtenbacher=np.loadtxt("Other_Cp/O2_Furtenbacher.txt",usecols=(0,1),unpack=True) 
            T_Furtenbacher_di_O2,di_Furtenbacher_O2=get_difference(Cp_dict,T_Furtenbacher,Cp_Furtenbacher,species)
            T_di_Furtenbacher_O2=np.vstack((T_Furtenbacher_di_O2,di_Furtenbacher_O2))
            plt.plot(T_Furtenbacher,Cp_Furtenbacher,color='purple',label="Furtenbacher et.al")
            T_Gurvich,Cp_Gurvich=np.loadtxt("Other_Cp/O2_Gurvich.txt",usecols=(0,1),unpack=True)
            plt.scatter(T_Gurvich,Cp_Gurvich,color="hotpink",marker="x",label="Gurvich")
            T_Jaffe,Cp_Jaffe=np.loadtxt("Other_Cp/O2_Jaffe.txt",usecols=(0,1),unpack=True)
            plt.scatter(T_Jaffe,Cp_Jaffe,color="gray",marker="x",label="Jaffe")
        if species=='N2':
            T_Jaffe,Cp_Jaffe=np.loadtxt("Other_Cp/N2_Jaffe.txt",usecols=(0,1),unpack=True)
            plt.scatter(T_Jaffe,Cp_Jaffe,color="gray",marker="x",label="Jaffe")
        if species=='NO':
            T_Jaffe,Cp_Jaffe=np.loadtxt("Other_Cp/NO_Jaffe.txt",usecols=(0,1),unpack=True)
            plt.scatter(T_Jaffe,Cp_Jaffe,color="gray",marker="x",label="Jaffe")
        if species=='H2O2':
            T_Chao,Cp_Chao=np.loadtxt("Other_Cp/H2O2_Chao.txt",usecols=(0,1),unpack=True)
            plt.plot(T_Chao,Cp_Chao,color="orange",label="Chao et.al")

        if species=='CH3OH':
            T_Chao,Cp_Chao=np.loadtxt("Other_Cp/CH3OH_Chao.txt",usecols=(0,1),unpack=True)
            plt.plot(T_Chao,Cp_Chao,color="orange",label="Chao et.al")

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

        # Plot JANAF specific heat
        if species in Cp_theorets:

            other_data=True
            th_temp=Cp_theorets[species]
            T_th=[]
            Cp_th=[]

            for t in th_temp:
                if (float(t) <=Tmax) & (float(t) >=Tmin):
                    T_th.append(float(t))
                    Cp_th.append(float(th_temp[t])) 
            # T_th_di,di_th=get_difference(Cp_dict,T_th,Cp_th,species)
            # T_di_th=np.vstack((T_th_di,di_th))
            # if species not in di_th_dict:
            #     di_th_dict[species]=T_di_th
            plt.plot(T_th,Cp_th,color="blue",label="Theorets")   

        # Plot specific heat results in this work 

        T=[]
        Cp=[]
        for t in cp_Temp:
            if (float(t) <=Tmax) & (float(t) >=Tmin):
                T.append(float(t))
                Cp.append(float(cp_Temp[t]))
        plt.plot(T,Cp,color="red",label="This work")
     

        if Tmax_this[species]!=Tmax:
            tmax=float(Tmax_this[species])
            plt.axvline(x=tmax,linestyle='--')  

        storepath="Cp_pictures"
        storename=storepath+"/"+species


        if species=='N2O':
            plt.legend(loc='lower left')
        else:
            plt.legend()


        plt.ylabel('$C_{p}$')
        plt.xlabel('T(K)')
        plt.title('Specific Heat Fit For '+species)
        # plt.savefig(storename)
        plt.show() 

        if other_data:
            plt.figure(figsize=(8,6))
            if species in Cp_JANAF:
                plt.plot(T_J_di,di_J,color="blue",label="JANAF")  
            if species in coe_nasa:
                plt.plot(T_nasa_di,di_nasa,color="black",label="NASA Glenn") 
            if species in coe_Capitelli:
                plt.plot(T_Ca_di,di_Ca,color="green",label="Capitelli et.al") 
            if species in coe_Burcat:
                plt.plot(T_Burcat_di,di_Burcat,color='green',label="Burcat")                 
            if species=='H2O':
                plt.plot(T_VT_di,di_VT,color='orange',label="Vidler & Tennyson")
                plt.plot(T_Harris_di,di_Harris,color='aqua',label="Harris et.al")
                plt.plot(T_Furtenbacher_di,di_Furtenbacher,color='purple',label="Furtenbacher et.al")
            if species=='O2':
                plt.plot(T_Furtenbacher_di_O2,di_Furtenbacher_O2,color='purple',label="Furtenbacher et.al")
            if species=='NH3':
                plt.plot(T_SS_di,di_SS_N,color='deepskyblue',label="C. Sousa-Silva et al.")
            if species=='PH3':
                plt.plot(T_SS_di,di_SS_P,color='deepskyblue',label="C. Sousa-Silva et al.")


            # plt.axhline(y=1.0,c='gray')
            # plt.axhline(y=2.0,c='gray')
            # plt.axhline(y=3.0,c='gray')
            # plt.axhline(y=5.0,c='gray')
            # storepath="Cp_Difference"
            storepath="Difference(%)"
            storename=storepath+"/"+species        
            plt.legend()
            plt.ylabel('Difference(%)')
            plt.xlabel('T(K)')
            plt.title('Specific Heat Difference(%) For '+species)
            # plt.savefig(storename)
            plt.show() 

    return di_J_dict, di_nasa_dict,di_Burcat_dict, di_Ca_dict, T_di_Furtenbacher, T_di_Furtenbacher_O2,T_di_SS_N, T_di_SS_P

# Define the nasa polynomial format using for fitting
def objective(x,a1,a2,a3,a4,a5,a6,a7):
    return R*(a1*x**(-2)+a2*x**(-1)+a3+a4*x+a5*x**2+a6*x**3+a7*x**4)

def get_coefficients(Cp_dict,Tmax_this):
    fit_dictionary=dict()
    for species in Cp_dict:
        Cpdict_temp=Cp_dict[species]            
        list_sorted=sorted(Cpdict_temp.items(),key=lambda item:item[0])
        T_this=[]
        Cp_this=[]
        for i in range(len(list_sorted)):
            T_this.append(list_sorted[i][0])
            Cp_this.append(list_sorted[i][1])

        Tmax=Tmax_this[species]

        Tmin_index=0
        if species=='NO+':
            Tmin_index=98
        
        # Fit the specific heat
        fit_Tinterval=[]
        fit_Cpinterval=[]
        if Tmax<=1000:
            temp_intervals=1
            fit_Tinterval.append(T_this[Tmin_index:int(Tmax-199)])
            fit_Cpinterval.append(Cp_this[Tmin_index:int(Tmax-199)])    
        else:
            temp_intervals=2
            fit_Tinterval.append(T_this[Tmin_index:801])
            fit_Cpinterval.append(Cp_this[Tmin_index:801])
            fit_Tinterval.append(T_this[800:int(Tmax-199)])
            fit_Cpinterval.append(Cp_this[800:int(Tmax-199)])

        coefficients = np.zeros([temp_intervals, 7])
        for ii in range(len(fit_Tinterval)):
            x=np.array(fit_Tinterval[ii])
            y=np.array(fit_Cpinterval[ii])
            popt, _ = curve_fit(objective, x, y)
            for i in range(7):
                coefficients[ii,i]=popt[i]

        fit_dictionary[species]=coefficients
    return fit_dictionary

def store_coefficients(Tmax_dict,fit_dictionary):
    molecule_list=[]
    Tmin_list=[]
    Tmax_list=[]
    a1_list=[]
    a2_list=[]
    a3_list=[]
    a4_list=[]
    a5_list=[]
    a6_list=[]
    a7_list=[]

    for molecule in fit_dictionary:
        Tmin=200.
        if molecule=='NO+':
            Tmin=298.15
        coefficient_temp=fit_dictionary[molecule]
        Tmax=Tmax_dict[molecule]
        if Tmax<=1000:
            temp_intervals=1
        else:
            temp_intervals=2
        
        for i in range(temp_intervals):
            if i ==0:
                Tmin_list.append(Tmin)
                if temp_intervals==1:
                    Tmax_list.append(Tmax)
                else:
                    Tmax_list.append(1000)
            if i ==1:
                Tmin_list.append(1000)
                Tmax_list.append(Tmax)
            molecule_list.append(molecule)
            a1_list.append(coefficient_temp[i][0])
            a2_list.append(coefficient_temp[i][1])
            a3_list.append(coefficient_temp[i][2])
            a4_list.append(coefficient_temp[i][3])
            a5_list.append(coefficient_temp[i][4])
            a6_list.append(coefficient_temp[i][5])
            a7_list.append(coefficient_temp[i][6])

        
    # fit_data = pd.DataFrame({'Molecule': molecule_list, 'Tmin':Tmin_list,'Tmax': Tmax_list,'a1':a1_list,'a2': a2_list,'a3': a3_list,'a4': a4_list,'a5': a5_list,'a6': a6_list,'a7': a7_list})
    # fit_data.to_csv("fit_coefficients.csv",index=False)

def calculate_and_plot_residuals(Tmax_this,Cp_dict,fit_dictionary):
    residual_dict=dict()
    for molecule in tqdm(Cp_dict):
        Tmax=Tmax_this[molecule]
        Tmin=200.
        if molecule=='NO+':
            Tmin=298
        # plot specific heat in this work after fitting
        cp_Temp=Cp_dict[molecule]
        T=[]
        Cp=[]
        for t in cp_Temp:
            if (float(t) <=Tmax) & (float(t) >=Tmin):
                T.append(float(t))
                Cp.append(float(cp_Temp[t]))


        #plot specific heat in this work after fitting
        if Tmax<=1000:
            x=np.arange(Tmin,Tmax+1., 1.)
            coefficient_this=fit_dictionary[molecule][0]
            Cp_fit_this=R*(coefficient_this[0]*x**(-2)+coefficient_this[1]*x**(-1)+coefficient_this[2]+coefficient_this[3]*x+coefficient_this[4]*x**2+coefficient_this[5]*x**3+coefficient_this[6]*x**4)
        
        else:
            x1=np.arange(Tmin,1000.,1.)
            coefficient_this1=fit_dictionary[molecule][0]
            Cp_fit_this1=R*(coefficient_this1[0]*x1**(-2)+coefficient_this1[1]*x1**(-1)+coefficient_this1[2]+coefficient_this1[3]*x1+coefficient_this1[4]*x1**2+coefficient_this1[5]*x1**3+coefficient_this1[6]*x1**4)
            x2=np.arange(1000.,Tmax+1., 1.)
            coefficient_this2=fit_dictionary[molecule][1]
            Cp_fit_this2=R*(coefficient_this2[0]*x2**(-2)+coefficient_this2[1]*x2**(-1)+coefficient_this2[2]+coefficient_this2[3]*x2+coefficient_this2[4]*x2**2+coefficient_this2[5]*x2**3+coefficient_this2[6]*x2**4)
            x=np.hstack((x1,x2))
            Cp_fit_this=np.hstack(( Cp_fit_this1, Cp_fit_this2))

        residual= Cp - Cp_fit_this
        array_temp=np.zeros(3)
        array_temp[0]=np.min(np.abs(residual))
        array_temp[1]=np.max(np.abs(residual))
        array_temp[2]=np.mean(np.abs(residual))

        if molecule not in residual_dict:
            residual_dict[molecule]=array_temp

        plt.figure(figsize=(8,6))
        fig,ax1 = plt.subplots()
        ax2 = ax1.twinx()          
        ax1.plot(x,Cp,color="blue",label="Original Cp")
        ax1.plot(x,Cp_fit_this,color="red",label="Fitted Cp")
        ax2.plot(x,residual,color="black")
        
        ax1.set_xlabel('T(k)')    
        ax1.set_ylabel('Cp',color = 'red')   
        ax2.set_ylabel('Residuals',color = 'black')

        plt.axhline(y=0.0,c="gray")
        # plt.plot(x,residual)    

        storepath="Fit_residuals"
        storename=storepath+"/"+molecule

        ax1.legend()
        # plt.ylabel('Residuals')
        # plt.xlabel('T(K)')
        plt.title('Residuals Plot For '+molecule)
        # plt.savefig(storename)
        plt.show()
    return residual_dict
#%%
if __name__ == '__main__':

    # Calculate the specific heat based on the partition functions from HITRAN database
    R=8.314 # Universal gas constant
    split_dict=create_split_dict()
    Cp_dict,Tmax_hitran=get_Cp(R,split_dict)

    # Read the specific heat results from other sources
    split_dict_theorets=create_split_dict_theorets()
    Cp_theorets=get_Cp_theorets(R,split_dict_theorets)
    Cp_JANAF=get_JANAF()
    coe_Capitelli=get_Capitelli()
    coe_Burcat=get_Burcat()
    coe_nasa=get_nasa(Cp_dict)
    coe_Barklem=get_Barklem()    
    Tmax_this=get_Tmax_this(Tmax_hitran)
    di_J_dict, di_nasa_dict, di_Burcat_dict,di_Ca_dict, T_di_Furtenbacher, T_di_Furtenbacher_O2, T_di_SS_N, T_di_SS_P=compare_and_plot_Cp(Cp_dict,Tmax_hitran,Cp_JANAF,coe_Capitelli,coe_nasa,coe_Burcat,coe_Barklem,Tmax_this,Cp_theorets)
    diff_dict=compare_difference(di_J_dict, di_nasa_dict,di_Burcat_dict, di_Ca_dict, T_di_Furtenbacher, T_di_Furtenbacher_O2,T_di_SS_N, T_di_SS_P,Cp_dict,Tmax_this)
    fit_dictionary=get_coefficients(Cp_dict,Tmax_this)
    store_coefficients(Tmax_this,fit_dictionary)
    residual_dict=calculate_and_plot_residuals(Tmax_this,Cp_dict,fit_dictionary)

#%%
diff_dict=compare_difference(di_J_dict, di_nasa_dict,di_Burcat_dict, di_Ca_dict, T_di_Furtenbacher, T_di_Furtenbacher_O2,T_di_SS_N, T_di_SS_P,Cp_dict,Tmax_this)
coe_Burcat.keys()
# %%
coe_Burcat['GeH4']
# %%
diff_dict['HBr']
# %%
diff_dict['CS']
# %%
diff_dict['HNO3']
# %%
coe_Capitelli=get_Capitelli()
coe_Capitelli['CO']
# %%
coe_Barklem['HI']
# %%


# %%
R=8.314 # Universal gas constant
split_dict_theorets=create_split_dict_theorets()
Cp_theorets=get_Cp_theorets(R,split_dict_theorets)
# %%
Cp_theorets
# %%
Cp_JANAF['H2O2']
# %%
