#%%
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pandas as pd
import json
import warnings
warnings.simplefilter('ignore', np.RankWarning)
# %%
# Add the split points for molecules in HITRAN  
split_list=[['CH4'],['H2O2'],['CLO'],['HCN','CO','HF','HI','N2','NH3','OCS','OH','PH3','C2H2','C2N2','C2H4','CH3CL','CH3F','OH','O2'],['SF6'],['CS','HO2','NO','NO+','SO'],['H2','HBr','HCL'],['SO3'],['C2H6'],['C4H2']]
split=[[200,500,1300,1500],[200,1500],[200,4000],[200,1000],[200,1000,2000,3000,4000],[200,1000,4000],[200,1000,5000],[200,500,650],[200,1000,2000,3000],[200,1000,2000]]
split_dict=dict()
for i in range(len(split_list)):
    for j in range(len(split_list[i])):
        if split_list[i][j] not in split_dict:
            split_dict[split_list[i][j]]=split[i]
# %%
# Calculate specific heat
Cp_dict=dict()
Tmax_dict=dict()
fit_dictionary=dict()
path = "data"
files= os.listdir(path) 
for file in tqdm(files): 
     if ".txt" in file:
        if not os.path.isdir(file): 
            filename=path+"/"+file
            molecule=file.split(".")[0]
            print(molecule)
            storepath="Cp_results"
            storename=storepath+"/"+molecule+".csv"
            # Read partition functions from database
            T,Q=np.loadtxt(filename,usecols=(0,1),unpack=True)

            #get interval 
            if molecule in split_dict:
                split=split_dict[molecule]
            else:
                split=[200]

            # devide the original data
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

            # Calculate Cp
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
                R=8.314 # Universal gas constant
                Cp_temp=R*((Q2/y)-(Q1/y)**2)+(5/2)*R
                for i in range(len(x)):
                    t=float(x[i])    
                    cp=float(Cp_temp[i])  
                    if molecule in Cp_dict:
                        Cp_dict[molecule][t]=cp
                    else:
                        Cp_dict[molecule]=dict()
                        Cp_dict[molecule][t]=cp 
            Cpdict_temp=Cp_dict[molecule]            
            list_sorted=sorted(Cpdict_temp.items(),key=lambda item:item[0])
            T_this=[]
            Cp_this=[]
            for i in range(len(list_sorted)):
                T_this.append(list_sorted[i][0])
                Cp_this.append(list_sorted[i][1])
            Cp_df=pd.DataFrame({'T(K)': T_this, 'Cp_Thiswork': Cp_this})
            Cp_df.to_csv(storename,index=False)

            # Store the maximum temperature in HITRAN database
            Tmax_this=T_this[-1]
            Tmax=min(Tmax_this,6000)
            Tmax_dict[molecule]=Tmax
            Tmin_index=0
            if molecule=='NO+':
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
                Cp_interval=np.array(fit_Cpinterval[ii])
                y=(Cp_interval/R)*x*x
                z1 = np.polyfit(x,y,6) 
                p1= np.poly1d(z1)
                for i in range(7):
                    coefficients[ii,i]=p1[i]
            fit_dictionary[molecule]=coefficients




# filename='Cp_thiswork.json'
# with open(filename,'w') as file_obj:
#         json.dump(Cp_dict,file_obj)

# filename1='Tmax_dict.json'
# with open(filename1,'w') as file_obj1:
#         json.dump(Tmax_dict,file_obj1)

#%%
# Store the fit coefficients in a csv file
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

    
fit_data = pd.DataFrame({'Molecule': molecule_list, 'Tmin':Tmin_list,'Tmax': Tmax_list,'a1':a1_list,'a2': a2_list,'a3': a3_list,'a4': a4_list,'a5': a5_list,'a6': a6_list,'a7': a7_list})
fit_data.to_csv("coefficients.csv",index=False)


# %%
# Generate a JANAF dictionary 
path = "JANAF"
files= os.listdir(path) 
Cp_JANAF=dict()
for file in tqdm(files): 
     if ".txt" in file:
        if not os.path.isdir(file): 
            filename=path+"/"+file
            molecule=file.split(".")[0]
            print(molecule)
            # Read partition functions from database
            T,Cp=np.loadtxt(filename,usecols=(0,1),unpack=True)

            for i in range(len(T)):
                t=T[i]    
                cp=float(Cp[i])  
                if molecule in Cp_JANAF:
                    Cp_JANAF[molecule][t]=cp
                else:
                    Cp_JANAF[molecule]=dict()
                    Cp_JANAF[molecule][t]=cp  


# filename2='Cp_JANAFdict.json'
# with open(filename2,'w') as file_obj2:
#         json.dump(Cp_JANAF,file_obj2)
# %%
# Generate a Capitelli dictionary 
path = "Capitelli"
files= os.listdir(path) 
Cp_Capitelli=dict()
for file in tqdm(files): 
     if ".txt" in file:
        if not os.path.isdir(file): 
            filename=path+"/"+file
            molecule=file.split(".")[0]
            print(molecule)
            # Read partition functions from database
            T,Cp=np.loadtxt(filename,usecols=(0,1),unpack=True)

            for i in range(len(T)):
                t=T[i]    
                cp=float(Cp[i])  
                if molecule in Cp_Capitelli:
                    Cp_Capitelli[molecule][t]=cp
                else:
                    Cp_Capitelli[molecule]=dict()
                    Cp_Capitelli[molecule][t]=cp
# %%
# Generate a Capitelli dictionary 
# Capitelli=np.loadtxt("Capitelli_fit.csv",delimiter=",",usecols=(3,4,5,6,7,8,9))
# Capitelli_specials=['CO','CO2','N2','NO','NO2','N2O','O2','O3']
# Capitelli_fit=dict()
# for i in range(len(Capitelli_specials)):
#     n=i*3
#     molecule=Capitelli_specials[i]
#     coefficients_Ca = np.zeros([3, 7])
#     for ii in range(3):
#         coefficients_Ca[ii,]=Capitelli[n]
#         n+=1
#     if molecule not in Capitelli_fit:
#         Capitelli_fit[molecule]=coefficients_Ca
# %%
# Read coefficients from original nasa glenn polynomials
temprange_nasa=dict()
fit_nasa=dict()
with open('nasa9.dat', 'r') as nasa_file:
    keep_searching = True

    temp_range = None
    temp_intervals = None
    coefficients = None
    


    while keep_searching:
        line = nasa_file.readline()
        
        if line=='':
            keep_searching = False
            


        elif (line.split()[0] in Cp_dict) or (line.split()[0]=='C2H2,acetylene'):
            print(line.split()[0])
            if line.split()[0]=='C2H2,acetylene':
                species='C2H2'
            else:
                species=line.split()[0]

            # Get reference for the species data:
            data = nasa_file.readline()
            temp_intervals = int(data[0:2])

            temp_range = np.zeros([temp_intervals, 2])
            coefficients = np.zeros([temp_intervals, 7])
            integration_cts = np.zeros([temp_intervals, 2])
            
            for ii in range(temp_intervals):
                temprange_data = nasa_file.readline()
                coeffs_line_1 = nasa_file.readline()
                coeffs_line_2 = nasa_file.readline()

                temp_range[ii][:]=[1,2]
                temp_range[ii,:] = [float(temprange_data[0:11]), float(temprange_data[11:22])]
                coefficients[ii,:] = [float(coeffs_line_1[0:16].replace('D','E')), float(coeffs_line_1[16:32].replace('D','E')),
                                float(coeffs_line_1[32:48].replace('D','E')), float(coeffs_line_1[48:64].replace('D','E')),
                                float(coeffs_line_1[64:80].replace('D','E')), float(coeffs_line_2[0:16].replace('D','E')),
                                float(coeffs_line_2[16:32].replace('D','E'))]
            
            temprange_nasa[species]=temp_range
            fit_nasa[species]=coefficients
                       
        else:
            pass

# %%
# Compare the specific heat results with existing sources
for molecule in tqdm(Cp_dict):
    Tmax=Tmax_dict[molecule]
    Tmin=200.
    if molecule=='NO+':
        Tmin=298.15
    # plot specific heat in this work after fitting
    cp_Temp=Cp_dict[molecule]
    T=[]
    Cp=[]
    for t in cp_Temp:
        if (float(t) <=Tmax) & (float(t) >=Tmin):
            T.append(float(t))
            Cp.append(float(cp_Temp[t]))
    plt.plot(T,Cp,color="violet",label="This work")

    # plot JANAF specific heat
    if molecule in Cp_JANAF:

        J_temp=Cp_JANAF[molecule]
        T_J=[]
        Cp_J=[]
        for t in J_temp:
            if (float(t) <=Tmax) & (float(t) >=Tmin):
                T_J.append(float(t))
                Cp_J.append(float(J_temp[t])) 
        plt.plot(T_J,Cp_J,color="blue",label="JANAF")  

    if molecule in Cp_Capitelli:

        Ca_temp=Cp_Capitelli[molecule]
        T_Ca=[]
        Cp_Ca=[]
        for t in Ca_temp:
            if (float(t) <=Tmax) & (float(t) >=Tmin):
                T_Ca.append(float(t))
                Cp_Ca.append(float(Ca_temp[t])) 
        plt.scatter(T_Ca,Cp_Ca,color="green",marker="x",label="Capitelli et.al")

    #plot specific heat in this work after fitting
    if Tmax<=1000:
        x=np.arange(Tmin, Tmax, 0.01)
        coefficient_this=fit_dictionary[molecule][0]
        Cp_fit_this=R*(coefficient_this[0]*x**(-2)+coefficient_this[1]*x**(-1)+coefficient_this[2]+coefficient_this[3]*x+coefficient_this[4]*x**2+coefficient_this[5]*x**3+coefficient_this[6]*x**4)
    
    else:
        x1=np.arange(Tmin, 1000, 0.01)
        coefficient_this1=fit_dictionary[molecule][0]
        Cp_fit_this1=R*(coefficient_this1[0]*x1**(-2)+coefficient_this1[1]*x1**(-1)+coefficient_this1[2]+coefficient_this1[3]*x1+coefficient_this1[4]*x1**2+coefficient_this1[5]*x1**3+coefficient_this1[6]*x1**4)
        if molecule in fit_nasa:
            coefficient_nasa1=fit_nasa[molecule][0]
            Cp_fit_nasa1=R*(coefficient_nasa1[0]*x**(-2)+coefficient_nasa1[1]*x**(-1)+coefficient_nasa1[2]+coefficient_nasa1[3]*x+coefficient_nasa1[4]*x**2+coefficient_nasa1[5]*x**3+coefficient_nasa1[6]*x**4)
        x2=np.arange(1000, Tmax, 0.01)
        coefficient_this2=fit_dictionary[molecule][1]
        Cp_fit_this2=R*(coefficient_this2[0]*x2**(-2)+coefficient_this2[1]*x2**(-1)+coefficient_this2[2]+coefficient_this2[3]*x2+coefficient_this2[4]*x2**2+coefficient_this2[5]*x2**3+coefficient_this2[6]*x2**4)
        if molecule in fit_nasa:
            coefficient_nasa2=fit_nasa[molecule][1]
            Cp_fit_nasa2=R*(coefficient_nasa2[0]*x**(-2)+coefficient_nasa2[1]*x**(-1)+coefficient_nasa2[2]+coefficient_nasa2[3]*x+coefficient_nasa2[4]*x**2+coefficient_nasa2[5]*x**3+coefficient_nasa2[6]*x**4)
        x=np.hstack((x1,x2))
        Cp_fit_this=np.hstack(( Cp_fit_this1, Cp_fit_this2))

    plt.plot(x,Cp_fit_this,color="red",label="This work_Fit")    

    # if molecule in Capitelli_fit:
    #     if Tmax<=1000:
    #         x_Ca=np.arange(298.15, Tmax, 0.01)
    #         coefficient_Ca=Capitelli_fit[molecule][0]
    #         Cp_fit_Ca=R*(coefficient_Ca[0]*x_Ca**(-2)+coefficient_Ca[1]*x_Ca**(-1)+coefficient_Ca[2]+coefficient_Ca[3]*x_Ca+coefficient_Ca[4]*x_Ca**2+coefficient_Ca[5]*x_Ca**3+coefficient_Ca[6]*x_Ca**4)
    #     else:
    #         x_Ca1=np.arange(298.15, 1000, 0.01)
    #         coefficient_Ca1=Capitelli_fit[molecule][0]
    #         Cp_fit_Ca1=R*(coefficient_Ca1[0]*x_Ca1**(-2)+coefficient_Ca1[1]*x_Ca1**(-1)+coefficient_Ca1[2]+coefficient_Ca1[3]*x_Ca1+coefficient_Ca1[4]*x_Ca1**2+coefficient_Ca1[5]*x_Ca1**3+coefficient_Ca1[6]*x_Ca1**4)
    #         x_Ca2=np.arange(1000, 3000, 0.01)
    #         coefficient_Ca2=Capitelli_fit[molecule][1]
    #         Cp_fit_Ca2=R*(coefficient_Ca2[0]*x_Ca2**(-2)+coefficient_Ca2[1]*x_Ca2**(-1)+coefficient_Ca2[2]+coefficient_Ca2[3]*x_Ca2+coefficient_Ca2[4]*x_Ca2**2+coefficient_Ca2[5]*x_Ca2**3+coefficient_Ca2[6]*x_Ca2**4)
    #         x_Ca3=np.arange(3000, Tmax, 0.01)
    #         coefficient_Ca3=Capitelli_fit[molecule][2]
    #         Cp_fit_Ca3=R*(coefficient_Ca3[0]*x_Ca3**(-2)+coefficient_Ca3[1]*x_Ca3**(-1)+coefficient_Ca3[2]+coefficient_Ca3[3]*x_Ca3+coefficient_Ca3[4]*x_Ca3**2+coefficient_Ca3[5]*x_Ca3**3+coefficient_Ca3[6]*x_Ca3**4)
    #         x_Ca=np.hstack((x_Ca1,x_Ca2,x_Ca3))
    #         Cp_fit_Ca=np.hstack(( Cp_fit_Ca1, Cp_fit_Ca2,Cp_fit_Ca3))
    #     print(molecule)
    #     print(Cp_fit_Ca)
    #     plt.plot(x_Ca,Cp_fit_Ca,color="aqua",label="Capitelli et.al")   

    # plot specific heat using nasa glenn polynomials
    if molecule in fit_nasa:
        if Tmax<=1000:
            x_nasa=np.arange(Tmin, Tmax, 0.01)
            coefficient_nasa=fit_nasa[molecule][0]
            Cp_fit_nasa=R*(coefficient_nasa[0]*x_nasa**(-2)+coefficient_nasa[1]*x_nasa**(-1)+coefficient_nasa[2]+coefficient_nasa[3]*x_nasa+coefficient_nasa[4]*x_nasa**2+coefficient_nasa[5]*x_nasa**3+coefficient_nasa[6]*x_nasa**4)
        else:
            x_nasa1=np.arange(Tmin, 1000, 0.01)
            coefficient_nasa1=fit_nasa[molecule][0]
            Cp_fit_nasa1=R*(coefficient_nasa1[0]*x_nasa1**(-2)+coefficient_nasa1[1]*x_nasa1**(-1)+coefficient_nasa1[2]+coefficient_nasa1[3]*x_nasa1+coefficient_nasa1[4]*x_nasa1**2+coefficient_nasa1[5]*x_nasa1**3+coefficient_nasa1[6]*x_nasa1**4)
            x_nasa2=np.arange(1000, Tmax, 0.01)
            coefficient_nasa2=fit_nasa[molecule][1]
            Cp_fit_nasa2=R*(coefficient_nasa2[0]*x_nasa2**(-2)+coefficient_nasa2[1]*x_nasa2**(-1)+coefficient_nasa2[2]+coefficient_nasa2[3]*x_nasa2+coefficient_nasa2[4]*x_nasa2**2+coefficient_nasa2[5]*x_nasa2**3+coefficient_nasa2[6]*x_nasa2**4)
            x_nasa=np.hstack((x_nasa1,x_nasa2))
            Cp_fit_nasa=np.hstack(( Cp_fit_nasa1, Cp_fit_nasa2))

        plt.plot(x_nasa,Cp_fit_nasa,color="black",label="NASA Glenn")   

    storepath="pictures"
    storename=storepath+"/"+molecule


    if molecule=='N2O':
        plt.legend(loc='lower left')
    else:
        plt.legend()
    plt.ylabel('$C_{p}$')
    plt.xlabel('T(K)')
    plt.title('Specific Heat Fit For '+molecule)
    plt.savefig(storename)
    plt.show()
# %%
# Draw the residuals versus fits plot
for molecule in tqdm(Cp_dict):
    Tmax=Tmax_dict[molecule]
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
        if molecule in fit_nasa:
            coefficient_nasa1=fit_nasa[molecule][0]
            Cp_fit_nasa1=R*(coefficient_nasa1[0]*x**(-2)+coefficient_nasa1[1]*x**(-1)+coefficient_nasa1[2]+coefficient_nasa1[3]*x+coefficient_nasa1[4]*x**2+coefficient_nasa1[5]*x**3+coefficient_nasa1[6]*x**4)
        x2=np.arange(1000.,Tmax+1., 1.)
        coefficient_this2=fit_dictionary[molecule][1]
        Cp_fit_this2=R*(coefficient_this2[0]*x2**(-2)+coefficient_this2[1]*x2**(-1)+coefficient_this2[2]+coefficient_this2[3]*x2+coefficient_this2[4]*x2**2+coefficient_this2[5]*x2**3+coefficient_this2[6]*x2**4)
        if molecule in fit_nasa:
            coefficient_nasa2=fit_nasa[molecule][1]
            Cp_fit_nasa2=R*(coefficient_nasa2[0]*x**(-2)+coefficient_nasa2[1]*x**(-1)+coefficient_nasa2[2]+coefficient_nasa2[3]*x+coefficient_nasa2[4]*x**2+coefficient_nasa2[5]*x**3+coefficient_nasa2[6]*x**4)
        x=np.hstack((x1,x2))
        Cp_fit_this=np.hstack(( Cp_fit_this1, Cp_fit_this2))

    resudual= Cp - Cp_fit_this

    plt.axhline(y=0.0,c="red")
    plt.scatter(x,resudual)    

    storepath="residual_plot"
    storename=storepath+"/"+molecule

    plt.legend()
    plt.ylabel('Residuals')
    plt.xlabel('T(K)')
    plt.title('Residuals Plot For '+molecule)
    plt.savefig(storename)
    plt.show()
# %%
