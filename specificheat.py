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
list=[['CH4'],['H2O2'],['CLO'],['HCN','CO','HF','HI','N2','NH3','OCS','OH','PH3','C2H2','C2N2','C2H4','CH3CL','CH3F','OH','O2'],['SF6'],['CS','HO2','NO','NO+','SO'],['H2','HBr','HCL'],['SO3'],['C2H6'],['C4H2']]
split=[[200,500,1300,1500],[200,1500],[200,4000],[200,1000],[200,1000,2000,3000,4000],[200,1000,4000],[200,1000,5000],[200,3500,200,500,650],[200,1000,2000,3000],[200,1000,2000]]
split_dict=dict()
for i in range(len(list)):
    for j in range(len(list[i])):
        if list[i][j] not in split_dict:
            split_dict[list[i][j]]=split[i]
# %%
# Calculate specific heat
Cp_dict=dict()
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

filename='Cp_thiswork.json'
with open(filename,'w') as file_obj:
        json.dump(Cp_dict,file_obj)
# %%
