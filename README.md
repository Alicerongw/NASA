# NASA polynomials

NASA polynomials are convenient empirical equations used to represent the thermodynamic properties of species. They are extensively used in various applications, including plasma modelling. However, recent work working on these polynomials are lacking and the NASA Glenn thermodynamic database hasn’t been updated for 20 years, meaning there are no polynomials available for many species. The aim of this project is to generate up-to-date NASA polynomials for a wider range of molecules. 

# Data Source

The HITRAN database (https://hitran.org/) has recently provided an extensive compilations of partition functions which can be used to compute specific heats. Values of the total internal partition functions can be directly accessed by downloading the text files in the 8th column of the Isotopologue Metadata table from HITRONonline (https://hitran.org/docs/iso-meta/) shown in the following picture. A registered account is needed for online access. 

![HITRAN](/Picture_used_in_readme/HITRAN_example.png) 
In these files, the total internal partition functions Q are given in intervals of 1K from 1 to a given maximum temperature T_max, which is generally 5000K. 

# Structure
The main strusture of the whole repository is shown in the figure.
![HITRAN](/Picture_used_in_readme/Structure.png) 

# Output
The most important NASA polynomials are stored in the Fit_coefficients.csv, where the coefficients for each intervals of every species are listed. Besides, the computed specific heats are stored in the folder Cp_results. Each species has its own csv files where the specific heats are reported in the main text and in 1 K increments from 200K to Tmax.

Besides, the plots comparing the present specific heat with existing results are stored in Cp_results. The plots about residuals for curve-fitting are available in the folder Fit_residuals_pictures while the mean Residuals,max Residuals and standard deviation of Residuals are saved in Fit_residuals.csv.

# Run Instructions
All the python codes were available in the file nasa_polynomial_HITRAN.py. All the output can be produced by running this file.

This program can be further used to compute NASA polynomials from other partition functions. For instance, it is still suitable if HITRAN database updated. The further use is straightforward that just download the new partitions function from HITRAN, rename it by the formula of species, store it in the folder Partition_functions and rerun it.