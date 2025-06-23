
load_pars <- function(){
  ## This is a function to compile the parameters, generated the needed ones and
  ## store them in a file that can easily be called for diverse applications
  
  # "food N" "food C" Stucture "N reserves" "C reserves" products
  #Conversion coefficients, organics (n = matrix of chemical indices)
  #        X_N  X_C     V    E_N    E_C   P  
  n_O <- matrix(
    + c(  0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
        + 0.00, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
        + 3.00, 2.50, 1.00, 2.50, 1.00, 0.50,  #O/C
        + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C, (N/N for E_N)
  
  #V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
  #E_N is N03- and N02- averaged
  #E_C is glucose C6H12O6 (Laminarin: c18h32o16 and mannitol c6h14o6)
  #We aren't using the X_N, X_C, or P collumn here
  
  #Molecular weights
  #t() is a matrix transpose function
  #organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
  w_O_step <- t(n_O)*matrix(c(12.011, 1.008, 15.999, 14.007), nrow=6, ncol=4, byrow= TRUE) #g/mol, molecular weights for organics
  w_O <- rowSums(w_O_step) #this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C)
  
  
  params <- data.frame(estimate=0, name='string', value=0, units='string', description='string')
  # name = notation of the parameter
  # value
  # units
  # description
  # estimate : boolean to point out which parameters if any to calibrate
  
  ## list parameters and attributes
  params[1,]=c(1, 'JENAM', 1.5e-4, 'molN/molV/h', 'maximum volume-specific assimilation rate of N before temperature correction') #or value 1.4e-4 #3e-4
  params[2,]=c(0, 'K_N', 2.5e-6, 'molNO3-NO2/L', 'half saturation constant of N uptake') #2.5e-6
  params[3,]=c(0, 'JCO2M', 0.0075, 'molDIC/molV/h', 'max volume-specific carbon dioxide assimilation rate')
  params[4,]=c(0, 'K_C', 4e-7, 'molDIC/L', 'half saturation constant of C uptake')
  params[5,]=c(1, 'JECAM', 0.282 , 'molC/molV/h', 'maximum volume-specific carbon assimilation rate') #or value 0.282
  params[6,]=c(0, 'rho_PSU', 0.5, 'molPSU/molV', 'Photosynthetic unit density')
  params[7,]=c(0, 'b_I', 0.5, '-', 'binding probability of photons to a Light SU')
  params[8,]=c(0, 'alpha', 1, 'm2/molPSU', 'Specific photon arrival cross section')
  params[9,]=c(0, 'k_I', 0.075, 'molY/molPSU/h', 'dissociation rate')
  params[10,]=c(0, 'y_I_C', 10, 'molY/molC', 'Yield factor of C reserve to photon')
  params[11,]=c(0, 'y_CO2_C', 1, 'molDIC/molC', 'Yield factor of C reserve to DIC')
  params[12,]=c(0, 'y_LO2', 0.125, 'molO2/molY', 'Yield factor of photon to O2')
  params[13,]=c(0, 'kE_C', 0.02, '1/h', 'C reserve turnover')
  params[14,]=c(0, 'kE_N', 0.04, '1/h', 'N reserve turnover')
  params[15,]=c(0, 'kappa_Ei', 0.9, '-', 'fraction of rejection flux from growth SU incorporated back into i-reserve')
  params[16,]=c(0, 'y_EN_V', 0.04, 'molN/molV', 'yield of structure on N reserve (percent of N in structure')
  params[17,]=c(0, 'y_EC_V', 1, 'molC/molV', 'yield of structure on C reserve (percent of C in structure')
  params[18,]=c(0, 'JENM', 4e-6, 'molN/molV/h', 'specific maintenance costs requiring N before temp correction') #4e-66
  params[19,]=c(0, 'JECM', 1e-6, 'molC/molV/h', 'specific maintenance costs requiring C before temp correction') #1e-6
  params[20,]=c(1, 'T_A', 3400, 'K', 'Arrhenius temperature') #Venolia: 6314.3 , Krasnow: 6403. 
  params[21,]=c(0, 'T_H', 11.65+273.15, 'K', 'Higher boundary of temperature tolerance') #Venolia: 13.39, Krasnow: 11.65 
  params[22,]=c(0, 'T_L', 273.15, 'K', 'Lower boundary of temperature tolerance') 
  params[23,]=c(0, 'T_AH', 16328, 'K', 'Arrhenius temperature outside T_H') #Venolia: 18702, Krasnow: 16328
  params[24,]=c(0, 'T_AL', 4391.9, 'K', 'Arrhenius temperature outside T_L')
  params[25,]=c(0, 'T_0', 20+273.15, 'K', 'reference temperature')
  params[26,]=c(0, 'w_V', w_O[3], 'g/mol', 'molecular weight of structure')
  params[27,]=c(0, 'w_EN', w_O[4], 'g/mol', 'molecular weight of N reserve')
  params[28,]=c(0, 'w_EC', w_O[5], 'g/mol', 'molecular weight of C reserve')
  params[29,]=c(0, 'w_O2', 32, 'g/mol', 'molecular weight of dioxygen')
  params[30,]=c(0, 'a_dw', 0.7, 'g DW', 'gram DW until growth is accelerated') #0.05 change 
  params[31,]=c(0, 'a_f', 25, '-', 'accelerator factor for early growth') #5? change 
  params[32,]=c(0, 'DM', 0.1342, '-', 'Dry matter %') 
  
  params$estimate <- as.numeric(params$estimate)
  params$value <- as.numeric(params$value)
  
  # create named vectors for easy use
  params_Lo <- setNames(params$value,params$name)
  cal <- setNames(params$estimate,params$name)
  
  out <- list(p=params_Lo, params=params, n_O = n_O, w_O = w_O)
  
  return(out)
  }


