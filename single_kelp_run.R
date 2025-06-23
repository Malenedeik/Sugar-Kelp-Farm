single_kelp_run <- function(forcings, deployment_date_time, harvest_date_time, params, initial_state) {
  
  time_steps <- as.numeric(difftime(harvest_date_time, deployment_date_time, units = "hours")) -1
  
  output <- list(
    W = numeric(time_steps),                  # DW of an individual kelp in the box
    K_length = numeric(time_steps),           # cm
    K_area = numeric(time_steps),             # cm^2
    structure = numeric(time_steps),          # M_V mass of structure (mol M_V)
    reserves_m_EC = numeric(time_steps),      # m_EC (mol C/mol M_V) reserve density of C reserve
    reserves_m_EN = numeric(time_steps),      # m_EN (mol N/mol M_V) reserve density of N reserve
    datetime = as.POSIXct(time_steps),        # saving the date and time
    info = numeric(time_steps),           #For checking/testing parameters and other outputs 
    acc_factor = numeric(time_steps)      #Checking the acceleration factor 
  )
  
  #Filtering out only forcings for specific time run 
  forcings <- forcings %>% 
    filter(datetime >= deployment_date_time & datetime <= harvest_date_time)
  
  #Setting up forcings 
  forcing_t <- list(
    T_field = forcings$temperature, #Kelvin
    I_field = forcings$I,           #PAR
    N_field = forcings$N,           #mol N /L 
    CO2_field = forcings$CO_2)      #mol CO2 / L
  
  #First row of output /initial conditions 
  output$W[1]=initial_state$W                   #g DW
  output$structure[1]=initial_state$M_V         #molM_V
  output$reserves_m_EC[1]=initial_state$m_EC    #mol C/molM_V
  output$reserves_m_EN[1]=initial_state$m_EN    #mol N/molM_V
  output$K_length[1] = (initial_state$W/0.00387)^(1/1.469)
  output$K_area[1] = 39.3790 * (initial_state$W/0.1342)^0.912530
  output$datetime[1] = forcings$datetime[1]
  #output$info[1] = params$JENAM
  
  
  for (t in 2:time_steps) {
  
    state <- c(m_EC = output$reserves_m_EC[t-1], m_EN = output$reserves_m_EN[t-1] , M_V = output$structure[t-1], W = output$W[t-1] ) #getting the current state of the kelp from previous 
    
    sol_t <- rates_Lo(t=t, state=state, parameters=params_Lo, forcing = forcing_t)
    
    #updating the output 
    output$W[t] = sol_t$W                                     #DW of individual kelp in box 
    output$K_length[t] = sol_t$L_allometric                   #length of individual kelp in box (cm)
    output$K_area[t] = sol_t$Area_allometric                  #area f individual kelp in box (cm2)
    output$structure[t] = output$structure[t-1] + sol_t$M_V  #M_V mass of structure (mol M_V), individual
    output$reserves_m_EC[t] = output$reserves_m_EC[t-1] +sol_t$m_EC #m_EC (mol C/mol M_V) reserve density of C reserve, individual 
    output$reserves_m_EN[t] =  output$reserves_m_EN[t-1] + sol_t$m_EN #m_EN (mol N/mol M_V) reserve density of N reserve, individual
    output$datetime[t] = forcings$datetime[t]
    output$info[t] = sol_t$info #checking the parameters
    output$acc_factor[t] = sol_t$acc_factor #checking the parameters
  }
  return(output)
}