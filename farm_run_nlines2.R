#trying to make more efficient code for farm_run_nlines.R 

farm_run_nlines <- function(pars, farms_pars, forcings){
  time_steps <- dim(forcings)[1]
  
  #creating matrices to store the results and adding initial values
  output_final <- list(
    ingested_carbon = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),   #how much CO2 is ingested by the kelp box
    ingested_nitrogen = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)), #how much nitrogen is ingested by the kelp box
    poc_carbon = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),        #how much concentration of CO2 is left after absorption by kelp box
    poc_nitrogen = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),      #how much concentration of nitrogen is left after absorption by kelp box
    W = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),                 #DW of an individual kelp in the box 
    K_length = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),          #cm
    K_area = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),            #cm2
    BM = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),                #DW of biomass of all kelp in the box
    structure = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),         #M_V mass of structure (mol M_V)
    reserves_m_EC = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),     #m_EC (mol C/mol M_V) reserve density of C reserve
    reserves_m_EN = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),     #m_EN (mol N/mol M_V) reserve density of N reserv
    nkelp = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),             #ind kelp in box
    current = array(0, dim=c(farms_pars$nbox, farms_pars$nlines, time_steps)),           #saving the current 
    datetime = array(NA, dim = c(farms_pars$nbox, farms_pars$nlines, time_steps))        #datetime for each timestep for plotting 
  ) #output list [n, m, t] : box n, line m, time t 
  
  
  temp_nlines = 3 
  #Creating a temporary matice to store 3 lines in
  output <- list(
    ingested_carbon = array(0, dim=c(farms_pars$nbox, 3, time_steps)),   #how much CO2 is ingested by the kelp box
    ingested_nitrogen = array(0, dim=c(farms_pars$nbox, 3, time_steps)), #how much nitrogen is ingested by the kelp box
    poc_carbon = array(0, dim=c(farms_pars$nbox, 3, time_steps)),        #how much concentration of CO2 is left after absorption by kelp box
    poc_nitrogen = array(0, dim=c(farms_pars$nbox, 3, time_steps)),      #how much concentration of nitrogen is left after absorption by kelp box
    W = array(0, dim=c(farms_pars$nbox, 3, time_steps)),                 #DW of an individual kelp in the box 
    K_length = array(0, dim=c(farms_pars$nbox, 3, time_steps)),          #cm
    K_area = array(0, dim=c(farms_pars$nbox, 3, time_steps)),            #cm2
    BM = array(0, dim=c(farms_pars$nbox, 3, time_steps)),                #DW of biomass of all kelp in the box
    structure = array(0, dim=c(farms_pars$nbox, 3, time_steps)),         #M_V mass of structure (mol M_V)
    reserves_m_EC = array(0, dim=c(farms_pars$nbox, 3, time_steps)),     #m_EC (mol C/mol M_V) reserve density of C reserve
    reserves_m_EN = array(0, dim=c(farms_pars$nbox, 3, time_steps)),     #m_EN (mol N/mol M_V) reserve density of N reserv
    nkelp = array(0, dim=c(farms_pars$nbox, 3, time_steps)),             #ind kelp in box
    current = array(0, dim=c(farms_pars$nbox, 3, time_steps)),           #saving the current 
    datetime = array(NA, dim = c(farms_pars$nbox, 3, time_steps))        #datetime for each timestep for plotting 
  ) #output list [n, m, t] : box n, line m, time t 
  
  
  ### Initial forcing pr timestep ### 
  #Check that temperatures are in kelvin
  if (mean(forcings$temperature<50)){                   # check if temperature in celsius
    forcings$temperature = forcings$temperature+273.15  # convert to kelvin
  }
  #Setting up the forcing for each timestep and longline, will change after each box run 
  forcing_t <- vector("list", temp_nlines)             # A list for each longline
  for (m in 1:temp_nlines) {
    forcing_t[[m]] <- list(
      T_field = forcings$temperature,   # Kelvin
      I_field = forcings$I,             # PAR
      N_field = forcings$N,             # mol N /L 
      CO2_field = forcings$CO_2,        # mol CO2 / L
      speed = forcings$speed )          # m/h
  }
  
  
  ### Box parameters ###
  boxheight = farms_pars$height                   # m , height of box, can change as kelp grows bigger
  boxlength = farms_pars$longline/farms_pars$nbox # m, length of box along longline  
  box_s_area = farms_pars$width * boxlength       # m2 , area of box surface 
  volumebox = box_s_area * boxheight * 1000       # L
  box_o_area = boxheight * farms_pars$width       # Area of the box opening 
  
  ### Farm design ###
  #Adjusting the number of kelp in the box based off the farm design 
  if(farms_pars$farm_design == 1) {            # Longline system 
    NKELP = farms_pars$density * boxlength
  } 
  if(farms_pars$farm_design == 2) {            # Vertical line system 
    NKELP = farms_pars$density * boxheight
  } 
  if(farms_pars$farm_design == 3) {            # Horizontal 5-line system 
    NKELP = farms_pars$density * boxlength * 5
  }
  if(farms_pars$farm_design == 4) {            # Net cultivation system 
    NKELP = farms_pars$density * (boxheight * boxlength)
  }
  
  ### Initial state ###
  for (m in 1:temp_nlines){
    for (n in 1:farms_pars$nbox){
      #Output for the initial time step (initial state/conditions)
      output$W[n,m,1]=farms_pars$Winit                                      # g DW 
      output$K_length[n,m,1] = (farms_pars$Winit/0.00387)^(1/1.469)         # cm
      output$K_area[n,m,1] = 39.3790 * (farms_pars$Winit/0.1342)^0.912530   # cm2
      output$structure[n,m,1]=farms_pars$M_V_init                           # molM_V
      output$reserves_m_EC[n,m,1]=farms_pars$m_EC_init                      # mol C/molM_V
      output$reserves_m_EN[n,m,1]=farms_pars$m_EN_init                      # mol N/molM_V
      output$nkelp[n,m,1]= NKELP                                            # ind / box 
      output$poc_carbon[n,m,1] = forcing_t[[m]]$CO2_field[1]
      output$poc_nitrogen[n,m,1] = forcing_t[[m]]$N_field[1]
      output$current[n,m,1] = forcing_t[[m]]$speed[1] 
      output$datetime[n,m,1] = forcings$datetime[1]
      output$ingested_carbon[n,m,1] = 0
      output$ingested_nitrogen[n,m,1] = 0
    }
  }
  
  delta = 2*farms_pars$fk/farms_pars$width # variation of flow speed along the longline
  spatial_resolution = farms_pars$longline/farms_pars$nbox
  
  
  ###
  
  ### Looping through each timestep ### 
  for (t in 2:time_steps){ #starting at timestep 2 , timestep 1 is initial state 
    
    #print(paste("Starting time step:", t))
    
    ### FLIP if current has changed direction 
    #Flip the matix when the flow direction changes 
    flip = 0 #This is to flip the matrix when the flow direction changes
    #current = forcings$speed[t]
    if (sign(forcings$speed[t]) != sign(forcings$speed[t-1])) { #only flips the matrix if the current is a different direction then previous timestep 
      flip=1
      for (m in 1:temp_nlines) {
        for (i in seq_along(output)) {
          output[[i]] = output[[i]][farms_pars$nbox:1, , ]
        } } }
    
    
    ### Looping through lines ###
    for (m in 1:temp_nlines){
      #print(paste("Lines, Time step:", t, "Line:", m))
      V1 = forcing_t[[m]]$speed[t] #saving the current speed entering the farm (first box)
      

        ### Looping through boxes ###
        for (n in 1:farms_pars$nbox){
          
          #TODO: current and friction effect on dissolved nutrients in water for each new box 
          # calculate current entering the box
          if (n == 1) { 
            output$current[n, m, t] = forcing_t[[m]]$speed[t]  
          } else {
            if (output$nkelp[n,m,t]>0){
              if (delta==0){ #no friction
                output$current[n,m,t] = output$current[n-1, m, t]
              }else{ #friction
                output$current[n,m,t]=output$current[n-1, m, t]*exp(-delta * spatial_resolution)
              }
              if (abs(output$current[n,m,t])<farms_pars$min_current){
                output$current[n,m,t]=output$current[n-1, m, t];
              } } }
          
          
          '#N of kelp in second timestep same as inital state 
          if(t == 2){
            output$nkelp[n,m,2] = output$nkelp[n,m,1]}'
          
          #Density 
          output$nkelp[n, m, t] = output$nkelp[n,m,t-1] #making the number of kelp in this timestep same as previous, here we can add mortality etc
          
          ### State of kelp in the box (from previous timestep)
          state <- c(
            m_EC = output$reserves_m_EC[n, m, t-1], 
            m_EN = output$reserves_m_EN[n, m, t-1] , 
            M_V = output$structure[n, m, t-1], 
            W = output$W[n, m, t-1] )
          
          
          #Solving for 1 kelp in the box at timestep t (DEB - model)
          sol_t <- rates_Lo(t=t, state=state, parameters=params_Lo, forcing = forcing_t[[m]])
          
          
          ### updating the output 
          output$W[n,m,t] = sol_t$W                                                        #DW of individual kelp in box 
          output$K_length[n,m,t] = sol_t$L_allometric                                      #length of individual kelp in box (cm)
          output$K_area[n,m,t] = sol_t$Area_allometric                                     #area f individual kelp in box (cm2)
          output$structure[n,m,t] = output$structure[n, m, t-1] + sol_t$M_V                #M_V mass of structure (mol M_V), individual
          output$reserves_m_EC[n,m,t] = output$reserves_m_EC[n, m, t-1] +sol_t$m_EC        #m_EC (mol C/mol M_V) reserve density of C reserve, individual 
          output$reserves_m_EN[n,m,t] =  output$reserves_m_EN[n, m, t-1] + sol_t$m_EN      #m_EN (mol N/mol M_V) reserve density of N reserve, individual 
          
          output$BM[n,m,t] = output$W[n, m, t] * output$nkelp[n,m,t]                       #Biomass of kelp in box in DW 
          output$ingested_carbon[n,m,t] = sol_t$J_CO2 * sol_t$M_V * output$nkelp[n,m,t]    #how much CO2 is ingested by the kelp box: mol c 
          output$ingested_nitrogen[n,m,t] = sol_t$J_EN_A * sol_t$M_V * output$nkelp[n,m,t] #how much nitrogen is ingested by the kelp box: mol N 
          
          #print(paste("Output updated, Time step:", t, "Line:", m))
          
          
          #current update , equation from Rosland et al., 2011
          Vn <- forcing_t[[m]]$speed[t] #current speed entrance of box 
          forcing_t[[m]]$speed[t] <- V1 * ( (1-((farms_pars$fk*boxlength)/farms_pars$width)) / (1+((farms_pars$fk*boxlength)/farms_pars$width)) )^n #current speed exit of box (entering the next box)
          
          
          #updating the output, what concentration of nitrogen and carbon does the kelp have available  
          output$poc_carbon[n,m,t] = forcing_t[[m]]$CO2_field[t]
          output$poc_nitrogen[n,m,t] = forcing_t[[m]]$N_field[t]
          
          #Concentration of carbon and nitrogen after box run (with current speed )
          #the nutrient concentration is same as initial forcing at timestep if the next box is a perimiter box 
          if( m == 1 || m ==temp_nlines) { #its a perimeter line
            forcing_t[[m]]$CO2_field[t] <- forcings$CO_2[t]
            #forcing_t[[m]]$CO2_field[t] <- forcings$CO_2[t] * ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) / ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) + output$ingested_carbon[n,m,t]))
            forcing_t[[m]]$N_field[t] <- forcings$N[t]
            #forcing_t[[m]]$N_field[t] <- forcings$N[t] * ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) / ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) + output$ingested_nitrogen[n,m,t]))
          } else{ #depletion of nutrients for other boxes 
            forcing_t[[m]]$CO2_field[t] <- forcing_t[[m]]$CO2_field[t] * ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) / 
                                                                            ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) + output$ingested_carbon[n,m,t]))
            forcing_t[[m]]$N_field[t] <- forcing_t[[m]]$N_field[t] * ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) / 
                                                                        ((box_o_area*(abs(Vn) + abs(forcing_t[[m]]$speed[t]))) + output$ingested_nitrogen[n,m,t]))
            
          }
          
          
          #Adding date and time for plotting 
          output$datetime[n,m,t] <- forcings$datetime[t]
          
          
          #Changing the height of the box if relevant 
          if(farms_pars$depth_change == TRUE){
            #change the depth and volume of box as the kelp grows and gets longer 
            if(!is.na(output$K_length[n,m,t]) && (output$K_length[n,m,t]/100) > boxheight) {
              boxheight = (output$K_length[n,m,t]/100) # + some depth for sinking of the rope between the buoys ? 
              volumebox = box_s_area * boxheight * 1000 #L
              box_o_area = boxheight * farms_pars$width #Area of the box opening
            }
          }
          
          
          #For last timestep Boxes back in original position
          if(t == time_steps) {
            if (sign(forcings$speed[t]) != sign(forcings$speed[1])) { #only flips the matrix if the current is a different direction then first timestep
              for (m in 1:temp_nlines) {
                for (i in seq_along(output)) {
                  output[[i]] = output[[i]][farms_pars$nbox:1, , ]
                } } }
          }
        }
    } 
    #print(paste("Time step:", t, "Out of:", time_steps))
  }
  
  #Fill in output_final 
  
  for (m in 1:farms_pars$nlines) {
    
    if (m == 1){
      output_final$ingested_carbon[ ,m, ] <- output$ingested_carbon[ ,1, ]
      output_final$ingested_nitrogen[ ,m, ] <- output$ingested_nitrogen[ ,1, ]
      output_final$poc_carbon[ ,m, ] <- output$poc_carbon[ ,1, ]
      output_final$poc_nitrogen[ ,m, ] <- output$poc_nitrogen[ ,1, ]
      output_final$W[ ,m, ] <- output$W[ ,1, ]
      output_final$K_length[ ,m, ] <- output$K_length[ ,1, ]
      output_final$K_area[ ,m, ] <- output$K_area[ ,1, ]
      output_final$BM[ ,m, ] <- output$BM[ ,1, ]
      output_final$structure[ ,m, ] <- output$structure[ ,1, ]
      output_final$reserves_m_EC[ ,m, ] <- output$reserves_m_EC[ ,1, ]
      output_final$reserves_m_EN[ ,m, ] <- output$reserves_m_EN[ ,1, ]
      output_final$nkelp[ ,m, ] <- output$nkelp[ ,1, ]
      output_final$current[ ,m, ] <- output$current[ ,1, ]
      output_final$datetime[ ,m, ] <- output$datetime[ ,1, ]
    }
    if(m > 1 & m < farms_pars$nlines) {
      output_final$ingested_carbon[ ,m, ] <- output$ingested_carbon[ ,2, ]
      output_final$ingested_nitrogen[ ,m, ] <- output$ingested_nitrogen[ ,2, ]
      output_final$poc_carbon[ ,m, ] <- output$poc_carbon[ ,2, ]
      output_final$poc_nitrogen[ ,m, ] <- output$poc_nitrogen[ ,2, ]
      output_final$W[ ,m, ] <- output$W[ ,2, ]
      output_final$K_length[ ,m, ] <- output$K_length[ ,2, ]
      output_final$K_area[ ,m, ] <- output$K_area[ ,2, ]
      output_final$BM[ ,m, ] <- output$BM[ ,2, ]
      output_final$structure[ ,m, ] <- output$structure[ ,2, ]
      output_final$reserves_m_EC[ ,m, ] <- output$reserves_m_EC[ ,2, ]
      output_final$reserves_m_EN[ ,m, ] <- output$reserves_m_EN[ ,2, ]
      output_final$nkelp[ ,m, ] <- output$nkelp[ ,2, ]
      output_final$current[ ,m, ] <- output$current[ ,2, ]
      output_final$datetime[ ,m, ] <- output$datetime[ ,2, ] 
    } 
    if (m == -1){
      output_final$ingested_carbon[ ,m, ] <- output$ingested_carbon[ ,3, ]
      output_final$ingested_nitrogen[ ,m, ] <- output$ingested_nitrogen[ ,3, ]
      output_final$poc_carbon[ ,m, ] <- output$poc_carbon[ ,3, ]
      output_final$poc_nitrogen[ ,m, ] <- output$poc_nitrogen[ ,3, ]
      output_final$W[ ,m, ] <- output$W[ ,3, ]
      output_final$K_length[ ,m, ] <- output$K_length[ ,3, ]
      output_final$K_area[ ,m, ] <- output$K_area[ ,3, ]
      output_final$BM[ ,m, ] <- output$BM[ ,3, ]
      output_final$structure[ ,m, ] <- output$structure[ ,3, ]
      output_final$reserves_m_EC[ ,m, ] <- output$reserves_m_EC[ ,3, ]
      output_final$reserves_m_EN[ ,m, ] <- output$reserves_m_EN[ ,3, ]
      output_final$nkelp[ ,m, ] <- output$nkelp[ ,3, ]
      output_final$current[ ,m, ] <- output$current[ ,3, ]
      output_final$datetime[ ,m, ] <- output$datetime[ ,3, ]
    }
    
  }
  
  #
  
  # 
   return(output_final)
}