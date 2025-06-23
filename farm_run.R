farm_run <- function(pars, farms_pars, forcings){
  time_steps <- dim(forcings)[1]
  #creating matrices to store the results and adding initial values
  output <- list(
    ingested_carbon = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),    #how much CO2 is ingested by the kelp box 
    ingested_nitrogen = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),  #how much nitrogen is ingested by the kelp box 
    poc_carbon = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),         #how much concentration of CO2 is left after absorption by kelp box 
    poc_nitrogen = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),       #how much concentration of nitrogen is left after absorption by kelp box
    current = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps), 
    W = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),                  #DW of an individual kelp in the box 
    K_length = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),           #cm
    K_area = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),             #cm2 
    BM = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),                 #DW of biomass of all kelp in the box
    structure = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),          #M_V mass of structure (mol M_V)
    reserves_m_EC = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),      #m_EC (mol C/mol M_V) reserve density of C reserve
    reserves_m_EN = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps),      #m_EN (mol N/mol M_V) reserve density of N reserve 
    maturity = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps), #delete ? 
    nkelp = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps)               #ind kelp in box 
    #boxvolume = matrix(data=0, nrow=farms_pars$nbox, ncol=time_steps) # L , changes as kelp grows gets longer
  )
  
  #Setting up the forcing for each timestep, will change after each box run 
  forcing_t <- list(
    T_field = forcings$temperature, #Kelvin
    I_field = forcings$I, #PAR
    N_field = forcings$N, #mol N /L 
    CO2_field = forcings$CO_2,
    speed = forcings$speed) #mol CO2 / L 
  
  #Setting up box parameters 
  boxheight = farms_pars$height #m , height of box, will change as kelp grows bigger
  boxlength = farms_pars$longline/farms_pars$nbox # m, length of box along longline  
  box_s_area = farms_pars$width * boxlength #m2 , area of box surface 
  volumebox = box_s_area * boxheight * 1000 #L
  box_o_area = boxheight * farms_pars$width #Area of the box opening 
  
  #Output for the initial time step (initial state/conditions)
  output$W[,1]=farms_pars$Winit #g DW 
  output$K_length[,1] = (farms_pars$Winit/0.00387)^(1/1.469)
  output$K_area[,1] = 39.3790 * (farms_pars$Winit/0.1342)^0.912530
  output$structure[,1]=farms_pars$M_V_init #molM_V
  output$reserves_m_EC[,1]=farms_pars$m_EC_init #mol C/molM_V
  output$reserves_m_EN[,1]=farms_pars$m_EN_init #mol N/molM_V
  #output$repro_buff[,1]=farms_pars$ERinit
  #output$maturity[,1]=farms_pars$Hinit
  output$nkelp[,1]= farms_pars$density * boxlength
  #output$boxvolume[,1] = volumebox
  output$poc_carbon[,1] = forcing_t$CO2_field[1]
  output$poc_nitrogen[,1] = forcing_t$N_field[1]
  
  
  
  delta = 2*farms_pars$fk/farms_pars$width # variation of flow speed along the longline
  spatial_resolution = farms_pars$longline/farms_pars$nbox
  
  #Check that temperatures are in kelvin
  if (mean(forcings$temperature<50)){ # temperature in celsius
    forcings$temperature = forcings$temperature+273.15  #convert to kelvin
    }

  #flip = 0 # for the first timestep, as a reference to the second time step 
  
  for (t in 2:time_steps){ #starting at timestep 2 , timestep 1 is initial state 
    
    #Flip the matix when the flow direction changes 
    flip = 0 #This is to flip the matrix when the flow direction changes
    current = forcings$speed[t]
    #if (forcings$speed[t]<0){
    if (sign(forcings$speed[t]) != sign(forcings$speed[t-1])) { #only flips the matrix if the current is a different direction then previous timestep 
      flip=1
      for(i in 1:length(output)){
        output[[i]]=output[[i]][farms_pars$nbox:1,]
      } }
    
    V1 = forcing_t$speed[t] #saving the current speed entering the farm (first box)

    
    ### Defining other boxes ###
    for (n in 1:farms_pars$nbox){
      
      #For the first box check that there is minimum current 
      if(n== 1){
        if (abs(forcings$speed[t])<farms_pars$min_current){ #maintain a minimum current at the boundary
          if (flip==1){
            forcings$speed[t] = -1 * farms_pars$min_current
          }else{
            forcings$speed[t] = farms_pars$min_current
          }
        }
        output$current[1,t]=forcings$speed[t]
      }
      else{
        #TODO: current and friction effect on dissolved nutrients in water for each new box 
        # calculate current entering the box
        if (output$nkelp[n,t]>0){
          if (delta==0){ #no friction
            output$current[n,t]=output$current[n-1,t]
          }else{ #friction
            output$current[n,t]=output$current[n-1,t]*exp(-delta * spatial_resolution)
          }
          if (abs(output$current[n,t])<farms_pars$min_current){
            output$current[n,t]=output$current[n-1,t];
          }
        }else{
          output$current[n,t]=output$current[n-1,t]
        }
      }
      
      #N of kelp in second timestep same as inital state 
      if(t == 2){
        output$nkelp[n,2] = output$nkelp[n,1]}
      
      #individual pumping rate based on DEB
      #DEB FOR KELP , loop of time step. for each time step run model 
      
      #State of kelp in the box
      if(t == 1){
        state <- c(m_EC = output$reserves_m_EC[n,t], m_EN = output$reserves_m_EN[n,t] , M_V = output$structure[n,t], W = output$W[n,t]) #for the first timestep
      }else{
        state <- c(m_EC = output$reserves_m_EC[n,t-1], m_EN = output$reserves_m_EN[n,t-1] , M_V = output$structure[n,t-1], W = output$W[n,t-1] ) #getting the current state of the kelp from previous 
      }
      
      #Solving for 1 kelp in the box at timestep t 
      sol_t <- rates_Lo(t=t, state=state, parameters=params_Lo, forcing = forcing_t)
      
      #updating the output 
      output$W[n,t] = sol_t$W #DW of individual kelp in box 
      output$K_length[n,t] = sol_t$L_allometric #length of individual kelp in box (cm)
      output$K_area[n,t] = sol_t$Area_allometric #area f individual kelp in box (cm2)
      output$structure[n,t] = output$structure[n,t-1] + sol_t$M_V  #M_V mass of structure (mol M_V), individual
      output$reserves_m_EC[n,t] = output$reserves_m_EC[n,t-1] +sol_t$m_EC #m_EC (mol C/mol M_V) reserve density of C reserve, individual 
      output$reserves_m_EN[n,t] =  output$reserves_m_EN[n,t-1] + sol_t$m_EN #m_EN (mol N/mol M_V) reserve density of N reserve, individual 
      
      output$BM[n,t] = output$W[n,t] * output$nkelp[n,t] #Biomass of kelp in box in DW 
      output$ingested_carbon[n,t] = sol_t$J_CO2 * sol_t$M_V * output$nkelp[n,t]#how much CO2 is ingested by the kelp box: mol c 
      output$ingested_nitrogen[n,t] = sol_t$J_EN_A * sol_t$M_V * output$nkelp[n,t] #how much nitrogen is ingested by the kelp box: mol N 
      #output$poc_carbon[n,t] =  forcing_t$CO2_field[t] - (output$ingested_carbon[n,t] /volumebox ) #how much concentration of CO2 is left after absorption by kelp box 
      #output$poc_nitrogen[n,t] = forcing_t$N_field[t] - (output$ingested_nitrogen[n,t] /volumebox)#how much concentration of nitrogen is left after absorption by kelp box
       
      #density 
      if(t < time_steps) {
        output$nkelp[n, t+1] = output$nkelp[n, t]} #making the number of kelp in the next timestep, here we can add mortality etc. 
      
      
      #updating the environmental data after each box run , without current speed 
      #forcing_t$CO2_field[t] <- output$poc_carbon[n,t] #CO2 left after this box run
      #forcing_t$N_field[t] <- output$poc_nitrogen[n,t] #Nitrite Nitrate after this box run
      
      #current update , equation from Rosland et al., 2011
      Vn <- forcing_t$speed[t] #current speed entrance of box 
      forcing_t$speed[t] <- V1 * ( (1-((farms_pars$fk*boxlength)/farms_pars$width)) / (1+((farms_pars$fk*boxlength)/farms_pars$width)) )^n #current speed exit of box (entering the next box)
      
      #Concentration of carbon and nitrogen after box run (with current speed )
      forcing_t$CO2_field[t] <- forcing_t$CO2_field[t] * ((box_o_area * (abs(Vn) + abs(forcing_t$speed[t])))  / 
                                                           ((box_o_area*(abs(Vn) + abs(forcing_t$speed[t]))) +output$ingested_carbon[n,t]))
      forcing_t$N_field[t] <- forcing_t$N_field[t] * ((box_o_area * (abs(Vn) + abs(forcing_t$speed[t]))) / 
                                                       ((box_o_area*(abs(Vn) + abs(forcing_t$speed[t]))) + output$ingested_nitrogen[n,t] ))
      
      #updating the output , after calculaction of nitrogen and carbon concentraion 
      output$poc_carbon[n,t] = forcing_t$CO2_field[t]
      output$poc_nitrogen[n,t] = forcing_t$N_field[t]
      
      #change the depth and volume of box as the kelp grows and gets longer 
      if(!is.na(output$K_length[n, t]) && (output$K_length[n,t]/100) > boxheight) {
        boxheight = (output$K_length[n,t]/100) # + some depth for sinking of the rope between the buoys ? 
        volumebox = box_s_area * boxheight * 1000 #L
        box_o_area = boxheight * farms_pars$width #Area of the box opening
      } else {}
    }
  }
  
  
  
  return(output)
}
