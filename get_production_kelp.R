library(future)
library(future.apply)
#Need single_kelp_run
#
# Code to run a kelp DEB model in each grid cell in norwecom at a specified depth. 
#forcings (datetime, Temp, irradiance, CO_2), needs to be for the same depth as specified, function does not change irradiance etc. 

# Author : Malene Dekke Eik 
##

get_production_kelp <- function(file, forcings, grid, depth, initial_state) {
  
  # Setting up multiple sessions 
  options(future.globals.maxSize = 10 * 1024^3) # memory 
  #options(future.globals.maxSize = 4 * 1024^3)
  plan(multisession, workers = 4)
  
  nc <- nc_open(file)
  
  totPP = ncvar_get(nc, varid="Nit")
  C = totPP
  nx <- dim(C)[1]
  ny <- dim(C)[2]
  
  # Extract time values
  nc_time_vals <- nc$dim$T$vals
  time_units <- nc$dim$T$units  # 'hours since 1950-01-01 01:00:00'
  origin_time <- as.POSIXct("1950-01-01 01:00:00", tz = "UTC")
  nc_datetime <- origin_time + as.difftime(nc_time_vals, units = "hours")
  
  # Changing the months October to December to the year 2019 , to fill nitrogen values in forcings (may change based on data set-up)  
  nc_datetime <- as.POSIXct(nc_datetime)  
  month_vals <- as.integer(format(nc_datetime, "%m"))
  year_vals <- as.integer(format(nc_datetime, "%Y"))
  year_vals[year_vals == 2020 & month_vals %in% 10:12] <- 2019
  
  # Replace with new datetime
  nc_datetime <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                                    year_vals,
                                    as.integer(format(nc_datetime, "%m")),
                                    as.integer(format(nc_datetime, "%d")),
                                    as.integer(format(nc_datetime, "%H")),
                                    as.integer(format(nc_datetime, "%M")),
                                    as.integer(format(nc_datetime, "%S"))),
                            tz = "UTC")
  
  # Depth 
  h <- ncvar_get(nc, varid='Topo') #depth at each x,y cordinate 
  zEdges <- nc$dim$Z$vals
  S.new <- array(NA, dim=c(dim(C)[1], dim(C)[2], dim(C)[3]))
  S.new[,,]=rep(h, dim(C)[3])
  S.new = apply(S.new, 1:2, function(x) x*zEdges)
  zm <- aperm(S.new, c(2,3,1))
  S.new <- zm-depth # subtracting the the depth from the calculated depth of each layer at each x,y 
  
  S.min <- apply(S.new, 1:2, function(x) which.min(abs(x)))
  
  #FOR TESTING smal grid
  #nx= 4
  #ny = 4
  
  # Dataframe with the time and date of forcings, and location of grid cells (EMPTY)
  nit_df <- expand.grid(
    x = 1:nx,
    y = 1:ny,
    datetime = forcings$datetime
  ) %>%
    mutate(nit_mol = NA_real_)
  
  
  # Function to process one grid column (i) to run parallell sessions 
  process_column <- function(i) {
    res_list <- vector("list", ny)
    
    for (j in 1:ny) {
      k <- S.min[i, j]  # correct depth layer
      orig_vals <- totPP[i, j, k, ]
      
      if (all(is.na(orig_vals))) {
        res_list[[j]] <- NULL
        next
      }
      
      # Interpolate values for each forcing datetime
      matched_vals <- approx(x = as.numeric(nc_datetime),
                             y = orig_vals,
                             xout = as.numeric(forcings$datetime),
                             method = "linear",
                             rule = 2)$y
      
      df_part <- data.frame(
        x = i,
        y = j,
        datetime = forcings$datetime,
        nit_mol = matched_vals * 1e-6  # B5mol to mol (mol / L )
      )
      
      res_list[[j]] <- df_part
    }
    
    do.call(rbind, res_list)
  }
  
  # Run parallel across columns
  nit_df_list <- future_lapply(1:nx, process_column)
  
  # Combiningg all grid cell data
  nit_df <- do.call(rbind, nit_df_list)
  
  #Write csv file if we want to save nit_df
  'library(here)
  write.csv(nit_df, here::here("data", "nit_df.csv"), row.names = FALSE)'
  
  
  
  # Empty dataframes to store last weight of singel kelp (DW g)
  single_kelp <- matrix(NA, nrow = nx, ncol = ny) 
  
  #Deployment and harvest time 
  #deployment_date_time = as.POSIXct("2019-11-01 00:00:00")
  deployment_date_time <- min(forcings$datetime)
  harvest_date_time <- max(forcings$datetime)
  
  ###
  #split nit_df in order to bring it into the function, without exceeded memory limit .. 
  nit_df_split <- split(nit_df, nit_df$x)
  
  process_single_kelp_column <- function(i, df_i, forcings_local) {
    sapply(1:ny, function(j) {
      cell_forcings <- forcings_local %>%
        mutate(N = df_i$nit_mol[df_i$y == j])
      
      if (all(is.na(cell_forcings$N))) return(NA)
      
      out <- tryCatch({
        single_kelp_run(
          forcings = cell_forcings,
          deployment_date_time = deployment_date_time,
          harvest_date_time = harvest_date_time,
          params = params,
          initial_state = initial_state
        )
      }, error = function(e) NULL)
      
      if (!is.null(out)) {
        final_dw <- tail(out$W, 1)
        return(final_dw)
      } else {
        return(NA)
      }
    })
  }
  
  # Run the first half 
  ix1 <- 1:floor(nx / 2)  #Indices for first half 
  nit_df_list_subset <- nit_df_split[as.character(ix1)]
  single_kelp_list_1 <- future_mapply( 
    FUN = process_single_kelp_column,
    i = ix1,
    df_i = nit_df_list_subset,
    MoreArgs = list(forcings_local = forcings),
    SIMPLIFY = FALSE
  )
  
  
  # Second half 
  ix2 <- (floor(nx / 2) + 1):nx
  nit_df_list_subset_2 <- nit_df_split[as.character(ix2)]
  single_kelp_list_2 <- future_mapply( 
    FUN = process_single_kelp_column,
    i = ix2,
    df_i = nit_df_list_subset_2,
    MoreArgs = list(forcings_local = forcings),
    SIMPLIFY = FALSE
  )
  
  # Fill the single_kelp list 
  for (idx in seq_along(ix1)) {
    i <- ix1[idx]
    single_kelp[i, ] <- single_kelp_list_1[[idx]]
  }
  for (idx in seq_along(ix2)) {
    i <- ix2[idx]
    single_kelp[i, ] <- single_kelp_list_2[[idx]]
  }
  
  
  
  lat <- ncvar_get(nc, varid="Latt")
  lon <- ncvar_get(nc, varid="Long")
  
  nc_close(nc)
  
  #TODO, comment out when not testing 
  #lat = lat[1:nx, 1:ny]
  #lon = lon[1:nx, 1:ny]
  
  #Add production per m2 
  nkelp = 100 #How many kelp pr square meter 
  production <- single_kelp * nkelp
  
  # 
  #prod <- list(single_kelp=single_kelp, production=production, lat=lat, lon=lon, topo=h)
  
  return(list(production=production, 
              single_kelp=single_kelp,
              topo=h, 
              lat=lat, 
              lon=lon))
  
  }