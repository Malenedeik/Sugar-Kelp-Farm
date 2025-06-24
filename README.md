This reposotory contains code associated with the master thesis of Malene Dekke Eik (June, 2025) : "Modelling production potential of sugar kelp in multi-use offshore area: optimal location and farm design". \
Any questions about the code or data can be directed towards Malene Dekke Eik (malene.deik@gmail.com) and Antonio Aguera Garcia (antonio.aguera@hi.no). 

The code for the sugar kelp DEB model is from Venolia et. al (2020) https://doi.org/10.1016/j.ecolmodel.2020.109151 ; the original code can be found at : https://github.com/CVenolia/SugarKelpDEB \
Data from the Austevoll farm was collected for the master thesis of Skaar (2019) https://bora.uib.no/bora-xmlui/handle/1956/21026

**Code file keys:** \
DEB_params.R: Function code for parameters for the DEB model \
KelpDEB_model_Venolia.R: Function code for the Kelp DEB model by Venolia using the R package deSolve, modified to add early growth parameter \
SolveR_R.R: Function code, The Newton_raphson solver for the net specific growth rate (r) in the Kelp DEB model, by Venolia \
Austevoll_Compare.qmd: The run file to compare the modelled growth of sugar kelp in the Austevoll farm to the measured growth \
area_W.qmd: The runfile to calculate the leaf area to weight allometric relationship and FW to DW ratio, using data from Austevoll \
 \
farm_model_kelp.qmd: The runfile for the kelp farm model, contains the plot code for the kelp farm growth \
farm_run_R: Function code for the single central kelp longline model \
farm_run_nlines2.R: Function code for the kelp farm model with multiple farm designs \
utsira.env.qmd: The runfile to get the production potential and nutrient availability in Utsira Nord from a NORWECOM.E2E nc file \
get_production_kelp.R: Function code for running a single kelp in each grid cell of a nc file \
single_kelp_run.R: Function code to model the growth of a single kelp using the DEB model \
 \
Forcing_data_plots.qmd: The runfile to get the forcings from Utsira Nord and contains the plots \
Plots_maps_current_wave-qmd: Runfile to plot maps, currents, and waves at Utsira Nord \
NOAA_KdPAR.qmd: The runfile to get the average kd PAR value from NOAA data, and plot  \
 
**Data file keys:** \
data/forcings_austevoll.csv: Forcing data for Austevoll farm location (hourly) \
data/Mastersheet_october2019_area-weight-fouling.xlsx: Measured growth data from Austevoll farm \
data/utsira_hourly_vars.txt: Forcing data for central point in Utsira Nord at surface for 2020 from NORWECOM.E2E (hourly) \
data/utsira_env_best_location_2020.txt: Forcing data for farm location at surface for 2020 from NORWECOM.E2E (hourly) \
data/forcings_farmloc_3m.csv: Forcing data for Utsira Nord farm location at 3 m depth (hourly) \
data/forcings_farmloc_6m.csv: Forcing data for Utsira Nord farm location at 6 m depth (hourly) \
data/forcings_farmloc_9m.csv: Forcing data for Utsira Nord farm location at 9 m depth (hourly) \
data/cmems_mod_glo_wav_anfc_0.083deg_PT3H-i_1740646698248.csv: Data from Copernicus for Sea surface primary swell wave significant height at central location of Utsira Nord (4.531135 59.274227) \
data/cmems_mod_glo_wav_anfc_0.083deg_PT3H-i_1740646731856.csv: Data from Copernicus for Sea surface wave maximum height at central location of Utsira Nord (4.531135 59.274227) \
data/cmems_mod_glo_wav_anfc_0.083deg_PT3H-i_1740646749295.csv: Data from Copernicus for Sea surface wave significant height at central location of Utsira Nord (4.531135 59.274227) \
data/noaacwNPPN20VIIRSkdparDaily_5fe0_5bba_1516_U1745571424213.nc: Data from NOAA for kd Par \
