# load data
setwd(here::here())
sd_sim_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_sim_data.csv") |> 
  dplyr::filter(!(!constant_Tau & Tau == 0.5)) ## exclude setting with Tau = 0.5, non-constant rate (not discussed)

# count the number of replicates where assay data had to be re-simulated 
## note: this is the same for all methods 
sd_sim_data |> 
  dplyr::filter(Method == "MLE_woUDSA") |> 
  dplyr::pull(Message_Detailed) |> 
  table()
