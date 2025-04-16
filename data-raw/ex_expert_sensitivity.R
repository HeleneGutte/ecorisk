# Code to prepare `sens_ac_ex` dataset goes here

# ex_expert_sensitivity <- create_template_sensitivity(
#   indicators = c("phytoplankton", "herring", "cod", "seabirds"),
#   pressures = c("temperature", "salinity", "oxygen", "nutrient", "fishing"),
#   type = c("direct", "direct_indirect"),
#   n_sensitivity_traits = 5,
#   mode_adaptive_capacity = "trait",
#   mode_uncertainty = "trait"
# )
#
# names(ex_expert_sensitivity) <- c("indicator", "pressure", "type",
#   "sens_feeding", "sens_behaviour", "sens_reproduction", "sens_habitat", "sens_general",
#   "ac_feeding", "ac_behaviour", "ac_reproduction", "ac_habitat", "ac_general",
#   "uncertainty_sens_feeding", "uncertainty_sens_behaviour", "uncertainty_sens_reproduction",
#   "uncertainty_sens_habitat", "uncertainty_sens_general",
#   "uncertainty_ac_feeding", "uncertainty_ac_behaviour", "uncertainty_ac_reproduction",
#   "uncertainty_ac_habitat", "uncertainty_ac_general")
# write.csv(ex_expert_sensitivity, file = "ex_expert_sensitivity.csv")

# Re-import file after filling out scores
ex_expert_sensitivity <- read.csv("data-raw/ex_expert_sensitivity.csv", sep = ";")

usethis::use_data(ex_expert_sensitivity, overwrite = TRUE)
