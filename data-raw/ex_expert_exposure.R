## Code to prepare `exposure_ex` dataset goes here
ex_expert_exposure <- create_template_exposure(
  pressures = c("temperature", "salinity", "oxygen", "nutrient", "fishing"),
  which_components = 1:4,
  mode_uncertainty = "component"
)

for(i in 2:5){
  set.seed(99+i)
  ex_expert_exposure[, i] <- sample(c(1:5), 5, replace = TRUE)
}
for(i in 6:9){
  set.seed(99+i)
  ex_expert_exposure[, i] <- sample(c(1:3), 5, replace = TRUE)
}

usethis::use_data(ex_expert_exposure, overwrite = TRUE)
