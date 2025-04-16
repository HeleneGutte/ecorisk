## Code to prepare `status_ex` dataset goes here
ex_expert_status <- data.frame(
  indicator = c("phytoplankton", "herring", "cod", "seabirds"),
  status = c("good", "undesired", "undesired", "good"),
  score = c(1, -1, -1, 1)
)

usethis::use_data(ex_expert_status, overwrite = TRUE)
