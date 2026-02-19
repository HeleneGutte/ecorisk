test_that("plotting functions give a known output", {
  plot_radar_test <- plot_radar(
    risk_scores = ex_output_risk_expert,
    aggregated_scores = ex_output_aggregate_risk_expert
    )
  vdiffr::expect_doppelganger("radar plot", plot_radar_test,
                              cran = FALSE)  
  
  plot_heatmap_test <- plot_heatmap(
    risk_scores = ex_output_risk_expert,
    aggregated_scores = ex_output_aggregate_risk_expert
  )
  
  vdiffr::expect_doppelganger("heatmap plot", plot_heatmap_test,
                              cran = FALSE) 
    
})
