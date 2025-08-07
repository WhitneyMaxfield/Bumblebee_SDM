library(dplyr)

# Fix potential conflicts with select function
select <- dplyr::select

# Get original and processed variables
original_vars <- names(env_df)
selected_vars <- vif_result@results$Variables
removed_vars <- setdiff(original_vars, selected_vars)

# Get final VIF values for retained variables
final_vif <- vif_result@results
names(final_vif)[names(final_vif) == "Variables"] <- "Variable"

# Create comprehensive results table
vif_summary <- data.frame(
  Variable = original_vars,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Status = case_when(
      Variable %in% selected_vars ~ "RETAINED",
      Variable %in% removed_vars ~ "REMOVED",
      TRUE ~ "UNKNOWN"
    )
  ) %>%
  left_join(final_vif, by = "Variable") %>%
  mutate(
    Final_VIF = round(VIF, 3),
    VIF_Category = case_when(
      is.na(VIF) ~ "Removed (High VIF)",
      VIF < 2.5 ~ "Low collinearity",
      VIF >= 2.5 & VIF < 5 ~ "Low-moderate collinearity", 
      VIF >= 5 & VIF < 10 ~ "Moderate collinearity",
      VIF >= 10 ~ "High collinearity",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Variable, Status, Final_VIF, VIF_Category) %>%
  arrange(Status, Final_VIF)

# Print the table
print(vif_summary)

# Count summary
cat("Total variables analyzed:", nrow(vif_summary), "\n")
cat("Variables retained:", sum(vif_summary$Status == "RETAINED"), "\n")
cat("Variables removed:", sum(vif_summary$Status == "REMOVED"), "\n")

# Save to CSV
write.csv(vif_summary, "vif_analysis_results.csv", row.names = FALSE)
cat("\nTable saved as 'vif_analysis_results.csv'\n")

# Optional: Create a more detailed table with removal order
if(length(vif_result@excluded) > 0) {
  removal_order <- data.frame(
    Removal_Order = 1:length(vif_result@excluded),
    Variable = vif_result@excluded,
    Reason = "High VIF (removed sequentially)"
  )
  
  print(removal_order)
  
  write.csv(removal_order, "vif_removal_order.csv", row.names = FALSE)
  cat("Removal order saved as 'vif_removal_order.csv'\n")
}


# Create a publication-ready table
pub_table <- vif_summary %>%
  mutate(
    Variable_Clean = case_when(
      grepl("wc2.1_2.5m_bio_1", Variable) ~ "Annual Mean Temperature",
      grepl("wc2.1_2.5m_bio_2", Variable) ~ "Mean Diurnal Range", 
      grepl("wc2.1_2.5m_bio_3", Variable) ~ "Isothermality",
      grepl("wc2.1_2.5m_bio_8", Variable) ~ "Mean Temp Wettest Quarter",
      grepl("wc2.1_2.5m_bio_9", Variable) ~ "Mean Temp Driest Quarter",
      grepl("wc2.1_2.5m_bio_13", Variable) ~ "Precipitation Wettest Month",
      grepl("wc2.1_2.5m_bio_14", Variable) ~ "Precipitation Driest Month", 
      grepl("wc2.1_2.5m_bio_15", Variable) ~ "Precipitation Seasonality",
      grepl("NLCD", Variable) ~ "Land Cover Class",
      TRUE ~ Variable
    )
  ) %>%
  select(Variable_Clean, Status, Final_VIF, VIF_Category)

print(pub_table)

write.csv(pub_table, "vif_results_publication.csv", row.names = FALSE)
cat("Publication-ready table saved as 'vif_results_publication.csv'\n")
