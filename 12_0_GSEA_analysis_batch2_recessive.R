# This code is going to do Gene set enrichment analysis

# Created 16-April-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clustering
rm(list = ls())

# Set the working directory to the path where files are located
# setwd("/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined")
setwd("P:/workspaces/lg-sit/working_data/MengYun/data/combined")

# Load libraries
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(forcats)
library(readxl)
library(openxlsx)


#########################################
# Step 1: Load gene list from a .txt file
#########################################
genes_in_solved <- readLines("genes_in_solved_cases_control_genes_excluded_recessive.txt")  
# Remove empty lines and trim whitespace
genes_in_solved <- trimws(genes_in_solved)
genes_in_solved <- genes_in_solved[genes_in_solved != ""] %>%
  unique(.)

genes_in_unsolved <- readLines("genes_in_unsolved_cases_control_genes_excluded_recessive.txt")
# Remove empty lines and trim whitespace
genes_in_unsolved <- trimws(genes_in_unsolved)
genes_in_unsolved <- genes_in_unsolved[genes_in_unsolved != ""] %>%
  unique(.)

genes_in_controls <- readLines("genes_in_control_case_genes_excluded_recessive.txt")
# Remove empty lines and trim whitespace
genes_in_controls <- trimws(genes_in_controls)
genes_in_controls <- genes_in_controls[genes_in_controls != ""] %>%
  unique(.)


####################################################
# Step 2: Convert gene identifiers (using g:Convert)
####################################################
solved_cases <- gconvert(genes_in_solved, organism = "hsapiens")$target
unsolved_cases <- gconvert(genes_in_unsolved, organism = "hsapiens")$target
controls <- gconvert(genes_in_controls, organism = "hsapiens")$target


#######################################
# Step3: do the main analysis
######################################
gsea_analysis <- function(gene_list_converted) {

  # Run enrichment analysis using g:Profiler
  gsea_result <- gost(
    query = gene_list_converted,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    correction_method = "gSCS",
    sources = c("GO:BP", "GO:MF", "GO:CC"),
    user_threshold = 0.01,
    domain_scope = "annotated",
    evcodes = TRUE
  )

  # Filter results by intersection size and gene set size
  filtered_results <- gsea_result$result %>%
    filter(intersection_size >= 2,
           term_size >= 15,
           term_size <= 500) 

  # Plot the enriched gene sets
  plot_data <- filtered_results %>%
    arrange(p_value) %>%
    mutate(term_name = fct_reorder(term_name,-log10(p_value)))
  
  gsea_plot <-
    ggplot(plot_data, aes(x = -log10(p_value), y = term_name)) +
    geom_point(aes(size = intersection_size, color = source)) +
    scale_size_continuous(name = "Intersecting Genes") +
    scale_color_brewer(palette = "Set1", name = "Source") +
    labs(title = " ",
         x = "-log10(Adjusted P-value)",
         y = "Gene Set") +
    theme_minimal() +
    theme(text = element_text(size = 12),
          legend.position = "right")

  # Interactive plot using gostplot
  gostplot(gsea_result, capped = TRUE, interactive = TRUE)

  # Convert the intersection column to a readable format
  gsea_clean <- filtered_results
  
  # convert Ensembl IDs into HGNC names
  ensembl_lists <- strsplit(gsea_clean$intersection, ",")
  all_ids <- unique(unlist(ensembl_lists))
  
  if (length(all_ids) > 0) {
  conversion_df <- gconvert(all_ids, organism = "hsapiens", target = "HGNC")
  id_to_symbol <- setNames(conversion_df$name, conversion_df$input)
  gsea_clean$intersection <- lapply(ensembl_lists, function(ids) id_to_symbol[ids])
  }
  
  # Identify list columns
  list_cols <- sapply(gsea_clean, is.list)
  
  # Convert lists to comma-separated character strings
  gsea_clean[list_cols] <- lapply(gsea_clean[list_cols], function(x) sapply(x, paste, collapse=", "))
  
  condition=deparse(substitute(gene_list_converted))
  
  if (condition=="solved_cases") {
    
    write.xlsx(gsea_clean, file = "Gene_set_enrichment_recessive.xlsx", sheetName = condition, rowNames=TRUE, colNames=TRUE)
    ggsave("Gene_set_enrichment_solved_cases_recessive.png", plot = gsea_plot, width = 8, height = 7, units = 'in', dpi = 300)
    
  } else {
    wb <- loadWorkbook("Gene_set_enrichment_recessive.xlsx")
    if (condition %in% names(wb)) {
      removeWorksheet(wb, condition)  # Remove the existing sheet
    }
    addWorksheet(wb, condition)
    writeData(wb, condition, gsea_clean, colNames = TRUE)
    saveWorkbook(wb, "Gene_set_enrichment_recessive.xlsx", overwrite=TRUE)
  }  
}

gsea_analysis(solved_cases)
gsea_analysis(unsolved_cases)
gsea_analysis(controls)

get_version_info(organism = "hsapiens")

