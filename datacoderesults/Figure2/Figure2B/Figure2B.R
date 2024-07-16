#load up the necessary packages and libraries and functions
installed.packages("barcodetrackR")
installed.packages("magrittr")
installed.packages("SummarizedExperiment")
library(dplyr)
library(barcodetrackR)
library (rlang)
library(ggplot2)
library(RColorBrewer)
require("magrittr")
require("barcodetrackR")
require("SummarizedExperiment")

setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/Figure2/Figure2B_Data_Code/Data")

#read main data
ZH33barcounts = read.delim("ZH33_combined_20231019_counts.txt", row.names = 1)
ZH19barcounts = read.delim("ZH19_combined_20240620_counts.txt", row.names = 1)
ZG66barcounts = read.delim("ZG66_combined_20230124_counts.txt", row.names = 1)
ZJ31barcounts = read.delim("ZJ31_combined_20220512_counts.txt", row.names = 1)
ZJ38barcounts = read.delim("ZJ38_combined_20240620_counts.txt", row.names = 1)

#read key file containing proper names and sampleIDs
ZH33barcountsmeta = read.delim("ZH33_combined_20231019_metadata_fig2B.txt")
ZH19barcountsmeta = read.delim("ZH19_combined_20240620_metadata_fig2B.txt")
ZG66barcountsmeta = read.delim("ZG66_combined_20230124_metadata_fig2B.txt")
ZJ31barcountsmeta = read.delim("ZJ31_combined_20220512_metadata_fig2B.txt")
ZJ38barcountsmeta = read.delim("ZJ38_combined_20240620_metadata_fig2B.txt")

#create the SEs
ZH33_SE <- create_SE(your_data = ZH33barcounts,
                   meta_data = ZH33barcountsmeta,
                   threshold = 0)

ZH19_SE <- create_SE(your_data = ZH19barcounts,
                     meta_data = ZH19barcountsmeta,
                     threshold = 0)

ZG66_SE <- create_SE(your_data = ZG66barcounts,
                     meta_data = ZG66barcountsmeta,
                     threshold = 0)

ZJ31_SE <- create_SE(your_data = ZJ31barcounts,
                     meta_data = ZJ31barcountsmeta,
                     threshold = 0)

ZJ38_SE <- create_SE(your_data = ZJ38barcounts,
                     meta_data = ZJ38barcountsmeta,
                     threshold = 0)

#create the function that will create the bar graph
clonal_count_unique_bar <- function(your_SE,
                         percent_threshold = 0,
                         plot_over,
                         plot_over_display_choices = NULL,
                         keep_numeric = TRUE,
                         group_by,
                         group_by_choices = NULL,
                         cumulative = FALSE,
                         point_size = 3,
                         line_size = 0.75,
                         text_size = 15,
                         your_title = NULL,
                         return_table = FALSE) {
  
  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  
  if (!(plot_over %in% coldata_names)) {
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if (!(group_by %in% coldata_names)) {
    stop("group_by must match a column name in colData(your_SE)")
  }
  
  if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% sort(unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
    plot_over_display_choices <- as.numeric(as.character(plot_over_display_choices))
  } else if (is.factor(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]], levels = levels(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  } else {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]], levels = unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  }
  
  group_by_choices <- group_by_choices %||% levels(as.factor(SummarizedExperiment::colData(your_SE)[[group_by]]))
  
  # More error handling
  if (!all(plot_over_display_choices %in% SummarizedExperiment::colData(your_SE)[, plot_over])) {
    stop("All elements of plot_over_display_choices must match cvalues in plot_over column")
  }
  if (!all(group_by_choices %in% SummarizedExperiment::colData(your_SE)[, group_by])) {
    stop("All elements of group_by_choices must match values in group_by column")
  }
  
  # extract metadata
  temp_subset <- your_SE[, (your_SE[[plot_over]] %in% plot_over_display_choices)]
  # Keep only the data included in group_by_choices
  temp_subset <- temp_subset[, (temp_subset[[group_by]] %in% group_by_choices)]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% tibble::as_tibble()
  your_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  your_data <- your_data[rowSums(your_data) > 0, , drop = FALSE]
  
  # calculate measure for each sample
  if (!cumulative) {
    calculated_index <- colSums(your_data > percent_threshold) %>%
      tibble::enframe(name = "SAMPLENAME", value = "index") %>%
      dplyr::mutate(index_type = "unique barcode count")
    ylabel <- "Unique Barcode Count"
  } else {
    
    # Get count data into tidy format (and remove zeros)
    tidy_counts <- your_data %>%
      dplyr::mutate(barcode_seq = rownames(your_data)) %>%
      dplyr::select(.data$barcode_seq, everything()) %>%
      tidyr::pivot_longer(cols = seq_len(ncol(your_data)) + 1, names_to = "sample", values_to = "count") %>%
      dplyr::filter(count > 0)
    
    # Add group by variable into tidy counts
    tidy_counts$group_by <- temp_subset_coldata[match(tidy_counts$sample, temp_subset_coldata$SAMPLENAME), ][[group_by]]
    
    # Add plot over variable into tidy counts
    tidy_counts$plot_over <- temp_subset_coldata[match(tidy_counts$sample, temp_subset_coldata$SAMPLENAME), ][[plot_over]]
    # tidy_counts$plot_over <- factor(tidy_counts$plot_over, levels = levels(plot_over_display_choices))
    
    # Order tidy counts by desired order
    if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
      tidy_counts_ordered <- tidy_counts %>%
        dplyr::arrange(group_by, plot_over)
    } else {
      tidy_counts_ordered <- tidy_counts %>%
        dplyr::arrange(group_by, factor(plot_over, levels = levels(plot_over_display_choices)))
    }
    
    # Only keep the first occurence of each barcode within each group_by category
    tidy_counts_filtered <- tidy_counts_ordered %>%
      dplyr::group_by(group_by) %>%
      dplyr::distinct(.data$barcode_seq, .keep_all = TRUE)
    
    # Summarize number of new barcodes and cumulative barcodes at each timepoint
    
    if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
      summarized_data <- tidy_counts_filtered %>%
        dplyr::group_by(.data$group_by, plot_over, .data$sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count))
    } else {
      summarized_data <- tidy_counts_filtered %>%
        dplyr::group_by(group_by, factor(plot_over, levels = levels(plot_over_display_choices)), sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count))
      colnames(summarized_data)[2] <- "plot_over"
      print(summarized_data)
    }
    
    # Put into proper structure
    calculated_index <- as.data.frame(summarized_data) %>%
      dplyr::select(.data$sample, .data$cumulative_count) %>%
      dplyr::rename(SAMPLENAME = .data$sample, index = .data$cumulative_count) %>%
      dplyr::mutate(index_type = "unique barcodes") %>%
      dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME))
    
    # Fix the fact that certain samples will be dropped if they have 0 new clones
    # Hacky way to do it. Need to go back and make it better later.
    if (length(unique(tidy_counts_filtered$sample)) < length(unique(tidy_counts_ordered$sample))) {
      samples_not_found <- unique(tidy_counts_ordered$sample)[unique(tidy_counts_ordered$sample) %in% unique(tidy_counts_filtered$sample) == FALSE]
      calculated_index <- rbind(calculated_index, data.frame(
        SAMPLENAME = samples_not_found,
        index = rep(0, length(samples_not_found)),
        index_type = rep("unique barcodes", length(samples_not_found))
      ))
      # Reorder
      calculated_index <- calculated_index %>%
        dplyr::arrange(factor(.data$SAMPLENAME), levels = unique(tidy_counts_ordered$sample))
      
      # Give the missing samples the same cumulative count as above samples
      for (i in seq_along(nrow(calculated_index))) {
        if (calculated_index[i, "SAMPLENAME"] %in% samples_not_found) {
          calculated_index[i, "index"] <- calculated_index[i - 1, "index"]
        }
      }
    }
    ylabel <- "cumulative barcode count"
  }
  
  # merge measures with colData
  plotting_data <- temp_subset_coldata %>%
    dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME)) %>%
    dplyr::left_join(calculated_index, by = "SAMPLENAME")
  
  # Make sure plot over is a factor if not numeric or specified to not keep numeric.
  if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric) {
  } else if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric == FALSE) {
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = unique(plot_over_display_choices))
  } else {
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = levels(plot_over_display_choices))
  }
  
  if (return_table) {
    return(calculated_index)
  }
  
  plotting_data$x_value <- plotting_data[[plot_over]]
  plotting_data$group_by <- plotting_data[[group_by]]

  custom_colors <- c("B_real" = rgb(207,93,109, maxColorValue=255), "Gr_real" = rgb(71,150,224, maxColorValue = 255), "T_real" = "black")
  
  # Create ggplot
  g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$x_value, y = .data$index, group = .data$group_by, colour = .data$group_by)) +
    ggplot2::geom_line(size = line_size) +
    ggplot2::geom_point(size = point_size, shape = 15) +
    ggplot2::labs(x = "Months Post-Transplant", col = group_by, y = ylabel) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   axis.text.x = element_text(size = 20, color = "black"), 
                   axis.text.y = element_text(size = 20, color = "black", angle = 90, vjust = 0.5, hjust = 0.5), 
                   axis.title.x = element_text(size = 20), 
                   axis.title.y = element_text(size = 20),
                   panel.border = element_rect(color = "black", fill = NA, size = 1),
                   axis.line = element_blank(),
                   axis.ticks.length=unit(.25, "cm"),
                   plot.title = element_text(size = 24, hjust = 0.5)) +
    ggplot2::ggtitle(your_title) +
    #ggplot2::guides(colour = FALSE) + 
    ggplot2::scale_color_manual(values = custom_colors)
  
  if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric) {
    g + ggplot2::scale_x_continuous(paste0("Months Post-Transplant"), breaks = (scales::breaks_pretty()), limits = c(0,140))
  } else {
    g + ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  }
  
  # g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$x_value, y = .data$index, fill = .data$group_by)) +
  #   ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.8) +
  #   ggplot2::labs(x = plot_over, y = ylabel, fill = group_by) +
  #   ggplot2::theme_classic() +
  #   ggplot2::theme(text = ggplot2::element_text(size = text_size)) +
  #   ggplot2::ggtitle(your_title)
  # 
  # if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric) {
  #   g + ggplot2::scale_x_continuous(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  # } else {
  #   g + ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  # }
  
  
  
}

#run the function to produce the graphs
clonal_count_unique_bar(your_SE = ZH33_SE, your_title = "ZH33",cumulative = FALSE, keep_numeric = TRUE, plot_over = "Months",group_by = "Cell.Type", group_by_choices = c("T_real","B_real","Gr_real"), text_size = 10)
clonal_count_unique_bar(your_SE = ZH19_SE, your_title = "ZH19",cumulative = FALSE, keep_numeric = TRUE, plot_over = "Months",group_by = "Cell.Type", group_by_choices = c("T_real","B_real","Gr_real"), text_size = 10)
clonal_count_unique_bar(your_SE = ZG66_SE, your_title = "ZG66",cumulative = FALSE, keep_numeric = TRUE, plot_over = "Months",group_by = "Cell.Type", group_by_choices = c("T_real","B_real","Gr_real"), text_size = 10)
clonal_count_unique_bar(your_SE = ZJ31_SE, your_title = "ZJ31",cumulative = FALSE, keep_numeric = TRUE, plot_over = "Months",group_by = "Cell.Type", group_by_choices = c("T_real","B_real","Gr_real"), text_size = 10)
clonal_count_unique_bar(your_SE = ZJ38_SE, your_title = "ZJ38",cumulative = FALSE, keep_numeric = TRUE, plot_over = "Months",group_by = "Cell.Type", group_by_choices = c("T_real","B_real","Gr_real"), text_size = 10)

