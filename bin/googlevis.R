#!/usr/bin/env Rscript

library(googleVis)


args = commandArgs(trailingOnly = TRUE)
report_file_for_sankey = args[1]
sankey_height = args[2]
sankey_width = args[3]
plot_file_for_sankey = args[4]


# get sankey plot
my_data = read.csv(report_file_for_sankey, header = TRUE)

if (is.null(sankey_height)) {
  sankey_height = length(row.names(refinement_df)) * 12
}

if (is.null(sankey_width)) {
  sankey_width = 3 * 200
}

if (sankey_width < 400) {sankey_width = 400}

sankey_color_setting = "{node: {colorMode: 'unique', labelPadding: 10 }, link:{colorMode: 'source'}}"

Sankey_plot_my_data <- gvisSankey(
    my_data,
    options = list(sankey = sankey_color_setting,
                   height = sankey_height, width = sankey_width)
    )

sink(file = plot_file_for_sankey, append = FALSE,
     type = c("output", "message"), split = FALSE)
print(Sankey_plot_my_data)
sink()