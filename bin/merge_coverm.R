#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript merge_tables.R <output.tsv> <input1.tsv> <input2.tsv> ...")
}

output_file <- args[1]
input_files <- args[-1]

message("Reading ", length(input_files), " files...")
dt_list <- lapply(input_files, fread)

row_names <- dt_list[[1]][[1]]

all_match <- sapply(dt_list, function(dt) identical(dt[[1]], row_names))

if (all(all_match)) {
    message("Rownames are identical, fast merging through `cbind`")
    data_cols <- lapply(dt_list, function(dt) dt[, -1, with = FALSE])
    merged_data <- do.call(cbind, data_cols)
    merged_dt <- cbind(data.table(Genome = row_names), merged_data)
} else {
    message("Rownames are not identical, safe merging through `merge`")
    merged_dt <- Reduce(function(x, y) merge(x, y, by = "Genome", all = TRUE), dt_list)
}

message("Write to output file: ", output_file)
fwrite(merged_dt, output_file, sep = "\t", quote = FALSE, na = "NA")