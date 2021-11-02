#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)
library(bit64)
library(ggplot2)
library(scales)

###########
# GLOBALS #
###########

hist_file <- snakemake@input[["hist"]]
plot_file <- snakemake@output[["plot"]]

# dev
# hist_file <- "output/030_merqury/guppy237.merqury/illumina.hist"

########
# MAIN #
########


# read data
hist_data <- fread(hist_file)

# plot
my_pal <- viridisLite::viridis(9)
kmer_plot <- ggplot(hist_data, aes(x = `V1`, y = V2)) +
    theme_minimal() +
    theme(legend.position = c(5/6, 2/4)) +
    geom_vline(xintercept = c(5, 10),
               linetype = 2,
               colour = my_pal[[4]]) +
    geom_vline(xintercept = 32,
               linetype = 3,
               colour = my_pal[[8]]) +
    geom_path(alpha = 0.75, colour = my_pal[[1]]) +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("21-mer depth") + ylab("Number of unique 21-mers")

# write output
ggsave(filename = plot_file,
       plot = kmer_plot,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
