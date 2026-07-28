#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(showtext)



# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line("grey90"),
  text = element_text(family="Helvetica", face="plain"),
  axis.title = element_text(colour="black"),
  title = element_text(colour = "#595A5C", face = "plain")
)


# Define functions


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  XGBoost_ROC <- fread("code/CompareEqtlsRelevance/XGBoost_ROC.txt")
  RF_ROC <- fread("code/CompareEqtlsRelevance/RandomForest_ROC.txt")

  ggplot(XGBoost_ROC, aes(FPR, TPR)) + geom_line(color="blue") +
    geom_line(data=RF_ROC, aes(FPR, TPR), color="salmon") +
    geom_abline(linetype="dashed") +
    guides(
      x = guide_axis(cap = "both"), # Cap both ends
      y = guide_axis(cap = "both") # Cap the upper end
    ) + theme(legend.position = "none", axis.text=element_text(size=7)) +
    coord_fixed()

  ggsave("aucplot.pdf", width=74, height=74, units="mm")


  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}