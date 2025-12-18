#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2021
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(data.table)
library(tidyverse)
library("rworldmap")
library("rworldxtra")
library("sf")
library(ggrepel)
library(extrafont)
library(readxl)
library(giscoR)


# Declare constants
MAPS <- tibble(map_name = c("i", "ii"), lon_min = c(-12, -128), lon_max = c(33, -60), lat_min = c(36, 28), lat_max = c(64, 54))
pop_colors <- c("African (Tanzania)" = "red",
                "African (African American)" = "tomato",
                "European (Ashkenazi)" = "darkturquoise",
                "European (Finnish)" = "steelblue4",
                "European (Estonian)" = "steelblue1",
                "European (West)" = "lightblue",
                "European (West and South)" = "lightblue3",
                "East Asian (Japanese)" = "forestgreen",
                "African (Moroccan)" = "deeppink3",
                "Middle-Eastern (Qatari)" = "pink4",
                "European (Toscani)" = "lightblue4",
                "European (Spanish)" = "deepskyblue3",
                "South Asian (Bangladeshi)" = "orange",
                "South Asian (Indian)" = "tan2",
                "East Asian (Malaysia)" = "darkolivegreen2",
                "East Asian (Chinese)" = "springgreen3")

# Declare constants
loadfonts()
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "plain"),
)


# Declare function definitions
plot_to_pdf <- function(plot, path, width = 6, height = 6) {
  pdf(path,
      width = width, height = height, useDingbats = F)
  par(xpd = NA)

  print(plot)

  dev.off()
  embed_fonts(path, outfile=path)
}

# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  #place_data <- fread("~/Desktop/map_for_poster.csv")
  #out <- geocode(place_data$Place)
  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  master_table <- fread(variables["master_table"]) %>% select(c("cohort_new_name", "N" = "dataset_n"))
  cohort_info <- fread("cohort_info.txt")
  genetic_axes <- fread("GeneticAxes_2024-05-21.txt.gz")
  lat_lon <- read_xlsx("mastertable_with_lat_long.xlsx", sheet="mastertable_empirical_2024-05-1")

  data_qc_metrics <- fread("cohort_qc_info.tsv")

  combined <- inner_join(cohort_info, master_table, by=c("cohort"="cohort_new_name")) %>%
    inner_join(lat_lon, by=c("cohort"="cohort_new_name")) %>%
    inner_join(genetic_axes, by=c("cohort"="Cohort"))

  fwrite(combined, "cohort_supplementary_table.txt", sep=";", col.names=T, row.names=F)

  world_st <- getMap(resolution = "low")
  world_sf <- st_as_sf(world_st)

  world1 <- giscoR::gisco_get_countries() %>%
    dplyr::filter(NAME_ENGL != "Antarctica")

  turbine_sf <- combined %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326)

  turbine_sf_t <- st_transform(turbine_sf, crs=2163)

  world <- ggplot(data = turbine_sf_t, label = cohort) +
    geom_sf(data = world1, inherit.aes = F, colour = "white", fill = "#e6e6e6") +
    geom_sf(aes(size = dataset_n, fill=ancestry),
            shape = 21, colour="white", show.legend = T, linewidth=(0.25 / (ggplot2::.pt * 72.27/96))) +
    scale_size_continuous(breaks = c(2000, 4000, 6000),
                          limits = c(0, 8000),
                          range = c(0, 8), name="Sample size") +
    scale_fill_manual(values=pop_colors, na.value="grey50") +
    guides(size= guide_legend(position ="bottom"),
           fill="none") +
    coord_sf(crs = st_crs("EPSG:3857")) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0),
                         "inches"),
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text = element_blank())

  world

  ggsave(file.path(".", sprintf(
    "out.world.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
         world, width=80, height=80,units="mm")

  mds <- ggplot(data = combined,
                aes(x = -MDS1, y = -MDS2, fill=ancestry)) +
    geom_point(aes(size = dataset_n, fill=ancestry),
               shape = 21, colour="white", linewidth=(0.25 / (ggplot2::.pt * 72.27/96))) +
    scale_fill_manual(values=pop_colors) +
    scale_x_continuous(expand = c(0.2,0.2)) +
    scale_y_continuous(expand = c(0.2,0.2)) +
    scale_size_continuous(breaks = c(2000, 4000, 6000), limits = c(0, 8000), range = c(0, 8), name="Sample size") +
    guides(size="none",
           fill="none") +
    coord_equal() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())

  mds

  ggsave(file.path(".", sprintf(
    "out.mds.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
         mds, width=42, height=42,units="mm")

  combined_ordered <- combined %>%
    mutate(ancestry_major = str_split_i(ancestry, " \\(", 1),
           ancestry_minor = str_extract(ancestry, "\\((.+)\\)", group=1)) %>%
    group_by(ancestry, ancestry_minor, ancestry_major) %>%
    summarise(N = sum(N), N_cohorts = n()) %>%
    ungroup() %>%
    arrange(N) %>% mutate(ancestry_minor = ordered(ancestry_minor, levels=rev(unique(ancestry_minor)))) %>%
    group_by(ancestry_major) %>%
    mutate(N_ancestry = sum(N)) %>%
    ungroup() %>%
    arrange(N_ancestry) %>% mutate(ancestry_major = ordered(ancestry_major, levels=rev(unique(ancestry_major))))

  n_cohorts <- ggplot(data = combined_ordered,
                aes(y = ancestry_minor, x=1, label=sprintf("N = %s, datasets = %s", N, N_cohorts), fill=ancestry)) +
    geom_text(colour="black") +
    scale_fill_manual(values=pop_colors) +
    guides(size= guide_legend(position = "bottom"),
           fill="none") +
    scale_x_discrete(name="Number of datasets") +
    theme() +
    facet_grid(rows = vars(ancestry_major), scales="free_y", space="free_y")

  n_cohorts

  sample_sizes <- ggplot(data = combined_ordered,
                         aes(y= ancestry_minor, x=N, fill=ancestry)) +
    geom_col(colour="white") +
    scale_fill_manual(values=pop_colors) +
    guides(size= guide_legend(position = "bottom"),
           fill="none") +
    scale_x_log10(name="Sample size per ancestry", breaks = c(1, 10, 100, 1000, 10000), limits=c(1,NA), labels=c(1, 10, 100, 1000, 10000)) +
    annotation_logticks(sides='l', outside=TRUE, linewidth=(0.25 / (ggplot2::.pt * 72.27/96)), short = unit(0.4, "mm"), mid = unit(0.5, "mm"), long = unit(0, "cm")) +
    coord_flip(clip = "off") +
    facet_grid(cols = vars(ancestry_major), scales="free", space="free", axes="all") +
    theme(strip.placement = "outside",
          strip.clip = "off",
          axis.text.x=element_text(angle=45, hjust=1, size=6),
          axis.title.x=element_blank(),
          strip.text.x.top=element_text(hjust=0))

  sample_sizes

  ggsave(file.path(".", sprintf(
    "out.sample_sizes.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
         sample_sizes, width=120, height=45,units="mm")

  sample_sizes <- ggplot(data = combined_ordered,
                         aes(y= N, x=ancestry_minor, fill=ancestry)) +
    geom_col(colour="white") +
    scale_fill_manual(values=pop_colors) +
    guides(size= guide_legend(position = "bottom"),
           fill="none") +
    scale_y_log10(name="Sample size", breaks = c(1, 10, 100, 1000, 10000), limits=c(1,NA), labels=c(1, 10, 100, 1000, 10000)) +
    coord_cartesian(clip = "off") +
    facet_grid(cols = vars(ancestry_major), scales="free", space="free", switch="x") +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle=90), strip.clip = "off",
          axis.text.x=element_text(angle=45, hjust=1), axis.title.x =element_blank())

  sample_sizes

  ggsave(file.path(".", sprintf(
    "out.sample_sizes_unflipped.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
         sample_sizes, width=80, height=40,units="mm")

  sample_sizes <- ggplot(data = combined_ordered,
                         aes(y= ancestry_minor, x=N, fill=ancestry)) +
    geom_col(colour="white") +
    scale_fill_manual(values=pop_colors) +
    guides(size= guide_legend(position = "bottom"),
           fill="none") +
    scale_x_log10(name="Sample size per ancestry", breaks = c(1, 10, 100, 1000, 10000), limits=c(1,NA), labels=c(1, 10, 100, 1000, 10000)) +
    annotation_logticks(sides='b', outside=TRUE, linewidth=(0.5 / (ggplot2::.pt * 72.27/96)), short = unit(0.4, "mm"), mid = unit(0.5, "mm"), long = unit(0, "cm")) +
    coord_cartesian(clip = "off") +
    facet_grid(rows = vars(ancestry_major), scales="free", space="free", switch="y") +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle=90), strip.clip = "off",
          axis.text.x=element_text(angle=45, hjust=1), axis.title.x =element_blank())

  sample_sizes

  ggsave(file.path(".", sprintf(
    "out.sample_sizes_flipped.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
         sample_sizes, width=50, height=80,units="mm")

  plot_to_pdf(
    world, file.path(".", sprintf(
      "out.world.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 10, height = 8)

  plot_to_pdf(
    mds, file.path(".", sprintf(
      "out.mds.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 8)

  world <- ggplot(data = world_sf) +
    geom_sf(colour = "white", fill = "#e6e6e6") +
    coord_sf(expand = FALSE, xlim=c(-175, 180), ylim=c(-60, 85)) +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          legend.position = NULL)

  plot_to_pdf(
    world, file.path(".", sprintf(
      "out.world.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 4)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}