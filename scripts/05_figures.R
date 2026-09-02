##============================================================================##
## 05_figures.R
##============================================================================##

## !! Requires output from all other scripts to run !! ##

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(foreach)
library(rnaturalearth)
library(countrycode)
library(sf)
library(maps)
library(readr)
library(cowplot)
library(readxl)
library(stringr)
library(patchwork)
library(grid)
library(forcats)
library(rphylopic)
library(metafor)


##============================================================================##
## Load richness and abundance data -----------
##============================================================================##

# Combine abundance with mining data
abundance <- readxl::read_excel('data/Abundance.xlsx') %>%
  left_join( readxl::read_excel('data/Mines.xlsx')) %>%
  as_tibble %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>% # Replace literal "NA" strings with actual NA values
  mutate(Country_std = countrycode(Country, origin = "country.name", destination = "country.name")) # Standardize country names


# Combine richness with mining data
richness <- readxl::read_excel('data/Richness.xlsx') %>%
  left_join( readxl::read_excel('data/Mines.xlsx')) %>%
  as_tibble %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%  # Replace literal "NA" strings with actual NA values
  mutate(Country_std = countrycode(Country, origin = "country.name", destination = "country.name")) # Standardize country names

##============================================================================##


##============================================================================##
## Figure 1 -----------
##============================================================================##
# Load country polygons
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# Standardize country names for joining
world_sf$Country_std <- countrycode(world_sf$name, origin = "country.name", destination = "country.name")
world_sf <- world_sf[!is.na(world_sf$Country_std), ]

richness  <- richness  %>% mutate(lat = readr::parse_number(lat), lon = readr::parse_number(lon))
abundance <- abundance %>% mutate(lat = readr::parse_number(lat), lon = readr::parse_number(lon))

# Group richness per unique mine
richness_complete <- richness %>%
  group_by(Study_ID, Mine_ID) %>%
  summarise(
    Country_std     = first(Country_std),
    realm           = first(realm),
    material_group  = first(material_group),
    n_obs           = n(),
    has_exact       = any(!is.na(lat) & !is.na(lon)),
    final_lat       = if_else(has_exact, first(lat[!is.na(lat) & !is.na(lon)]), NA_real_),
    final_lon       = if_else(has_exact, first(lon[!is.na(lat) & !is.na(lon)]), NA_real_),
    coord_source    = if_else(has_exact, "exact", "random"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_obs))

# Group abundance per unique mine
abundance_complete <- abundance %>%
  group_by(Study_ID, Mine_ID) %>%
  summarise(
    Country_std     = first(Country_std),
    realm           = first(realm),
    material_group  = first(material_group),
    n_obs           = n(),
    has_exact       = any(!is.na(lat) & !is.na(lon)),
    final_lat       = if_else(has_exact, first(lat[!is.na(lat) & !is.na(lon)]), NA_real_),
    final_lon       = if_else(has_exact, first(lon[!is.na(lat) & !is.na(lon)]), NA_real_),
    coord_source    = if_else(has_exact, "exact", "random"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_obs))



# Uncertain mine locations are assigned to random positions within their country's polygon 
set.seed(27) 

sample_points_in_country <- function(country_geom, n) {
  pts <- st_sample(country_geom, size = n, type = "random")
  st_coordinates(pts)
}

## Richness 
richness_geo <- richness_complete %>%
  left_join(world_sf %>% dplyr::select(Country_std, geometry), by = "Country_std") %>%
  st_as_sf()

rich_unknown <- richness_geo %>% dplyr::filter(coord_source == "random")
rich_exact   <- richness_geo %>% dplyr::filter(coord_source == "exact")

rich_unknown_new <- rich_unknown %>%
  st_set_geometry("geometry") %>%
  dplyr::group_by(Country_std) %>%
  dplyr::group_modify(function(df, key) {
    n <- nrow(df)
    coords <- sample_points_in_country(df$geometry[1], n = n)
    df$final_lon <- coords[,1]
    df$final_lat <- coords[,2]
    df
  }) %>%
  dplyr::ungroup() %>%
  st_drop_geometry()

richness_complete <- bind_rows(
  rich_exact %>% st_drop_geometry(),
  rich_unknown_new
) %>%
  filter(!is.na(final_lat), !is.na(final_lon)) %>%  
  arrange(desc(n_obs))


## Abundance 
abundance_geo <- abundance_complete %>%
  left_join(world_sf %>% dplyr::select(Country_std, geometry), by = "Country_std") %>%
  st_as_sf()

abun_unknown <- abundance_geo %>% dplyr::filter(coord_source == "random")
abun_exact   <- abundance_geo %>% dplyr::filter(coord_source == "exact")

abun_unknown_new <- abun_unknown %>%
  st_set_geometry("geometry") %>%
  dplyr::group_by(Country_std) %>%
  dplyr::group_modify(function(df, key) {
    n <- nrow(df)
    coords <- sample_points_in_country(df$geometry[1], n = n)
    df$final_lon <- coords[,1]
    df$final_lat <- coords[,2]
    df
  }) %>%
  dplyr::ungroup() %>%
  st_drop_geometry()

abundance_complete <- bind_rows(
  abun_exact %>% st_drop_geometry(),
  abun_unknown_new
) %>%
  filter(!is.na(final_lat), !is.na(final_lon)) %>% 
  arrange(desc(n_obs))

# Ensure plotting order (larger symbols underneath smaller ones)
richness_complete  <- richness_complete  |> arrange(desc(n_obs))
abundance_complete <- abundance_complete |> arrange(desc(n_obs))


### Panel (a) - richness ####

# Base world map
world <- map_data("world") %>%
  dplyr::filter(region != "Antarctica")

# Define colors and shapes
realm_colors <- c("terrestrial" = "#1b9e77", "freshwater" = "#7570b3")

material_shapes <- c(
  "metals"      = 21,  
  "coal"        = 22,   
  "others"      = 24,  
  "aggregates"  = 23,  
  "gold"        = 25    
)



richness_complete <- richness_complete %>%
  mutate(obs_bin = case_when(
    n_obs <= 5              ~ "1–5",
    n_obs > 5  & n_obs <=10 ~ "6–10",
    n_obs > 10 & n_obs <=15 ~ "11–15",
    n_obs > 15              ~ ">15"
  )) %>%
  mutate(obs_bin = factor(obs_bin, levels = c("1–5", "6–10", "11–15", ">15")))



# Plot
richness_panel <- ggplot() +
  # Base world map
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray80", linewidth = 0.1) +
  
  # Richness observation points
  geom_point(data = richness_complete,
             aes(x = final_lon, y = final_lat,
                 size = obs_bin,
                 shape = material_group,
                 fill = realm,       
                 alpha = coord_source),
             color = "black",        
             stroke = 0.3) +          
  scale_fill_manual(
    values = realm_colors,
    labels = c("freshwater" = "Freshwater", "terrestrial" = "Terrestrial")
  ) +
  scale_shape_manual(
    values = material_shapes,
    labels = c(
      "metals"     = "Metals",
      "coal"       = "Coal",
      "others"     = "Others",
      "aggregates" = "Aggregates",
      "gold"       = "Gold"
    )
  )+
  scale_size_manual(
    name = "No. of observations",
    values = c("1–5" = 2, "6–10" = 4, "11–15" = 6, ">15" = 8)
  ) +
  scale_alpha_manual(values = c("exact" = 1, "random" = 0.5), guide = "none") +
  guides(
    fill  = guide_legend(
      override.aes = list(
        shape  = 21,     
        colour = "black",
        stroke = 0.3,
        size   = 5,
        alpha  = 1
      )
    ),
    shape = guide_legend(override.aes = list(size = 3.6, stroke = 1.1))
  ) +
  # Coordinate and theme
  coord_fixed(xlim = c(-180, 180), ylim = c(-60, 90)) +
  theme_minimal() +
  labs(
    title = "",
    fill  = "Realm",
    shape = "Commodity group"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


### Panel (b) - abundance ####

# Base world map
world <- map_data("world") %>%
  dplyr::filter(region != "Antarctica")

# Define colors and shapes
realm_colors <- c("terrestrial" = "#1b9e77", "freshwater" = "#7570b3")
material_shapes <- c(
  "metals"      = 21,  
  "coal"        = 22,   
  "others"      = 24,  
  "aggregates"  = 23,  
  "gold"        = 25    
)




abundance_complete <- abundance_complete %>%
  mutate(obs_bin = case_when(
    n_obs <= 25               ~ "1–25",
    n_obs > 25  & n_obs <= 50 ~ "26–50",
    n_obs > 50  & n_obs <= 75 ~ "51–75",
    n_obs > 75                ~ ">75"
  )) %>%
  mutate(obs_bin = factor(obs_bin, levels = c("1–25", "26–50", "51–75", ">75")))




# Plot
abundance_panel <- ggplot() +
  # Base world map
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray80", linewidth = 0.1) +
  
  # abundance observation points
  geom_point(data = abundance_complete,
             aes(x = final_lon, y = final_lat,
                 size = obs_bin,
                 shape = material_group,
                 fill = realm,     
                 alpha = coord_source),
             color = "black",        
             stroke = 0.3)  +         # outline thickness
  scale_fill_manual(
    values = realm_colors,
    labels = c("freshwater" = "Freshwater", "terrestrial" = "Terrestrial")
  ) +
  scale_shape_manual(
    values = material_shapes,
    labels = c(
      "metals"     = "Metals",
      "coal"       = "Coal",
      "others"     = "Others",
      "aggregates" = "Aggregates",
      "gold"       = "Gold"
    )
  )+
  scale_size_manual(
    name = "No. of observations",
    values = c("1–25" = 2, "26–50" = 4, "51–75" = 6, ">75" = 8)
  ) +
  scale_alpha_manual(values = c("exact" = 1, "random" = 0.5), guide = "none") +
  guides(
    fill  = guide_legend(
      override.aes = list(
        shape  = 21,     
        colour = "black",
        stroke = 0.3,
        size   = 5,
        alpha  = 1
      )
    ),
    shape = guide_legend(override.aes = list(size = 3.6, stroke = 1.1))
  ) +
  # Coordinate and theme
  coord_fixed(xlim = c(-180, 180), ylim = c(-60, 90)) +
  theme_minimal() +
  labs(
    title = "",
    fill  = "Realm",
    shape = "Material group"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


### Final combined plot ####

# Remove repeated legends
richness_panel_a <- richness_panel +
  guides(fill = "none", shape = "none", alpha = "none")+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.1),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "gray80")
  )


abundance_panel_b <- abundance_panel +
  guides(fill = "none", shape = "none", alpha = "none")+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.1),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "gray80")
  )


shared_legend <- get_legend(
  richness_panel +
    guides(size = "none", alpha = "none") +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.5, "cm")
    )
)


panels <- plot_grid(
  richness_panel_a, 
  abundance_panel_b,
  ncol = 1, 
  labels = c("(a) Richness", "(b) Abundance"), 
  label_size = 14
)

final_plot <- plot_grid(
  panels,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)


ggsave(
  "./figures/Figure 1.png",
  plot = final_plot,
  width = 12,  
  height = 8,
  dpi = 1000,
  bg = "white"
)


ggsave(
  "./figures/Figure 1.tif",
  plot = final_plot,
  width = 12,  
  height = 8,
  dpi = 1000,
  bg = "white"
)

##============================================================================##


##============================================================================##
## Figure 2 --------------
##============================================================================##
dir_rich_base <- "./results/richness/Study_Mine_VCV/Bracken"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken"

## Input files
rich_taxo_file  <- file.path(
  dir_rich_base,
  "taxo_group_realm",
  "richness_taxo_group_realm_summary.csv"
)
rich_realm_file <- file.path(
  dir_rich_base,
  "realm",
  "richness_realm_summary.csv"
)

abun_taxo_file  <- file.path(
  dir_abun_base,
  "taxo_group_realm",
  "abundance_taxo_group_realm_summary.csv"
)
abun_realm_file <- file.path(
  dir_abun_base,
  "realm",
  "abundance_realm_summary.csv"
)

## Read data
rich_taxo_raw  <- read.csv(rich_taxo_file,  check.names = FALSE)
rich_realm_raw <- read.csv(rich_realm_file, check.names = FALSE)

abun_taxo_raw  <- read.csv(abun_taxo_file,  check.names = FALSE)
abun_realm_raw <- read.csv(abun_realm_file, check.names = FALSE)


## Settings ##
cols <- c(
  freshwater  = "#0072B2",
  terrestrial = "#009E73"
)

# Silhouettes 
uuid_map <- c(
  terrestrial_plants        = "43afe2df-ab6c-47e6-a105-c0a82b8af1c5",
  terrestrial_vertebrates   = "fe09db98-924b-4ad4-bcfc-2b660732ff9d",
  terrestrial_invertebrates = "fbb67694-0fd2-4e90-b88d-efab9cbac37c",
  terrestrial_others        = "18b8f0a9-b280-4fb9-b750-495ceec4ef87",
  freshwater_plants         = "4e396beb-1555-47ee-9997-47fc970accf6",
  freshwater_vertebrates    = "63028840-34fa-44fb-8883-41b2866e61b8",
  freshwater_invertebrates  = "2757cb60-acda-422d-a318-85b812a658c2",
  freshwater_others         = "a0753559-d626-4b3a-9af7-39c79409b8f4"
)

## Taxon order used by every panel
taxon_order <- c("plants", "vertebrates", "invertebrates", "others")

## Shared left margin for richness panels
richness_left_margin <- 90

## Consistent row order for both richness and abundance panels
wanted_groups <- c(
  "freshwater_plants",
  "freshwater_vertebrates",
  "freshwater_invertebrates",
  "freshwater_others",
  "terrestrial_plants",
  "terrestrial_vertebrates",
  "terrestrial_invertebrates",
  "terrestrial_others"
)

order_map <- c(
  freshwater_plants         = 1,
  freshwater_vertebrates    = 2,
  freshwater_invertebrates  = 3,
  freshwater_others         = 4,
  terrestrial_plants        = 1,
  terrestrial_vertebrates   = 2,
  terrestrial_invertebrates = 3,
  terrestrial_others        = 4
)

## Icon nudges
icon_dx_map <- c(
  freshwater_plants         =  0.07 * 2,
  freshwater_vertebrates    = -0.04 * 2,
  freshwater_invertebrates  = -0.01 * 2,
  freshwater_others         =  0.03 * 2,
  terrestrial_plants        =  0.03 * 2,
  terrestrial_vertebrates   = -0.02 * 2,
  terrestrial_invertebrates =  0.07 * 2,
  terrestrial_others        =  0.04 * 2
)

## Icon heights
icon_height_map <- c(
  terrestrial_plants        = 0.37 * 1.15,
  terrestrial_vertebrates   = 0.26 * 1.15,
  terrestrial_invertebrates = 0.44 * 1.15,
  terrestrial_others        = 0.32 * 1.15,
  freshwater_plants         = 0.46 * 1.15,
  freshwater_vertebrates    = 0.23 * 1.15,
  freshwater_invertebrates  = 0.32 * 1.15,
  freshwater_others         = 0.47 * 1.15
)

N_dx_map <- c(
  freshwater_plants         = 0,
  freshwater_vertebrates    = 0,
  freshwater_invertebrates  = 0,
  freshwater_others         = 0,
  terrestrial_plants        = 0,
  terrestrial_vertebrates   = 0,
  terrestrial_invertebrates = 0,
  terrestrial_others        = 0
)

## Fetch silhouettes
img_keys <- unique(unname(uuid_map))
img_list <- setNames(
  lapply(img_keys, rphylopic::get_phylopic),
  img_keys
)


get_cell <- function(df, row, col, default = NA_real_) {
  match_idx <- which(tolower(df$stratum) == row)
  cn <- colnames(df)
  
  if (length(match_idx) == 1 && col %in% cn) {
    return(df[[col]][match_idx])
  } else {
    return(default)
  }
}

prepare_realm_df <- function(df_raw) {
  df <- df_raw %>%
    mutate(realm = tolower(stratum)) %>%
    rename(p = pval)
  
  df$p <- suppressWarnings(as.numeric(df$p))
  
  df %>%
    mutate(
      p_label = ifelse(
        is.na(p),
        "",
        ifelse(
          p < 0.01,
          "p < 0.01",
          paste0("p = ", sprintf("%.2f", p))
        )
      )
    )
}

finalize_taxo_df <- function(df) {
  realms_present <- unique(df$realm)
  
  scaffold <- expand.grid(
    realm       = realms_present,
    taxon_group = taxon_order,
    stringsAsFactors = FALSE
  )
  
  df <- scaffold %>%
    left_join(df, by = c("realm", "taxon_group")) %>%
    mutate(
      group_key      = paste0(realm, "_", taxon_group),
      uuid           = uuid_map[group_key],
      icon_dx        = unname(icon_dx_map[group_key]),
      N_dx           = unname(N_dx_map[group_key]),
      order_in_realm = unname(order_map[group_key]),
      plot_row       = ifelse(is.na(plot_row), FALSE, plot_row),
      k              = ifelse(is.na(k), 0, k),
      n_studies      = ifelse(is.na(n_studies), 0, n_studies),
      taxon_group    = factor(taxon_group, levels = taxon_order)
    ) %>%
    arrange(realm, taxon_group) %>%
    group_by(realm) %>%
    mutate(
      label_text = paste0(as.character(taxon_group), " (", k, ")"),
      label = factor(label_text, levels = rev(label_text))
    ) %>%
    ungroup()
  
  df
}

prepare_richness_taxo <- function(df_raw) {
  taxo_names <- strsplit(df_raw$stratum, "_")
  
  df <- df_raw %>%
    mutate(
      realm = tolower(sapply(taxo_names, `[`, 1)),
      taxon = tolower(
        sapply(taxo_names, function(x) paste(x[-1], collapse = "_"))
      )
    ) %>%
    rename(p = pval)
  
  df$p <- suppressWarnings(as.numeric(df$p))
  
  df <- df %>%
    mutate(
      taxon_group = case_when(
        taxon == "plants"        ~ "plants",
        taxon == "vertebrates"   ~ "vertebrates",
        taxon == "invertebrates" ~ "invertebrates",
        TRUE                      ~ "others"
      ),
      plot_row =
        is.finite(log_RR) &
        is.finite(ci_lower) &
        is.finite(ci_upper),
      p_label = ifelse(
        is.na(p),
        "",
        ifelse(
          p < 0.01,
          "p < 0.01",
          paste0("p = ", sprintf("%.2f", p))
        )
      ),
      N_label = paste0("(", k, ")")
    )
  
  finalize_taxo_df(df)
}

prepare_abundance_taxo <- function(df_raw) {
  df <- data.frame(
    group_key = wanted_groups,
    stringsAsFactors = FALSE
  )
  
  df$k <- sapply(df$group_key, function(g) {
    val <- get_cell(df_raw, g, "k", default = NA)
    if (is.na(val)) 0 else as.numeric(val)
  })
  
  df$n_studies <- sapply(df$group_key, function(g) {
    val <- get_cell(df_raw, g, "n_studies", default = NA)
    if (is.na(val)) 0 else as.numeric(val)
  })
  
  df$log_RR <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "log_RR", default = NA))
  })
  
  df$ci_lower <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "ci_lower", default = NA))
  })
  
  df$ci_upper <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "ci_upper", default = NA))
  })
  
  df$p <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "pval", default = NA))
  })
  
  split_keys <- strsplit(df$group_key, "_")
  df$realm <- sapply(split_keys, `[`, 1)
  df$taxon_group <- sapply(
    split_keys,
    function(x) paste(x[-1], collapse = "_")
  )
  
  df <- df %>%
    mutate(
      plot_row =
        is.finite(log_RR) &
        is.finite(ci_lower) &
        is.finite(ci_upper),
      p_label = ifelse(
        !plot_row,
        "",
        ifelse(
          is.na(p),
          "",
          ifelse(
            p < 0.01,
            "p < 0.01",
            paste0("p = ", sprintf("%.2f", p))
          )
        )
      ),
      N_label = paste0("(", k, ")")
    ) %>%
    select(-group_key)
  
  finalize_taxo_df(df)
}


build_core_plot <- function(
    df_taxo,
    df_realm,
    x_limits,
    x_breaks,
    metric_name,
    size_limits,
    size_breaks,
    show_x_axis = TRUE) {
  
  panel_width <- diff(x_limits)
  
  if (metric_name == "richness") {
    x_left <- x_limits[1] - 0.17 * panel_width
    left_plot_margin  <- richness_left_margin
    right_plot_margin <- 10
    x_N <- x_limits[1] - 0.06 * panel_width
  } else {
    x_left <- x_limits[1] - 0.07 * panel_width
    left_plot_margin  <- 40
    right_plot_margin <- 10
    x_N <- x_limits[1] - 0.06 * panel_width
  }
  
  ## Percent-change labels for significant taxonomic-group effects only
  ## Labels are placed at a fixed small offset to the right of x = 0
  pct_lab_df <- df_taxo %>%
    filter(plot_row, !is.na(p), p < 0.05) %>%
    mutate(
      pct_mean = (exp(log_RR)   - 1) * 100,
      pct_lo   = (exp(ci_lower) - 1) * 100,
      pct_hi   = (exp(ci_upper) - 1) * 100,
      pct_label = sprintf(
        "%s%.0f%% [%.0f, %.0f]%%",
        ifelse(pct_mean >= 0, "+", ""),
        pct_mean,
        pct_lo,
        pct_hi
      ),
      x_whisker = 0 + 0.02 * panel_width
    )
  
  p <- ggplot(df_taxo, aes(y = label)) +
    ## Realm-level CI band
    geom_rect(
      data = df_realm,
      aes(xmin = ci_lower, xmax = ci_upper, fill = realm),
      ymin = -Inf,
      ymax = Inf,
      inherit.aes = FALSE,
      alpha = 0.05,
      show.legend = FALSE
    ) +
    ## Realm-level CI limits and mean
    geom_vline(
      data = df_realm,
      aes(xintercept = ci_lower, colour = realm),
      linetype = "22",
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    geom_vline(
      data = df_realm,
      aes(xintercept = ci_upper, colour = realm),
      linetype = "22",
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    geom_vline(
      data = df_realm,
      aes(xintercept = log_RR, colour = realm),
      linewidth = 1.4,
      show.legend = FALSE
    ) +
    ## Zero line
    geom_vline(
      xintercept = 0,
      linetype = 2,
      colour = "black",
      linewidth = 1.0
    ) +
    ## Taxonomic-group CIs and means
    geom_segment(
      data = df_taxo %>% filter(plot_row),
      aes(
        x = ci_lower,
        xend = ci_upper,
        y = label,
        yend = label
      ),
      linewidth = 1,
      colour = "grey25",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df_taxo %>% filter(plot_row),
      aes(x = log_RR, size = n_studies),
      colour = "grey25",
      stroke = 0.6
    ) +
    scale_size_area(
      max_size = 7,
      limits   = size_limits,
      breaks   = size_breaks,
      guide    = "none"
    ) +
    ## Percent-change labels for significant taxonomic-group effects
    geom_text(
      data = pct_lab_df,
      aes(x = x_whisker, y = label, label = pct_label),
      hjust = -0.05,
      size = 4.0,
      fontface = "plain",
      colour = "black",
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    ## Number of observations shown beside the taxon silhouettes/rows
    geom_text(
      data = df_taxo,
      aes(
        x = x_N,
        y = label,
        label = N_label
      ),
      hjust = 1,
      vjust = 0.5,
      size = 4.2,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = cols,
      name = NULL,
      breaks = names(cols),
      guide = "none"
    ) +
    scale_fill_manual(
      values = cols,
      name = NULL,
      breaks = names(cols),
      guide = "none"
    ) +
    labs(
      x = if (metric_name == "richness") {
        "Change in species richness (lnRR)"
      } else {
        "Change in species abundance (lnRR)"
      },
      y = NULL
    ) +
    theme_bw(base_size = 20) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey90",
        linewidth = 0.35
      ),
      strip.background  = element_blank(),
      strip.text.y.left = element_blank(),
      strip.placement   = "outside",
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.title.x      = element_text(
        size = 19,
        margin = margin(t = 10)
      ),
      axis.text.x       = element_text(size = 16),
      legend.position   = "none",
      plot.margin       = margin(
        8,
        right_plot_margin,
        14,
        left_plot_margin
      )
    ) +
    coord_cartesian(
      xlim = x_limits,
      clip = "off"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_discrete(drop = FALSE)
  
  if (!show_x_axis) {
    p <- p +
      theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  ## Add silhouettes only to richness panels
  if (metric_name == "richness") {
    rows_icon <- df_taxo %>%
      distinct(label, uuid, group_key, icon_dx) %>%
      filter(!is.na(uuid)) %>%
      mutate(
        icon_height = unname(icon_height_map[group_key])
      )
    
    for (i in seq_len(nrow(rows_icon))) {
      p <- p + rphylopic::add_phylopic(
        img    = img_list[[rows_icon$uuid[i]]],
        x      = x_left + rows_icon$icon_dx[i],
        y      = rows_icon$label[i],
        height = rows_icon$icon_height[i]
      )
    }
  }
  
  p
}

## Prepare data
rich_taxo_df  <- prepare_richness_taxo(rich_taxo_raw)
rich_realm_df <- prepare_realm_df(rich_realm_raw)

abun_taxo_df  <- prepare_abundance_taxo(abun_taxo_raw)
abun_realm_df <- prepare_realm_df(abun_realm_raw)

x_limits_rich <- c(-2.8, 2.8)
x_breaks_rich <- seq(-2, 2, by = 2)

x_limits_abun <- c(-8.36, 8.36)
x_breaks_abun <- c(-6, -3, 0, 3, 6)


## Split by realm
rich_fw <- rich_taxo_df %>%
  filter(realm == "freshwater") %>%
  mutate(label = fct_drop(label))

rich_tr <- rich_taxo_df %>%
  filter(realm == "terrestrial") %>%
  mutate(label = fct_drop(label))

rrealm_fw <- rich_realm_df %>%
  filter(realm == "freshwater")

rrealm_tr <- rich_realm_df %>%
  filter(realm == "terrestrial")

abun_fw <- abun_taxo_df %>%
  filter(realm == "freshwater") %>%
  mutate(label = fct_drop(label))

abun_tr <- abun_taxo_df %>%
  filter(realm == "terrestrial") %>%
  mutate(label = fct_drop(label))

arealm_fw <- abun_realm_df %>%
  filter(realm == "freshwater")

arealm_tr <- abun_realm_df %>%
  filter(realm == "terrestrial")


all_plotted_n <- c(
  rich_fw$n_studies[rich_fw$plot_row],
  rich_tr$n_studies[rich_tr$plot_row],
  abun_fw$n_studies[abun_fw$plot_row],
  abun_tr$n_studies[abun_tr$plot_row]
)

all_plotted_n <- all_plotted_n[is.finite(all_plotted_n)]

if (length(all_plotted_n) == 0) {
  stop("No finite n_studies values were found for plotted taxonomic groups.")
}

size_limits_global <- c(0, max(all_plotted_n))
size_breaks_global <- pretty(all_plotted_n, n = 4)
size_breaks_global <- size_breaks_global[
  size_breaks_global >= size_limits_global[1] &
    size_breaks_global <= size_limits_global[2]
]


## Richness panels
p_rich_fw <- build_core_plot(
  df_taxo     = rich_fw,
  df_realm    = rrealm_fw,
  x_limits    = x_limits_rich,
  x_breaks    = x_breaks_rich,
  metric_name = "richness",
  size_limits = size_limits_global,
  size_breaks = size_breaks_global,
  show_x_axis = FALSE
)

p_rich_tr <- build_core_plot(
  df_taxo     = rich_tr,
  df_realm    = rrealm_tr,
  x_limits    = x_limits_rich,
  x_breaks    = x_breaks_rich,
  metric_name = "richness",
  size_limits = size_limits_global,
  size_breaks = size_breaks_global,
  show_x_axis = TRUE
)


## Abundance panels
p_abun_fw <- build_core_plot(
  df_taxo     = abun_fw,
  df_realm    = arealm_fw,
  x_limits    = x_limits_abun,
  x_breaks    = x_breaks_abun,
  metric_name = "abundance",
  size_limits = size_limits_global,
  size_breaks = size_breaks_global,
  show_x_axis = FALSE
)

p_abun_tr <- build_core_plot(
  df_taxo     = abun_tr,
  df_realm    = arealm_tr,
  x_limits    = x_limits_abun,
  x_breaks    = x_breaks_abun,
  metric_name = "abundance",
  size_limits = size_limits_global,
  size_breaks = size_breaks_global,
  show_x_axis = TRUE
)


## Panel titles
## Only the (a)-(d) prefixes are bold
p_rich_fw <- p_rich_fw +
  labs(title = "<b>&#40;a&#41;</b> Freshwater species richness") +
  theme(plot.title = ggtext::element_markdown(
    hjust = 0, size = 19, colour = "grey20",
    margin = margin(b = 7)
  ))

p_abun_fw <- p_abun_fw +
  labs(title = "<b>&#40;b&#41;</b> Freshwater species abundance") +
  theme(plot.title = ggtext::element_markdown(
    hjust = 0, size = 19, colour = "grey20",
    margin = margin(b = 7)
  ))

p_rich_tr <- p_rich_tr +
  labs(title = "<b>&#40;c&#41;</b> Terrestrial species richness") +
  theme(plot.title = ggtext::element_markdown(
    hjust = 0, size = 19, colour = "grey20",
    margin = margin(b = 7)
  ))

p_abun_tr <- p_abun_tr +
  labs(title = "<b>&#40;d&#41;</b> Terrestrial species abundance") +
  theme(plot.title = ggtext::element_markdown(
    hjust = 0, size = 19, colour = "grey20",
    margin = margin(b = 7)
  ))

## Combine 2 x 2 panel grid
p_combined <-
  (p_rich_fw | p_abun_fw) /
  (p_rich_tr | p_abun_tr)


## Save
ggsave(
  "./figures/Figure 2.png",
  p_combined,
  width  = 360,
  height = 250,
  units  = "mm",
  bg     = "white",
  dpi    = 1000
)

ggsave(
  "./figures/Figure 2.tif",
  p_combined,
  width  = 360,
  height = 250,
  units  = "mm",
  bg     = "white",
  dpi    = 1000,
  compression = "lzw"
)
##============================================================================##


##============================================================================##
## Figure 3  --------
##============================================================================##
dir_rich_base <- "./results/richness/Study_Mine_VCV/Bracken/material_group_moderator"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/material_group_moderator"

rich_means_file <- file.path(dir_rich_base, "richness_material_moderator_means_all.csv")
abun_means_file <- file.path(dir_abun_base, "abundance_material_moderator_means_all.csv")

cols <- c(
  "Freshwater"  = "#0072B2",
  "Terrestrial" = "#009E73"
)

material_order <- c("Coal", "Gold", "Metals", "Aggregates")

kept_materials <- c("metals", "coal", "gold", "aggregates")


## Number of studies
abundance <- readxl::read_excel("data/Abundance.xlsx") %>%
  left_join(readxl::read_excel("data/Mines.xlsx")) %>%
  as_tibble() %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA")))

richness <- readxl::read_excel("data/Richness.xlsx") %>%
  left_join(readxl::read_excel("data/Mines.xlsx")) %>%
  as_tibble() %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA")))

rich_study_counts <- richness %>%
  mutate(
    realm = str_to_lower(realm),
    material_group = str_to_lower(material_group)
  ) %>%
  filter(material_group %in% kept_materials) %>%
  distinct(Study_ID, realm, material_group) %>%
  count(realm, material_group, name = "n_studies")

abun_study_counts <- abundance %>%
  mutate(
    realm = str_to_lower(realm),
    material_group = str_to_lower(material_group)
  ) %>%
  filter(material_group %in% kept_materials) %>%
  distinct(Study_ID, realm, material_group) %>%
  count(realm, material_group, name = "n_studies")


pretty_realm <- function(x) {
  case_when(
    str_to_lower(x) == "freshwater"  ~ "Freshwater",
    str_to_lower(x) == "terrestrial" ~ "Terrestrial",
    TRUE ~ str_to_title(x)
  )
}

pretty_material <- function(x) {
  x <- str_to_lower(x)
  case_when(
    x == "coal"       ~ "Coal",
    x == "gold"       ~ "Gold",
    x == "metals"     ~ "Metals",
    x == "aggregates" ~ "Aggregates",
    TRUE              ~ str_to_title(x)
  )
}

sig_stars <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 1e-4   ~ "****",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ ""
  )
}

read_realm_material <- function(path, outcome_name, study_counts) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(stratum_col == "realm") %>%
    mutate(
      outcome        = outcome_name,
      realm          = str_to_lower(stratum),
      realm_label    = pretty_realm(realm),
      material_group = str_to_lower(material_group),
      material_label = pretty_material(material_group),
      modelled       = is.finite(log_RR),
      star_label     = sig_stars(pval),
      k_label        = paste0("(k=", k, ")"),
      pct_label      = ifelse(
        modelled & !is.na(pval) & pval < 0.05,
        sprintf("%s%.0f%% [%.0f, %.0f]%%",
                ifelse(pct_change >= 0, "+", ""),
                pct_change, pct_change_lower, pct_change_upper),
        ""
      )
    ) %>%
    filter(material_group %in% kept_materials, modelled) %>%
    left_join(
      study_counts,
      by = c("realm", "material_group")
    )
}


## Panel function
format_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p))))
}

make_panel <- function(df, outcome_name, x_title, show_legend = FALSE) {
  
  keep_materials  <- df %>% distinct(material_label) %>% pull(material_label)
  material_levels <- material_order[material_order %in% keep_materials]
  
  y_base <- tibble(
    material_label = factor(material_levels, levels = material_levels),
    y = rev(seq_along(material_levels))
  )
  
  df2 <- df %>%
    mutate(material_label = factor(material_label, levels = material_levels)) %>%
    inner_join(y_base, by = "material_label") %>%
    mutate(
      y_pos = case_when(
        realm_label == "Freshwater"  ~ y + 0.14,
        realm_label == "Terrestrial" ~ y - 0.14
      ),
      k_y    = y_pos,
      star_y = y_pos
    )
  
  if (outcome_name == "Abundance") {
    max_abs    <- 12
    x_breaks   <- seq(-10, 10, by = 5)
  } else {
    max_abs <- max(abs(c(df2$ci_lower, df2$ci_upper)), na.rm = TRUE)
    max_abs <- ceiling(max_abs * 1.05 * 2) / 2
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    x_breaks <- waiver()
  }
  
  panel_width <- 2 * max_abs
  df2 <- df2 %>%
    mutate(
      text_x = ci_upper + 0.015 * panel_width,
      star_x = text_x + nchar(k_label) * 0.018 * panel_width + 0.045 * panel_width
    )
  
  qm_df <- df2 %>%
    distinct(realm_label, QM_pval) %>%
    arrange(realm_label) %>%
    mutate(
      qm_label = paste0("Q<sub>M</sub> ", format_p(QM_pval)),
      vjust_val = 1.3 + (row_number() - 1) * 1.6
    )
  
  p <- ggplot(df2, aes(x = log_RR, y = y_pos, colour = realm_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey45") +
    geom_segment(
      aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos),
      linewidth = 0.9
    ) +
    geom_point(aes(size = n_studies)) +
    geom_text(
      aes(x = text_x, y = k_y, label = k_label),
      colour = "black", hjust = 0, size = 2.8, show.legend = FALSE
    ) +
    geom_text(
      aes(x = star_x, y = star_y, label = star_label),
      colour = "black", hjust = 0, size = 4.2, show.legend = FALSE
    ) +
    geom_richtext(
      data = qm_df,
      aes(x = Inf, y = Inf, label = qm_label, colour = realm_label, vjust = vjust_val),
      hjust = 1.05, size = 2.6, fontface = "italic",
      inherit.aes = FALSE, show.legend = FALSE,
      label.color = NA, fill = NA,
      label.padding = unit(c(0, 0, 0, 0), "pt"),
      label.margin  = unit(c(0, 0, 0, 0), "pt")
    ) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_size_area(
      max_size = 4,
      limits = size_limits_global,
      guide = "none"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(
      breaks = y_base$y,
      labels = as.character(y_base$material_label),
      limits = c(0.5, max(y_base$y) + 0.45)
    ) +
    coord_cartesian(xlim = c(-max_abs, max_abs), clip = "off") +
    labs(x = x_title, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      plot.margin = margin(2, 15, 2, 6)
    )
  
  p
}


## Read data
rich <- read_realm_material(rich_means_file, "Richness", rich_study_counts)
abun <- read_realm_material(abun_means_file, "Abundance", abun_study_counts)

## Common point-size scale across richness and abundance
size_limits_global <- c(
  0,
  max(c(rich$n_studies, abun$n_studies), na.rm = TRUE)
)

## Build panels
p_rich <- make_panel(
  rich,
  "Richness",
  "Change in species richness (lnRR)",
  show_legend = TRUE
) +
  labs(title = "<b>&#40;a&#41;</b> Species richness") +
  theme(
    plot.title = ggtext::element_markdown(
      hjust = 0, size = 10, colour = "black",
      margin = margin(b = 4)
    )
  )

p_abun <- make_panel(
  abun,
  "Abundance",
  "Change in species abundance (lnRR)",
  show_legend = FALSE
) +
  labs(title = "<b>&#40;b&#41;</b> Species abundance") +
  theme(
    plot.title = ggtext::element_markdown(
      hjust = 0, size = 10, colour = "black",
      margin = margin(b = 4)
    )
  )

p_final <- p_rich + p_abun +
  plot_layout(guides = "collect") &
  theme(
    legend.position    = "bottom",
    legend.margin      = margin(t = 0, b = 0),
    legend.box.margin  = margin(t = -4, b = 0)
  )


ggsave(
  "./figures/Figure 3.png",
  p_final, width = 7, height = 4, dpi = 600
)
ggsave(
  "./figures/Figure 3.tif",
  p_final, width = 7, height = 4, dpi = 600, compression = "lzw"
)

##============================================================================##

##============================================================================##
## Figure 4 -----------------
##============================================================================##
TRAIT_RESULTS <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/trait_moderator/abundance_trait_univariate_means_all.csv"
TRAIT_TAXA <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/trait_moderator/taxa_trait_table.csv"

OUTDIR <- "./figures/Final"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FOCAL_TRAITS <- c("Disp", "Drft", "Exit")

point_col <- "#0072B2"  

trait_titles <- c(
  Disp = "Female dispersal",
  Drft = "Occurrence in drift",
  Exit = "Ability to exit water"
)

panel_letters <- c(
  Disp = "a",
  Drft = "b",
  Exit = "c"
)

level_lookup <- list(
  Disp = c("1" = "Low", "2" = "High"),
  Exit = c("1" = "Absent", "2" = "Present"),
  Drft = c("1" = "Rare", "2" = "Common", "3" = "Abundant")
)

## Display order
level_order_by_trait <- list(
  Disp = c("Low", "High"),
  Drft = c("Rare", "Common", "Abundant"),
  Exit = c("Absent", "Present")
)

format_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p))))
}

sig_stars <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 1e-4   ~ "****",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ ""
  )
}


## Number of studies
trait_taxa <- read_csv(TRAIT_TAXA, show_col_types = FALSE) %>%
  select(genus, all_of(FOCAL_TRAITS)) %>%
  distinct()

study_counts <- abundance %>%
  filter(!is.na(genus), !is.na(Study_ID)) %>%
  select(Study_ID, genus) %>%
  distinct() %>%
  inner_join(trait_taxa, by = "genus") %>%
  pivot_longer(
    cols = all_of(FOCAL_TRAITS),
    names_to = "trait",
    values_to = "level"
  ) %>%
  filter(!is.na(level)) %>%
  mutate(level = as.character(level)) %>%
  distinct(Study_ID, trait, level) %>%
  count(trait, level, name = "n_studies")


## Read + prepare
df <- read_csv(TRAIT_RESULTS, show_col_types = FALSE) %>%
  filter(trait %in% FOCAL_TRAITS) %>%
  mutate(level = as.character(level)) %>%
  left_join(study_counts, by = c("trait", "level")) %>%
  rowwise() %>%
  mutate(level_label = unname(level_lookup[[trait]][level])) %>%
  ungroup() %>%
  mutate(
    modelled   = is.finite(log_RR),
    star_label = ifelse(modelled, sig_stars(pval), "")
  ) %>%
  filter(modelled)

has_qm <- "QM_pval" %in% names(df)

## Shared point-size scale across all 3 traits
size_limits <- c(0, max(df$n_studies, na.rm = TRUE))

## Shared x-range across all 3 traits, so panels are directly comparable
max_abs <- max(abs(c(df$ci_lower, df$ci_upper)), na.rm = TRUE)
max_abs <- ceiling(max_abs * 1.15 * 2) / 2
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
x_limits <- c(-max_abs, max_abs)
x_breaks <- seq(-floor(max_abs), floor(max_abs), by = 2)


build_trait_panel <- function(trait_code, show_x_axis) {
  
  df_tr <- df %>%
    filter(trait == trait_code) %>%
    mutate(
      level_label = factor(level_label, levels = level_order_by_trait[[trait_code]]),
      y_label     = paste0(as.character(level_label), " (k=", k, ")")
    ) %>%
    arrange(level_label) %>%
    mutate(y_label = factor(y_label, levels = rev(y_label)))
  
  df_tr <- df_tr %>% mutate(star_x = ci_upper + 0.06 * max_abs)
  
  qm_label <- if (has_qm) paste0("Q<sub>M</sub> ", format_p(unique(df_tr$QM_pval)[1])) else ""
  
  p <- ggplot(df_tr, aes(x = log_RR, y = y_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey45") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0,
                   linewidth = 0.9, colour = point_col) +
    geom_point(aes(size = n_studies), colour = point_col) +
    geom_text(
      aes(x = star_x, y = y_label, label = star_label),
      colour = "black", fontface = "bold", hjust = 0, size = 4.0
    ) +
    scale_size_area(
      max_size = 2.5,
      limits = size_limits,
      guide = "none"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(
      title = paste0(
        "<b>&#40;", panel_letters[trait_code], "&#41;</b> ",
        trait_titles[trait_code]
      ),
      x = if (show_x_axis) "Change in species abundance (lnRR)" else NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = ggtext::element_markdown(
        size = 8.5, hjust = 0, margin = margin(b = 4)
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.text.y  = element_text(size = 8),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8),
      plot.margin  = margin(4, 20, 4, 6)
    )
  
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (has_qm) {
    p <- p + annotate(
      "richtext", x = Inf, y = Inf, label = qm_label,
      hjust = 1, vjust = 1, size = 2.5, fontface = "italic", colour = "grey30",
      label.color = NA, fill = NA
    )
  }
  
  p
}


## Build pabels + combine
p_disp <- build_trait_panel("Disp", show_x_axis = FALSE)
p_drft <- build_trait_panel("Drft", show_x_axis = FALSE)
p_exit <- build_trait_panel("Exit", show_x_axis = TRUE)

## Panel heights proportional to number of levels in each trait (2, 3, 2),
## so row spacing stays consistent across panels
n_levels <- sapply(FOCAL_TRAITS, function(tr) length(level_order_by_trait[[tr]]))

p_combined <- p_disp / p_drft / p_exit +
  plot_layout(heights = n_levels[c("Disp", "Drft", "Exit")])


## Save
ggsave(
  file.path(OUTDIR, "Figure 4.png"),
  p_combined, width = 90, height = 100, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  file.path(OUTDIR, "Figure 4.tif"),
  p_combined, width = 90, height = 100, units = "mm", dpi = 600, bg = "white", compression = "lzw"
)

##============================================================================##


##============================================================================##
## Figure S1 -------------
##============================================================================##
dir_rich_base <- "./results/richness/Study_Mine_VCV/Bracken/realm/richness_realm_models.rds"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/realm/abundance_realm_models.rds"


get_model_by_realm <- function(model_list, realm_name) {
  nms <- tolower(names(model_list))
  idx <- match(tolower(realm_name), nms)
  
  if (is.na(idx)) {
    stop(
      "Could not find realm '", realm_name, "' in model list. Available names are: ",
      paste(names(model_list), collapse = ", ")
    )
  }
  
  model_list[[idx]]
}

extract_sampling_variance <- function(m) {
  if (!is.null(m$vi)) {
    return(as.numeric(m$vi))
  }
  
  if (!is.null(m$V)) {
    V <- m$V
    
    if (is.matrix(V) || inherits(V, "Matrix")) {
      return(as.numeric(diag(as.matrix(V))))
    } else {
      return(as.numeric(V))
    }
  }
  
}

extract_effects_from_model <- function(m, realm_name, metric_name) {
  yi <- as.numeric(m$yi)
  vi <- extract_sampling_variance(m)
  
  if (length(yi) != length(vi)) {
    stop(
      "Length mismatch for ", metric_name, " / ", realm_name,
      ": length(yi) = ", length(yi),
      ", length(vi) = ", length(vi)
    )
  }
  
  tibble(
    metric = metric_name,
    realm = tolower(realm_name),
    RR_log = yi,
    SamplingVariance = vi
  ) %>%
    filter(
      is.finite(RR_log),
      is.finite(SamplingVariance),
      SamplingVariance >= 0
    ) %>%
    mutate(
      se = sqrt(SamplingVariance),
      lowerCI = RR_log - 1.96 * se,
      upperCI = RR_log + 1.96 * se
    ) %>%
    arrange(RR_log) %>%
    mutate(plottingID = row_number())
}

extract_model_summary <- function(m, realm_name, metric_name) {
  pr <- predict(m)
  
  tibble(
    metric = metric_name,
    realm = tolower(realm_name),
    plottingID = -1,
    logRR = as.numeric(pr$pred)[1],
    lowerCI = as.numeric(pr$ci.lb)[1],
    upperCI = as.numeric(pr$ci.ub)[1]
  )
}

create_forest_plot <- function(data, model_data, realm_name, metric_name, color, xlim_vals = NULL) {
  
  dat <- data %>%
    filter(
      realm == tolower(realm_name),
      metric == metric_name
    ) %>%
    arrange(RR_log) %>%
    mutate(plottingID = row_number())
  
  md <- model_data %>%
    filter(
      realm == tolower(realm_name),
      metric == metric_name
    )
  
  p <- ggplot(
    data = dat,
    aes(y = plottingID, x = RR_log, xmin = lowerCI, xmax = upperCI)
  ) +

    geom_linerange(linewidth = 0.45, color = "grey60") +
    geom_point(
      size = 0.8,
      shape = 21,
      fill = color,
      colour = "black",
      stroke = 0.2
    ) +
    
    geom_rect(
      data = md,
      aes(xmin = lowerCI, xmax = upperCI),
      ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE,
      fill = color,
      alpha = 0.10
    ) +
    geom_vline(
      data = md,
      aes(xintercept = logRR),
      color = color,
      linewidth = 1
    ) +
    geom_vline(
      data = md,
      aes(xintercept = lowerCI),
      color = color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_vline(
      data = md,
      aes(xintercept = upperCI),
      color = color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    
    
    # Model estimate shown as diamond below rows
    geom_pointrange(
      data = md,
      aes(y = plottingID, x = logRR, xmin = lowerCI, xmax = upperCI),
      inherit.aes = FALSE,
      color = color,
      shape = 18,
      linewidth = 0.7,
      size = 1.1
    ) +
    geom_point(
      data = md,
      aes(y = plottingID, x = logRR),
      inherit.aes = FALSE,
      color = color,
      shape = 18,
      size = 5
    ) +
    
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.04))) +
    
    labs(
      x = expression(paste("Effect size (", lnRR, ")")),
      y = NULL
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(10, 10, 10, 10),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length.x = grid::unit(0.15, "cm")
    )
  
  if (!is.null(xlim_vals)) {
    p <- p + coord_cartesian(xlim = xlim_vals, clip = "off")
  } else {
    p <- p + coord_cartesian(clip = "off")
  }
  
  return(p)
}

add_panel_title <- function(p, panel_label, title_text) {
  
  title_expr <- bquote(
    bold(.(panel_label)) ~ .(title_text)
  )
  
  cowplot::ggdraw() +
    cowplot::draw_plot(
      p,
      x = 0,
      y = 0,
      width = 1,
      height = 0.94
    ) +
    cowplot::draw_label(
      title_expr,
      x = 0.03,
      y = 0.97,
      hjust = 0,
      vjust = 1,
      size = 15
    )
}

## Load RDS model objects
rich_models <- readRDS(dir_rich_base)
abun_models <- readRDS(dir_abun_base)

print(names(rich_models))
print(names(abun_models))

## Extract individual effect sizes and model summaries
rich_terrestrial <- get_model_by_realm(rich_models, "terrestrial")
rich_freshwater  <- get_model_by_realm(rich_models, "freshwater")

abun_terrestrial <- get_model_by_realm(abun_models, "terrestrial")
abun_freshwater  <- get_model_by_realm(abun_models, "freshwater")

forest_df <- bind_rows(
  extract_effects_from_model(rich_terrestrial, "terrestrial", "Richness"),
  extract_effects_from_model(rich_freshwater,  "freshwater",  "Richness"),
  extract_effects_from_model(abun_terrestrial, "terrestrial", "Abundance"),
  extract_effects_from_model(abun_freshwater,  "freshwater",  "Abundance")
)

## Remove the extreme freshwater-abundance effect size from the plotted rows only.
## This record has an implausibly large negative lnRR, around -100, which stretches
## the x-axis and makes all other individual effect sizes unreadable (will be 
## mentioned in the fig caption instead). The model summary shown in the panel 
## is still taken from the fitted model
forest_df <- forest_df %>%
  filter(!(metric == "Abundance" & realm == "freshwater" & lowerCI < -99))

model_summary_df <- bind_rows(
  extract_model_summary(rich_terrestrial, "terrestrial", "Richness"),
  extract_model_summary(rich_freshwater,  "freshwater",  "Richness"),
  extract_model_summary(abun_terrestrial, "terrestrial", "Abundance"),
  extract_model_summary(abun_freshwater,  "freshwater",  "Abundance")
)

## Create four panels
col_terrestrial <- "#009E73"
col_freshwater  <- "#0072B2"

rich_terr_plot <- create_forest_plot(
  data = forest_df,
  model_data = model_summary_df,
  realm_name = "terrestrial",
  metric_name = "Richness",
  color = col_terrestrial
)

rich_fresh_plot <- create_forest_plot(
  data = forest_df,
  model_data = model_summary_df,
  realm_name = "freshwater",
  metric_name = "Richness",
  color = col_freshwater
)

abun_terr_plot <- create_forest_plot(
  data = forest_df,
  model_data = model_summary_df,
  realm_name = "terrestrial",
  metric_name = "Abundance",
  color = col_terrestrial
)

abun_fresh_plot <- create_forest_plot(
  data = forest_df,
  model_data = model_summary_df,
  realm_name = "freshwater",
  metric_name = "Abundance",
  color = col_freshwater
)

rich_terr_plot  <- add_panel_title(rich_terr_plot,  "(a)", "Terrestrial species richness")
rich_fresh_plot <- add_panel_title(rich_fresh_plot, "(b)", "Freshwater species richness")
abun_terr_plot  <- add_panel_title(abun_terr_plot,  "(c)", "Terrestrial species abundance")
abun_fresh_plot <- add_panel_title(abun_fresh_plot, "(d)", "Freshwater species abundance")


## Combine into one final 4-panel figure
final_forest_plot <- plot_grid(
  rich_terr_plot, rich_fresh_plot,
  abun_terr_plot, abun_fresh_plot,
  nrow = 2,
  align = "hv"
)


## Save
ggsave(
  "./figures/Figure S1.png",
  plot = final_forest_plot,
  width = 12,
  height = 10,
  dpi = 1000,
  bg = "white"
)

ggsave(
  "./figures/Figure S1.tif",
  plot = final_forest_plot,
  width = 12,
  height = 10,
  dpi = 1000,
  bg = "white",
  compression = "lzw"
)
##============================================================================##


##============================================================================##
## Figures S2-S3 ---------------
##============================================================================##
RE_LABEL_RICH <- "Study_Mine_VCV"
RE_LABEL_ABUN <- "Study_Mine_Binomial_VCV"
RUN_DIR       <- "Bracken"   

dir_mod <- "./figures/"

## Settings (identical to Figure 2)
cols <- c(
  freshwater  = "#0072B2",
  terrestrial = "#009E73"
)

uuid_map <- c(
  terrestrial_plants        = "43afe2df-ab6c-47e6-a105-c0a82b8af1c5",
  terrestrial_vertebrates   = "fe09db98-924b-4ad4-bcfc-2b660732ff9d",
  terrestrial_invertebrates = "fbb67694-0fd2-4e90-b88d-efab9cbac37c",
  terrestrial_others        = "18b8f0a9-b280-4fb9-b750-495ceec4ef87",
  freshwater_plants         = "4e396beb-1555-47ee-9997-47fc970accf6",
  freshwater_vertebrates    = "63028840-34fa-44fb-8883-41b2866e61b8",
  freshwater_invertebrates  = "2757cb60-acda-422d-a318-85b812a658c2",
  freshwater_others         = "a0753559-d626-4b3a-9af7-39c79409b8f4"
)

taxon_order <- c("plants", "vertebrates", "invertebrates", "others")

richness_left_margin <- 90

wanted_groups <- c(
  "freshwater_plants",
  "freshwater_vertebrates",
  "freshwater_invertebrates",
  "freshwater_others",
  "terrestrial_plants",
  "terrestrial_vertebrates",
  "terrestrial_invertebrates",
  "terrestrial_others"
)

order_map <- c(
  freshwater_plants         = 1,
  freshwater_vertebrates    = 2,
  freshwater_invertebrates  = 3,
  freshwater_others         = 4,
  terrestrial_plants        = 1,
  terrestrial_vertebrates   = 2,
  terrestrial_invertebrates = 3,
  terrestrial_others        = 4
)

icon_dx_map <- c(
  freshwater_plants         =  0.07 * 2,
  freshwater_vertebrates    = -0.04 * 2,
  freshwater_invertebrates  = -0.01 * 2,
  freshwater_others         =  0.03 * 2,
  terrestrial_plants        =  0.03 * 2,
  terrestrial_vertebrates   = -0.02 * 2,
  terrestrial_invertebrates =  0.07 * 2,
  terrestrial_others        =  0.04 * 2
)

icon_height_map <- c(
  terrestrial_plants        = 0.37 * 1.15,
  terrestrial_vertebrates   = 0.26 * 1.15,
  terrestrial_invertebrates = 0.44 * 1.15,
  terrestrial_others        = 0.32 * 1.15,
  freshwater_plants         = 0.46 * 1.15,
  freshwater_vertebrates    = 0.23 * 1.15,
  freshwater_invertebrates  = 0.32 * 1.15,
  freshwater_others         = 0.47 * 1.15
)

N_dx_map <- c(
  freshwater_plants         = 0,
  freshwater_vertebrates    = 0,
  freshwater_invertebrates  = 0,
  freshwater_others         = 0,
  terrestrial_plants        = 0,
  terrestrial_vertebrates   = 0,
  terrestrial_invertebrates = 0,
  terrestrial_others        = 0
)

## Fetch silhouettes 
img_keys <- unique(unname(uuid_map))
img_list <- setNames(lapply(img_keys, rphylopic::get_phylopic), img_keys)


get_cell <- function(df, row, col, default = NA_real_) {
  match_idx <- which(tolower(df$stratum) == row)
  cn <- colnames(df)
  if (length(match_idx) == 1 && col %in% cn) {
    return(df[[col]][match_idx])
  } else {
    return(default)
  }
}

prepare_realm_df <- function(df_raw) {
  df <- df_raw %>%
    mutate(realm = tolower(stratum)) %>%
    rename(p = pval)
  
  df$p <- suppressWarnings(as.numeric(df$p))
  
  df %>%
    mutate(
      p_label = ifelse(
        is.na(p), "",
        ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p)))
      )
    )
}

finalize_taxo_df <- function(df) {
  realms_present <- unique(df$realm)
  
  scaffold <- expand.grid(
    realm       = realms_present,
    taxon_group = taxon_order,
    stringsAsFactors = FALSE
  )
  
  df <- scaffold %>%
    left_join(df, by = c("realm", "taxon_group")) %>%
    mutate(
      group_key      = paste0(realm, "_", taxon_group),
      uuid           = uuid_map[group_key],
      icon_dx        = unname(icon_dx_map[group_key]),
      N_dx           = unname(N_dx_map[group_key]),
      order_in_realm = unname(order_map[group_key]),
      plot_row      = ifelse(is.na(plot_row), FALSE, plot_row),
      k             = ifelse(is.na(k), 0, k),
      n_studies     = ifelse(is.na(n_studies), 0, n_studies),
      taxon_group   = factor(taxon_group, levels = taxon_order)
    ) %>%
    arrange(realm, taxon_group) %>%
    group_by(realm) %>%
    mutate(
      label_text = paste0(as.character(taxon_group), " (", k, ")"),
      label = factor(label_text, levels = rev(label_text))
    ) %>%
    ungroup()
  
  df
}

prepare_richness_taxo <- function(df_raw) {
  taxo_names <- strsplit(df_raw$stratum, "_")
  
  df <- df_raw %>%
    mutate(
      realm = tolower(sapply(taxo_names, `[`, 1)),
      taxon = tolower(sapply(taxo_names, function(x) paste(x[-1], collapse = "_")))
    ) %>%
    rename(p = pval)
  
  df$p <- suppressWarnings(as.numeric(df$p))
  
  df <- df %>%
    mutate(
      taxon_group = case_when(
        taxon == "plants"        ~ "plants",
        taxon == "vertebrates"   ~ "vertebrates",
        taxon == "invertebrates" ~ "invertebrates",
        TRUE                     ~ "others"
      ),
      plot_row  = is.finite(log_RR) & is.finite(ci_lower) & is.finite(ci_upper),
      p_label = ifelse(
        is.na(p), "",
        ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p)))
      ),
      N_label = paste0("(", k, ")")
    )
  
  finalize_taxo_df(df)
}

prepare_abundance_taxo <- function(df_raw) {
  df <- data.frame(
    group_key = wanted_groups,
    stringsAsFactors = FALSE
  )
  
  df$k <- sapply(df$group_key, function(g) {
    val <- get_cell(df_raw, g, "k", default = NA)
    if (is.na(val)) 0 else as.numeric(val)
  })
  
  df$n_studies <- sapply(df$group_key, function(g) {
    val <- get_cell(df_raw, g, "n_studies", default = NA)
    if (is.na(val)) 0 else as.numeric(val)
  })
  
  df$log_RR <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "log_RR", default = NA))
  })
  
  df$ci_lower <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "ci_lower", default = NA))
  })
  
  df$ci_upper <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "ci_upper", default = NA))
  })
  
  df$p <- sapply(df$group_key, function(g) {
    as.numeric(get_cell(df_raw, g, "pval", default = NA))
  })
  
  split_keys <- strsplit(df$group_key, "_")
  df$realm <- sapply(split_keys, `[`, 1)
  df$taxon_group <- sapply(split_keys, function(x) paste(x[-1], collapse = "_"))
  
  df <- df %>%
    mutate(
      plot_row = is.finite(log_RR) &
        is.finite(ci_lower) &
        is.finite(ci_upper),
      p_label = ifelse(
        !plot_row, "",
        ifelse(is.na(p), "",
               ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p))))
      ),
      N_label = paste0("(", k, ")")
    ) %>%
    select(-group_key)
  
  finalize_taxo_df(df)
}


auto_x_limits <- function(df_taxo, df_realm, metric_name,
                          fallback_limits = c(-3, 3)) {
  vals <- c(
    df_taxo$ci_lower[df_taxo$plot_row], df_taxo$ci_upper[df_taxo$plot_row],
    df_taxo$log_RR[df_taxo$plot_row],
    df_realm$ci_lower, df_realm$ci_upper, df_realm$log_RR
  )
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) return(fallback_limits)
  
  max_abs <- max(abs(vals))
  
  margin_mult <- if (metric_name == "richness") 1.30 else 1.20
  
  bound <- max_abs * margin_mult
  
  step <- if (bound <= 3) 0.5 else if (bound <= 6) 1 else if (bound <= 12) 2 else 5
  bound <- ceiling(bound / step) * step
  
  c(-bound, bound)
}

auto_x_breaks <- function(x_limits) {
  bound <- x_limits[2]
  step <- if (bound <= 3) 1 else if (bound <= 6) 2 else if (bound <= 12) 3 else 5
  brks <- seq(-bound, bound, by = step)
  brks[abs(brks) <= bound]
}

build_core_plot <- function(
    df_taxo,
    df_realm,
    x_limits,
    x_breaks,
    metric_name,
    size_limits,
    size_breaks,
    show_x_axis = TRUE) {
  
  panel_width <- diff(x_limits)
  
  if (metric_name == "richness") {
    x_left <- x_limits[1] - 0.17 * panel_width
    left_plot_margin  <- richness_left_margin
    right_plot_margin <- 10
    x_N <- x_limits[1] - 0.06 * panel_width
  } else {
    x_left <- x_limits[1] - 0.07 * panel_width
    left_plot_margin  <- 40
    right_plot_margin <- 10
    x_N <- x_limits[1] - 0.06 * panel_width
  }
  
  ## Percent-change labels for significant taxonomic-group effects only
  pct_reliable_lnRR_cap <- 3
  
  pct_lab_df <- df_taxo %>%
    filter(
      plot_row, !is.na(p), p < 0.05,
      abs(log_RR) <= pct_reliable_lnRR_cap,
      abs(ci_lower) <= pct_reliable_lnRR_cap,
      abs(ci_upper) <= pct_reliable_lnRR_cap
    ) %>%
    mutate(
      pct_mean = (exp(log_RR)   - 1) * 100,
      pct_lo   = (exp(ci_lower) - 1) * 100,
      pct_hi   = (exp(ci_upper) - 1) * 100,
      pct_label = sprintf(
        "%s%.0f%% [%.0f, %.0f]%%",
        ifelse(pct_mean >= 0, "+", ""),
        pct_mean,
        pct_lo,
        pct_hi
      ),
      x_whisker = 0 + 0.02 * panel_width
    )
  
  bar_df <- df_taxo %>%
    filter(plot_row) %>%
    mutate(
      trunc_lo = ci_lower < x_limits[1],
      trunc_hi = ci_upper > x_limits[2],
      ci_lower_plot = pmax(ci_lower, x_limits[1]),
      ci_upper_plot = pmin(ci_upper, x_limits[2])
    )
  
  seg_line_solid  <- bar_df %>% filter(!trunc_lo & !trunc_hi)
  seg_line_dashed <- bar_df %>% filter(trunc_lo | trunc_hi)
  seg_lo_arrow    <- bar_df %>% filter(trunc_lo)
  seg_hi_arrow    <- bar_df %>% filter(trunc_hi)
  
  p <- ggplot(df_taxo, aes(y = label)) +
    geom_rect(
      data = df_realm,
      aes(xmin = ci_lower, xmax = ci_upper, fill = realm),
      ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE,
      alpha = 0.05,
      show.legend = FALSE
    ) +
    geom_vline(
      data = df_realm,
      aes(xintercept = ci_lower, colour = realm),
      linetype = "22", linewidth = 1.1, show.legend = FALSE
    ) +
    geom_vline(
      data = df_realm,
      aes(xintercept = ci_upper, colour = realm),
      linetype = "22", linewidth = 1.1, show.legend = FALSE
    ) +
    geom_vline(
      data = df_realm,
      aes(xintercept = log_RR, colour = realm),
      linewidth = 1.4, show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black", linewidth = 1.0) +
    geom_segment(
      data = seg_line_solid,
      aes(x = ci_lower_plot, xend = ci_upper_plot, y = label, yend = label),
      linewidth = 1, colour = "grey25", inherit.aes = FALSE
    ) +
    geom_segment(
      data = seg_line_dashed,
      aes(x = ci_lower_plot, xend = ci_upper_plot, y = label, yend = label),
      linewidth = 1, colour = "grey25", linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = seg_lo_arrow,
      aes(x = ci_lower_plot, xend = ci_lower_plot, y = label, yend = label),
      arrow = arrow(length = unit(0.09, "in"), type = "open", ends = "first"),
      linewidth = 1, colour = "grey25", inherit.aes = FALSE
    ) +
    geom_segment(
      data = seg_hi_arrow,
      aes(x = ci_upper_plot, xend = ci_upper_plot + 0.001 * diff(x_limits),
          y = label, yend = label),
      arrow = arrow(length = unit(0.09, "in"), type = "open", ends = "last"),
      linewidth = 1, colour = "grey25", inherit.aes = FALSE
    ) +
    geom_point(
      data = bar_df,
      aes(x = log_RR, size = n_studies),
      colour = "grey25",
      stroke = 0.6
    ) +
    scale_size_area(
      max_size = 7,
      limits   = size_limits,
      breaks   = size_breaks,
      guide    = "none"
    ) +
    geom_text(
      data = pct_lab_df,
      aes(x = x_whisker, y = label, label = pct_label),
      hjust = -0.05,
      size = 4.0,
      fontface = "plain",
      colour = "black",
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = df_taxo,
      aes(
        x = x_N,
        y = label,
        label = N_label
      ),
      hjust = 1,
      vjust = 0.5,
      size = 4.2,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = cols, name = NULL, breaks = names(cols), guide = "none") +
    scale_fill_manual(values = cols, name = NULL, breaks = names(cols), guide = "none") +
    labs(
      x = if (metric_name == "richness") {
        "Change in species richness (lnRR)"
      } else {
        "Change in species abundance (lnRR)"
      },
      y = NULL
    ) +
    theme_bw(base_size = 20) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
      strip.background   = element_blank(),
      strip.text.y.left  = element_blank(),
      strip.placement    = "outside",
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.title.x = element_text(size = 19, margin = margin(t = 10)),
      axis.text.x        = element_text(size = 16),
      legend.position    = "none",
      plot.margin        = margin(8, right_plot_margin, 14, left_plot_margin)
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_discrete(drop = FALSE)
  
  if (!show_x_axis) {
    p <- p +
      theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  if (metric_name == "richness") {
    rows_icon <- df_taxo %>%
      distinct(label, uuid, group_key, icon_dx) %>%
      filter(!is.na(uuid)) %>%
      mutate(icon_height = unname(icon_height_map[group_key]))
    
    for (i in seq_len(nrow(rows_icon))) {
      p <- p + rphylopic::add_phylopic(
        img    = img_list[[rows_icon$uuid[i]]],
        x      = x_left + rows_icon$icon_dx[i],
        y      = rows_icon$label[i],
        height = rows_icon$icon_height[i]
      )
    }
  }
  
  p
}


robustness_paths <- function(outcome, re_label, sensitivity, run_dir = RUN_DIR) {
  check_dir <- switch(sensitivity,
                      geary    = "geary",
                      outliers = "outliers",
                      stop("sensitivity must be 'geary' or 'outliers'"))
  file_name <- switch(sensitivity,
                      geary    = "summary_geary.csv",
                      outliers = "summary_no_outliers.csv")
  
  base <- file.path("./results/robustness_checks", outcome, re_label, run_dir)
  
  list(
    taxo  = file.path(base, "taxo_group_realm", check_dir, file_name),
    realm = file.path(base, "realm",            check_dir, file_name)
  )
}


## Build one combined 2x2 SI figure 
build_SI_figure <- function(sensitivity = c("geary", "outliers"),
                            out_name) {
  
  sensitivity <- match.arg(sensitivity)
  
  rich_paths <- robustness_paths("richness",  RE_LABEL_RICH, sensitivity)
  abun_paths <- robustness_paths("abundance", RE_LABEL_ABUN, sensitivity)
  
  rich_taxo_raw  <- read.csv(rich_paths$taxo,  check.names = FALSE)
  rich_realm_raw <- read.csv(rich_paths$realm, check.names = FALSE)
  abun_taxo_raw  <- read.csv(abun_paths$taxo,  check.names = FALSE)
  abun_realm_raw <- read.csv(abun_paths$realm, check.names = FALSE)
  
  rich_taxo_df  <- prepare_richness_taxo(rich_taxo_raw)
  rich_realm_df <- prepare_realm_df(rich_realm_raw)
  abun_taxo_df  <- prepare_abundance_taxo(abun_taxo_raw)
  abun_realm_df <- prepare_realm_df(abun_realm_raw)
  
  x_limits_rich <- auto_x_limits(rich_taxo_df, rich_realm_df, "richness")
  x_breaks_rich <- auto_x_breaks(x_limits_rich)
  
  x_limits_abun <- switch(sensitivity,
                          outliers = c(-5.4, 5.4),
                          geary    = c(-6, 6))
  x_breaks_abun <- switch(sensitivity,
                          outliers = c(-5, -2.5, 0, 2.5, 5),
                          geary    = auto_x_breaks(x_limits_abun))
  
  rich_fw  <- filter(rich_taxo_df,  realm == "freshwater")  %>% mutate(label = fct_drop(label))
  rich_tr  <- filter(rich_taxo_df,  realm == "terrestrial") %>% mutate(label = fct_drop(label))
  rrealm_fw <- filter(rich_realm_df, realm == "freshwater")
  rrealm_tr <- filter(rich_realm_df, realm == "terrestrial")
  
  abun_fw  <- filter(abun_taxo_df,  realm == "freshwater")  %>% mutate(label = fct_drop(label))
  abun_tr  <- filter(abun_taxo_df,  realm == "terrestrial") %>% mutate(label = fct_drop(label))
  arealm_fw <- filter(abun_realm_df, realm == "freshwater")
  arealm_tr <- filter(abun_realm_df, realm == "terrestrial")
  
  all_plotted_n <- c(
    rich_fw$n_studies[rich_fw$plot_row],
    rich_tr$n_studies[rich_tr$plot_row],
    abun_fw$n_studies[abun_fw$plot_row],
    abun_tr$n_studies[abun_tr$plot_row]
  )
  
  all_plotted_n <- all_plotted_n[is.finite(all_plotted_n)]
  
  if (length(all_plotted_n) == 0) {
    stop("No finite n_studies values were found for plotted taxonomic groups.")
  }
  
  size_limits_global <- c(0, max(all_plotted_n))
  size_breaks_global <- pretty(all_plotted_n, n = 4)
  size_breaks_global <- size_breaks_global[
    size_breaks_global >= size_limits_global[1] &
      size_breaks_global <= size_limits_global[2]
  ]
  
  p_rich_fw <- build_core_plot(
    rich_fw, rrealm_fw, x_limits_rich, x_breaks_rich, "richness",
    size_limits_global, size_breaks_global, show_x_axis = FALSE
  )
  p_rich_tr <- build_core_plot(
    rich_tr, rrealm_tr, x_limits_rich, x_breaks_rich, "richness",
    size_limits_global, size_breaks_global, show_x_axis = TRUE
  )
  p_abun_fw <- build_core_plot(
    abun_fw, arealm_fw, x_limits_abun, x_breaks_abun, "abundance",
    size_limits_global, size_breaks_global, show_x_axis = FALSE
  )
  p_abun_tr <- build_core_plot(
    abun_tr, arealm_tr, x_limits_abun, x_breaks_abun, "abundance",
    size_limits_global, size_breaks_global, show_x_axis = TRUE
  )
  
  p_rich_fw <- p_rich_fw +
    labs(title = "<b>&#40;a&#41;</b> Freshwater species richness") +
    theme(plot.title = ggtext::element_markdown(
      hjust = 0, size = 19, colour = "grey20", margin = margin(b = 7)
    ))
  
  p_abun_fw <- p_abun_fw +
    labs(title = "<b>&#40;b&#41;</b> Freshwater species abundance") +
    theme(plot.title = ggtext::element_markdown(
      hjust = 0, size = 19, colour = "grey20", margin = margin(b = 7)
    ))
  
  p_rich_tr <- p_rich_tr +
    labs(title = "<b>&#40;c&#41;</b> Terrestrial species richness") +
    theme(plot.title = ggtext::element_markdown(
      hjust = 0, size = 19, colour = "grey20", margin = margin(b = 7)
    ))
  
  p_abun_tr <- p_abun_tr +
    labs(title = "<b>&#40;d&#41;</b> Terrestrial species abundance") +
    theme(plot.title = ggtext::element_markdown(
      hjust = 0, size = 19, colour = "grey20", margin = margin(b = 7)
    ))
  
  p_combined <- (p_rich_fw | p_abun_fw) / (p_rich_tr | p_abun_tr)
  
  out_png <- file.path(dir_mod, paste0(out_name, ".png"))
  out_tif <- file.path(dir_mod, paste0(out_name, ".tif"))
  
  ggsave(out_png, p_combined, width = 360, height = 250, units = "mm",
         bg = "white", dpi = 1000)
  ggsave(out_tif, p_combined, width = 360, height = 250, units = "mm",
         bg = "white", dpi = 1000, compression = "lzw")
  
  invisible(p_combined)
}


## Build the two SI figures

fig_outliers <- build_SI_figure(
  sensitivity = "outliers",
  out_name    = "Figure S2 (outliers)"
)


fig_geary <- build_SI_figure(
  sensitivity = "geary",
  out_name    = "Figure S3 (geary)"
)


##============================================================================##



##============================================================================##
## Figure S4 ------
##============================================================================##

dir_rich_base <- "./results/richness/Study_Mine_VCV"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV"

methods <- c("Bracken", "Median", "HotDeck", "Poisson", "unweighted")

method_display <- c(
  Bracken     = "Bracken",
  Median      = "Median",
  HotDeck     = "HotDeck",
  Poisson     = "Poisson",
  unweighted  = "Unweighted"
)

cols_methods <- c(
  Bracken    = "#0072B2",
  Median     = "#E69F00",
  HotDeck    = "#009E73",
  Poisson    = "#D55E00",
  Unweighted = "#CC79A7"
)

RICH_X_LIMITS <- c(-4, 4)
RICH_X_BREAKS <- c(-4, -2, 0, 2, 4)

ABUN_X_LIMITS <- c(-10, 10)
ABUN_X_BREAKS <- c(-10, -5, 0, 5, 10)

## PhyloPic UUID map (same as other figures)
uuid_map <- c(
  terrestrial_plants        = "43afe2df-ab6c-47e6-a105-c0a82b8af1c5",
  terrestrial_vertebrates   = "fe09db98-924b-4ad4-bcfc-2b660732ff9d",
  terrestrial_invertebrates = "fbb67694-0fd2-4e90-b88d-efab9cbac37c",
  terrestrial_others        = "18b8f0a9-b280-4fb9-b750-495ceec4ef87",
  freshwater_plants         = "4e396beb-1555-47ee-9997-47fc970accf6",
  freshwater_vertebrates    = "63028840-34fa-44fb-8883-41b2866e61b8",
  freshwater_invertebrates  = "2757cb60-acda-422d-a318-85b812a658c2",
  freshwater_others         = "a0753559-d626-4b3a-9af7-39c79409b8f4"
)

## Icon positions
icon_dx_map <- c(
  freshwater_plants         =  0.06 * 2,
  freshwater_vertebrates    = -0.15 * 2,
  freshwater_invertebrates  = -0.03 * 2,
  freshwater_others         = -0.01 * 2,
  terrestrial_plants        = -0.03 * 2,
  terrestrial_vertebrates   = -0.12 * 2,
  terrestrial_invertebrates =  0.06 * 2,
  terrestrial_others        =  0.01 * 2
)
icon_height_map <- c(
  terrestrial_plants        = 0.37 * 1.15,
  terrestrial_vertebrates   = 0.26 * 1.15,
  terrestrial_invertebrates = 0.44 * 1.15,
  terrestrial_others        = 0.32 * 1.15,
  freshwater_plants         = 0.46 * 1.15,
  freshwater_vertebrates    = 0.23 * 1.15,
  freshwater_invertebrates  = 0.32 * 1.15,
  freshwater_others         = 0.47 * 1.15
)

ICON_X_FRAC_RICH <- 0.17
ICON_X_FRAC_ABUN <- 0.07
N_X_FRAC         <- 0.06
LEFT_MARGIN_RICH <- 90
LEFT_MARGIN_ABUN <- 40

n_methods   <- length(methods)
dodge_span  <- 0.36
method_offsets <- setNames(
  seq(-dodge_span / 2, dodge_span / 2, length.out = n_methods),
  unname(method_display[methods])
)

tick_half_height <- 0.03   

taxon_order <- c("plants", "vertebrates", "invertebrates", "others")


read_richness_method <- function(method) {
  df <- read.csv(
    file.path(dir_rich_base, method, "taxo_group_realm",
              "richness_taxo_group_realm_summary.csv"),
    check.names = FALSE
  )
  
  taxo_names <- strsplit(df$stratum, "_")
  
  df %>%
    mutate(
      realm       = tolower(sapply(taxo_names, `[`, 1)),
      taxon       = tolower(sapply(taxo_names, function(x) paste(x[-1], collapse = "_"))),
      taxon_group = case_when(
        taxon == "plants"        ~ "plants",
        taxon == "vertebrates"   ~ "vertebrates",
        taxon == "invertebrates" ~ "invertebrates",
        TRUE                     ~ "others"
      ),
      group_key = paste0(realm, "_", taxon_group),
      uuid      = unname(uuid_map[group_key]),
      method    = unname(method_display[method])
    ) %>%
    filter(is.finite(log_RR), is.finite(ci_lower), is.finite(ci_upper)) %>%
    select(realm, taxon, taxon_group, group_key, uuid, method, k, log_RR, ci_lower, ci_upper)
}

read_abundance_method <- function(method) {
  df <- read.csv(
    file.path(dir_abun_base, method, "taxo_group_realm",
              "abundance_taxo_group_realm_summary.csv"),
    check.names = FALSE
  )
  
  taxo_names <- strsplit(df$stratum, "_")
  
  df %>%
    mutate(
      realm       = tolower(sapply(taxo_names, `[`, 1)),
      taxon       = tolower(sapply(taxo_names, function(x) paste(x[-1], collapse = "_"))),
      taxon_group = case_when(
        taxon == "plants"        ~ "plants",
        taxon == "vertebrates"   ~ "vertebrates",
        taxon == "invertebrates" ~ "invertebrates",
        TRUE                     ~ "others"
      ),
      group_key = paste0(realm, "_", taxon_group),
      uuid      = unname(uuid_map[group_key]),
      method    = unname(method_display[method]),
      tax_group_order = case_when(
        realm == "freshwater"  & taxon_group == "vertebrates"   ~ 1,
        realm == "freshwater"  & taxon_group == "invertebrates" ~ 2,
        realm == "freshwater"  & taxon_group == "others"        ~ 3,
        realm == "terrestrial" & taxon_group == "plants"        ~ 1,
        realm == "terrestrial" & taxon_group == "vertebrates"   ~ 2,
        realm == "terrestrial" & taxon_group == "invertebrates" ~ 3,
        TRUE ~ 99
      )
    ) %>%
    select(realm, taxon, taxon_group, group_key, uuid, method, k, log_RR, ci_lower, ci_upper,
           tax_group_order)
}


## Read and combine all methods
df_rich_all <- bind_rows(lapply(methods, read_richness_method))
df_abun_all <- bind_rows(lapply(methods, read_abundance_method))

rich_order <- df_rich_all %>%
  filter(method == "Bracken") %>%
  arrange(realm, log_RR) %>%
  pull(group_key) %>%
  unique()
df_rich_all$group_key <- factor(df_rich_all$group_key, levels = rev(rich_order))

canonical_abun <- expand.grid(
  realm       = c("freshwater", "terrestrial"),
  taxon_group = factor(taxon_order, levels = taxon_order),
  stringsAsFactors = FALSE
) %>%
  arrange(realm, taxon_group) %>%
  mutate(group_key = paste0(realm, "_", taxon_group))

abun_order <- canonical_abun$group_key
df_abun_all$group_key <- factor(df_abun_all$group_key, levels = rev(abun_order))

icon_df_rich <- df_rich_all %>%
  filter(method == "Bracken") %>%
  distinct(realm, group_key, uuid, k) %>%
  mutate(
    N_label     = paste0("(", k, ")"),
    icon_dx     = unname(icon_dx_map[as.character(group_key)]),
    icon_height = unname(icon_height_map[as.character(group_key)])
  )

icon_df_abun <- canonical_abun %>%
  left_join(
    df_abun_all %>% filter(method == "Bracken") %>% distinct(realm, group_key, k),
    by = c("realm", "group_key")
  ) %>%
  mutate(
    k           = ifelse(is.na(k), 0, k),
    N_label     = paste0("(", k, ")"),
    uuid        = unname(uuid_map[group_key]),
    icon_dx     = unname(icon_dx_map[group_key]),
    icon_height = unname(icon_height_map[group_key]),
    group_key   = factor(group_key, levels = rev(abun_order))
  )

img_keys <- unique(na.omit(c(icon_df_rich$uuid, icon_df_abun$uuid)))
img_list <- setNames(lapply(img_keys, rphylopic::get_phylopic), img_keys)


build_panel <- function(df, icon_df, x_limits, x_breaks, x_title,
                        show_x_axis, icon_x_frac, N_x_frac = N_X_FRAC,
                        left_margin) {
  
  panel_width <- diff(x_limits)
  x_left      <- x_limits[1] - icon_x_frac * panel_width
  x_N         <- x_limits[1] - N_x_frac * panel_width
  
  df <- df %>%
    mutate(
      y_num         = as.numeric(group_key) + unname(method_offsets[method]),
      trunc_lo      = ci_lower < x_limits[1],
      trunc_hi      = ci_upper > x_limits[2],
      ci_lower_plot = pmax(ci_lower, x_limits[1]),
      ci_upper_plot = pmin(ci_upper, x_limits[2]),
      log_RR_plot   = pmin(pmax(log_RR, x_limits[1]), x_limits[2])
    )
  
  seg_line_solid  <- df %>% filter(!trunc_lo & !trunc_hi)
  seg_line_dashed <- df %>% filter(trunc_lo | trunc_hi)
  seg_lo_flat     <- df %>% filter(!trunc_lo)
  seg_lo_arrow    <- df %>% filter(trunc_lo)
  seg_hi_flat     <- df %>% filter(!trunc_hi)
  seg_hi_arrow    <- df %>% filter(trunc_hi)
  
  icon_df <- icon_df %>%
    mutate(y_num = as.numeric(group_key))
  
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1, colour = "grey40") +
    geom_segment(
      data = seg_line_solid,
      aes(x = ci_lower_plot, xend = ci_upper_plot, y = y_num, yend = y_num, colour = method),
      linewidth = 0.9
    ) +
    geom_segment(
      data = seg_line_dashed,
      aes(x = ci_lower_plot, xend = ci_upper_plot, y = y_num, yend = y_num, colour = method),
      linewidth = 0.9, linetype = "dashed", show.legend = FALSE
    ) +
    geom_segment(
      data = seg_lo_flat,
      aes(x = ci_lower_plot, xend = ci_lower_plot,
          y = y_num - tick_half_height, yend = y_num + tick_half_height, colour = method),
      linewidth = 0.9, show.legend = FALSE
    ) +
    geom_segment(
      data = seg_hi_flat,
      aes(x = ci_upper_plot, xend = ci_upper_plot,
          y = y_num - tick_half_height, yend = y_num + tick_half_height, colour = method),
      linewidth = 0.9, show.legend = FALSE
    ) +
    geom_segment(
      data = seg_lo_arrow,
      aes(x = ci_lower_plot, xend = ci_lower_plot - 0.001 * diff(x_limits),
          y = y_num, yend = y_num, colour = method),
      arrow = arrow(length = unit(0.07, "in"), type = "open", ends = "last"),
      linewidth = 0.9, show.legend = FALSE
    ) +
    geom_segment(
      data = seg_hi_arrow,
      aes(x = ci_upper_plot, xend = ci_upper_plot + 0.001 * diff(x_limits),
          y = y_num, yend = y_num, colour = method),
      arrow = arrow(length = unit(0.07, "in"), type = "open", ends = "last"),
      linewidth = 0.9, show.legend = FALSE
    ) +
    geom_point(
      data = df,
      aes(x = log_RR_plot, y = y_num, colour = method),
      size = 1.8
    ) +
    geom_text(
      data = icon_df,
      aes(x = x_N, y = y_num, label = N_label),
      hjust = 1, vjust = 0.5, size = 3.6, inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = cols_methods, name = NULL,
      breaks = names(cols_methods), limits = names(cols_methods), drop = FALSE
    ) +
    labs(x = if (show_x_axis) x_title else NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title    = element_blank(),
      plot.margin = margin(10, 15, 10, left_margin)
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    scale_x_continuous(breaks = x_breaks)
  
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  attr(p, "x_left") <- x_left
  p
}

add_icons <- function(p, icon_df_realm) {
  rows <- icon_df_realm %>% filter(!is.na(uuid))
  x_left <- attr(p, "x_left")
  
  for (i in seq_len(nrow(rows))) {
    p <- p + rphylopic::add_phylopic(
      img    = img_list[[rows$uuid[i]]],
      x      = x_left + rows$icon_dx[i],
      y      = as.numeric(rows$group_key)[i],
      height = rows$icon_height[i]
    )
  }
  p
}


## Build the 4 panels
p_rich_fw <- build_panel(
  df_rich_all %>% filter(realm == "freshwater"),
  icon_df_rich %>% filter(realm == "freshwater"),
  RICH_X_LIMITS, RICH_X_BREAKS,
  "Change in species richness (lnRR)", show_x_axis = FALSE,
  icon_x_frac = ICON_X_FRAC_RICH, left_margin = LEFT_MARGIN_RICH
)
p_rich_tr <- build_panel(
  df_rich_all %>% filter(realm == "terrestrial"),
  icon_df_rich %>% filter(realm == "terrestrial"),
  RICH_X_LIMITS, RICH_X_BREAKS,
  "Change in species richness (lnRR)", show_x_axis = TRUE,
  icon_x_frac = ICON_X_FRAC_RICH, left_margin = LEFT_MARGIN_RICH
)

p_abun_fw <- build_panel(
  df_abun_all %>% filter(realm == "freshwater"),
  icon_df_abun %>% filter(realm == "freshwater"),
  ABUN_X_LIMITS, ABUN_X_BREAKS,
  "Change in species abundance (lnRR)", show_x_axis = FALSE,
  icon_x_frac = ICON_X_FRAC_ABUN, left_margin = LEFT_MARGIN_ABUN
)
p_abun_tr <- build_panel(
  df_abun_all %>% filter(realm == "terrestrial"),
  icon_df_abun %>% filter(realm == "terrestrial"),
  ABUN_X_LIMITS, ABUN_X_BREAKS,
  "Change in species abundance (lnRR)", show_x_axis = TRUE,
  icon_x_frac = ICON_X_FRAC_ABUN, left_margin = LEFT_MARGIN_ABUN
)

p_rich_fw <- add_icons(p_rich_fw, icon_df_rich %>% filter(realm == "freshwater"))
p_rich_tr <- add_icons(p_rich_tr, icon_df_rich %>% filter(realm == "terrestrial"))


p_combined <- (p_rich_fw + p_abun_fw + p_rich_tr + p_abun_tr + guide_area()) +
  plot_layout(
    design = "
    AB
    CD
    EE
    ",
    heights = c(1, 1, 0.12),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    legend.key.size = unit(1.4, "lines")
  )


## Save
ggsave(
  "./figures/Figure S4.png",
  p_combined, width = 320, height = 260, units = "mm", dpi = 600, bg = "white"
)

ggsave(
  "./figures/Figure S4.tif",
  p_combined, width = 320, height = 260, units = "mm", dpi = 600, bg = "white",
  compression = "lzw"
)
##============================================================================##


##============================================================================##
## Figure S5 ----
##============================================================================##
col_terrestrial <- "#009E73"
col_freshwater  <- "#0072B2"

format_pval <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) "P < 0.001" else paste0("P = ", sprintf("%.3f", p))
}


get_y_lim <- function(funnel_data, tick_max) {
  
  max_abs_resid <- max(
    abs(funnel_data$resid[is.finite(funnel_data$resid)]),
    na.rm = TRUE
  )
  
  y_lim <- max(
    tick_max * 1.02,
    max_abs_resid * 1.03
  )
  
  y_lim
}


build_funnel_panel <- function(funnel_data, egger_summary, realm_val,
                               colour_val, outcome_label, panel_label,
                               y_lim, tick_max,
                               show_y_title = TRUE,
                               show_x_axis = TRUE) {
  
  fd <- filter(funnel_data, stratum == realm_val)
  er <- filter(egger_summary, stratum == realm_val)
  
  slope_txt <- if (nrow(er) == 1) {
    paste0(
      "Slope = ", sprintf("%.2f", er$estimate),
      " (95% CI: ", sprintf("%.2f", er$ci_lower), ", ",
      sprintf("%.2f", er$ci_upper), ")"
    )
  } else {
    "Small-study test unavailable"
  }
  
  p_txt <- if (nrow(er) == 1) {
    format_pval(er$pval)
  } else {
    ""
  }
  # ensure that (c) doesn't turn into the copyright symbol
  panel_title <- paste0(
    "<b>",
    sub("\\(", "&#40;", sub("\\)", "&#41;", panel_label)),
    "</b> ",
    outcome_label, " (", realm_val, ")"
  )
  
  ## annotation position
  x_range <- range(fd$effective_n, na.rm = TRUE)
  x_ann   <- x_range[2] 
  
  y_slope <- y_lim * 0.90
  y_pval  <- y_lim * 0.81
  
  p <- ggplot(fd, aes(x = effective_n, y = resid)) +
    
    geom_point(
      alpha = 0.6,
      size = 2,
      colour = colour_val
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    
    annotate(
      "text",
      x = x_ann,
      y = y_slope,
      label = slope_txt,
      hjust = 1,
      vjust = 1,
      size = 4.2
    ) +
    
    annotate(
      "text",
      x = x_ann,
      y = y_pval,
      label = p_txt,
      hjust = 1,
      vjust = 1,
      size = 4.2
    ) +
    
    scale_y_continuous(
      breaks = seq(-tick_max, tick_max, by = 5),
      limits = c(-y_lim, y_lim),
      expand = expansion(mult = 0)
    ) +
    
    labs(
      title = panel_title,
      x = if (show_x_axis) "Effective sample size" else NULL,
      y = if (show_y_title) "Meta-analytic residuals" else NULL
    ) +
    
    theme_minimal(base_size = 15) +
    
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      ),
      
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      
      plot.title = ggtext::element_markdown(
        hjust = 0,
        size = 15,
        margin = margin(b = 6)
      ),
      
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!show_x_axis) {
    p <- p +
      theme(
        axis.title.x = element_blank()
      )
  }
  
  p
}


make_pair_funnel_figure <- function(funnel_data, egger_summary, out_path,
                                    outcome_label, tick_max,
                                    panel_labels = c("(a)", "(b)")) {
  
  y_lim <- get_y_lim(funnel_data, tick_max)
  
  p_terr <- build_funnel_panel(
    funnel_data   = funnel_data,
    egger_summary = egger_summary,
    realm_val     = "terrestrial",
    colour_val    = col_terrestrial,
    outcome_label = outcome_label,
    panel_label   = panel_labels[1],
    y_lim         = y_lim,
    tick_max      = tick_max,
    show_y_title  = TRUE,
    show_x_axis   = TRUE
  )
  
  p_fresh <- build_funnel_panel(
    funnel_data   = funnel_data,
    egger_summary = egger_summary,
    realm_val     = "freshwater",
    colour_val    = col_freshwater,
    outcome_label = outcome_label,
    panel_label   = panel_labels[2],
    y_lim         = y_lim,
    tick_max      = tick_max,
    show_y_title  = FALSE,
    show_x_axis   = TRUE
  )
  
  combined <- p_terr + p_fresh + patchwork::plot_layout(ncol = 2)
  
  ggsave(
    out_path,
    combined,
    width = 11,
    height = 5,
    dpi = 600
  )
  
  combined
}


make_four_panel_funnel_figure <- function(funnel_data_r, egger_summary_r,
                                          funnel_data_a, egger_summary_a,
                                          out_path) {
  
  y_lim_r <- get_y_lim(funnel_data_r, tick_max = 10)
  y_lim_a <- get_y_lim(funnel_data_a, tick_max = 15)
  
  p_a <- build_funnel_panel(
    funnel_data   = funnel_data_r,
    egger_summary = egger_summary_r,
    realm_val     = "terrestrial",
    colour_val    = col_terrestrial,
    outcome_label = "Species richness",
    panel_label   = "(a)",
    y_lim         = y_lim_r,
    tick_max      = 10,
    show_y_title  = TRUE,
    show_x_axis   = FALSE
  )
  
  p_b <- build_funnel_panel(
    funnel_data   = funnel_data_r,
    egger_summary = egger_summary_r,
    realm_val     = "freshwater",
    colour_val    = col_freshwater,
    outcome_label = "Species richness",
    panel_label   = "(b)",
    y_lim         = y_lim_r,
    tick_max      = 10,
    show_y_title  = FALSE,
    show_x_axis   = FALSE
  )
  
  p_c <- build_funnel_panel(
    funnel_data   = funnel_data_a,
    egger_summary = egger_summary_a,
    realm_val     = "terrestrial",
    colour_val    = col_terrestrial,
    outcome_label = "Species abundance",
    panel_label   = "(c)",
    y_lim         = y_lim_a,
    tick_max      = 15,
    show_y_title  = TRUE,
    show_x_axis   = TRUE
  )
  
  p_d <- build_funnel_panel(
    funnel_data   = funnel_data_a,
    egger_summary = egger_summary_a,
    realm_val     = "freshwater",
    colour_val    = col_freshwater,
    outcome_label = "Species abundance",
    panel_label   = "(d)",
    y_lim         = y_lim_a,
    tick_max      = 15,
    show_y_title  = FALSE,
    show_x_axis   = TRUE
  )
  
  combined_all <- (p_a + p_b) / (p_c + p_d)
  
  ggsave(
    out_path,
    combined_all,
    width = 11,
    height = 10,
    dpi = 600
  )
  
  combined_all
}



## Read data 
out_dir_r <- file.path(
  "results", "richness", "Study_Mine_VCV", "Bracken", "realm"
)

funnel_data_r <- read_csv(
  file.path(out_dir_r, "richness_realm_funnel_data.csv"),
  show_col_types = FALSE
)

egger_summary_r <- read_csv(
  file.path(out_dir_r, "richness_realm_egger_summary.csv"),
  show_col_types = FALSE
)


out_dir_a <- file.path(
  "results", "abundance", "Study_Mine_Binomial_VCV", "Bracken", "realm"
)

funnel_data_a <- read_csv(
  file.path(out_dir_a, "abundance_realm_funnel_data.csv"),
  show_col_types = FALSE
)

egger_summary_a <- read_csv(
  file.path(out_dir_a, "abundance_realm_egger_summary.csv"),
  show_col_types = FALSE
)


## Four-panel figure
make_four_panel_funnel_figure(
  funnel_data_r   = funnel_data_r,
  egger_summary_r = egger_summary_r,
  funnel_data_a   = funnel_data_a,
  egger_summary_a = egger_summary_a,
  out_path        = "./figures/Figure S5.png"
)


make_four_panel_funnel_figure(
  funnel_data_r   = funnel_data_r,
  egger_summary_r = egger_summary_r,
  funnel_data_a   = funnel_data_a,
  egger_summary_a = egger_summary_a,
  out_path        = "./figures/Figure S5.tif"
)

##============================================================================##




##============================================================================##
## Figures S6-S7 -----
##============================================================================##
dir_rich_base <- "./results/richness/Study_Mine_VCV/Bracken/material_group_moderator"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/material_group_moderator"

rich_means_file <- file.path(dir_rich_base, "richness_material_moderator_means_all.csv")
abun_means_file <- file.path(dir_abun_base, "abundance_material_moderator_means_all.csv")

material_order <- c("Coal", "Gold", "Metals", "Aggregates")

canonical_group_order <- c(
  "Freshwater - Plants", "Freshwater - Vertebrates",
  "Freshwater - Invertebrates", "Freshwater - Microorganisms and algae",
  "Terrestrial - Plants", "Terrestrial - Vertebrates",
  "Terrestrial - Invertebrates", "Terrestrial - Microorganisms and algae"
)

group_cols <- c(
  "Freshwater - Invertebrates"              = "#08306B",
  "Freshwater - Vertebrates"                = "#2171B5",
  "Freshwater - Microorganisms and algae"   = "#6BAED6",
  "Freshwater - Plants"                     = "#BDD7E7",
  "Terrestrial - Invertebrates"             = "#00441B",
  "Terrestrial - Vertebrates"               = "#238B45",
  "Terrestrial - Plants"                    = "#66C2A4",
  "Terrestrial - Microorganisms and algae"  = "#C7E9C0"
)

pretty_realm <- function(x) {
  case_when(
    str_to_lower(x) == "freshwater"  ~ "Freshwater",
    str_to_lower(x) == "terrestrial" ~ "Terrestrial",
    TRUE ~ str_to_title(x)
  )
}

pretty_material <- function(x) {
  x <- str_to_lower(x)
  case_when(
    x == "coal"       ~ "Coal",
    x == "gold"       ~ "Gold",
    x == "metals"     ~ "Metals",
    x == "aggregates" ~ "Aggregates",
    TRUE              ~ str_to_title(x)
  )
}

pretty_taxo <- function(x) {
  x <- str_to_lower(x)
  case_when(
    x == "vertebrates"   ~ "Vertebrates",
    x == "invertebrates" ~ "Invertebrates",
    x == "plants"        ~ "Plants",
    x == "others"        ~ "Microorganisms and algae",
    TRUE                 ~ str_to_title(x)
  )
}

sig_stars <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 1e-4   ~ "****",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ ""
  )
}

read_taxo_material <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(stratum_col == "taxo_group_realm") %>%
    mutate(stratum = str_to_lower(stratum)) %>%
    tidyr::separate(
      stratum, into = c("realm", "taxo_group"), sep = "_",
      extra = "merge", remove = FALSE
    ) %>%
    mutate(
      realm_label    = pretty_realm(realm),
      taxo_label     = pretty_taxo(taxo_group),
      material_group = str_to_lower(material_group),
      material_label = pretty_material(material_group),
      modelled       = is.finite(log_RR),
      star_label     = sig_stars(pval),
      k_label        = paste0("(k=", k, ")"),
      group_label    = paste(realm_label, taxo_label, sep = " - ")
    ) %>%
    filter(material_group != "", modelled)
}

build_commodity_blocks <- function(df, legend_levels, row_spacing = 1, inner_pad = 0.55) {
  
  df <- df %>%
    mutate(
      material_label = factor(material_label, levels = material_order),
      group_label    = factor(group_label, levels = legend_levels),
      commodity_order = match(as.character(material_label), material_order),
      group_order     = match(as.character(group_label), legend_levels)
    ) %>%
    filter(!is.na(commodity_order), !is.na(group_order)) %>%
    arrange(commodity_order, group_order) %>%
    group_by(material_label) %>%
    arrange(group_order, .by_group = TRUE) %>%
    mutate(row_in_block = row_number(), n_block = n()) %>%
    ungroup()
  
  y_vals      <- numeric(nrow(df))
  block_info  <- list()
  current_top <- 0
  
  for (mat in unique(df$material_label)) {
    block_rows <- which(df$material_label == mat)
    n_block    <- length(block_rows)
    
    block_height <- if (n_block == 1) 2 * inner_pad else (n_block - 1) * row_spacing + 2 * inner_pad
    y_top    <- current_top
    y_bottom <- current_top - block_height
    y_mid    <- (y_top + y_bottom) / 2
    
    y_vals[block_rows] <- if (n_block == 1) {
      y_mid
    } else {
      seq(from = y_top - inner_pad, to = y_bottom + inner_pad, length.out = n_block)
    }
    
    block_info[[as.character(mat)]] <- data.frame(
      material_label = mat, y_top = y_top, y_bottom = y_bottom, y_mid = y_mid
    )
    current_top <- y_bottom
  }
  
  df$y      <- y_vals
  df$k_y    <- df$y
  df$star_y <- df$y
  
  list(df = df, blocks = bind_rows(block_info))
}

format_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "p < 0.01", paste0("p = ", sprintf("%.2f", p))))
}


build_taxo_commodity_plot <- function(df_raw, x_title, label_gap_frac = 0.02,
                                      legend_text_size = 8, legend_key_size = 0.8,
                                      left_margin = 85) {
  
  legend_levels <- canonical_group_order[canonical_group_order %in% unique(df_raw$group_label)]
  
  built   <- build_commodity_blocks(df_raw, legend_levels)
  df      <- built$df
  blocks  <- built$blocks
  
  max_abs <- max(abs(c(df$ci_lower, df$ci_upper)), na.rm = TRUE)
  max_abs <- ceiling(max_abs * 1.15 * 2) / 2
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  panel_width <- 2 * max_abs
  

  df <- df %>%
    mutate(
      text_x     = ci_upper + 0.02 * max_abs,
      k_digits   = nchar(as.character(k)),
      star_x     = text_x + 0.09 * max_abs + pmax(k_digits - 1, 0) * 0.045 * max_abs
    )
  
  label_x <- -max_abs * (1 + label_gap_frac)
  
  ggplot(df, aes(x = log_RR, y = y, colour = group_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey45") +
    geom_hline(data = blocks %>% slice(1), aes(yintercept = y_top),
               inherit.aes = FALSE, colour = "black", linewidth = 0.35) +
    geom_hline(data = blocks, aes(yintercept = y_bottom),
               inherit.aes = FALSE, colour = "black", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0, linewidth = 0.9) +
    geom_point(size = 3) +
    geom_text(
      aes(x = text_x, y = k_y, label = k_label),
      colour = "black", hjust = 0, size = 2.8, show.legend = FALSE
    ) +
    geom_text(
      aes(x = star_x, y = star_y, label = star_label),
      colour = "black", fontface = "bold", hjust = 0, size = 4.0, show.legend = FALSE
    ) +
    geom_text(
      data = blocks, aes(x = label_x, y = y_mid, label = material_label),
      inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = 3.4
    ) +
    scale_colour_manual(
      values = group_cols, breaks = legend_levels, name = NULL,
      guide = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    scale_y_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(-max_abs, max_abs), clip = "off") +
    labs(x = x_title, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.title.x      = element_text(size = 10),
      legend.position   = "bottom",
      legend.text       = element_text(size = legend_text_size),
      legend.key.size   = unit(legend_key_size, "lines"),
      legend.margin     = margin(t = 2, b = 0),
      legend.box.margin = margin(t = 4, b = 0),
      plot.margin = margin(6, 20, 6, left_margin)
    )
}


## Read data
si_rich_raw <- read_taxo_material(rich_means_file)
si_abun_raw <- read_taxo_material(abun_means_file)


## Build plots
p_S1     <- build_taxo_commodity_plot(si_rich_raw, "Change in species richness (lnRR)",
                                      legend_text_size = 8, legend_key_size = 0.8)
p_S_abun <- build_taxo_commodity_plot(si_abun_raw, "Change in species abundance (lnRR)",
                                      legend_text_size = 9, legend_key_size = 0.9)


## Save
ggsave("./figures/Figure S6 (richness).png",
       p_S1, width = 200, height = 148, units = "mm", dpi = 600, bg = "white")
ggsave("./figures/Figure S6 (richness).tif",
       p_S1, width = 200, height = 148, units = "mm", dpi = 600, bg = "white", compression = "lzw")

ggsave("./figures/Figure S7 (abundance).png",
       p_S_abun, width = 200, height = 148, units = "mm", dpi = 600, bg = "white")
ggsave("./figures/Figure S7 (abundance).tif",
       p_S_abun, width = 200, height = 148, units = "mm", dpi = 600, bg = "white", compression = "lzw")

##============================================================================##




##============================================================================##
## Figure S9 -------------
##============================================================================##
dir_rich_base <- "./results/richness/Study_Mine_VCV/Bracken/realm/richness_realm_models.rds"
dir_abun_base <- "./results/abundance/Study_Mine_Binomial_VCV/Bracken/realm/abundance_realm_models.rds"
mines_path    <- "data/Mines.xlsx"

get_model_by_realm <- function(model_list, realm_name) {
  nms <- tolower(names(model_list))
  idx <- match(tolower(realm_name), nms)
  
  if (is.na(idx)) {
    stop(
      "Could not find realm '", realm_name, "' in model list. Available names are: ",
      paste(names(model_list), collapse = ", ")
    )
  }
  
  model_list[[idx]]
}

extract_sampling_variance <- function(m) {
  if (!is.null(m$vi)) {
    return(as.numeric(m$vi))
  }
  
  if (!is.null(m$V)) {
    V <- m$V
    
    if (is.matrix(V) || inherits(V, "Matrix")) {
      return(as.numeric(diag(as.matrix(V))))
    } else {
      return(as.numeric(V))
    }
  }
}


extract_effects_with_distance <- function(m, realm_name, metric_name, mines_df) {
  yi <- as.numeric(m$yi)
  vi <- extract_sampling_variance(m)
  
  df_used <- attr(m, "._df_used")
  df_used <- as.data.frame(df_used)
  
  tibble(
    metric               = metric_name,
    realm                = tolower(realm_name),
    RR_log               = yi,
    SamplingVariance     = vi,
    distance_disturbed_in_meters = as.numeric(
      gsub("[^0-9.]", "", df_used$distance_disturbed_in_meters)
    )
  ) %>%
    filter(
      is.finite(RR_log),
      is.finite(SamplingVariance), SamplingVariance >= 0,
      is.finite(distance_disturbed_in_meters)
    ) %>%
    mutate(
      distance_km = distance_disturbed_in_meters / 1000,
      weight      = 1 / SamplingVariance
    )
}

rich_models <- readRDS(dir_rich_base)
abun_models <- readRDS(dir_abun_base)

mines_df <- readxl::read_excel(mines_path)

rich_terrestrial <- get_model_by_realm(rich_models, "terrestrial")
rich_freshwater  <- get_model_by_realm(rich_models, "freshwater")

abun_terrestrial <- get_model_by_realm(abun_models, "terrestrial")
abun_freshwater  <- get_model_by_realm(abun_models, "freshwater")

dist_df <- bind_rows(
  extract_effects_with_distance(rich_terrestrial, "terrestrial", "Richness",  mines_df),
  extract_effects_with_distance(rich_freshwater,  "freshwater",  "Richness",  mines_df),
  extract_effects_with_distance(abun_terrestrial, "terrestrial", "Abundance", mines_df),
  extract_effects_with_distance(abun_freshwater,  "freshwater",  "Abundance", mines_df)
)

## Quick sanity check
dist_df %>% count(metric, realm)

## Scale point size by inverse-variance weight
dist_df <- dist_df %>%
  group_by(metric) %>%
  mutate(
    wt_range = diff(range(weight, na.rm = TRUE)),
    wt_size  = ifelse(
      wt_range > 0,
      1.5 + 4 * (weight - min(weight, na.rm = TRUE)) / wt_range,
      3
    )
  ) %>%
  ungroup()

cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
dist_df$realm_lbl <- cap_first(dist_df$realm)
dist_df$metric    <- factor(dist_df$metric, levels = c("Richness", "Abundance"))

col_terrestrial <- "#009E73"
col_freshwater  <- "#0072B2"
pal_realm <- c(Terrestrial = col_terrestrial, Freshwater = col_freshwater)

base_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(2.5, "mm"),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin       = margin(10, 10, 10, 10)
  )

## Richness panel
rich_dd   <- dplyr::filter(dist_df, metric == "Richness")
rich_ymax <- max(abs(rich_dd$RR_log), na.rm = TRUE)

richness_plot <- ggplot(rich_dd, aes(x = distance_km, y = RR_log, colour = realm_lbl)) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.5) +
  geom_point(aes(size = wt_size), alpha = 0.6, show.legend = c(size = FALSE)) +
  scale_size_identity() +
  scale_colour_manual(values = pal_realm, name = NULL) +
  scale_y_continuous(limits = c(-rich_ymax, rich_ymax)) +
  labs(
    x = "Distance from mine (km)",
    y = expression(paste("Effect size (", lnRR, ")")),
    title = "<b>(a)</b> Species richness"
  ) +
  base_theme +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(
    plot.title = ggtext::element_markdown(size = 14, colour = "black"),
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = "white", colour = NA)
  )

## Abundance panel 
abun_dd   <- dplyr::filter(dist_df, metric == "Abundance")
abun_ymax <- max(abs(abun_dd$RR_log), na.rm = TRUE)

abundance_plot <- ggplot(abun_dd, aes(x = distance_km, y = RR_log, colour = realm_lbl)) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.5) +
  geom_point(aes(size = wt_size), alpha = 0.85, show.legend = FALSE) +
  scale_size_identity() +
  scale_colour_manual(values = pal_realm, name = NULL) +
  scale_x_continuous(breaks = seq(0, 90, by = 30)) +
  scale_y_continuous(limits = c(-abun_ymax, abun_ymax)) +
  labs(
    x = "Distance from mine (km)",
    y = NULL,
    title = "<b>(b)</b> Species abundance"
  ) +
  base_theme +
  theme(
    plot.title = ggtext::element_markdown(size = 14, colour = "black"),
    legend.position = "none"
  )

distance_plot <- cowplot::plot_grid(richness_plot, abundance_plot, nrow = 1, align = "hv")
distance_plot

ggsave(
  "./figures/Figure S9.png",
  plot = distance_plot,
  width = 12, height = 6, dpi = 1000, bg = "white"
)
ggsave(
  "./figures/Figure S9.tif",
  plot = distance_plot,
  width = 12, height = 6, dpi = 1000, bg = "white", compression = "lzw"
)
##============================================================================##




