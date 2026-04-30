#!/usr/bin/env Rscript

local_lib <- file.path(getwd(), ".Rlib")
dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))

ensure_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
    repos <- "https://cloud.r-project.org"
    # Su macOS preferiamo binari per evitare compilation/toolchain.
    if (is_macos) {
      try(
        install.packages(
          missing,
          repos = repos,
          lib = local_lib,
          type = "binary",
          dependencies = TRUE
        ),
        silent = TRUE
      )
    }
    missing2 <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing2) > 0) {
      install.packages(missing2, repos = repos, lib = local_lib, dependencies = TRUE)
    }
  }
  invisible(TRUE)
}

ensure_pkgs(c("sf", "dplyr", "readr", "ggplot2", "viridis"))
ensure_pkgs(c("stringr"))

suppressPackageStartupMessages({
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop(
      "Pacchetto 'sf' non disponibile.\n",
      "Suggerimento macOS: aggiorna R e installa binari CRAN, oppure installa dipendenze via Homebrew:\n",
      "  brew install gdal proj geos udunits openssl@3\n",
      "Poi reinstalla 'sf'."
    )
  }
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(viridis)
})

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

url_communes_geojson <- "https://raw.githubusercontent.com/gregoiredavid/france-geojson/master/communes.geojson"
url_restaurants_csv <- "https://raw.githubusercontent.com/holtzy/R-graph-gallery/master/DATA/data_on_french_states.csv"

path_communes_geojson <- file.path("data", "communes.geojson")
path_restaurants_csv <- file.path("data", "ristoranti_francia.csv")

if (!file.exists(path_communes_geojson)) {
  message("Scarico GeoJSON comuni...")
  download.file(url_communes_geojson, destfile = path_communes_geojson, mode = "wb", quiet = TRUE)
}

if (!file.exists(path_restaurants_csv)) {
  message("Scarico CSV ristoranti...")
  download.file(url_restaurants_csv, destfile = path_restaurants_csv, mode = "wb", quiet = TRUE)
}

message("Leggo dati...")
communes <- st_read(path_communes_geojson, quiet = TRUE)

restaurants_raw <- read_delim(
  path_restaurants_csv,
  delim = ";",
  quote = "\"",
  skip = 1,
  col_names = c("id", "reg", "dep", "depcom", "dciris", "an", "typequ", "nb_equip"),
  col_types = cols(
    id = col_character(),
    reg = col_character(),
    dep = col_character(),
    depcom = col_character(),
    dciris = col_character(),
    an = col_character(),
    typequ = col_character(),
    nb_equip = col_double()
  )
)

restaurants <- restaurants_raw %>%
  mutate(
    an = as.integer(an),
    nb_equip = as.numeric(nb_equip)
  ) %>%
  filter(an == 2016, typequ == "A504") %>%
  group_by(depcom) %>%
  summarise(nb_ristoranti = sum(nb_equip, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(depcom), depcom != "", nb_ristoranti > 0)

code_candidates <- c("code", "code_insee", "insee", "INSEE_COM", "code_commune", "codeCommune", "CODGEO")
code_col <- intersect(names(communes), code_candidates)
if (length(code_col) == 0) {
  stop(
    "Non trovo nel GeoJSON una colonna codice comune (INSEE). Colonne disponibili: ",
    paste(names(communes), collapse = ", ")
  )
}
code_col <- code_col[1]

communes <- communes %>%
  mutate(
    depcom = str_pad(as.character(.data[[code_col]]), width = 5, side = "left", pad = "0")
  )

communes_w <- communes %>%
  left_join(restaurants, by = "depcom") %>%
  filter(!is.na(nb_ristoranti), nb_ristoranti > 0)

communes_w_4326 <- st_transform(communes_w, 4326)
centroids_4326 <- st_point_on_surface(communes_w_4326)
xy_4326 <- st_coordinates(centroids_4326)

pts <- communes_w_4326 %>%
  mutate(lon = xy_4326[, 1], lat = xy_4326[, 2]) %>%
  st_drop_geometry()

# Sud della Francia: filtro con bounding box (approssimazione pratica e riproducibile)
south_bbox <- list(lon_min = -5.5, lon_max = 8.7, lat_min = 42.0, lat_max = 45.6)
pts_south <- pts %>%
  filter(
    lon >= south_bbox$lon_min, lon <= south_bbox$lon_max,
    lat >= south_bbox$lat_min, lat <= south_bbox$lat_max
  )

communes_south <- communes_w_4326 %>%
  mutate(lon = xy_4326[, 1], lat = xy_4326[, 2]) %>%
  filter(
    lon >= south_bbox$lon_min, lon <= south_bbox$lon_max,
    lat >= south_bbox$lat_min, lat <= south_bbox$lat_max
  )

south_outline <- st_union(st_geometry(communes_south))
south_outline <- st_sfc(south_outline, crs = st_crs(communes_south))

# KDE in metri (Lambert-93) per una densità più sensata
pts_sf <- st_as_sf(pts_south, coords = c("lon", "lat"), crs = 4326)
pts_2154 <- st_transform(pts_sf, 2154)
outline_2154 <- st_transform(south_outline, 2154)

xy_2154 <- st_coordinates(pts_2154)
df_kde <- pts_south %>%
  mutate(x = xy_2154[, 1], y = xy_2154[, 2])

# stat_density_2d_filled non gestisce "weight": approssimiamo replicando punti,
# con cap per mantenere tempi/memoria ragionevoli.
rep_w <- pmax(1, round(df_kde$nb_ristoranti / max(df_kde$nb_ristoranti, na.rm = TRUE) * 25))
df_kde_rep <- df_kde[rep(seq_len(nrow(df_kde)), rep_w), , drop = FALSE]

out_png <- file.path("output", "mappa_densita_ristoranti_sud_francia.png")

message("Creo mappa densità (salvo in output/)...")
p <- ggplot() +
  stat_density_2d_filled(
    data = df_kde_rep,
    aes(x = x, y = y, fill = after_stat(level)),
    alpha = 0.85,
    contour_var = "ndensity",
    bins = 20
  ) +
  geom_sf(data = outline_2154, fill = NA, color = "grey15", linewidth = 0.4) +
  scale_fill_viridis_d(option = "magma", direction = -1, name = "Densità\n(relativa)") +
  coord_sf(crs = 2154, lims_method = "geometry_bbox") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Densità dei ristoranti nel Sud della Francia",
    subtitle = "KDE pesata (A504, 2016) su centroidi comunali",
    x = NULL,
    y = NULL,
    caption = "Dati: R Graph Gallery (CSV) + france-geojson (communes)"
  )

ggsave(out_png, p, width = 10, height = 8, dpi = 200)
message("Fatto: ", out_png)

