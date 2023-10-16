my_path <- "~/GIT_repos/R_Geo_workshop"

# load packages -------------------------------------------------
library(dplyr)
library(tectonicr)
library(sf)
library(ggplot2)
library(patchwork)


# import wsm2016 data and manipulate -------------------------------------------------
wsm2016 <- read.csv(file.path(my_path, "Data/wsm2016_donwload.csv")) |>
  rename_all(tolower) |>
  filter(azi != 999, depth <= 40, quality != "E") |>
  mutate(
    unc = ifelse(is.na(sd), quantise_wsm_quality(quality), sd),
    unc = ifelse(unc == 0, 1, unc),
    regime = ifelse(regime == "SS", "S", regime),
    regime = ifelse(regime == "TF", "T", regime),
    regime = ifelse(regime == "NF", "N", regime),
    regime = ifelse(regime == "U", NA, regime)
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = "WGS84", remove = FALSE)

# save wsm2016 into shp file
# saveRDS(wsm2016, "~/GIT_repos/R_Geo_workshop/Data/wsm2016.rds")
write_sf(wsm2016, file.path(my_path, "Data/wsm2016.shp"))

# load crop shape file -------------------------------------------------
crop_shp <- read_sf(file.path(my_path, "Data/crop.shp"))
stress_df <- sf::st_intersection(wsm2016, crop_shp) |> # crop wsm2016 with crop shape file
  mutate(radius = scales::rescale(unc, from = c(40, 0), to = c(0.1, 1)))
# import countries and plates sf object -------------------------------------------------
countries <- rnaturalearth::ne_countries(50, returnclass = "sf")
data("plates")


# DEM ---------------------------------------------------------------------
etopo2022 <- terra::rast("C:/Users/tobis/Downloads/ETOPO_2022_v1_60s_N90W180_bed.tif")
crop_extent <- st_bbox(stress_df) + c(-2, -2, 2, 2)
dem <- terra::crop(etopo2022, crop_extent)
elevation <- terra::values(dem)

slope <- terra::terrain(dem, "slope", unit = "radians")
aspect <- terra::terrain(dem, "aspect", unit = "radians")

hill <- terra::shade(slope, aspect,
  angle = c(45, 45, 45, 80),
  direction = c(225, 270, 315, 135)
)
hill <- Reduce(mean, hill)

# initialize map -------------------------------------------------
map <- ggplot() +
  tidyterra::geom_spatraster(data = hill, show.legend = FALSE) +
  scico::scale_fill_scico(palette = "grayC", begin = 0, end = .5) +
  ggnewscale::new_scale_fill() +
  tidyterra::geom_spatraster(data = dem, alpha = .25, show.legend = FALSE) +
  scico::scale_fill_scico(
    name = "Topography (m a.s.l.)",
    palette = "grayC",
    limits = range(elevation),
    oob = scales::squish,
    breaks = seq(-8000, 10000, 2000),
    na.value = NA
  ) +
  geom_sf(data = countries, color = "grey80", fill = NA, alpha = .5) +
  geom_sf(
    data = plates,
    color = "red",
    lwd = 1,
    alpha = .5
  ) +
  scale_alpha_discrete(name = "Quality rank", range = c(1, 0.4)) +
  theme_minimal() +
  theme_bw() +
  labs(x = "", y = "")

data_map <- map +
  geom_spoke(
    data = stress_df,
    aes(
      x = lon,
      y = lat,
      angle = deg2rad(90 - azi),
      color = regime,
      alpha = quality,
      radius = radius,
    ),
    # radius = .5,
    position = "center_spoke",
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Tectonic regime", values = stress_colors(),
    breaks = names(stress_colors())
  ) +
  coord_sf(
    xlim = range(stress_df$lon),
    ylim = range(stress_df$lat)
  )
#data_map

# load plate motions -------------------------------------------------
data("cpm_models")
models <- unique(cpm_models$model)


# Loop through all cpm models
s <- data.frame()
for (i in models) {
  motion <- filter(cpm_models, model == i)
  rel_motion <- equivalent_rotation(motion, fixed = "sa", rot = "nz")

  pb_pair <- ifelse(
    rel_motion$plate.rot < rel_motion$plate.fix,
    paste0(rel_motion$plate.rot, "-", rel_motion$plate.fix),
    paste0(rel_motion$plate.fix, "-", rel_motion$plate.rot)
  )
  file_ending <- paste0(pb_pair, "_", i, ".png")

  # anatolia_eu <- data.frame(
  #   plate.rot = "AN",
  #   lat = 30.7, lon = 32.6, angle = 1.2,
  #   plate.fix = "EU"
  # ) # McClusky, S., et al. (2000)
  #
  # tarim_eu <- data.frame(plate.rot = "AN", lat = 30.7, lon = 32.6, angle = 1.2, plate.fix = "EU") #


  # azimuth transformation -------------------------------------------------
  stress_df_trn <- PoR_shmax(stress_df, rel_motion, type = "in")
  stress_df2 <- cbind(stress_df, stress_df_trn)




  #circular_var(stress_df2$azi.PoR, w = 1 / stress_df2$unc)
  muci <- confidence_interval(stress_df2$azi.PoR, w = 1 / stress_df2$unc)

  #rayleigh_test(stress_df2$azi.PoR, mu = 90)

  # Rose diagram -------------------------------------------------
  png(
    file.path(my_path, paste0("Analysis_results/rose_", file_ending)),
    width = 5, height = 5, units = "in", res = 300
  )
  rose(
    stress_df2$azi.PoR,
    weights = 1 / stress_df2$unc,
    sub = paste0("Mean ± sde (°): ", round(muci$mu, 1), " ± ", round(muci$conf.angle, 1)),
    mtext = motion$model[1], main = pb_pair
  )
  rose_line(90, col = "dodgerblue", lty = 3, radius = 1.1)
  dev.off()


  trajectories <- eulerpole_smallcircles(rel_motion, 40)

  plot_deviation <- map +
    geom_sf(
      data = trajectories,
      lty = 2
    ) +
    geom_spoke(
      data = stress_df2,
      aes(
        x = lon,
        y = lat,
        angle = deg2rad(90 - azi),
        color = deviation_norm(dev),
        alpha = quality,
        radius = radius
      ),
      # radius = .5,
      position = "center_spoke",
      na.rm = TRUE
    ) +
    scale_color_continuous(
      type = "viridis",
      limits = c(0, 90),
      name = "|Deviation| in (\u00B0)",
      breaks = seq(0, 90, 22.5)
    ) +
    coord_sf(
      xlim = range(stress_df2$lon),
      ylim = range(stress_df2$lat)
    )



  # Plate boundary distance -------------------------------------------------
  plate_boundary <- filter(plates, pair %in% c(pb_pair, "ap-nz", "nd-nz"))


  stress_df2$distance <- distance_from_pb(
    x = stress_df2,
    PoR = rel_motion,
    pb = plate_boundary,
    tangential = FALSE
  )


  # Transformed azimuth vs. distance ----------------------------------------
  azi_plot <- ggplot(stress_df2, aes(x = distance, y = azi.PoR)) +
    coord_cartesian(ylim = c(0, 180)) +
    labs(x = "Distance from plate boundary (\u00B0)", y = "Azimuth in PoR (\u00B0)") +
    geom_hline(yintercept = c(0, 45, 90, 135, 180), lty = 3) +
    geom_pointrange(
      aes(
        ymin = azi.PoR - unc, ymax = azi.PoR + unc,
        color = stress_df2$regime, alpha = stress_df2$quality
      ),
      size = .25
    ) +
    scale_y_continuous(
      breaks = seq(-180, 360, 45),
      sec.axis = sec_axis(
        ~.,
        name = NULL,
        breaks = c(0, 45, 90, 135, 180),
        labels = c("Outward", "Tan (L)", "Inward", "Tan (R)", "Outward")
      )
    ) +
    scale_alpha_discrete(name = "Quality rank", range = c(1, 0.1)) +
    scale_color_manual(
      name = "Tectonic regime", values = stress_colors(), breaks = names(stress_colors())
    ) +
    theme_bw()


  stress_roll <- dplyr::arrange(stress_df2, distance)

  stress_roll$r_mean <- roll_circstats(
    stress_roll$azi.PoR,
    w = 1 / stress_roll$unc,
    FUN = circular_mean,
    width = 51
  )

  stress_roll$r_conf95 <- roll_confidence(
    stress_roll$azi.PoR,
    w = 1 / stress_roll$unc,
    width = 51
  ) / 2


  dist_plot1 <- azi_plot +
    pammtools::geom_stepribbon(
      data = stress_roll,
      aes(x = distance, ymin = r_mean - r_conf95, ymax = r_mean + r_conf95),
      color = NA, fill = "grey", alpha = 1 / 3
    ) +
    geom_step(
      data = stress_roll,
      aes(distance, r_mean),
      color = "red"
    )
  #dist_plot1

  # Dispersion vs distance --------------------------------------------------
  stress_roll$roll_disp <- roll_dispersion(
    x = stress_roll$azi.PoR,
    y = stress_roll$prd,
    w = 1 / stress_roll$unc,
    width = 51
  )

  dist_plot2 <-
    ggplot(stress_roll, aes(x = distance, y = cdist)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Distance from plate boundary (\u00B0)",
      y = expression(D[w](sigma["Hmax"], sigma["PB"]))
    ) +
    geom_hline(yintercept = c(0.15, .33, .7), lty = 3) +
    geom_point(aes(color = regime)) +
    scale_color_manual(name = "Tectonic regime", values = stress_colors(), breaks = names(stress_colors())) +
    geom_step(
      data = stress_roll,
      aes(distance, roll_disp),
      color = "red"
    ) +
    labs(title = "Dispersion") +
    theme_bw()





  # Interpolation -----------------------------------------------------------
  mean_SH_PoR <- PoR_stress2grid(
    stress_df,
    PoR = rel_motion, gridsize = .5, R_range = seq(50, 350, 100)
  )



  plot_interp_mdr <- map +
    geom_sf(data = trajectories, lty = 2) +
    geom_spoke(
      data = stress_df,
      aes(lon, lat, angle = deg2rad(90 - azi)),
      radius = .5, color = "grey30", position = "center_spoke"
    ) +
    geom_spoke(
      data = mean_SH_PoR,
      aes(lon, lat, angle = deg2rad(90 - azi), alpha = sd, color = mdr),
      radius = 1, position = "center_spoke", size = 1
    ) +
    scale_alpha(name = "Standard deviation", range = c(1, .25)) +
    scale_color_viridis_c(
      limits = c(0, 1),
      name = "Wavelength\n(R-normalized mean distance)",
      breaks = seq(0, 1, .25)
    ) +
    coord_sf(
      xlim = range(stress_df$lon),
      ylim = range(stress_df$lat)
    ) +
    facet_wrap(~R)




  mean_SH_PoR_compact <- mean_SH_PoR |>
    compact_grid() |>
    mutate(cdist = circular_distance(azi.PoR, 135))

  interp_map <- ggplot() +
    tidyterra::geom_spatraster(data = hill, show.legend = FALSE) +
    scico::scale_fill_scico(palette = "grayC", begin = 0, end = .5) +
    ggnewscale::new_scale_fill() +
    # tidyterra::geom_spatraster(data = dem_cropped, alpha = .8) +
    # tidyterra::scale_fill_wiki_c(
    #   name = "Topography (m a.s.l.)",
    #   limits = range(elevation),
    #   oob = scales::squish,
    #   breaks = seq(-5000, 5000, 1000),
    #   na.value = NA
    # )+
    geom_sf(data = countries, color = "grey80", fill = NA) +
    ggforce::geom_voronoi_tile(
      data = mean_SH_PoR_compact,
      aes(lon, lat, fill = cdist),
      max.radius = .7, normalize = FALSE, alpha = .8
    ) +
    geom_sf(
      data = trajectories,
      lty = 2
    ) +
    geom_spoke(
      data = mean_SH_PoR_compact,
      aes(lon, lat, angle = deg2rad(90 - azi)),
      radius = .5, color = "white", alpha = .2
    ) +
    scale_fill_viridis_c("Angular distance", limits = c(0, 1)) +
    geom_sf(data = countries, fill = NA) +
    geom_sf(
      data = plates,
      color = "red",
      lwd = 2,
      alpha = .5
    ) +
    coord_sf(
      xlim = range(stress_df$lon),
      ylim = range(stress_df$lat)
    ) +
    theme_bw() +
    labs(
      title = "Stress field analysis of the Andes", subtitle = paste0("Stress tested against ", pb_pair, " (", i, ")"),
      x = "", y = ""
    )

  (data_map | interp_map) / (dist_plot1 / dist_plot2) +
    plot_layout(guides = "collect")
  ggsave(
    file.path(my_path, paste0("Analysis_results/stress_", file_ending)),
    width = 8, height = 6, scale = 1.5
  )


  si <- stress_df2 |> mutate(
    model = i

    )
  si$ep_lat  <-  rep(rel_motion$lat[1], nrow(stress_df2))
  si$ep_lon <- rep(rel_motion$lon[1], nrow(stress_df2))
  si$ep_angle <- rep(rel_motion$angle[1], nrow(stress_df2))

  s <- plyr::rbind.fill(s, si)
}

# final summary of stress analysis
group_by(s, model) |>
  summarise(
    ep_lat = unique(ep_lat),
    ep_lon = unique(ep_lon),
    ep_angle = unique(ep_angle),
    cmean = circular_mean(azi.PoR, w = 1/unc) |> round(1),
    ci = confidence_angle(azi.PoR, w = 1/unc) |> round(1),
    cmedian = circular_median(azi.PoR, w = 1/unc) |> round(1),
    cvar = circular_var(azi.PoR, w = 1/unc) |> round(2),
    cdisp = circular_dispersion(azi.PoR, prd, w = 1/unc) |> round(2)
  ) |> gt::gt()

# width of plate boundarz zone:
group_by(s, model) |>
  filter(cdist <= 0.25) |>
  summarise(
    pbz_width = round(max(abs(distance)), 1)
  )|> gt::gt()
