testthat::test_that("the packaged constants are readable", {
  constants <- hectorconfig::hector_constants()
  testthat::expect_equal(constants$HECTOR_plate_radius, 226)
})

testthat::test_that("the configuration engine can initialise its geometry", {
  parameters <- hectorconfig::hector_engine_parameters()
  testthat::expect_equal(parameters$fov, 452)
  testthat::expect_equal(parameters$dngalprobes, 19)
  testthat::expect_true(parameters$excl_radius > 0)
  testthat::expect_equal(parameters$delta_poly, 5)
})

testthat::test_that("example tables can be read", {
  tile <- system.file("extdata", "Hexas_G23_tile_263_NOT_CONFIGURED.csv",
                      package = "hectorconfig")
  guide <- system.file("extdata", "Guides_G23_tile_263_NOT_CONFIGURED.csv",
                       package = "hectorconfig")
  tables <- hectorconfig::hector_read_tables(tile, guide)
  testthat::expect_true(all(c("ID", "MagnetX", "MagnetY", "type") %in%
                              names(tables$tile)))
  testthat::expect_true(nrow(tables$tile) > 0)
})

testthat::test_that("configured output handles pre-existing geometry columns", {
  tile <- system.file("extdata", "Hexas_G23_tile_263_NOT_CONFIGURED.csv",
                      package = "hectorconfig")
  guide <- system.file("extdata", "Guides_G23_tile_263_NOT_CONFIGURED.csv",
                       package = "hectorconfig")
  tables <- hectorconfig::hector_read_tables(tile, guide)

  target_ids <- tables$tile$ID[!grepl("Sky", tables$tile$ID)][1:2]
  guide_ids <- tables$guide$ID[1:2]
  target_pos <- data.frame(x = c(1, 2), y = c(3, 4),
                           row.names = target_ids)
  guide_pos <- data.frame(x = c(5, 6), y = c(7, 8),
                          row.names = guide_ids)
  final_config <- list(
    pos = target_pos,
    rads = c(1, 1),
    angs = c(0, 0),
    azAngs = c(0, 0),
    angs_azAng = c(0, 0),
    gpos = guide_pos,
    grads = c(1, 1),
    gangs = c(0, 0),
    gazAngs = c(0, 0),
    gangs_gazAng = c(0, 0)
  )

  output <- hectorconfig:::.hector_output_tables(
    tile = tables$tile,
    guide = tables$guide,
    sky_data = tables$tile[grepl("Sky", tables$tile$ID), , drop = FALSE],
    final_config = final_config,
    fdata = NULL,
    engine = NULL
  )

  testthat::expect_equal(anyDuplicated(names(output$hexas)), 0L)
  testthat::expect_equal(anyDuplicated(names(output$guides)), 0L)
  testthat::expect_equal(nrow(output$guides), 2)
})

testthat::test_that("configured tables write and reload as ordinary CSV files", {
  tile <- system.file("extdata", "Hexas_G23_tile_263_NOT_CONFIGURED.csv",
                      package = "hectorconfig")
  guide <- system.file("extdata", "Guides_G23_tile_263_NOT_CONFIGURED.csv",
                       package = "hectorconfig")
  tables <- hectorconfig::hector_read_tables(tile, guide)

  target_ids <- tables$tile$ID[!grepl("Sky", tables$tile$ID)][1:2]
  guide_ids <- tables$guide$ID[1:2]
  final_config <- list(
    pos = data.frame(x = c(1, 2), y = c(3, 4), row.names = target_ids),
    rads = c(1, 1), angs = c(0, 0), azAngs = c(0, 0),
    angs_azAng = c(0, 0),
    gpos = data.frame(x = c(5, 6), y = c(7, 8), row.names = guide_ids),
    grads = c(1, 1), gangs = c(0, 0), gazAngs = c(0, 0),
    gangs_gazAng = c(0, 0)
  )
  output <- hectorconfig:::.hector_output_tables(
    tile = tables$tile,
    guide = tables$guide,
    sky_data = tables$tile[grepl("Sky", tables$tile$ID), , drop = FALSE],
    final_config = final_config,
    fdata = NULL,
    engine = NULL
  )

  hexafile <- tempfile(fileext = ".csv")
  guidefile <- tempfile(fileext = ".csv")
  hectorconfig::hector_write_tables(output, hexafile, guidefile)
  reloaded <- hectorconfig::hector_read_tables(hexafile, guidefile)

  testthat::expect_false(anyDuplicated(names(reloaded$tile)))
  testthat::expect_false(anyDuplicated(names(reloaded$guide)))
  testthat::expect_true(all(c("x", "y", "angs") %in% names(reloaded$tile)))
  testthat::expect_true(all(c("x", "y", "angs") %in% names(reloaded$guide)))
  testthat::expect_equal(nrow(reloaded$tile), nrow(output$hexas))
  testthat::expect_equal(nrow(reloaded$guide), nrow(output$guides))
})
