# Public configuration workflow -----------------------------------------------

#' Configure a HECTOR field
#'
#' This is the package equivalent of `HECTOR_ClusterFieldsTest.R`. It reads a
#' distortion-corrected tile and guide table, configures galaxy, standard, and
#' guide probes, and optionally writes robot tables and a diagnostic plot.
#'
#' @param tile_file Path to the distortion-corrected hexabundle CSV.
#' @param guide_file Path to the distortion-corrected guide-star CSV.
#' @param hexafile_out Optional output path for the configured hexabundle CSV.
#' @param guidefile_out Optional output path for the configured guide CSV.
#' @param plot_file Optional output path for a PDF diagnostic plot.
#' @param constants Optional named list of constants. If supplied,
#'   `constants_path` is ignored.
#' @param constants_path Optional YAML constants path.
#' @param visualise If `TRUE`, draw intermediate plots during configuration.
#' @param plot_pause Seconds to pause after each visual redraw. This is used
#'   only when `visualise = TRUE`; set it to zero for the fastest redraws.
#' @param plot_device Graphics device for intermediate plots. `"current"`
#'   uses the active device; `"quartz"` opens a separate macOS Quartz window,
#'   reproducing the device used by the original driver script.
#' @return A list containing `config`, `hexas`, `guides`, `inputs`, and the
#'   configuration engine environment. When automatic guide placement leaves
#'   unresolved clashes, `config$manual_intervention` is `TRUE` and the output
#'   tables are retained for repair in the Shiny checker.
#' @export
configure_hector <- function(tile_file, guide_file, hexafile_out = NULL,
                             guidefile_out = NULL, plot_file = NULL,
                             constants = NULL, constants_path = NULL,
                             visualise = FALSE, plot_pause = 0.1,
                             plot_device = c("current", "quartz")) {
  plot_device <- match.arg(plot_device)
  engine <- .hector_new_engine(
    constants = constants,
    constants_path = constants_path,
    visualise = visualise,
    plot_pause = plot_pause
  )
  tables <- hector_read_tables(tile_file, guide_file)
  prepared <- .hector_prepare_field(tables$tile, tables$guide, engine)
  engine$pos_master <- prepared$pos_master
  engine$visualise <- isTRUE(visualise)
  if (isTRUE(visualise)) {
    .hector_open_plot_device(plot_device)
  }

  message("Configuring field ", normalizePath(tile_file, mustWork = FALSE))
  config <- engine$configure_HECTOR(pos_master = prepared$pos_master)
  if (isTRUE(config$manual_intervention)) {
    message("Configuration returned with flags requiring manual intervention: ",
            config$flags)
  }
  outputs <- .hector_output_tables(
    tile = tables$tile,
    guide = tables$guide,
    sky_data = prepared$sky_data,
    final_config = config,
    fdata = prepared$fdata,
    engine = engine
  )

  if (!is.null(hexafile_out) || !is.null(guidefile_out)) {
    if (is.null(hexafile_out) || is.null(guidefile_out)) {
      stop("hexafile_out and guidefile_out must be supplied together.",
           call. = FALSE)
    }
    hector_write_tables(outputs, hexafile_out, guidefile_out)
  }
  if (!is.null(plot_file)) {
    engine$plot_configured_field(
      filename = plot_file,
      pos = config$pos,
      angs = config$angs,
      gpos = config$gpos,
      gangs = config$gangs,
      fieldflags = config$flags,
      aspdf = TRUE
    )
  }

  list(
    config = config,
    hexas = outputs$hexas,
    guides = outputs$guides,
    inputs = prepared,
    engine = engine
  )
}
