# Input/output handling -------------------------------------------------------

.hector_required_tile_columns <- c("ID", "MagnetX", "MagnetY", "type")
.hector_required_guide_columns <- c("ID", "MagnetX", "MagnetY", "type")

.hector_check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(data)
}

#' Read HECTOR tile and guide tables
#'
#' @param tile_file Path to the distortion-corrected hexabundle CSV.
#' @param guide_file Path to the distortion-corrected guide-star CSV.
#' @return A list with `tile` and `guide` data frames.
#' @export
hector_read_tables <- function(tile_file, guide_file) {
  if (!file.exists(tile_file)) stop("Tile file does not exist: ", tile_file,
                                    call. = FALSE)
  if (!file.exists(guide_file)) stop("Guide file does not exist: ", guide_file,
                                     call. = FALSE)

  tile <- read.table(
    file = tile_file, header = TRUE, sep = ",", comment.char = "#",
    check.names = TRUE, stringsAsFactors = FALSE,
    colClasses = c(ID = "character")
  )
  guide <- read.table(
    file = guide_file, header = TRUE, sep = ",", comment.char = "#",
    check.names = TRUE, stringsAsFactors = FALSE,
    colClasses = c(ID = "character")
  )
  .hector_check_columns(tile, .hector_required_tile_columns, "Tile table")
  .hector_check_columns(guide, .hector_required_guide_columns, "Guide table")
  list(tile = tile, guide = guide)
}

.hector_prepare_field <- function(tile, guide, engine) {
  sky <- grepl("Sky", tile$ID)
  sky_data <- tile[sky, , drop = FALSE]
  science <- tile[!sky, c("ID", "MagnetX", "MagnetY", "r_mag", "type"),
                  drop = FALSE]
  guide_input <- guide[, c("ID", "MagnetX", "MagnetY", "r_mag", "type"),
                       drop = FALSE]
  fdata <- rbind(science, guide_input)
  fdata$x <- fdata$MagnetX / 1000
  fdata$y <- fdata$MagnetY / 1000
  fdata$r <- sqrt(fdata$x^2 + fdata$y^2)

  keep <- fdata$r < (engine$fov / 2 - engine$excl_radius)
  fdata <- fdata[keep & !is.na(keep), , drop = FALSE]
  if (nrow(fdata) < engine$dngalprobes) {
    stop("The field contains fewer than ", engine$dngalprobes,
         " usable galaxy targets.", call. = FALSE)
  }

  # IDs are used as row names because the legacy engine propagates row names
  # through its position tables. Enforce uniqueness so output joins remain
  # deterministic.
  if (anyDuplicated(fdata$ID)) {
    stop("Input IDs must be unique across target and guide tables.", call. = FALSE)
  }
  rownames(fdata) <- fdata$ID

  type <- as.numeric(fdata$type)
  pos <- fdata[type == 1, c("x", "y"), drop = FALSE]
  stdpos <- fdata[type == 0, c("x", "y"), drop = FALSE]
  guidepos <- fdata[type == 2, c("x", "y"), drop = FALSE]
  if (nrow(pos) < engine$dngalprobes) {
    stop("The field contains fewer than ", engine$dngalprobes,
         " type-1 galaxy targets.", call. = FALSE)
  }
  list(
    fdata = fdata,
    sky_data = sky_data,
    pos_master = list(pos = pos, stdpos = stdpos, guidepos = guidepos)
  )
}

.hector_output_tables <- function(tile, guide, sky_data, final_config, fdata,
                                  engine) {
  configured_ids <- rownames(final_config$pos)
  selected_idx <- match(configured_ids, tile$ID)
  if (anyNA(selected_idx)) {
    stop("Some configured target IDs could not be found in the tile table.",
         call. = FALSE)
  }
  selected <- tile[selected_idx, , drop = FALSE]

  # The input tile may already contain geometry columns from an earlier
  # configuration pass. Keep the newly calculated geometry authoritative and
  # append only the non-derived input metadata. This also keeps target and
  # sky-fibre tables column-compatible for rbind().
  derived_columns <- c("x", "y", "rads", "angs", "azAngs", "angs_azAng")
  # `probe` is generated below. Exclude it when the input is already a
  # configured table so repeated configure/save cycles cannot create
  # duplicate probe columns.
  input_columns <- setdiff(names(tile), c("ID", "probe", derived_columns))

  properties <- data.frame(
    probe = seq_len(nrow(final_config$pos)),
    ID = configured_ids,
    final_config$pos,
    rads = final_config$rads,
    angs = final_config$angs,
    azAngs = final_config$azAngs,
    angs_azAng = final_config$angs_azAng,
    row.names = NULL,
    check.names = FALSE
  )
  target_table <- cbind(
    properties,
    selected[, input_columns, drop = FALSE]
  )
  rownames(target_table) <- NULL

  if (nrow(sky_data)) {
    sky_ids <- sky_data$ID
    sky_properties <- data.frame(
      probe = rep(-99, nrow(sky_data)),
      ID = sky_ids,
      x = sky_data$MagnetX / 1000,
      y = sky_data$MagnetY / 1000,
      rads = rep(-99, nrow(sky_data)),
      angs = rep(-99, nrow(sky_data)),
      azAngs = rep(-99, nrow(sky_data)),
      angs_azAng = rep(-99, nrow(sky_data)),
      row.names = NULL,
      check.names = FALSE
    )
    sky_table <- cbind(
      sky_properties,
      sky_data[, input_columns, drop = FALSE]
    )
    rownames(sky_table) <- NULL
    final_table <- rbind(target_table, sky_table)
  } else {
    final_table <- target_table
  }

  guide_ids <- rownames(final_config$gpos)
  guide_idx <- match(guide_ids, guide$ID)
  if (anyNA(guide_idx)) {
    stop("Some configured guide IDs could not be found in the guide table.",
         call. = FALSE)
  }
  guide_rows <- guide[guide_idx, , drop = FALSE]
  guide_input_columns <- setdiff(names(guide_rows), c("probe", derived_columns))
  guide_table <- cbind(
    guide_rows[, guide_input_columns, drop = FALSE],
    x = final_config$gpos$x,
    y = final_config$gpos$y,
    rads = final_config$grads,
    angs = final_config$gangs,
    azAngs = final_config$gazAngs,
    angs_azAng = final_config$gangs_gazAng
  )
  list(hexas = final_table, guides = guide_table)
}

#' Write configured HECTOR tables
#'
#' @param tables A list returned by `configure_hector()`.
#' @param hexafile_out Output path for the configured hexabundle table.
#' @param guidefile_out Output path for the configured guide table.
#' @export
hector_write_tables <- function(tables, hexafile_out, guidefile_out) {
  if (is.null(tables$hexas) || is.null(tables$guides)) {
    stop("tables must contain `hexas` and `guides` data frames.", call. = FALSE)
  }
  write.table(tables$hexas, file = hexafile_out, row.names = FALSE,
              sep = ",", quote = FALSE)
  write.table(tables$guides, file = guidefile_out, row.names = FALSE,
              sep = ",", quote = FALSE)
  invisible(tables)
}
