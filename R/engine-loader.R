# Internal engine construction -------------------------------------------------

# The original configuration code predates package-style state management and
# uses <<- for several quantities that are shared by its functions. Each call
# therefore gets a fresh environment. This keeps those assignments local to a
# configuration or Shiny session rather than leaking them into .GlobalEnv.

.hector_default_constants_path <- function() {
  path <- system.file("constants", "HECTOR_CONSTANTS.yaml", package = "hectorconfig")
  if (!nzchar(path)) {
    stop("The package constants file could not be located.", call. = FALSE)
  }
  path
}

.hector_bind_dependencies <- function(engine) {
  # The legacy source is evaluated in a child environment. Explicitly binding
  # dependencies here is more robust than relying on the package namespace's
  # import environment being visible from that child during sys.source().
  bindings <- list(
    graphics.off = grDevices::graphics.off,
    dev.off = grDevices::dev.off,
    dev.flush = grDevices::dev.flush,
    lines = graphics::lines,
    mtext = graphics::mtext,
    pdf = grDevices::pdf,
    points = graphics::points,
    polygon = graphics::polygon,
    text = graphics::text,
    magplot = magicaxis::magplot,
    as.PolySet = PBSmapping::as.PolySet,
    joinPolys = PBSmapping::joinPolys,
    draw.circle = plotrix::draw.circle,
    withTimeout = R.utils::withTimeout,
    st_as_sf = sf::st_as_sf,
    st_distance = sf::st_distance,
    actionButton = shiny::actionButton,
    downloadButton = shiny::downloadButton,
    downloadHandler = shiny::downloadHandler,
    fileInput = shiny::fileInput,
    fluidPage = shiny::fluidPage,
    h5 = shiny::h5,
    h6 = shiny::h6,
    hr = shiny::hr,
    mainPanel = shiny::mainPanel,
    observeEvent = shiny::observeEvent,
    plotOutput = shiny::plotOutput,
    renderPlot = shiny::renderPlot,
    renderText = shiny::renderText,
    req = shiny::req,
    runApp = shiny::runApp,
    shinyApp = shiny::shinyApp,
    sidebarLayout = shiny::sidebarLayout,
    sidebarPanel = shiny::sidebarPanel,
    tabPanel = shiny::tabPanel,
    tabsetPanel = shiny::tabsetPanel,
    textOutput = shiny::textOutput,
    useShinyjs = shinyjs::useShinyjs,
    Polygon = sp::Polygon,
    Polygons = sp::Polygons,
    SpatialPoints = sp::SpatialPoints,
    SpatialPolygons = sp::SpatialPolygons,
    read_yaml = yaml::read_yaml
  )
  list2env(bindings, envir = engine)
  invisible(engine)
}

.hector_open_plot_device <- function(plot_device) {
  if (identical(plot_device, "quartz")) {
    if (!grepl("darwin", R.version$os, ignore.case = TRUE)) {
      stop("plot_device = 'quartz' is only available on macOS.",
           call. = FALSE)
    }
    grDevices::quartz()
  }
  invisible(NULL)
}

.hector_new_engine <- function(constants = NULL, constants_path = NULL,
                               visualise = FALSE, plot_pause = 0.1) {
  if (is.null(constants)) {
    if (is.null(constants_path)) {
      constants_path <- .hector_default_constants_path()
    }
    if (!file.exists(constants_path)) {
      stop("constants_path does not exist: ", constants_path, call. = FALSE)
    }
    constants <- yaml::read_yaml(constants_path)
  }

  if (!is.list(constants) || is.null(constants$HECTOR_plate_radius)) {
    stop("constants must be a list containing HECTOR_plate_radius.", call. = FALSE)
  }

  engine <- new.env(parent = asNamespace("hectorconfig"))
  engine$hector_constants <- constants
  engine$visualise <- isTRUE(visualise)
  engine$plot_pause <- max(0, as.numeric(plot_pause)[1])
  .hector_bind_dependencies(engine)

  source_file <- system.file("engine", "HECTOR_Config_v3.5.1.R",
                             package = "hectorconfig")
  if (!nzchar(source_file)) {
    stop("The packaged HECTOR configuration engine could not be located.",
         call. = FALSE)
  }
  sys.source(source_file, envir = engine)
  engine$visualise <- isTRUE(visualise)
  engine
}

#' Read HECTOR constants
#'
#' @param path Optional path to a YAML constants file. If omitted, the
#'   constants shipped with the package are used.
#' @return A named list of constants.
#' @export
hector_constants <- function(path = NULL) {
  if (is.null(path)) path <- .hector_default_constants_path()
  if (!file.exists(path)) {
    stop("Constants file does not exist: ", path, call. = FALSE)
  }
  yaml::read_yaml(path)
}

#' Return the active geometric parameters of an engine
#'
#' @param engine An internal engine environment, normally obtained indirectly
#'   from `configure_hector()`.
#' @return A named list containing the main field and probe parameters.
#' @keywords internal
hector_engine_parameters <- function(engine = .hector_new_engine()) {
  names <- c(
    "fov", "skybuffer", "wceg", "fbr", "tip_w", "tip_l", "probe_w",
    "probe_l", "cable_w", "cable_l", "excl_radius", "dngalprobes",
    "nstdprobes", "ngprobesmin", "ngprobesmax", "gtip_w", "gtip_l",
    "gprobe_w", "gprobe_l", "gcable_w", "gcable_l", "gexcl_radius",
    "timeout", "delta_poly"
  )
  as.list(setNames(lapply(names, get, envir = engine), names))
}
