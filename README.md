# hectorconfig_R

Note from calofost: this package is based on code written pre-AI but the packaging and documentation was done by AI ChatGPT. I have tested what I could, but please use a healthy dose of scepticism and kindly report issues.

This repository contains the installable R package `hectorconfig`, which
configures the positions and orientations of HECTOR
hexabundles, standard-star probes, and guide probes after a field has been
selected and distortion-corrected. It is based on Caroline Foster's
`HECTOR_Config_v3.5.1.R` configuration engine and the accompanying cluster
field driver.

The package also includes the post-configuration Shiny checker used to inspect
a field, rotate individual probes, check conflicts, reset probe orientations,
and download updated hexabundle and guide tables.

## What the configuration does

The workflow is:

1. Read the distortion-corrected hexabundle and guide-star CSV files.
2. Convert `MagnetX` and `MagnetY` from micrometres to millimetres.
3. Remove sky fibres from the target-selection table while retaining them for
   the final hexabundle output.
4. Keep positions inside the HECTOR field of view.
5. Select cable-exit gaps by counting the positions in front of each probe.
6. Set initial probe angles towards the selected cable exits.
7. Search for non-overlapping orientations, trying wider angle ranges and
   cable-exit changes when necessary.
8. Add two non-conflicting standard-star probes.
9. Attempt to configure six guide probes. If automatic placement leaves
   unresolved guide clashes, all available guide candidates are tested before
   unresolved candidates are retained, the result is flagged for manual
   intervention, and the tables can still be opened in the Shiny checker.
   The guide solver tests discrete 180-degree flips, then performs a local
   stochastic angular search around each flipped state. For guide-guide
   clashes, both guide probes are included in that search. If a guide-target
   clash remains, the implicated target/standard and guide angles are searched
   jointly, while trials that create target-target clashes are rejected.
   Every trial also checks the complete target and guide probe polygons against
   the circular field boundary; out-of-field configurations are rejected and
   flagged as `FieldFail` if retained for manual intervention.
   The joint target-guide search is intentionally kept small because it is a
   low-dimensional search. Its 180-degree single/pair starts are attempted
   before the broader guide-only fallback; a larger joint search is used only
   as a last chance for unusually difficult conflicts.
10. Check for the wedge condition involving probe heads, fibre bend radius,
    and nearby probe tips.
11. Write configured robot tables and, optionally, a diagnostic PDF.

The legacy engine approximates probe shapes as polygons. The rounded probe
head is sampled in `delta_poly`-degree increments; the default is 5 degrees.
The conflict search is stochastic, so a later revision should expose and
record a random seed for exact reproducibility.

## Installation

### For a first-time R user

Install a current version of R from [CRAN](https://cran.r-project.org/).
RStudio is optional, but is a convenient interface on macOS and Windows.
Open RStudio or the R console and run the following block once:

```r
install.packages("remotes")
install.packages(c(
  "PBSmapping", "plotrix", "R.utils", "sf", "shiny", "shinyjs", "sp",
  "testthat", "yaml"
))
remotes::install_github("asgr/magicaxis", dependencies = TRUE)
```

`magicaxis` is installed from its development repository because it is not
assumed to be available in the same form on every machine. Its GitHub
dependencies are installed by `dependencies = TRUE`; the other package
dependencies above are installed from CRAN.

On Windows, R may ask to install Rtools when a package must be built from
source; accept the matching Rtools version. On Linux, `sf` may also require
the system libraries GDAL, GEOS, and PROJ. The CRAN `sf` documentation lists
the commands for the Linux distribution in use.

Now install the package directly from GitHub within R:

```r
remotes::install_github("calofost/hectorconfig_R", dependencies = TRUE)
library(hectorconfig)
packageVersion("hectorconfig")
```

The package also declares `asgr/magicaxis` in its `Remotes` field, so a fresh
`remotes::install_github()` should install it automatically. The explicit
commands above make the process easier to diagnose if a dependency fails.
The old machine-specific paths are not used: the default constants file is
shipped with the package, and a different constants file can be supplied
explicitly.

To check the installation before configuring a field:

```r
required <- c("magicaxis", "PBSmapping", "plotrix", "R.utils", "sf",
              "shiny", "shinyjs", "sp", "yaml")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Still missing: ", paste(missing, collapse = ", "))
hectorconfig::hector_engine_parameters()
```

## Example configuration

The package includes the supplied G23 tile and guide files as example data.
This complete example can be run immediately after installation:

```r
library(hectorconfig)

tile <- system.file(
  "extdata",
  "Hexas_G23_tile_263_NOT_CONFIGURED.csv",
  package = "hectorconfig"
)
guide <- system.file(
  "extdata",
  "Guides_G23_tile_263_NOT_CONFIGURED.csv",
  package = "hectorconfig"
)

out_dir <- file.path(tempdir(), "hectorconfig-example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result <- configure_hector(
  tile_file = tile,
  guide_file = guide,
  hexafile_out = file.path(out_dir, "hexas_configured.csv"),
  guidefile_out = file.path(out_dir, "guides_configured.csv"),
  plot_file = file.path(out_dir, "configuration.pdf"),
  visualise = FALSE
)

result$config$flags
result$config$swaps
result$config$guide_count
result$config$manual_intervention
list.files(out_dir)
```

The example writes its output to a temporary directory. Replace `out_dir`
with a permanent directory when configuring a real field.

To watch the intermediate configuration, use `visualise = TRUE`. The
original driver explicitly opened a separate macOS Quartz device; that
behavior is available with `plot_device = "quartz"`:

```r
result <- configure_hector(
  tile_file = tile,
  guide_file = guide,
  visualise = TRUE,
  plot_device = "quartz",
  plot_pause = 0.2
)
```

To use the active RStudio graphics device instead, leave
`plot_device = "current"` (the default). `plot_pause` is only
applied during visualisation and is ignored when `visualise = FALSE`;
it is a redraw aid, not the mechanism that creates the graphics device.

The bundled G23 example contains six guide candidates. The solver now tests
180-degree flips of conflicting guide probes before its wider angular search,
but the result can still depend on the candidate geometry and stochastic
search. If unresolved candidates remain, the package retains six guide rows,
marks the result with `GuideFail`, and writes the tables for manual repair in
the checker. The output must not be sent to the robot until
`result$config$manual_intervention` is `FALSE` and the checker reports no
conflicts.

For your own input files, use the same function with explicit paths:

```r
result <- configure_hector(
  tile_file = "input/hexas.csv",
  guide_file = "input/guides.csv",
  hexafile_out = "output/hexas_configured.csv",
  guidefile_out = "output/guides_configured.csv",
  plot_file = "output/configuration.pdf"
)
```

All paths are arguments to the function. The caller's working directory is
not changed. To use another constants file:

```r
result <- configure_hector(
  tile_file = "input/hexas.csv",
  guide_file = "input/guides.csv",
  constants_path = "configuration/HECTOR_CONSTANTS.yaml"
)
```

The returned object contains:

- `config`: positions, angles, selected guides, swaps, field flags, guide
  counts, and a `manual_intervention` status;
- `hexas`: the configured hexabundle table, including retained sky-fibre rows;
- `guides`: the configured guide table;
- `inputs`: the prepared position tables used by the engine;
- `engine`: the private engine environment, useful for advanced diagnostics.

## Checking and reloading files saved by the Shiny app

The two download buttons write ordinary comma-separated text files. Always
download both files together: the Hexa file contains target, standard, and
sky-fibre rows, while the Guide file contains the guide rows. After editing a
field, check the files with:

```r
hexas <- read.csv("output/hexas_configured.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
guides <- read.csv("output/guides_configured.csv",
                   stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(!anyDuplicated(names(hexas)))
stopifnot(!anyDuplicated(names(guides)))
stopifnot(all(c("ID", "x", "y", "angs", "fibre_type") %in% names(hexas)))
stopifnot(all(c("ID", "x", "y", "angs", "fibre_type") %in% names(guides)))

# The important round-trip test: the files should be readable by the package
# again, not merely parseable as text.
reloaded <- hector_read_tables(
  "output/hexas_configured.csv",
  "output/guides_configured.csv"
)
```

If reload fails, check that both files came from the same field and that their
headers have not been edited. A repeated Hexa download should contain one
`probe` column, never `probe.1` or a second `probe` column.

## Input tables

The tile and guide files are comma-separated tables. Both require at least:

```text
ID, MagnetX, MagnetY, type
```

The current tiling workflow uses:

- `type == 1`: galaxy targets;
- `type == 0`: standard stars;
- `type == 2`: guide stars;
- IDs containing `Sky`: sky-fibre rows retained in the final hexabundle file.

The input examples in `inst/extdata/` are the supplied G23 tile and guide
tables. Additional columns are carried through to the output where possible.

## Shiny checker

Create the app and launch it with:

```r
library(hectorconfig)
shiny::runApp(hector_checker_app())
```

Or launch it directly:

```r
hector_run_app()
```

The checker expects the configured hexabundle and guide CSV files. Double-click
near a probe head to select it, click to choose a new tail direction, and use
`Fix Probe` only after checking the reported conflict count. The original app
instructions intentionally warn that manual changes do not automatically
check whether a probe has been oriented outside the field of view; that remains
an important limitation of this initial package version.

## Main geometric defaults

The packaged YAML constants set a HECTOR plate radius of 226 mm. The R engine
then derives a 452 mm field diameter and uses three cable-exit directions at
`-90`, `30`, and `150` degrees. The current probe defaults are:

- 19 galaxy probes and 2 standard-star probes;
- 6 guide probes are required for a complete automatic configuration;
  unresolved guide clashes produce six flagged rows for manual repair;
- 45 mm probe and cable lengths;
- 45 mm fibre bend radius;
- 5-degree polygon sampling for the probe head;
- 600 seconds per conflict-minimisation attempt.

These values are code/configuration defaults, not assumptions hidden in the
public functions. A future version should expose the non-YAML geometry values
as a formal configuration object rather than relying on the legacy engine
environment.

## Development notes

This first packaging pass deliberately keeps the mature configuration engine
largely intact. Each configuration run receives a fresh private environment,
which contains the legacy engine's mutable state. This removes dependence on
`.GlobalEnv` and makes the constants and file paths explicit without changing
the geometry algorithm all at once.

The next validation steps are:

1. Run the package tests in an R installation with the dependencies available.
2. Compare the configured output for the supplied G23 example with the current
   known-good output.
3. Add deterministic seeds and structured failure conditions.
4. Replace the remaining polygon-geometry dependencies where this improves
   installation and portability.
5. Translate the tested, public workflow to Python and compare R/Python output
   tables field by field.

## Related Python repository

The Python translation is maintained separately in
[`calofost/hectorconfig_python`](https://github.com/calofost/hectorconfig_python).
It has the corresponding explicit interface:

```python
result = configure_hector(
    tile_file="Hexas_G23_tile_263_NOT_CONFIGURED.csv",
    guide_file="Guides_G23_tile_263_NOT_CONFIGURED.csv",
    hexafile_out="Hexas_G23_tile_263_CONFIGURED.csv",
    guidefile_out="Guides_G23_tile_263_CONFIGURED.csv",
)
```

The R output remains the reference implementation. The Python package is
useful for development and comparison, but its current solver has performed
substantially worse than this R implementation in testing and should not be
treated as a validated replacement for observing.
