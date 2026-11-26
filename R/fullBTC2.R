#' fullBTC2
#'
#' The land of the downscaling along the transitions in the land classes of the project Bending the Curve
#' 2.0 (\code{\link{calcLandReportBtC2}}), as grid cell products: \code{cell.land_btc2_0.5.mz} (Mha by grid
#' cell and class) and three NetCDF files named after the modelling protocol (MP 3.4), with shares of the total
#' grid cell area (\code{pixel_area}, Mha; the area of the 0.5 degree cell on a sphere of mean earth radius, or
#' the land of the cell where that is larger):
#' \itemize{
#'   \item \code{<projectName>-<scenario>-<model>-<runNr>.nc}: \code{LC_area_share} of the classes, coded 1 to
#'   12 in the order of the classes;
#'   \item \code{<projectName>-<scenario>-convyear-<model>-<runNr>.nc}: \code{Fraction}, the area of classes 4
#'   to 11 by period of conversion. Area that a class gains in the grid cell in a time step gets the period of
#'   that step (code of the last two digits of its start and end year, 2095 to 2100 = 9510); area it loses
#'   leaves all periods in proportion. Area of the first year has no period;
#'   \item \code{<projectName>-<scenario>-fate-aband-recon-<model>-<runNr>.nc}: classes 6 and 7 split into
#'   \code{LC_area_share_forest} and \code{LC_area_share_nonforest} (\code{calcLandReportBtC2} with
#'   \code{fate = TRUE}).
#' }
#' The scenario is the run title without the flag and without underscores (BTC2v12_INT_26: INT26).
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' @param projectName first part of the NetCDF file names
#' @param model name of the model in the NetCDF file names
#' @param runNr run number in the NetCDF file names
#' @param logFile name of the log of the consistency checks
#' @author Patrick v. Jeetze, Florian Humpenoeder
#' @seealso \code{\link{downscaleRun}}
fullBTC2 <- function(rev = numeric_version("0"), projectName = "BTC2", model = "MAgPIE", runNr = "run1",
                     logFile = "consistencyCheck_BtC2.log") {
  # Writes a NetCDF file of variables on the 0.5 degree grid with the global attributes of the protocol;
  # fill(nc, var) writes the values of var.
  .writeNC <- function(name, vars, fill) {
    nc <- ncdf4::nc_create(paste0(name, ".nc"), vars, force_v4 = TRUE)
    for (v in vars) {
      fill(nc, v)
    }
    ncdf4::ncatt_put(nc, 0, "title", name)
    ncdf4::ncatt_put(nc, 0, "scenario", scenario)
    ncdf4::ncatt_put(nc, 0, "institution", "PIK")
    ncdf4::ncatt_put(nc, 0, "contact",
                     "F. Humpenöder, humpenoeder@pik-potsdam.de; P. v. Jeetze, vjeetze@pik-potsdam.de")
    ncdf4::ncatt_put(nc, 0, "history", date())
    ncdf4::nc_close(nc)
  }

  # Values by grid cell (rows) as array on the grid, the other dimensions as given.
  .toGrid <- function(x, dims) {
    n <- prod(dims)
    out <- array(NA_real_, c(length(lon), length(lat), n))
    x <- matrix(x, ncol = n)
    for (k in seq_len(n)) {
      out[cbind(index, k)] <- x[, k]
    }
    dim(out) <- c(length(lon), length(lat), dims)
    return(out)
  }

  landHr <- calcOutput("LandReportBtC2", input = "magpie", aggregate = FALSE)
  write.magpie(landHr, "cell.land_btc2_0.5.mz", comment = "unit: Mha per grid-cell")
  fateHr <- calcOutput("LandReportBtC2", input = "magpie", fate = TRUE, aggregate = FALSE)

  cfg <- readSource("MagpieRunFiles", subtype = "config")
  scenario <- sub(paste0("^", cfg$info$flag, "_"), "", cfg$title)
  stem <- paste(projectName, gsub("_", "", scenario), sep = "-")
  classes <- getItems(landHr, dim = 3)
  years <- getYears(landHr, as.integer = TRUE)

  # grid cells on the 0.5 degree grid, the latitude from north to south
  lon <- seq(-179.75, 179.75, by = 0.5)
  lat <- -seq(-89.75, 89.75, by = 0.5)
  coords <- getCoords(landHr)
  index <- cbind(match(coords[["x"]], lon), match(coords[["y"]], lat))
  # total area of the grid cell (sphere of mean earth radius), at least its land
  landArea <- dimSums(landHr, dim = 3)
  lat1 <- (coords[["y"]] + 0.25) * pi / 180
  lat2 <- (coords[["y"]] - 0.25) * pi / 180
  pixelArea <- 0.5 * pi / 180 * 6371.0088^2 * abs(sin(lat1) - sin(lat2)) / 1e4
  pixelArea <- pmax(pixelArea, as.vector(landArea[, 1, ]))

  dimLon <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  dimLat <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  dimTime <- ncdf4::ncdim_def("time", "years", years, calendar = "standard", unlim = TRUE)
  .dimClass <- function(selected) {
    codes <- match(selected, classes)
    return(ncdf4::ncdim_def("lc_class", paste(codes, selected, sep = "=", collapse = "/"), codes,
                            longname = "land cover class"))
  }
  areaVar <- ncdf4::ncvar_def("pixel_area", units = "million ha", longname = "total area of pixel",
                              dim = list(dimLon, dimLat), missval = NaN, prec = "double", compression = 9)
  .fill <- function(values) {
    return(function(nc, v) {
      ncdf4::ncvar_put(nc, v, if (v$name == "pixel_area") .toGrid(pixelArea, 1)[, , 1] else values[[v$name]])
    })
  }

  # land cover classes; magclass order cell x year x class, NetCDF order class x time
  shareVar <- ncdf4::ncvar_def("LC_area_share", units = "share of total pixel area",
                               longname = "share of total pixel area occupied by various LC classes",
                               dim = list(dimLon, dimLat, .dimClass(classes), dimTime), missval = NaN,
                               prec = "double", compression = 9)
  share <- aperm(as.array(landHr / pixelArea), c(1, 3, 2))
  .writeNC(paste(stem, model, runNr, sep = "-"), list(shareVar, areaVar),
           .fill(list(LC_area_share = .toGrid(share, dim(share)[2:3]))))

  # area by period of conversion
  ageClasses <- classes[4:11]
  periods <- paste0(formatC((years[-length(years)]) %% 100, width = 2, flag = "0"),
                    formatC(years[-1] %% 100, width = 2, flag = "0"))
  periods[length(periods)] <- sub("00$", "10", periods[length(periods)])
  dimConv <- ncdf4::ncdim_def("convyear", paste(periods, collapse = ","), as.integer(periods),
                              longname = "conversion year")
  fracVar <- ncdf4::ncvar_def("Fraction", units = "share of total pixel area",
                              longname = paste("share of total pixel area that was restored in each 5-year",
                                               "(or 10-year, depending on IAM) period"),
                              dim = list(dimLon, dimLat, dimTime, .dimClass(ageClasses), dimConv), missval = NaN,
                              prec = "double", compression = 9)
  # one class at a time (the whole grid in double precision would take about 8 GB)
  .fillFraction <- function(nc, v) {
    if (v$name == "pixel_area") {
      return(.fill(NULL)(nc, v))
    }
    for (k in seq_along(ageClasses)) {
      area <- as.array(landHr[, , ageClasses[k]])[, , 1]
      untracked <- area[, 1]
      cohorts <- matrix(0, nrow(index), length(periods))
      frac <- array(0, c(nrow(index), length(years), length(periods)))
      for (i in seq_along(years)[-1]) {
        before <- untracked + rowSums(cohorts)
        keep <- pmin(area[, i] / replace(before, before <= 0, 1), 1)
        untracked <- untracked * keep
        cohorts <- cohorts * keep
        cohorts[, i - 1] <- cohorts[, i - 1] + pmax(area[, i] - before, 0)
        frac[, i, ] <- cohorts / pixelArea
      }
      ncdf4::ncvar_put(nc, v, .toGrid(frac, dim(frac)[2:3]), start = c(1, 1, 1, k, 1),
                       count = c(-1, -1, -1, 1, -1))
    }
  }
  .writeNC(paste(stem, "convyear", model, runNr, sep = "-"), list(fracVar, areaVar), .fillFraction)

  # fate of classes 6 and 7
  fateClasses <- getItems(fateHr, dim = 3.1)
  fateVars <- lapply(c(forest = "forest", nonforest = "non-forest"), function(f) {
    ncdf4::ncvar_def(paste0("LC_area_share_", sub("-", "", f)), units = "1",
                     longname = paste("share of total pixel area converted into", f),
                     dim = list(dimLon, dimLat, .dimClass(fateClasses), dimTime), missval = NaN,
                     prec = "double", compression = 9)
  })
  fateValues <- lapply(c(forest = "forest", nonforest = "nonforest"), function(f) {
    x <- aperm(as.array(collapseDim(fateHr[, , f], dim = 3.2) / pixelArea), c(1, 3, 2))
    return(.toGrid(x, dim(x)[2:3]))
  })
  names(fateValues) <- paste0("LC_area_share_", names(fateValues))
  .writeNC(paste(stem, "fate-aband-recon", model, runNr, sep = "-"), c(fateVars, list(areaVar)),
           .fill(fateValues))

  toolWriteMadratLog(logPath = logFile)
}
