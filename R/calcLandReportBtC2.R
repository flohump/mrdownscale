#' calcLandReportBtC2
#'
#' The land of the downscaling along the transitions (\code{\link{calcLandTransitionsHighRes}}) in the land
#' classes of the project Bending the Curve 2.0 (BtC2), by grid cell. The classes only reclassify the
#' downscaled land (\code{cell.land_0.5.mz}); the land itself is not changed.
#'
#' @details
#' Classes of the state of the land, from the land split (\code{\link{calcLandReport}} with
#' \code{outputFormat = "downscaledmagpieSplit"}):
#' \itemize{
#'   \item cropland_other and cropland_2Gbioen: cropland without and with the crop area of bioenergy crops;
#'   \item grassland: pasture and rangeland;
#'   \item nat_regen_forest: primary and secondary forest, NDC forest, and afforestation if it is not
#'   established as plantation (\code{s32_aff_plantation});
#'   \item forest_planted: timber plantations, other planted forest, and afforestation as plantation;
#'   \item built_up: urban land.
#' }
#' Classes of the history of natural land (other land and secondary forest), along the grid cell transitions
#' from \code{startYear} on (\code{\link{toolGriddedLandOrigin}}): area that came from cropland
#' (recon_abn_crop), pasture (recon_abn_grass), or timber plantations and other planted forest (recon_forest)
#' keeps its class while it stays natural, also when other land matures to secondary forest; the rest is other
#' or nat_regen_forest. Conversions take the area that entered last; maturation and the conversion to NDC forest
#' and to afforestation that is not a plantation take the area that entered first.
#'
#' Reconverted managed forest (recon_forest): secondary forest protected from harvest by the reconversion
#' targets of the run (\code{p22_reconv_forest}), placed in the grid cells of each cluster in proportion to
#' the secondary forest that is not recovering land; placed area is kept.
#'
#' Restoration of the optimisation (the lower bounds \code{p35_land_restoration}): in each cluster the bound
#' as share of the transitions into the restored type in the optimisation stage; the transitions of each
#' grid cell are split by this share. Conversions take restored land (and in runs with reconversion targets
#' all land from cropland and pasture) only after all other area of the type in the grid cell.
#'
#' Rehabilitated cropland and grassland (rehab_crop, rehab_grass) are not yet filled.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @param startYear last year before the history classes start
#' @return land by grid cell and class (Mha) in the years of the input
#' @author Patrick v. Jeetze, Florian Humpenoeder
#' @seealso \code{\link{fullBTC2}}, \code{\link{toolGriddedLandOrigin}}
calcLandReportBtC2 <- function(input = "magpie", startYear = 2020) {
  # Land report of the run in one of the MAgPIE formats, with the arguments of fullDOWNSCALEDMAGPIE (same cache).
  .report <- function(outputFormat) {
    x <- calcOutput("LandReport", outputFormat = outputFormat, input = input, harmonizationPeriod = c(2015, 2050),
                    yearsSubset = seq(1995, 2100, 5), harmonization = "none", downscaling = "transitions",
                    aggregate = FALSE)
    return(x[, "y1985", , invert = TRUE])
  }

  # Splits the transitions into secondary forest and other land into the restored types by the restoration
  # share of their cluster (restoration bound per transitions of the cluster); the restored part has the source
  # "<source>Restored". The bounds apply to all transitions of the model, which magpie4::landTransitions
  # splits into the stages snvRelocation, peatRewetting and optimisation.
  .splitRestoration <- function(transitionsHr, restorationLr, mapping) {
    stage <- getItems(transitionsHr, "stage", full = TRUE)
    fromMain <- getItems(transitionsHr, "from", full = TRUE)
    toMain <- getItems(transitionsHr, "to", full = TRUE)
    sources <- list(secdforest = c("crop", "past", "forestry"), other = c("crop", "past"))
    years <- getYears(transitionsHr)
    restored <- list()
    for (k in names(sources)) {
      select <- which(stage %in% c("snvRelocation", "peatRewetting", "optimisation") & toMain == k &
                        fromMain %in% sources[[k]])
      flowLr <- toolAggregate(dimSums(transitionsHr[, , select], dim = 3), mapping, from = "cell", to = "cluster")
      shareLr <- pmin(restorationLr[getItems(flowLr, dim = 1), years, k] / setNames(flowLr, k), 1)
      shareLr[is.na(shareLr) | is.infinite(shareLr)] <- 0
      shareHr <- toolAggregate(shareLr, mapping, from = "cluster", to = "cell")[getCells(transitionsHr), , ]
      restored[[k]] <- transitionsHr[, , select] * setNames(shareHr, NULL)
      getNames(restored[[k]]) <- paste(stage[select], paste0(fromMain[select], "Restored"),
                                       getItems(transitionsHr, "fromType", full = TRUE)[select], toMain[select],
                                       getItems(transitionsHr, "toType", full = TRUE)[select], sep = ".")
      transitionsHr[, , select] <- transitionsHr[, , select] * setNames(1 - shareHr, NULL)
    }
    return(mbind(transitionsHr, restored[["secdforest"]], restored[["other"]]))
  }

  # Places the reconverted managed forest of each cluster (cluster x years) in its grid cells, in proportion to
  # the free secondary forest of the grid cells (grid cells x years); the area placed in a year is kept.
  .placeForest <- function(forestLr, freeHr, cellCluster) {
    out <- freeHr * 0
    placed <- as.vector(freeHr[, 1, ]) * 0
    years <- getYears(freeHr)
    for (i in seq_along(years)) {
      free <- pmax(as.vector(freeHr[, i, ]) - placed, 0)
      freeLr <- tapply(free, cellCluster, sum)
      missing <- pmax(as.vector(forestLr[names(freeLr), years[i], ]) - tapply(placed, cellCluster, sum)[names(freeLr)],
                      0)
      share <- pmin(missing / replace(freeLr, freeLr <= 0, 1), 1)
      placed <- placed + free * share[cellCluster]
      out[, i, ] <- placed
    }
    return(out)
  }

  if (input != "magpie") {
    stop("Can only report for input = magpie")
  }
  landHr <- .report("downscaledmagpie")
  splitHr <- .report("downscaledmagpieSplit")
  transitionsHr <- calcOutput("LandTransitionsHighRes", input = input, aggregate = FALSE)[["transitionsHr"]]
  years <- getYears(landHr)
  trackYears <- years[getYears(landHr, as.integer = TRUE) > startYear]
  start <- paste0("y", startYear)
  transitionsHr <- transitionsHr[, trackYears, ]
  mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")
  cfg <- readSource("MagpieRunFiles", subtype = "config")
  reconversion <- isTRUE(as.numeric(cfg$gms$s22_reconversion) == 1)

  # =================================================================================================
  # Classes of the state of the land
  # =================================================================================================

  splitItems <- getItems(splitHr, dim = 3)
  classOf <- stats::setNames(rep(NA_character_, length(splitItems)), splitItems)
  classOf[startsWith(splitItems, "crop.")] <- "cropland_other"
  classOf[startsWith(splitItems, "crop.kbe.")] <- "cropland_2Gbioen"
  classOf[startsWith(splitItems, "past.")] <- "grassland"
  classOf[startsWith(splitItems, "primforest.") | startsWith(splitItems, "secdforest.")] <- "nat_regen_forest"
  classOf[startsWith(splitItems, "forestry.")] <- "forest_planted"
  classOf["forestry.ndc.total"] <- "nat_regen_forest"
  affPlantation <- as.logical(cfg$gms$s32_aff_plantation)
  classOf["forestry.aff.total"] <- if (affPlantation) "forest_planted" else "nat_regen_forest"
  classOf[startsWith(splitItems, "urban.")] <- "built_up"
  stateClasses <- c("cropland_other", "cropland_2Gbioen", "grassland", "nat_regen_forest", "forest_planted",
                    "built_up")
  stateHr <- mbind(lapply(stateClasses, function(class) {
    setNames(dimSums(splitHr[, , names(which(classOf == class))], dim = 3), class)
  }))

  # =================================================================================================
  # Classes of the history of natural land
  # =================================================================================================

  # restored land as source "<source>Restored"; labels "stage.from.to" with the sub-types of other land and
  # forestry as destination, and as source managed forest (timber plantations, other planted forest) as
  # "timber", NDC forest and afforestation as "forestry"
  transitionsHr <- .splitRestoration(transitionsHr, readSource("MagpieFulldataGdx", subtype = "restoration"),
                                     mapping)
  restored <- endsWith(getItems(transitionsHr, "from", full = TRUE), "Restored")
  fromLabel <- getItems(transitionsHr, "fromType", full = TRUE)
  fromLabel[fromLabel == "total"] <- getItems(transitionsHr, "from", full = TRUE)[fromLabel == "total"]
  fromLabel[fromLabel %in% c("plant", "other_planted")] <- "timber"
  fromLabel[fromLabel %in% c("aff", "ndc")] <- "forestry"
  fromLabel[restored] <- paste0(sub("Restored$", "", fromLabel[restored]), "Restored")
  toLabel <- getItems(transitionsHr, "toType", full = TRUE)
  toLabel[toLabel == "total"] <- getItems(transitionsHr, "to", full = TRUE)[toLabel == "total"]
  label <- paste(getItems(transitionsHr, "stage", full = TRUE), fromLabel, toLabel, sep = ".")
  transitionsHr <- mbind(lapply(unique(label), function(l) {
    setNames(dimSums(transitionsHr[, , which(label == l)], dim = 3), l)
  }))
  getSets(transitionsHr, fulldim = FALSE)[3] <- "stage.from.to"

  # origins of each history class; conversions take them only after all other area where they are protected
  originsOf <- list(recon_abn_crop = c("crop", "cropRestored"), recon_abn_grass = c("past", "pastRestored"),
                    recon_forest = c("timber", "timberRestored"))
  protectedOrigins <- c("cropRestored", "pastRestored", if (reconversion) c("crop", "past"))
  types <- c("othernat", "youngsecdf", "secdforest")
  naturalForestry <- c("ndc", if (!affPlantation) "aff")
  oldest <- c("youngsecdf.secdforest", paste(rep(c("othernat", "youngsecdf", "secdforest"),
                                                 each = length(naturalForestry)), naturalForestry, sep = "."))
  outflowOrder <- c(default = "newest", stats::setNames(rep("oldest", length(oldest)), oldest))
  last <- list(default = protectedOrigins, othernat.youngsecdf = character(0), youngsecdf.secdforest = character(0))
  startItems <- c("other.othernat.total", "other.youngsecdf.total")
  stopifnot(all(startItems %in% splitItems))
  startHr <- mbind(setNames(splitHr[, start, startItems], c("othernat", "youngsecdf")),
                   setNames(landHr[, start, "secdforest"], "secdforest"))
  originHr <- toolGriddedLandOrigin(transitionsHr, startHr, types, order = outflowOrder, last = last)

  # the same on the cluster sums: protected land that the model converts itself (MAgPIE protects it only as
  # area of its type in the cluster); the rest of the protected land taken comes from the placement
  transitionsLr <- toolAggregate(transitionsHr, mapping, from = "cell", to = "cluster")
  originLr <- toolGriddedLandOrigin(transitionsLr, toolAggregate(startHr, mapping, from = "cell", to = "cluster"),
                                    types, order = outflowOrder, last = last)
  takenByModel <- attr(originLr, "takenLast")
  takenByPlacement <- pmax(attr(originHr, "takenLast") - takenByModel, 0)

  # the history classes in all years (zero up to startYear), each of its origins in the tracked types
  .origins <- function(x, type, origins) {
    items <- intersect(paste(type, origins, sep = "."), getItems(x, dim = 3))
    area <- new.magpie(getCells(landHr), years, NULL, fill = 0)
    if (length(items) > 0) {
      area[, trackYears, ] <- dimSums(x[, trackYears, items], dim = 3)
    }
    return(area)
  }
  recoveringHr <- mbind(lapply(names(originsOf), function(class) {
    setNames(.origins(originHr, types, originsOf[[class]]), class)
  }))
  recoveringOther <- .origins(originHr, c("othernat", "youngsecdf"), unlist(originsOf))
  recoveringSecdforest <- .origins(originHr, "secdforest", unlist(originsOf))
  otherHr <- setNames(landHr[, , "other"] - recoveringOther, "other")
  stateHr[, , "nat_regen_forest"] <- stateHr[, , "nat_regen_forest"] - recoveringSecdforest

  # reconverted managed forest of the run, on the secondary forest that is not recovering land
  forestLr <- readSource("MagpieFulldataGdx", subtype = "reconversionForest")
  cellCluster <- mapping[["cluster"]][match(getCells(landHr), mapping[["cell"]])]
  freeHr <- setNames(landHr[, , "secdforest"] - recoveringSecdforest, NULL)
  placedHr <- .placeForest(forestLr[unique(cellCluster), years, ], freeHr, cellCluster)
  recoveringHr[, , "recon_forest"] <- recoveringHr[, , "recon_forest"] + placedHr
  stateHr[, , "nat_regen_forest"] <- stateHr[, , "nat_regen_forest"] - placedHr
  missingForest <- dimSums(forestLr[unique(cellCluster), years, ], dim = c(1, 3)) - dimSums(placedHr, dim = c(1, 3))
  if (max(missingForest) > 1e-6) {
    worst <- getYears(missingForest)[which.max(missingForest)]
    toolStatusMessage("warn", paste0("Reconverted managed forest not placed for lack of secondary forest, up to ",
                                     signif(max(missingForest), 3), " Mha (", worst, ")"))
  }

  if (max(takenByModel) > 1e-6) {
    toolStatusMessage("note", paste0("MAgPIE converted protected land (it protects it only as area of its type ",
                                     "in the cluster), up to ", signif(max(takenByModel), 3), " Mha in a time ",
                                     "step (", names(which.max(takenByModel)), ")"))
  }
  if (max(takenByPlacement) > 1e-6) {
    toolStatusMessage("warn", paste0("The grid cell placement converted protected land beyond the model, up to ",
                                     signif(max(takenByPlacement), 3), " Mha in a time step (",
                                     names(which.max(takenByPlacement)), ")"))
  }

  rehabHr <- setNames(stateHr[, , c("cropland_other", "grassland")] * 0, c("rehab_crop", "rehab_grass"))
  x <- mbind(stateHr[, , setdiff(stateClasses, "built_up")], recoveringHr, rehabHr, otherHr,
             stateHr[, , "built_up"])
  getSets(x, fulldim = FALSE) <- c("x.y.iso", "year", "class")
  toolExpectLessDiff(dimSums(x, dim = 3), dimSums(landHr, dim = 3), 1e-6,
                     "The classes cover all land of the grid cells")
  toolExpectTrue(min(x) >= -1e-6, "No class has negative area")
  return(list(x = x, isocountries = FALSE, unit = "Mha", cache = FALSE,
              description = "MAgPIE land by grid cell in the land classes of Bending the Curve 2.0"))
}
