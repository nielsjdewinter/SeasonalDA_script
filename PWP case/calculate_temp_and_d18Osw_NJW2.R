##' Calculate bootstrapped mean values for age, d18O, d13C, and D47 and calculate temperature and d18Ow
##'
apply_calibration_and_d18O_with_uncertainty <- function(data, calib, group) {
  stop("todo: write wrapper function that ties it all together")
}

##' Calculate bootstrapped means for each group in a dataframe
##'
##' @param data A dataframe or tibble that should have columns `outlier`,
##'   `broadid` ("other" for samples), `D47_final`, `d18O_PDB_vit`,
##'   `d13C_PDB_vit`, and the grouping column.
##' @param group The column name in `data` with a character or factor column
##'   that contains the binning information that you want to calculate
##'   bootstrapped averages for.
##' @param d13C The column name in `data` with the carbon isotope values.
##' @param d18O The column name in `data` with the oxygen isotope values.
##' @param D47 The column name in `data` with the clumped isotope values.
##' @param Nsim The number of bootstraps you want to run. Defaults to the
##'   number of rows in `calib`.
##' @author Ilja J. Kocken (minor modifications by Niels J. de Winter)
##' @examples
##' data |> apply_calibration_and_d18O_with_uncertainty(
##'   calib = calib,
##'   group = spec,
##'   d13C = d13C_PDB_vit,
##'   d18O = d18O_PDB_vit,
##'   D47 = D47_final,
##'   Nsim = 5000)
bootstrap_means <- function(data,
                            group,
                            d13C,
                            d18O,
                            D47,
                            Nsim = 1e5) {
  data |>
    # make sure that there are no NAs in group or in your d13C etc!
    group_by({{group}}) |>
    # subset only relevant columns for now unfortunately this gets rid of
    # potentially useful columns. You can left_join them back with
    # distinct(data, id_col1, id_col2)
    select({{group}}, {{d13C}}, {{d18O}}, {{D47}}) |>
    # R magic, search for nesting/unnesting if you want to understand what
    # happens here.
    nest() |>
    # create Nsim bootstrapped copies of the data
    mutate(boot = purrr::map(data,
                      ~ infer::rep_slice_sample(.x,
                                                # we resample using all data
                                                prop = 1,
                                                replace = TRUE,
                                                reps = Nsim))) |>
    # get rid of the raw data, leave only the bootstrapped values
    select(-data) |>
    # calculate summaries for the bootstrapped data, Nsim times
    mutate(summ = purrr::map(boot, ~ .x |>
                              summarize(
                                # here they get these new simpler names
                                d13C = mean({{d13C}}, na.rm = TRUE),
                                d18O = mean({{d18O}}, na.rm = TRUE),
                                D47 = mean({{D47}}, na.rm = TRUE)))) |>
    # get rid of the bootstrapped values
    select(-boot) |>
    # unfold the bootstraps, we're back to a simple tibble now
    unnest(summ)
}

##' Calculate bootstrapped means for each group in a dataframe
##'
##' @param data A dataframe or tibble that should have columns `outlier`,
##'   `broadid` ("other" for samples), `D47_final`, `d18O_PDB_vit`,
##'   `d13C_PDB_vit`, and the grouping column.
##' @param group The column name in `data` with a character or factor column
##'   that contains the binning information that you want to calculate
##'   bootstrapped averages for.
##' @param d13C The column name in `data` with the carbon isotope values.
##' @param d18O The column name in `data` with the oxygen isotope values.
##' @param D47 The column name in `data` with the clumped isotope values.
##' @param Nsim The number of bootstraps you want to run. Defaults to the
##'   number of rows in `calib`.
##' @author Ilja J. Kocken (minor modifications by Niels J. de Winter)
##' @examples
##' data |> apply_calibration_and_d18O_with_uncertainty(
##'   calib = calib,
##'   group = spec,
##'   d13C = d13C_PDB_vit,
##'   d18O = d18O_PDB_vit,
##'   D47 = D47_final,
##'   Nsim = 5000)
bootstrap_medians <- function(data,
                            group,
                            d13C,
                            d18O,
                            D47,
                            Nsim = 1e5) {
  data |>
    # make sure that there are no NAs in group or in your d13C etc!
    group_by({{group}}) |>
    # subset only relevant columns for now unfortunately this gets rid of
    # potentially useful columns. You can left_join them back with
    # distinct(data, id_col1, id_col2)
    select({{group}}, {{d13C}}, {{d18O}}, {{D47}}) |>
    # R magic, search for nesting/unnesting if you want to understand what
    # happens here.
    nest() |>
    # create Nsim bootstrapped copies of the data
    mutate(boot = purrr::map(data,
                      ~ infer::rep_slice_sample(.x,
                                                # we resample using all data
                                                prop = 1,
                                                replace = TRUE,
                                                reps = Nsim))) |>
    # get rid of the raw data, leave only the bootstrapped values
    select(-data) |>
    # calculate summaries for the bootstrapped data, Nsim times
    mutate(summ = purrr::map(boot, ~ .x |>
                              summarize(
                                # here they get these new simpler names
                                d13C = median({{d13C}}, na.rm = TRUE),
                                d18O = median({{d18O}}, na.rm = TRUE),
                                D47 = median({{D47}}, na.rm = TRUE)))) |>
    # get rid of the bootstrapped values
    select(-boot) |>
    # unfold the bootstraps, we're back to a simple tibble now
    unnest(summ)
}

##' Calculate the d18Ow from the d18Occ and temperature
##'
##' according to Kim & O'neil 1997 as modified by Bemis et al., 1998
##'
##' @param d18Occ The oxygen isotope composition of the calcite.
##' @param temperature  The formation temperature (in °C).
##' @return The oxygen isotope composition of the sea water.
##' @author Ilja J. Kocken (minor modifications by Niels J. de Winter)
d18Ow_from_d18O_T <- function(d18Oc, T, mineral = "calcite"){
    if(mineral == "calcite"){
        # Kim, S.-T. and O’Neil, J. R.: Equilibrium and nonequilibrium oxygen isotope effects in synthetic carbonates, Geochimica et Cosmochimica Acta, 61, 3461–3475, https://doi.org/10.1016/S0016-7037(97)00169-5, 1997.
        return(((d18Oc / 1000 + 1) / exp(((18.03 * 10 ^ 3) / T - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92)
    }else if(mineral == "aragonite"){
        # Grossman, E. L. and Ku, T.-L.: Oxygen and carbon isotope fractionation in biogenic aragonite: temperature effects, Chemical Geology: Isotope Geoscience section, 59, 59–74, 1986.
        return(d18Oc - (20.6 - (T - 273.15)) / 4.34 + 0.2)
    }
}

d18Oc_from_T_d18Ow <- function(T, d18Ow, mineral = "calcite"){
    if(mineral == "calcite"){
      # Kim, S.-T. and O’Neil, J. R.: Equilibrium and nonequilibrium oxygen isotope effects in synthetic carbonates, Geochimica et Cosmochimica Acta, 61, 3461–3475, https://doi.org/10.1016/S0016-7037(97)00169-5, 1997.
      return(((exp(((18.03 * 10 ^ 3) / T - 32.42) / 1000) * (((d18Ow - 30.92) / 1.03092) / 1000 + 1)) - 1) * 1000)
    }else if(mineral == "aragonite"){
      # Grossman, E. L. and Ku, T.-L.: Oxygen and carbon isotope fractionation in biogenic aragonite: temperature effects, Chemical Geology: Isotope Geoscience section, 59, 59–74, 1986.
      return(d18Oc + (20.6 - (T - 273.15)) / 4.34 - 0.2)
    }
}

T_from_d18O_d18Ow <- function(d18Oc, d18Ow, mineral = "calcite"){
    if(mineral == "calcite"){
      # Kim, S.-T. and O’Neil, J. R.: Equilibrium and nonequilibrium oxygen isotope effects in synthetic carbonates, Geochimica et Cosmochimica Acta, 61, 3461–3475, https://doi.org/10.1016/S0016-7037(97)00169-5, 1997.
      return((18.03 * 1000) / (log((d18Oc / 1000 + 1) / (((d18Ow - 30.92) / 1.03092) / 1000 + 1)) * 1000 + 32.42))
    }else if(mineral == "aragonite"){
      # Grossman, E. L. and Ku, T.-L.: Oxygen and carbon isotope fractionation in biogenic aragonite: temperature effects, Chemical Geology: Isotope Geoscience section, 59, 59–74, 1986.
      return(20.6 - 4.34 * (d18Oc - d18Ow + 0.2))
    }
}

##' @param calib A dataframe with draws from the bootstrapped (or Bayesian)
##'   temperature regression. Should have columns `slope` and `intercept`,
##'   which are related via `clumpedr::revcal()`.
##' @author Ilja J. Kocken (minor modifications by Niels J. de Winter)
calc_temp_d18Ow <- function(boot, calib, Nsim = NULL, mineral = "calcite") {
  if (is.null(Nsim)) {
    # we simulate the same number of bootstraps for easy combination
    Nsim <- nrow(calib)
  } else {
    # take a subset of the calibration of size Nsim
    calib <- calib[sample(nrow(calib), replace = TRUE, size = Nsim), ]
  }

  boot |>
    # append the slope/intercept pairs of the temperature calibration
    # this is why we made sure that they are Nsim long as well.
    mutate(slope = calib$slope,
           intercept = calib$intercept) |>
    # calculate temperature using the parameters
    # this relies on my clumpedr package
    # https://github.com/isoverse/clumpedr/
    # you can also just copy its revcal function from here:
    # https://github.com/isoverse/clumpedr/blob/master/R/calibration.R#L72
    mutate(temp = clumpedr::revcal(D47, slope = slope, intercept = intercept,
                                   # we have to use ignorecnf because the confidence calculations
                                   # in clumpedr are WRONG!
                                   ignorecnf = TRUE)) |>
    # get rid of calibration intercept and slope
    select(-slope, -intercept) |>
    # calculate d18Ow using the function above
    # we do not take into account potential uncertainty in these parameters,
    # but this is likely nothing.
    mutate(d18Ow = d18Ow_from_d18O_T(d18O, temp + 273.15, mineral))
}

##' Summarize the bootstrapped values into a mean, sd, and the 68% and 95% CIs
##'
##' @param boot Output of `apply_calibration_and_d18O_with_uncertainty()`
##' @param group The group to summarize by.
our_summary <- function(boot, group) {
  boot |>
    group_by({{group}}) |>
    ggdist::median_qi(.exclude = "replicate",
                      .width = c(.68, .95),
                      na.rm = TRUE)
}

Meinicke_D47_to_T <- function(D47){
    # Meinicke, N., Reimi, M. A., Ravelo, A. C., and Meckler, A. N.: Coupled Mg/Ca and Clumped Isotope Measurements Indicate Lack of Substantial Mixed Layer Cooling in the Western Pacific Warm Pool During the Last ∼5 Million Years, Paleoceanography and Paleoclimatology, 36, e2020PA004115, https://doi.org/10.1029/2020PA004115, 2021.
    return(sqrt(0.0397E6 / (D47 - 0.1518)))
}