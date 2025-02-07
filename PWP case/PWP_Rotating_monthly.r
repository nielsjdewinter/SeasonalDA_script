# Data preparation step
# Optimize monthly identification by rotating months
# Plot optimized monthly for all specimens
# Project "PWP case"

require(tidyverse)
require(RColorBrewer)

# Load functions for transferring between d18O, T, d18Ow and D47
# Source functions for bootstrapping temperature and d18Ow
source("calculate_temp_and_d18Osw_NJW2.r")

# Load calibration data
Meinicke_slopes <- read.csv("Meinicke_cal/slope_boot.csv", header = FALSE)
Meinicke_intercepts <- read.csv("Meinicke_cal/inter_boot.csv", header = FALSE)
calib <- data.frame(
    slope = Meinicke_slopes$V1,
    intercept = Meinicke_intercepts$V1
)

# Get data
sample_data_dated <- read.csv("6_PWP_sample_data_dated2_strict.csv")

# ------------------------------------------------------------------------------
# Use rotating seasons script on all specimens using all methods
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SG105
# ------------------------------------------------------------------------------
# Find temperature seasons for SG105
SG105 <- sample_data_dated |>
    filter(Specimen == "SG105", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG105_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG105_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG105 <- SG105 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG105 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG105_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG105_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG105_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG105_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG105_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG105$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG105_months_try$summerT - SG105_months_try$winterT)
SG105 <- SG105 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG105 <- bootstrap_medians(SG105,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG105 <- calc_temp_d18Ow(boot_SG105, calib, mineral = "aragonite")
Summary_data_stats_SG105_strict <- our_summary(boot_SG105, month_score)
Summary_data_stats_SG105_strict <- mutate(
    Summary_data_stats_SG105_strict,
    Specimen = first(SG105$Specimen),
    Species = first(SG105$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG105, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG105_summer <- first(Summary_data_stats_SG105_strict$month_score[which(Summary_data_stats_SG105_strict$temp == max(Summary_data_stats_SG105_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG105_summer - 1) %% 12 + 1
)
SG105$month <- vapply(SG105$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG105$month_score <- vapply(SG105$month_score,
    function(x) (x - SG105_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG105_strict$month <- vapply(Summary_data_stats_SG105_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG105_strict$month_score <- vapply(Summary_data_stats_SG105_strict$month_score,
    function(x) (x - SG105_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# Plot D47 data colored by month
# SG105_D47plot <- ggplot(SG105) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG105$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG105$Specimen), first(SG105$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()


# ------------------------------------------------------------------------------
# SG107
# ------------------------------------------------------------------------------
# Find temperature seasons for SG107
SG107 <- sample_data_dated |>
    filter(Specimen == "SG107", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG107_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG107_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG107 <- SG107 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG107 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG107_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG107_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG107_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG107_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG107_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG107$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG107_months_try$summerT - SG107_months_try$winterT)
SG107 <- SG107 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG107 <- bootstrap_medians(SG107,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG107 <- calc_temp_d18Ow(boot_SG107, calib, mineral = "aragonite")
Summary_data_stats_SG107_strict <- our_summary(boot_SG107, month_score)
Summary_data_stats_SG107_strict <- mutate(
    Summary_data_stats_SG107_strict,
    Specimen = first(SG107$Specimen),
    Species = first(SG107$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG107, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG107_summer <- first(Summary_data_stats_SG107_strict$month_score[which(Summary_data_stats_SG107_strict$temp == max(Summary_data_stats_SG107_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG107_summer - 1) %% 12 + 1
)
SG107$month <- vapply(SG107$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG107$month_score <- vapply(SG107$month_score,
    function(x) (x - SG107_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG107_strict$month <- vapply(Summary_data_stats_SG107_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG107_strict$month_score <- vapply(Summary_data_stats_SG107_strict$month_score,
    function(x) (x - SG107_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG107_D47plot <- ggplot(SG107) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG107$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG107$Specimen), first(SG107$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()


# ------------------------------------------------------------------------------
# SG113
# ------------------------------------------------------------------------------
# Find temperature seasons for SG113
SG113 <- sample_data_dated |>
    filter(Specimen == "SG113", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG113_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG113_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG113 <- SG113 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG113 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG113_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG113_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG113_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG113_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG113_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG113$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG113_months_try$summerT - SG113_months_try$winterT)
SG113 <- SG113 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG113 <- bootstrap_medians(SG113,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG113 <- calc_temp_d18Ow(boot_SG113, calib, mineral = "aragonite")
Summary_data_stats_SG113_strict <- our_summary(boot_SG113, month_score)
Summary_data_stats_SG113_strict <- mutate(
    Summary_data_stats_SG113_strict,
    Specimen = first(SG113$Specimen),
    Species = first(SG113$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG113, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG113_summer <- first(Summary_data_stats_SG113_strict$month_score[which(Summary_data_stats_SG113_strict$temp == max(Summary_data_stats_SG113_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG113_summer - 1) %% 12 + 1
)
SG113$month <- vapply(SG113$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG113$month_score <- vapply(SG113$month_score,
    function(x) (x - SG113_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG113_strict$month <- vapply(Summary_data_stats_SG113_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG113_strict$month_score <- vapply(Summary_data_stats_SG113_strict$month_score,
    function(x) (x - SG113_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG113_D47plot <- ggplot(SG113) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG113$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG113$Specimen), first(SG113$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()


# ------------------------------------------------------------------------------
# SG115
# ------------------------------------------------------------------------------
# Find temperature seasons for SG115
SG115 <- sample_data_dated |>
    filter(Specimen == "SG115", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG115_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG115_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG115 <- SG115 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG115 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG115_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG115_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG115_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG115_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG115_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG115$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG115_months_try$summerT - SG115_months_try$winterT)
SG115 <- SG115 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG115 <- bootstrap_medians(SG115,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG115 <- calc_temp_d18Ow(boot_SG115, calib, mineral = "aragonite")
Summary_data_stats_SG115_strict <- our_summary(boot_SG115, month_score)
Summary_data_stats_SG115_strict <- mutate(
    Summary_data_stats_SG115_strict,
    Specimen = first(SG115$Specimen),
    Species = first(SG115$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG115, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG115_summer <- first(Summary_data_stats_SG115_strict$month_score[which(Summary_data_stats_SG115_strict$temp == max(Summary_data_stats_SG115_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG115_summer - 1) %% 12 + 1
)
SG115$month <- vapply(SG115$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG115$month_score <- vapply(SG115$month_score,
    function(x) (x - SG115_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG115_strict$month <- vapply(Summary_data_stats_SG115_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG115_strict$month_score <- vapply(Summary_data_stats_SG115_strict$month_score,
    function(x) (x - SG115_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG115_D47plot <- ggplot(SG115) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG115$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG115$Specimen), first(SG115$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()


# ------------------------------------------------------------------------------
# SG116
# ------------------------------------------------------------------------------
# Find temperature seasons for SG116
SG116 <- sample_data_dated |>
    filter(Specimen == "SG116", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG116_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG116_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG116 <- SG116 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG116 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG116_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG116_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG116_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG116_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG116_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG116$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG116_months_try$summerT - SG116_months_try$winterT)
SG116 <- SG116 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG116 <- bootstrap_medians(SG116,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG116 <- calc_temp_d18Ow(boot_SG116, calib, mineral = "aragonite")
Summary_data_stats_SG116_strict <- our_summary(boot_SG116, month_score)
Summary_data_stats_SG116_strict <- mutate(
    Summary_data_stats_SG116_strict,
    Specimen = first(SG116$Specimen),
    Species = first(SG116$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG116, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG116_summer <- first(Summary_data_stats_SG116_strict$month_score[which(Summary_data_stats_SG116_strict$temp == max(Summary_data_stats_SG116_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG116_summer - 1) %% 12 + 1
)
SG116$month <- vapply(SG116$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG116$month_score <- vapply(SG116$month_score,
    function(x) (x - SG116_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG116_strict$month <- vapply(Summary_data_stats_SG116_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG116_strict$month_score <- vapply(Summary_data_stats_SG116_strict$month_score,
    function(x) (x - SG116_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG116_D47plot <- ggplot(SG116) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG116$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG116$Specimen), first(SG116$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()


# ------------------------------------------------------------------------------
# SG126
# ------------------------------------------------------------------------------
# Find temperature seasons for SG126
SG126 <- sample_data_dated |>
    filter(Specimen == "SG126", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG126_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG126_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG126 <- SG126 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG126 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG126_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG126_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG126_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG126_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG126_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG126$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG126_months_try$summerT - SG126_months_try$winterT)
SG126 <- SG126 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG126 <- bootstrap_medians(SG126,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG126 <- calc_temp_d18Ow(boot_SG126, calib, mineral = "aragonite")
Summary_data_stats_SG126_strict <- our_summary(boot_SG126, month_score)
Summary_data_stats_SG126_strict <- mutate(
    Summary_data_stats_SG126_strict,
    Specimen = first(SG126$Specimen),
    Species = first(SG126$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG126, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG126_summer <- first(Summary_data_stats_SG126_strict$month_score[which(Summary_data_stats_SG126_strict$temp == max(Summary_data_stats_SG126_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG126_summer - 1) %% 12 + 1
)
SG126$month <- vapply(SG126$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG126$month_score <- vapply(SG126$month_score,
    function(x) (x - SG126_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG126_strict$month <- vapply(Summary_data_stats_SG126_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG126_strict$month_score <- vapply(Summary_data_stats_SG126_strict$month_score,
    function(x) (x - SG126_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG126_D47plot <- ggplot(SG126) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG126$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG126$Specimen), first(SG126$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()

# ------------------------------------------------------------------------------
# SG127
# ------------------------------------------------------------------------------
# Find temperature seasons for SG127
SG127 <- sample_data_dated |>
    filter(Specimen == "SG127", OUTLIER == FALSE, Instrument != "GasBench")

# Find optimal seasons using each Intcorr method

# Keep track of mean (unbootstrapped) temperatures for each combination of seasons
SG127_months_try <- as.data.frame(matrix(ncol = 5 * 12 + 4, nrow = 31))
colnames(SG127_months_try) <- c(paste(rep(c("dc", "dw", "D", "N", "T"), each = 12), "_", rep(1:12, length(c("dw", "dc", "D", "N", "T"))), sep = ""), "winterT", "summerT", "winterdw", "summerdw") 

# Loop through roughly a quarter of a year, day by day
for(i in 1:31){
    # Assign all datapoints to one of monthlys
    SG127 <- SG127 |>
        mutate(month_score = ceiling(((ShellChron_DOY + i) %% 365) / (365 / 12)))
    
    # Calculate mean clumped and d18O values per season
    month_summary <- SG127 |>
        group_by(month_score) |>
        summarize(
            N = n(),
            D47_seasonmean = mean(D47_corrETF_strict, na.rm = TRUE),
            d18Oc_seasonmean = mean(Final_d18O, na.rm = TRUE)
        )

    # Add temperatures and d18Ow values
    month_summary <- month_summary |>
        mutate(
            T_seasonmean = Meinicke_D47_to_T(D47_seasonmean),
            d18Ow_seasonmean = d18Ow_from_d18O_T(d18Oc_seasonmean, T_seasonmean, mineral = "aragonite")
    )
    month_summary$T_seasonmean <- month_summary$T_seasonmean - 273.15 # Convert temperatures to Celcius

    # # VIRTUAL temperatures and d18Ow values
    # month_summary <- month_summary |>
    #     mutate(
    #         T_seasonmean = D47_seasonmean,
    #         d18Ow_seasonmean = d18Oc_seasonmean
    # )

    # Report NA's in case one or more seasons is not represented in the data
    if(nrow(month_summary) < 12){
        for(j in 1:12){
            if(!(j %in% month_summary$month_score)){ # Find missing seasons
                month_summary <- rbind(month_summary,
                    c(j, 0, rep(NA, 12)) # Add NA rows for missing seasons
                )
            }
        }
    }

    # Add means to matrix
    months_long <- pivot_longer(month_summary, cols = colnames(month_summary)[-1])
    months_long <- months_long[order(months_long$name), ]
    SG127_months_try[i, 1:(12 * 5)] <- months_long$value # Add values to store data frame

    # Add summer and winter extremes
    summer <- which.max(months_long$value[which(months_long$name == "T_seasonmean")])
    winter <- which.min(months_long$value[which(months_long$name == "T_seasonmean")])
    SG127_months_try$summerT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == summer)]
    SG127_months_try$winterT[i] <- months_long$value[which(months_long$name == "T_seasonmean" & months_long$month_score == winter)]
    SG127_months_try$summerdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == summer)]
    SG127_months_try$winterdw[i] <- months_long$value[which(months_long$name == "d18Ow_seasonmean" & months_long$month_score == winter)]

    # Remove season classification from data before moving on to the next step
    SG127$month_score <- NULL
}

# Find optimal monthly definition and assign to data
besti <- which.max(SG127_months_try$summerT - SG127_months_try$winterT)
SG127 <- SG127 |>
    mutate(month_score = ceiling(((ShellChron_DOY + besti) %% 365) / (365 / 12)))

# Bootstrap stats for optimal temperature monthlyity
boot_SG127 <- bootstrap_medians(SG127,
                            month_score,
                            Final_d13C,
                            Final_d18O,
                            D47_corrETF_strict,
                            Nsim = 1e5)

boot_SG127 <- calc_temp_d18Ow(boot_SG127, calib, mineral = "aragonite")
Summary_data_stats_SG127_strict <- our_summary(boot_SG127, month_score)
Summary_data_stats_SG127_strict <- mutate(
    Summary_data_stats_SG127_strict,
    Specimen = first(SG127$Specimen),
    Species = first(SG127$Species),
    N = rep(vapply(unique(month_score), function(x) length(which(filter(SG127, !is.na(D47_corrETF_strict))$month_score == x)), FUN.VALUE = 1), 2),
    method = "strict"
)

# Label months, defining the warmest month as August (8th month; in agreement with present-day monthly temperature in the Netherlands and Belgium)
SG127_summer <- first(Summary_data_stats_SG127_strict$month_score[which(Summary_data_stats_SG127_strict$temp == max(Summary_data_stats_SG127_strict$temp, na.rm = TRUE))]) # Define summer as the warmest month
month_labels <- data.frame(
    month_name = month.name,
    month_score = (c(1:12) - 8 + SG127_summer - 1) %% 12 + 1
)
SG127$month <- vapply(SG127$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
SG127$month_score <- vapply(SG127$month_score,
    function(x) (x - SG127_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

Summary_data_stats_SG127_strict$month <- vapply(Summary_data_stats_SG127_strict$month_score,
    function(x) month_labels$month_name[which(month_labels$month_score == x)],
    FUN.VALUE = "January"
) # Apply new month definition
Summary_data_stats_SG127_strict$month_score <- vapply(Summary_data_stats_SG127_strict$month_score,
    function(x) (x - SG127_summer + 8 - 1) %% 12 + 1,
    FUN.VALUE = 1
) # Now update the month score to align summer with August

# month colorscale based on month score
# monthcolors2 <- brewer.pal(12, "Set3")
# names(monthcolors2) <- month.name

# # Plot D47 data colored by month
# SG127_D47plot <- ggplot(SG127) +
#     geom_pointrange(
#         aes(
#             x = ShellChron_Age,
#             y = D47_corrETF_strict,
#             ymin = D47_corrETF_strict - 2 * D47_sd_ext,
#             ymax = D47_corrETF_strict + 2 * D47_sd_ext,
#             color = month
#         )
#     ) +
#     scale_x_continuous("Age [days]",
#         breaks = seq(0, max(SG127$ShellChron_Age), 300)) +
#     scale_y_continuous(Delta[47] ~ "(\u2030"~"I-CDES)") +
#     ggtitle(paste("D47 vs Age", first(SG127$Specimen), first(SG127$Species))) +
#     scale_color_manual(values = monthcolors2) +
#     theme_classic()

# ------------------------------------------------------------------------------
# Combine results from all specimens
Summary_data_stats_all <- rbind(
    Summary_data_stats_SG105_strict,
    Summary_data_stats_SG107_strict,
    Summary_data_stats_SG113_strict,
    Summary_data_stats_SG115_strict,
    Summary_data_stats_SG116_strict,
    Summary_data_stats_SG126_strict,
    Summary_data_stats_SG127_strict
)

# Combine updated monthly assignments for all datapoints
sample_data_dated2 <- rbind(
    SG105,
    SG107,
    SG113,
    SG115,
    SG116,
    SG126,
    SG127
) |>
mutate(
    temp_strict_SD = (Meinicke_D47_to_T(D47_corrETF_strict - D47_SD) - Meinicke_D47_to_T(D47_corrETF_strict + D47_SD)) / 2, # Add uncertainties on T and d18Ow based on SD on per-aliquot D47
    d18Ow_strict_SD = (d18Ow_from_d18O_T(Final_d18O + d18O_sd_ext, temp_strict + temp_strict_SD + 273.15) - d18Ow_from_d18O_T(Final_d18O - d18O_sd_ext, temp_strict - temp_strict_SD + 273.15)) / 2 # Add uncertainties on T and d18Ow based on SD on per-aliquot D47 and d18Oc
)

# Save results of seasonal optimization
write.csv(Summary_data_stats_all, "Rotating_monthly_data_summary.csv", row.names = FALSE)
write.csv(sample_data_dated2, "PWP_proxy_data_months.csv", row.names = FALSE)
