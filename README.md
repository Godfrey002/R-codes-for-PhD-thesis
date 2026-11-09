
**R script for antifungal assay summary statistics, ANOVA, post hoc and bar charts
**

# ============================================================
# Concentrations used: C1, C2, 25, 50, 75 only
# Manual labels only: (a) to (f) centered over 50 and near top edge
# No automatic panel tags anywhere
# Error bars represent STANDARD ERROR (SE), not SD
# Font-safe PDF export using cairo_pdf + sans
# ============================================================
# -----------------------------
# 1. Load packages
# -----------------------------
required <- c(
  "readxl", "dplyr", "ggplot2", "agricolae",
  "patchwork", "purrr", "tibble"
)

new <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl)
library(dplyr)
library(ggplot2)
library(agricolae)
library(patchwork)
library(purrr)
library(tibble)

# -----------------------------
# 2. Load data
# -----------------------------
data <- read_excel("file_name.xlsx")
names(data) <- trimws(names(data))
required_cols <- c("organism", "concentration", "rep", "inhibition_mm")
if (!all(required_cols %in% names(data))) {
  stop(
    "Your file must contain these columns: organism, concentration, rep, inhibition_mm.\n",
    "Current columns are: ", paste(names(data), collapse = ", ")
  )
}
# -----------------------------
# 3. Clean + enforce order
# -----------------------------
conc_levels <- c("C1","C2", "25", "50", "75")
data <- data %>%
  mutate(
    organism = trimws(as.character(organism)),
    concentration = as.character(concentration),
    concentration = trimws(concentration),
    concentration = gsub("µg/ml|ug/ml|mg/ml", "", concentration, ignore.case = TRUE),
    concentration = gsub("\\s+", "", concentration),

    inhibition_mm = as.numeric(inhibition_mm)
  )

# Keep only valid concentration levels
data <- data %>%
  filter(concentration %in% conc_levels)
# Convert to ordered factor
data$concentration <- factor(data$concentration, levels = conc_levels)
if (nrow(data) == 0) {
  stop(
    "No valid concentration values were found.\n",
    "Allowed values are: C1,C2, 25, 50, 75"
  )
}
# Keep organism order as it appears in the file
org_levels <- unique(data$organism)
# -----------------------------
# 4. Analysis function
# -----------------------------
analyze_org <- function(df, panel_letter) {
  df <- df %>% mutate(concentration = factor(concentration, levels = conc_levels))
  # ANOVA
  model <- aov(inhibition_mm ~ concentration, data = df)
  # Duncan's multiple range test
  duncan <- agricolae::duncan.test(model, "concentration", group = TRUE, console = FALSE)
  letters <- duncan$groups %>%
    as.data.frame() %>%
    rownames_to_column("concentration") %>%
    as_tibble() %>%
    mutate(concentration = factor(concentration, levels = conc_levels))

  # Summary stats using STANDARD ERROR (SE)
  summary_df <- df %>%
    group_by(concentration) %>%
    summarise(
      n = sum(!is.na(inhibition_mm)),
      mean = mean(inhibition_mm, na.rm = TRUE),
      sd = sd(inhibition_mm, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(concentration = factor(concentration, levels = conc_levels)) %>%
    left_join(letters, by = "concentration") %>%
    arrange(concentration)
  # Consistent y-limit across all panels
  y_limit <- max(data$inhibition_mm, na.rm = TRUE) * 1.25
  # Duncan letters above each bar
  # Panel label centered over 50 but positioned near the top edge
  letter_offset <- max(1.2, 0.05 * y_limit)
  panel_y <- y_limit * 0.93

  p <- ggplot(summary_df, aes(x = concentration, y = mean, fill = concentration)) +
    geom_col(
      width = 0.6,
      color = "grey30",
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.15,
      linewidth = 0.5,
      color = "black"
    ) +
    geom_text(
      aes(y = mean + se + letter_offset, label = groups),
      size = 4.0,
      fontface = "bold",
      family = "sans",
      vjust = 0
    ) +
    annotate(
      "text",
      x = factor("50", levels = conc_levels),
      y = panel_y,
      label = paste0("(", panel_letter, ")"),
      size = 4.6,
      fontface = "bold",
      family = "sans",
      vjust = 0
    ) +
   scale_fill_manual(
  values = c(
   
    "C1" =  "#000000",  # Black
    "C2" = "#999999",  # Grey
    "25"      = "#56B4E9",  # Sky Blue
    "50"      = "#009E73",  # Bluish Green
    "75"      = "#E69F00"   # Orange
  ),
  drop = FALSE
) +
    scale_x_discrete(limits = conc_levels, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = expression(paste("Concentration (", mu, "g/mL)")),
      y = "Mycelial radial growth inhibition (mm)"
    ) +
    theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 8)),
      axis.title.y = element_text(size = 10, margin = margin(r = 12)),
      axis.text = element_text(size = 9),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      axis.ticks.length = unit(2.5, "mm"),
      plot.margin = margin(8, 12, 18, 16)
    )

  list(plot = p, model = model, summary = summary_df, duncan = duncan)
}

# -----------------------------
# 5. Generate plots
# -----------------------------
results <- map2(org_levels, letters[seq_along(org_levels)], ~ analyze_org(filter(data, organism == .x), .y))
plots <- map(results, "plot")
# -----------------------------
# 6. Combine plots
#    No automatic tags. Only manual (a–f) centered over 50.
# -----------------------------
final_plot <- wrap_plots(plots, ncol = 2)
# -----------------------------
# 7. Export
# -----------------------------
ggsave(
  "Figure1_FINAL.tiff",
  final_plot,
  dpi = 600,
  width = 17,
  height = 22,
  units = "cm",
  compression = "lzw"
)
ggsave(
  "Figure1_FINAL.pdf",
  final_plot,
  width = 17,
  height = 22,
  units = "cm",
  device = cairo_pdf
)
# grDevices::embedFonts("Figure1_FINAL.pdf", outfile = "Figure1_FINAL_embedded.pdf")
# -----------------------------
# 8. Save tables
# -----------------------------
anova_tables <- map2_df(results, org_levels, function(res, org_name) {
  out <- as.data.frame(summary(res$model)[[1]])
  out$Source <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c(ncol(out), 1:(ncol(out) - 1))]
  out$Organism <- org_name
  out
})
summary_tables <- map2_df(results, org_levels, function(res, org_name) {
  out <- res$summary
  out$Organism <- org_name
  out
})
write.csv(anova_tables, "ANOVA_table_all_organisms.csv", row.names = FALSE)
write.csv(summary_tables, "Summary_stats_all_organisms.csv", row.names = FALSE)
cat("✔ Figure generated using only Control, 25, 50, 75\n")
cat("✔ Manual (a–f) labels centered over 50 and near the top edge\n")
cat("✔ Error bars represent standard error (SE)\n")
cat("✔ PDF exported with cairo_pdf and sans font family\n")



R script for germination assay summary statistics, ANOVA, post hoc and bar charts


# ========================================================================
# EFFECT OF NANOPARTICLES ON SEED GERMINATION PARAMETERS
#
# TREATMENT CODING:
# Control         = Control
# Zinc Oxide NP      = A
# Iron oxide NP      = B      
# Zinc/Iron Oxide NP = C
#
# PARAMETERS:
# 1. Dry weight (g)
# 2. Fresh weight (g)
# 3. Percentage germination (%)
# 4. Root length (mm)
# 5. Shoot length (mm)
# 6. Root-Shoot Ratio
# 7. Seedling length (mm)
# 8. Seedling vigor index (SVI)
#
# STATISTICS:
# - Mean ± SE
# - One-way ANOVA
# - Tukey HSD
# - Significance letters
# - Arcsine-square-root transformation for percentage germination
#
# FIGURES:
# - Individual publication-ready figures
# - Combined Figure 1(a-h)
# - CVD-friendly colours
# - theme_classic()
# - White background
# - No gridlines
# - Cairo PDF + sans font
# - TIFF 600 dpi
# - 1.25 × maximum-value y-axis limit
# - Same treatment order across all panels ========================================================================
# 1. INSTALL REQUIRED PACKAGES
========================================================================
required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "multcompView",
  "patchwork",
  "openxlsx",
  "grid"
)
installed_packages <- rownames(
  installed.packages()
)
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg)
  }
}
# ========================================================================
# 2. LOAD PACKAGES
#========================================================================
library(readxl)
library(dplyr)
library(ggplot2)
library(multcompView)
library(patchwork)
library(openxlsx)
library(grid)
# ========================================================================
# 3. SELECT EXCEL FILE
# ========================================================================
#
# Select:
# Germination raw.xlsx
========================================================================
file_path <- file.choose()
========================================================================
# 4. CREATE OUTPUT FOLDERS
# ========================================================================
output_dir <- "NP_Seed_Germination_Results"
dir.create(
  output_dir,
  showWarnings = FALSE
)
dir.create(
  file.path(
    output_dir,
    "Individual_Figures"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  file.path(
    output_dir,
    "Statistics"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)
# ========================================================================
# 5. READ EXCEL FILE
# ========================================================================
#
# Your workbook uses Sheet1.
#
# ========================================================================
raw_data <- read_excel(
  file_path,
  sheet = "Sheet1"
)
# ========================================================================
# 6. CHECK REQUIRED COLUMNS
# ========================================================================
required_columns <- c(
  "Replicate",
  "Treatment",
  "Dry_weight_g",
  "Fresh_weight_g",
  "Percentage_germination",
  "Root_length_mm",
  "Shoot_length_mm",
  "Root_shoot_ratio",
  "Seedling_length_mm",
  "SVI"

)
missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)
if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required columns are missing:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}

# ========================================================================
# 7. RECODE TREATMENT GROUPS
# ========================================================================
#
# Original Excel names:
#
# Control
# Zinc Oxide NP
# Iron oxide NP 
# Zinc/Iron Oxide NP
#
# Statistical/graphical labels:
#
# Control
# A
# B
# C
#
# ========================================================================
data <- raw_data %>%
  mutate(
    Treatment_original =
      trimws(
        as.character(
          Treatment
        )
      ),

    Treatment = case_when(
      Treatment_original %in% c(
        "Control",
        "control"
      ) ~ "Control",
      Treatment_original %in% c(
        "Zinc Oxide NP",
        "A"
      ) ~ "A",
      Treatment_original %in% c(
        "Iron oxide NP",
        "B"
      ) ~ "B",
      Treatment_original %in% c(
        "Zinc/Iron Oxide NP",
                "C"
      ) ~ "C",

      TRUE ~ NA_character_
    ),
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "A",
        "B",
        "C"
      )
    )

  ) %>%

  mutate(
    across(
      all_of(
        c(

          "Dry_weight_g",
          "Fresh_weight_g",
          "Percentage_germination",
          "Root_length_mm",
          "Shoot_length_mm",
          "Root_shoot_ratio",
          "Seedling_length_mm",
          "SVI"
        )
      ),
      as.numeric
    )
  )
# ========================================================================
# 8. VERIFY TREATMENT NAMES
# ========================================================================
if (any(is.na(data$Treatment))) {
  invalid_treatments <-
    unique(
      data$Treatment_original[
        is.na(data$Treatment)
      ]
    )
  stop(
    paste(
      "The following treatment names were not recognised:",
      paste(
        invalid_treatments,
        collapse = ", "
      ),
      "\n\nExpected Excel treatment names include:",
      "\nControl",
      "\nZinc Oxide NP",
      "\nIron oxide NP ",
      "\nZinc/Iron Oxide NP"
    )
  )

}
# ========================================================================
# 9. CHECK PERCENTAGE GERMINATION
# ========================================================================
invalid_percentage <- data %>%
  filter(
    !is.na(
      Percentage_germination
    ),
    Percentage_germination < 0 |
      Percentage_germination > 100
  )

if (
  nrow(invalid_percentage) > 0
) {
  stop(
    "Percentage_germination contains values outside 0-100. Please check the Excel data."
  )
}

# ========================================================================
# 10. PARAMETERS
# ========================================================================
parameters <- c(
  "Dry_weight_g",
  "Fresh_weight_g",
  "Percentage_germination",
  "Root_length_mm",
  "Shoot_length_mm",
  "Root_shoot_ratio",
  "Seedling_length_mm",
  "SVI"
)


# ========================================================================
# 11. EXACT Y-AXIS LABELS
# ========================================================================
y_labels <- c(
  Dry_weight_g =
    "Dry weight (g)",
  Fresh_weight_g =
    "Fresh weight (g)",
  Percentage_germination =
    "Percentage germination (%)",
  Root_length_mm =
    "Root length (mm)",
  Shoot_length_mm =
    "Shoot length (mm)",
  Root_shoot_ratio =
    "Root-Shoot Ratio",
  Seedling_length_mm =
    "Seedling length (mm)",
  SVI =
    "SVI"
)
# ========================================================================
# 12. PANEL LABELS
# ========================================================================
panel_labels <- c(
  Dry_weight_g = "(a)",
  Fresh_weight_g = "(b)",
  Percentage_germination = "(c)",
  Root_length_mm = "(d)",
  Shoot_length_mm = "(e)",
  Root_shoot_ratio = "(f)",
  Seedling_length_mm = "(g)",
  SVI = "(h)"
)


# ========================================================================
# 13. TREATMENT LABELS
# ========================================================================
#
# These are the labels that will appear on the x-axis and in
# statistical output.
#
# ========================================================================
treatment_labels <- c(
  "Control",
  "A",
  "B",
  "C"
)
# ========================================================================
# 14. CVD-FRIENDLY COLOUR PALETTE
# ========================================================================
treatment_colors <- c(
  "Control" = "#000000",
  "A" = "#0072B2",
  "B" = "#D55E00",
  "C" = "#009E73"
)
# ========================================================================
# 15. STANDARD ERROR FUNCTION
# ========================================================================
calculate_se <- function(x) {
  x <- x[
    !is.na(x)
  ]
  if (
    length(x) <= 1
  ) {
    return(NA_real_)
  }
  sd(x) /
    sqrt(
      length(x)
    )
}
# ========================================================================
# 16. DESCRIPTIVE STATISTICS FUNCTION
# ========================================================================
#
# Descriptive statistics are calculated from the ORIGINAL data.
#
# ========================================================================

calculate_descriptive <- function(
  df,
  parameter
) {
  df %>%
    group_by(
      Treatment
    ) %>%
    summarise(
      n =
        sum(
          !is.na(
            .data[[parameter]]
          )
        ),
      Mean =
        mean(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SD =
        sd(
          .data[[parameter]],
          na.rm = TRUE
        ),
      SE =
        calculate_se(
          .data[[parameter]]
        ),
      Minimum =
        min(
          .data[[parameter]],
          na.rm = TRUE
        ),
      Maximum =
        max(
          .data[[parameter]],
          na.rm = TRUE
        ),
      .groups = "drop"

    ) %>%
    mutate(
      Parameter =
        parameter
    ) %>%
    select(
      Parameter,
      Treatment,
      n,
      Mean,
      SD,
      SE,
      Minimum,

      Maximum
    )
}
# ========================================================================
# 17. ANOVA FUNCTION
# ========================================================================
#
# Percentage germination:
#
#   Percentage
#       ↓
#   proportion
#       ↓
#   arcsine-square-root transformation
#       ↓
#   ANOVA
#
# All other parameters:
#
#   Original values
#       ↓
#   ANOVA
#
# ========================================================================

perform_anova <- function(
  df,
  parameter
) {
  # ----------------------------------------------------------------------
  # PREPARE DATA
  # ----------------------------------------------------------------------
  analysis_data <- df %>%
    select(

      Treatment,
      Value =
        all_of(parameter)
    ) %>%
    filter(

      !is.na(Treatment),

      !is.na(Value)
    )

  # ----------------------------------------------------------------------
  # TRANSFORM ONLY PERCENTAGE GERMINATION
  # ----------------------------------------------------------------------
  transformation_used <- "None"
  if (
    parameter ==
    "Percentage_germination"
  ) {
    analysis_data <-
      analysis_data %>%
      mutate(
        Value =
          Value / 100,
        Value =
          asin(
            sqrt(Value)
          )

      )
    transformation_used <-
      "Arcsine-square-root transformation"

  }

  # ----------------------------------------------------------------------
  # ANOVA MODEL
  # ----------------------------------------------------------------------
  model <- aov(
    Value ~ Treatment,
    data = analysis_data
  )

  # ----------------------------------------------------------------------
  # ANOVA TABLE
  # ----------------------------------------------------------------------
  anova_table <-
    summary(model)[[1]]
  # ----------------------------------------------------------------------
  # ETA-SQUARED
  # ----------------------------------------------------------------------
  SS_treatment <-
    anova_table[
      "Treatment",
      "Sum Sq"
    ]


  SS_total <-

    sum(
      anova_table[
        ,
        "Sum Sq"
      ]
    )
  eta_squared <-

    SS_treatment /
    SS_total
  # ----------------------------------------------------------------------
  # PUBLICATION TABLE
  # ----------------------------------------------------------------------
  publication_table <- data.frame(
    Parameter =
      parameter,
    Analysis_scale =
      transformation_used,
    Source =
      rownames(
        anova_table
      ),

    df =
      anova_table[
        ,
        "Df"
      ],

    Sum_Squares =
      anova_table[
        ,
        "Sum Sq"
      ],

    Mean_Square =
      anova_table[
        ,
        "Mean Sq"
      ],

    F_value =
      anova_table[
        ,
        "F value"
      ],

    p_value =
      anova_table[
        ,
        "Pr(>F)"
      ],

    Eta_squared =
      NA_real_
  )

  publication_table$Eta_squared[

    publication_table$Source ==
      "Treatment"

  ] <- eta_squared

  list(

    model =
      model,

    table =
      publication_table
  )

}
# ========================================================================
# 18. TUKEY HSD + SIGNIFICANCE LETTERS
# ========================================================================

perform_tukey <- function(
  model,
  parameter
) {

  # ----------------------------------------------------------------------
  # TUKEY HSD
  # ----------------------------------------------------------------------

  tukey_result <-
  TukeyHSD(
    model,
    "Treatment",
    conf.level = 0.95
  )
  # ----------------------------------------------------------------------
  # TUKEY TABLE
  # ----------------------------------------------------------------------
  tukey_table <-
    as.data.frame(

      tukey_result$Treatment
    )

  tukey_table$Comparison <-
    rownames(
      tukey_table
    )

  rownames(
    tukey_table
  ) <- NULL


  tukey_table <-

    tukey_table %>%

    select(

      Comparison,

      everything()

    ) %>%

    mutate(

      Parameter =
        parameter

    ) %>%

    select(

      Parameter,

      everything()
    )

  # ----------------------------------------------------------------------
  # SIGNIFICANCE LETTERS
  # ----------------------------------------------------------------------

  letter_object <-

    multcompLetters4(

      model,

      tukey_result

    )

  letter_values <-

    letter_object$
    Treatment$
    Letters


  letters_table <-

    data.frame(

      Treatment =
        names(
          letter_values
        ),

      Letters =
        unname(
          letter_values
        ),

      stringsAsFactors =
        FALSE

    )


  letters_table$Treatment <-

    factor(

      letters_table$Treatment,

      levels = c(

        "Control",
        "A",
        "B",
        "C"

      )

    )


  letters_table <-

    letters_table %>%

    arrange(
      Treatment
    )


  list(

    tukey =
      tukey_table,

    letters =
      letters_table

  )

}


# ========================================================================
# 19. STORAGE OBJECTS
# ========================================================================

all_descriptive <- list()

all_anova <- list()

all_tukey <- list()

all_letters <- list()

all_models <- list()

all_plots <- list()


# ========================================================================
# 20. MAIN ANALYSIS LOOP
# ========================================================================

for (
  parameter in parameters
) {


  message(
    "Analysing: ",
    parameter
  )


  # ======================================================================
  # 20A. DESCRIPTIVE STATISTICS
  # ======================================================================

  descriptive <-

    calculate_descriptive(

      data,

      parameter

    )


  all_descriptive[[parameter]] <-
    descriptive


  # ======================================================================
  # 20B. ANOVA
  # ======================================================================

  anova_result <-

    perform_anova(

      data,

      parameter

    )


  model <-
    anova_result$model


  all_models[[parameter]] <-
    model


  all_anova[[parameter]] <-
    anova_result$table


  # ======================================================================
  # 20C. TUKEY + SIGNIFICANCE LETTERS
  # ======================================================================

  tukey_result <-

    perform_tukey(

      model,

      parameter

    )


  all_tukey[[parameter]] <-
    tukey_result$tukey


  all_letters[[parameter]] <-
    tukey_result$letters


  # ======================================================================
  # 20D. PREPARE PLOT DATA
  # ======================================================================

  plot_data <-

    descriptive %>%

    left_join(

      tukey_result$letters,

      by = "Treatment"

    )


  # ======================================================================
  # 20E. CALCULATE Y-AXIS LIMIT
  # ======================================================================
  #
  # The graph always uses ORIGINAL data.
  #
  # For percentage germination, the figure remains in %.
  #
  # ======================================================================

  max_value <-

    max(

      plot_data$Mean +
        plot_data$SE,

      na.rm = TRUE

    )


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <-

      max(

        plot_data$Mean,

        na.rm = TRUE

      )

  }


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <- 1

  }


  # ======================================================================
  # 1.25 × MAXIMUM VALUE
  # ======================================================================

  y_max <-

    1.25 *
    max_value


  # ======================================================================
  # 20F. SIGNIFICANCE LETTER POSITION
  # ======================================================================

  letter_offset <-

    0.035 *
    y_max


  plot_data <-

    plot_data %>%

    mutate(

      Letter_y =

        Mean +

        SE +

        letter_offset

    )


  # ======================================================================
  # 20G. KEEP LETTERS INSIDE PANEL
  # ======================================================================

  plot_data$Letter_y <-

    pmin(

      plot_data$Letter_y,

      0.95 *
        y_max

    )


  # ======================================================================
  # 20H. CREATE PLOT
  # ======================================================================

  p <-

    ggplot(

      plot_data,

      aes(

        x = Treatment,

        y = Mean,

        fill = Treatment

      )

    ) +


    # --------------------------------------------------------------------
    # BARS
    # --------------------------------------------------------------------

    geom_col(

      width = 0.68,

      colour = "black",

      linewidth = 0.5

    ) +


    # --------------------------------------------------------------------
    # ERROR BARS = MEAN ± SE
    # --------------------------------------------------------------------

    geom_errorbar(

      aes(

        ymin = pmax(
          0,
          Mean - SE
        ),

        ymax =
          Mean + SE

      ),

      width = 0.15,

      linewidth = 0.55,

      colour = "black"

    ) +


    # --------------------------------------------------------------------
    # SIGNIFICANCE LETTERS
    # --------------------------------------------------------------------

    geom_text(

      aes(

        x = Treatment,

        y = Letter_y,

        label = Letters

      ),

      family = "sans",

      size = 3.0,

      fontface = "bold",

      colour = "black",

      vjust = 0

    ) +


    # --------------------------------------------------------------------
    # PANEL LABEL
    # TOP-LEFT INSIDE PANEL
    # --------------------------------------------------------------------

    annotate(

      "text",

      x = 0.65,

      y = 0.985 * y_max,

      label =
        panel_labels[parameter],

      hjust = 0,

      vjust = 1,

      family = "sans",

      fontface = "bold",

      size = 3.2

    ) +


    # --------------------------------------------------------------------
    # CVD-FRIENDLY COLOURS
    # --------------------------------------------------------------------

    scale_fill_manual(

      values =
        treatment_colors,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # Y-AXIS
    # --------------------------------------------------------------------

    scale_y_continuous(

      limits = c(

        0,

        y_max

      ),

      expand =
        expansion(
          mult =
            c(
              0,
              0
            )
        ),

      n.breaks = 5

    ) +


    # --------------------------------------------------------------------
    # SAME X-AXIS ACROSS ALL PANELS
    # --------------------------------------------------------------------

    scale_x_discrete(

      limits = c(

        "Control",
        "A",
        "B",
        "C"

      ),

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # AXIS LABELS
    # --------------------------------------------------------------------

    labs(

      x = "Treatment",

      y =
        y_labels[parameter]

    ) +


    # --------------------------------------------------------------------
    # PUBLICATION THEME
    # --------------------------------------------------------------------

    theme_classic(

      base_size = 9,

      base_family = "sans"

    ) +


    theme(

      # ------------------------------------------------------------------
      # FONT
      # ------------------------------------------------------------------

      text =
        element_text(

          family = "sans",

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.x =
        element_text(

          size = 8,

          colour = "black",

          angle = 0,

          hjust = 0.5,
    margin = margin(t = 3)

      ),


      # ------------------------------------------------------------------
      # Y-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.y =
        element_text(

          size = 8,

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.x =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              t = 5
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.y =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              r = 5
            )

        ),


      # ------------------------------------------------------------------
      # NO GRIDLINES
      # ------------------------------------------------------------------

      panel.grid.major =
        element_blank(),

      panel.grid.minor =
        element_blank(),


      # ------------------------------------------------------------------
      # NO LEGEND
      # ------------------------------------------------------------------

      legend.position =
        "none",


      # ------------------------------------------------------------------
      # BREATHING ROOM
      # ------------------------------------------------------------------

      plot.margin =
        margin(

          t = 5,

          r = 5,

          b = 5,

          l = 5,

          unit = "mm"

        ),


      # ------------------------------------------------------------------
      # AXIS LINE WIDTH
      # ------------------------------------------------------------------

      axis.line =
        element_line(

          linewidth = 0.5,

          colour = "black"

        ),

      axis.ticks =
        element_line(

          linewidth = 0.5,

          colour = "black"

        )

    )


  # ======================================================================
  # 20I. STORE PLOT
  # ======================================================================

  all_plots[[parameter]] <-
    p


  # ======================================================================
  # 20J. SAVE INDIVIDUAL PDF
  # ======================================================================

  ggsave(

    filename = file.path(

      output_dir,

      "Individual_Figures",

      paste0(

        parameter,

        ".pdf"

      )

    ),

    plot = p,

    device = cairo_pdf,

    family = "sans",

    width = 105,

    height = 65,

    units = "mm"

  )


  # ======================================================================
  # 20K. SAVE INDIVIDUAL TIFF — 600 DPI
  # ======================================================================

  ggsave(

    filename = file.path(

      output_dir,

      "Individual_Figures",

      paste0(

        parameter,

        "_600dpi.tiff"

      )

    ),

    plot = p,

    device = "tiff",

    dpi = 600,

    width = 105,

    height = 65,

    units = "mm",

    compression = "lzw"

  )

}


# ========================================================================
# 21. COMBINE EIGHT PANELS
# ========================================================================

combined_figure <-

  all_plots[["Dry_weight_g"]] +

  all_plots[["Fresh_weight_g"]] +

  all_plots[["Percentage_germination"]] +

  all_plots[["Root_length_mm"]] +

  all_plots[["Shoot_length_mm"]] +

  all_plots[["Root_shoot_ratio"]] +

  all_plots[["Seedling_length_mm"]] +

  all_plots[["SVI"]] +

  plot_layout(

    ncol = 3,

    nrow = 3,

    guides = "keep"

  )


# ========================================================================
# 22. MULTI-PANEL SPACING
# ========================================================================

combined_figure <-

  combined_figure &

  theme(

    legend.position =
      "none",

    panel.spacing =
      unit(
        10,
        "mm"
      ),

    plot.margin =
      margin(

        t = 4,

        r = 4,

        b = 4,

        l = 4,

        unit = "mm"

      )

  )


# ========================================================================
# 23. SAVE COMBINED PDF
# ========================================================================

ggsave(

  filename = file.path(

    output_dir,

    "Figure_1_a_to_h.pdf"

  ),

  plot =
    combined_figure,

  device =
    cairo_pdf,

  family =
    "sans",

  width =
    180,

  height =
    190,

  units =
    "mm"

)


# ========================================================================
# 24. SAVE COMBINED TIFF — 600 DPI
# ========================================================================

ggsave(

  filename = file.path(

    output_dir,

    "Figure_1_a_to_h_600dpi.tiff"

  ),

  plot =
    combined_figure,

  device =
    "tiff",

  dpi = 600,

  width =
    180,

  height =
    190,

  units =
    "mm",

  compression =
    "lzw"

)


# ========================================================================
# 25. COMBINE DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  bind_rows(
    all_descriptive
  )


# ========================================================================
# 26. COMBINE ANOVA RESULTS
# ========================================================================

anova_final <-

  bind_rows(
    all_anova
  )


# ========================================================================
# 27. COMBINE TUKEY RESULTS
# ========================================================================

tukey_final <-

  bind_rows(
    all_tukey
  )


# ========================================================================
# 28. COMBINE SIGNIFICANCE LETTERS
# ========================================================================

letters_final <-

  bind_rows(

    lapply(

      names(all_letters),

      function(parameter_name) {

        all_letters[[parameter_name]] %>%

          mutate(

            Parameter =
              parameter_name

          ) %>%

          select(

            Parameter,

            everything()

          )

      }

    )

  )


# ========================================================================
# 29. ROUND DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  descriptive_final %>%

  mutate(

    Mean =
      round(
        Mean,
        4
      ),

    SD =
      round(
        SD,
        4
      ),

    SE =
      round(
        SE,
        4
      ),

    Minimum =
      round(
        Minimum,
        4
      ),

    Maximum =
      round(
        Maximum,
        4
      )

  )


# ========================================================================
# 30. ROUND ANOVA RESULTS
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Sum_Squares =
      round(
        Sum_Squares,
        4
      ),

    Mean_Square =
      round(
        Mean_Square,
        4
      ),

    F_value =
      round(
        F_value,
        4
      ),

    p_value =
      signif(
        p_value,
        4
      ),

    Eta_squared =
      round(
        Eta_squared,
        4
      )

  )


# ========================================================================
# 31. SIGNIFICANCE SYMBOLS
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

   Significance = case_when(

  is.na(p_value)
  ~ "",

  p_value < 0.05
  ~ "Significant",

  TRUE
  ~ "Not significant"

    )

  )


# ========================================================================
# 32. PUBLICATION-READY ANOVA TABLE
# ========================================================================

anova_publication <-

  anova_final %>%

  select(

    Parameter,

    Analysis_scale,

    Source,

    df,

    Sum_Squares,

    Mean_Square,

    F_value,

    p_value,

    Significance,

    Eta_squared

  )


# ========================================================================
# 33. WRITE STATISTICAL RESULTS TO EXCEL
# ========================================================================

write.xlsx(

  list(

    Descriptive_Statistics =
      descriptive_final,

    ANOVA_Publication =
      anova_publication,

    Tukey_HSD =
      tukey_final,

    Significance_Letters =
      letters_final

  ),

  file = file.path(

    output_dir,

    "Statistical_Results.xlsx"

  ),

  overwrite = TRUE

)


# ========================================================================
# 34. SAVE CSV TABLES
# ========================================================================

write.csv(

  descriptive_final,

  file.path(

    output_dir,

    "Statistics",

    "Descriptive_Statistics.csv"

  ),

  row.names = FALSE

)


write.csv(

  anova_publication,

  file.path(

    output_dir,

    "Statistics",

    "ANOVA_Publication_Table.csv"

  ),

  row.names = FALSE

)


write.csv(

  tukey_final,

  file.path(

    output_dir,

    "Statistics",

    "Tukey_HSD.csv"

  ),

  row.names = FALSE

)


write.csv(

  letters_final,

  file.path(

    output_dir,

    "Statistics",

    "Significance_Letters.csv"

  ),

  row.names = FALSE

)


# ========================================================================
# 35. PRINT COMPLETION MESSAGE
# ========================================================================

cat("\n")
cat("====================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("====================================================\n\n")


cat(

  "Output directory:\n",

  normalizePath(
    output_dir
  ),

  "\n\n"

)


cat(

  "Combined PDF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_h.pdf"

  ),

  "\n\n"

)


cat(

  "Combined 600-dpi TIFF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_h_600dpi.tiff"

  ),

  "\n\n"

)


cat(

  "Statistical workbook:\n",

  file.path(

    output_dir,

    "Statistical_Results.xlsx"

  ),

  "\n\n"

)


# ========================================================================
# 36. PRINT ANOVA TABLE
# ========================================================================

cat("\n")
cat("====================================================\n")
cat("PUBLICATION-READY ANOVA TABLE\n")
cat("====================================================\n\n")

print(
  anova_publication
)


# ========================================================================
# 37. PRINT DESCRIPTIVE STATISTICS
# ========================================================================

cat("\n")
cat("====================================================\n")
cat("DESCRIPTIVE STATISTICS\n")
cat("====================================================\n\n")

print(
  descriptive_final
)


# ========================================================================
# 38. PRINT TRANSFORMATION INFORMATION
# ========================================================================

cat("\n")
cat("====================================================\n")
cat("PERCENTAGE GERMINATION TRANSFORMATION\n")
cat("====================================================\n\n")

cat(
  "Percentage germination was transformed as:\n"
)

cat(
  "asin(sqrt(percentage / 100))\n\n"
)

cat(
  "ANOVA and Tukey HSD use the transformed values.\n"
)

cat(
  "The figure displays the original percentage scale.\n\n"
)


# 




R script for growth promotion assay summary statistics, ANOVA, post hoc and bar charts

========================================================================
# TREATMENT CODING
# ========================================================================
#
# Control 1     = C1
# Control 2     = C2
# Zinc Oxide NP    = A
# Iron oxide NP      = B
# Zinc/Iron Oxide NP   = C
#
# ========================================================================
# REPLICATIONS
# ========================================================================
#
# 8 replications per treatment
#
# ========================================================================
# PARAMETERS ANALYSED
# ========================================================================
#
# 1. Plant height (mm)
# 2. No. of leaves
# 3. Stem girth (mm)
# 4. Leaf length (mm)
# 5. Leaf width (mm)
# 6. Leaf area (mm²)
# 7. No. of branches
# 8. Shoot dry weight (g)
# 9. Shoot fresh weight (g)
# 10. Root dry weight (g)
# 11. Root fresh weight (g)
# 12. Root length (mm)

# ========================================================================
# 1. INSTALL REQUIRED PACKAGES
# ========================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "multcompView",
  "patchwork",
  "openxlsx",
  "grid"
)

installed_packages <- rownames(
  installed.packages()
)

for (pkg in required_packages) {

  if (!pkg %in% installed_packages) {

    install.packages(pkg)

  }

}

# ========================================================================
# 2. LOAD PACKAGES
# ========================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(multcompView)
library(patchwork)
library(openxlsx)
library(grid)


# ========================================================================
# 3. SELECT EXCEL FILE
# ========================================================================

file_path <- file.choose()


# ========================================================================
# 4. CREATE OUTPUT DIRECTORIES
# ========================================================================

output_dir <- "NP_Plant_Growth_Results"

dir.create(
  output_dir,
  showWarnings = FALSE
)

dir.create(
  file.path(
    output_dir,
    "Individual_Figures"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    output_dir,
    "Statistics"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)


# ========================================================================
# 5. READ EXCEL FILE
# ========================================================================
#
# Your worksheet is assumed to be named Sheet1.
#
# ========================================================================

raw_data <- read_excel(
  file_path,
  sheet = "Sheet1"
)


# ========================================================================
# 6. REQUIRED COLUMNS
# ========================================================================

required_columns <- c(

  "Replicate",

  "Treatment",

  "Plant_height_mm",

  "No_of_leaves",

  "Stem_girth_mm",

  "Leaf_length_mm",

  "Leaf_width_mm",

  "Leaf_area_mm2",

  "No_of_branches",

  "Shoot_dry_weight_g",

  "Shoot_fresh_weight_g",

  "Root_dry_weight_g",

  "Root_fresh_weight_g",

  "Root_length_mm"

)


missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)


if (
  length(missing_columns) > 0
) {

  stop(

    paste(

      "The following required columns are missing:",

      paste(
        missing_columns,
        collapse = ", "
      )

    )

  )

}


# ========================================================================
# 7. RECODE TREATMENTS
# ========================================================================
#
# Original Excel labels:
#
# Control 1
# Control 2
# Zinc Oxide NP
# Iron oxide NP
# Zinc/Iron Oxide NP

#
# Analysis/figure labels:
#
# C1
# C2
# A
# B
# C
#
# ========================================================================

data <- raw_data %>%

  mutate(

    Treatment_original =

      trimws(
        as.character(
          Treatment
        )
      ),

    Treatment = case_when(

      Treatment_original %in% c(
        "Control 1",
        "control 1",
        "Control1"
      )

      ~ "C1",


      Treatment_original %in% c(
        "Control 2",
        "control 2",
        "Control2"
      )

      ~ "C2",


      Treatment_original %in% c(
        "Zinc Oxide NP",
        "Zinc Oxide NP",
        "Zinc Oxide NP",
        
        "A"
      )

      ~ "A",


      Treatment_original %in% c(
        "Iron oxide NP",
        "Iron oxide NP",
        "Iron oxide NP",
        "B"
      )

      ~ "B",


      Treatment_original %in% c(
        "Zinc/Iron Oxide NP",
        "Zinc/Iron Oxide NP",
        "Zinc/Iron Oxide NP",
                "C"
      )

      ~ "C",


      TRUE ~ NA_character_

    ),


    Treatment = factor(

      Treatment,

      levels = c(

        "C1",

        "C2",

        "A",

        "B",

        "C"

      )

    )

  ) %>%

  mutate(

    across(

      all_of(

        c(

          "Plant_height_mm",

          "No_of_leaves",

          "Stem_girth_mm",

          "Leaf_length_mm",

          "Leaf_width_mm",

          "Leaf_area_mm2",

          "No_of_branches",

          "Shoot_dry_weight_g",

          "Shoot_fresh_weight_g",

          "Root_dry_weight_g",

          "Root_fresh_weight_g",

          "Root_length_mm"

        )

      ),

      as.numeric

    )

  )


# ========================================================================
# 8. VERIFY TREATMENT RECODING
# ========================================================================

if (
  any(
    is.na(
      data$Treatment
    )
  )
) {

  invalid_treatments <-

    unique(

      data$Treatment_original[
        is.na(data$Treatment)
      ]

    )


  stop(

    paste(

      "The following treatment names were not recognised:",

      paste(
        invalid_treatments,
        collapse = ", "
      ),

      "\n\nExpected Excel treatment names:",

      "\nControl 1",

      "\nControl 2",

      "\nZinc Oxide NP",

      "\nIron oxide NP",

      "\nZinc/Iron Oxide NP"

    )

  )

}


# ========================================================================
# 9. CHECK NUMBER OF REPLICATES
# ========================================================================

replication_check <-

  data %>%

  count(
    Treatment
  )

print(
  replication_check
)


# ========================================================================
# 10. PARAMETERS
# ========================================================================

parameters <- c(

  "Plant_height_mm",

  "No_of_leaves",

  "Stem_girth_mm",

  "Leaf_length_mm",

  "Leaf_width_mm",

  "Leaf_area_mm2",

  "No_of_branches",

  "Shoot_dry_weight_g",

  "Shoot_fresh_weight_g",

  "Root_dry_weight_g",

  "Root_fresh_weight_g",

  "Root_length_mm"

)


# ========================================================================
# 11. EXACT Y-AXIS LABELS
# ========================================================================

y_labels <- c(

  Plant_height_mm =
    "Plant height (mm)",

  No_of_leaves =
    "No. of leaves",

  Stem_girth_mm =
    "Stem girth (mm)",

  Leaf_length_mm =
    "Leaf length (mm)",

  Leaf_width_mm =
    "Leaf width (mm)",

  Leaf_area_mm2 =
    "Leaf area (mm²)",

  No_of_branches =
    "No. of branches",

  Shoot_dry_weight_g =
    "Shoot dry weight (g)",

  Shoot_fresh_weight_g =
    "Shoot fresh weight (g)",

  Root_dry_weight_g =
    "Root dry weight (g)",

  Root_fresh_weight_g =
    "Root fresh weight (g)",

  Root_length_mm =
    "Root length (mm)"

)


# ========================================================================
# 12. PANEL LABELS
# ========================================================================

panel_labels <- c(

  Plant_height_mm = "(a)",

  No_of_leaves = "(b)",

  Stem_girth_mm = "(c)",

  Leaf_length_mm = "(d)",

  Leaf_width_mm = "(e)",

  Leaf_area_mm2 = "(f)",

  No_of_branches = "(g)",

  Shoot_dry_weight_g = "(h)",

  Shoot_fresh_weight_g = "(i)",

  Root_dry_weight_g = "(j)",

  Root_fresh_weight_g = "(k)",

  Root_length_mm = "(l)"

)


# ========================================================================
# 13. FIVE-COLOUR CVD-FRIENDLY PALETTE
# ========================================================================
#
# Designed for five simultaneous treatment groups with strong
# separation in hue and lightness.
#
# C1 = black
# C2 = blue
# A         = orange
# B         = green
# C         = purple
#
# ========================================================================

treatment_colors <- c(

  "C1" = "#000000",

  "C2" = "#0072B2",

  "A" = "#E69F00",

  "B" = "#009E73",

  "C" = "#CC79A7"

)


# ========================================================================
# 14. TREATMENT ORDER
# ========================================================================

treatment_order <- c(

  "C1",

  "C2",

  "A",

  "B",

  "C"

)


# ========================================================================
# 15. STANDARD ERROR FUNCTION
# ========================================================================

calculate_se <- function(x) {

  x <- x[
    !is.na(x)
  ]


  if (
    length(x) <= 1
  ) {

    return(
      NA_real_
    )

  }


  sd(x) /
    sqrt(
      length(x)
    )

}


# ========================================================================
# 16. DESCRIPTIVE STATISTICS FUNCTION
# ========================================================================

calculate_descriptive <- function(

  df,

  parameter

) {

  df %>%

    group_by(
      Treatment
    ) %>%

    summarise(

      n =
        sum(
          !is.na(
            .data[[parameter]]
          )
        ),

      Mean =
        mean(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SD =
        sd(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SE =
        calculate_se(
          .data[[parameter]]
        ),

      Minimum =
        min(
          .data[[parameter]],
          na.rm = TRUE
        ),

      Maximum =
        max(
          .data[[parameter]],
          na.rm = TRUE
        ),

      .groups =
        "drop"

    ) %>%

    mutate(

      Parameter =
        parameter

    ) %>%

    select(

      Parameter,

      Treatment,

      n,

      Mean,

      SD,

      SE,

      Minimum,

      Maximum

    )

}


# ========================================================================
# 17. ANOVA FUNCTION
# ========================================================================

perform_anova <- function(

  df,

  parameter

) {


  analysis_data <-

    df %>%

    select(

      Treatment,

      Value =
        all_of(parameter)

    ) %>%

    filter(

      !is.na(Treatment),

      !is.na(Value)

    )


  # ----------------------------------------------------------------------
  # ONE-WAY ANOVA
  # ----------------------------------------------------------------------

  model <-

    aov(

      Value ~ Treatment,

      data = analysis_data

    )


  # ----------------------------------------------------------------------
  # ANOVA TABLE
  # ----------------------------------------------------------------------

  anova_table <-

    summary(
      model
    )[[1]]


  # ----------------------------------------------------------------------
  # ETA-SQUARED
  # ----------------------------------------------------------------------

  SS_treatment <-

    anova_table[
      "Treatment",
      "Sum Sq"
    ]


  SS_total <-

    sum(
      anova_table[
        ,
        "Sum Sq"
      ]
    )


  eta_squared <-

    SS_treatment /
    SS_total


  # ----------------------------------------------------------------------
  # PUBLICATION TABLE
  # ----------------------------------------------------------------------

  publication_table <-

    data.frame(

      Parameter =
        parameter,

      Source =
        rownames(
          anova_table
        ),

      df =
        anova_table[
          ,
          "Df"
        ],

      Sum_Squares =
        anova_table[
          ,
          "Sum Sq"
        ],

      Mean_Square =
        anova_table[
          ,
          "Mean Sq"
        ],

      F_value =
        anova_table[
          ,
          "F value"
        ],

      p_value =
        anova_table[
          ,
          "Pr(>F)"
        ],

      Eta_squared =
        NA_real_

    )


  publication_table$Eta_squared[

    publication_table$Source ==
      "Treatment"

  ] <-

    eta_squared


  list(

    model =
      model,

    table =
      publication_table

  )

}


# ========================================================================
# 18. TUKEY HSD + SIGNIFICANCE LETTERS
# ========================================================================
#
# α = 0.05
#
# ========================================================================

perform_tukey <- function(

  model,

  parameter

) {


  # ----------------------------------------------------------------------
  # TUKEY HSD
  # ----------------------------------------------------------------------

  tukey_result <-

    TukeyHSD(

      model,

      "Treatment",

      conf.level = 0.95

    )


  # ----------------------------------------------------------------------
  # TUKEY TABLE
  # ----------------------------------------------------------------------

  tukey_table <-

    as.data.frame(

      tukey_result$Treatment

    )


  tukey_table$Comparison <-

    rownames(
      tukey_table
    )


  rownames(
    tukey_table
  ) <- NULL


  tukey_table <-

    tukey_table %>%

    select(

      Comparison,

      everything()

    ) %>%

    mutate(

      Parameter =
        parameter,

      Alpha =
        0.05

    ) %>%

    select(

      Parameter,

      Alpha,

      everything()

    )


  # ----------------------------------------------------------------------
  # SIGNIFICANCE LETTERS
  # ----------------------------------------------------------------------

  letter_object <-

    multcompLetters4(

      model,

      tukey_result

    )


  letter_values <-

    letter_object$
    Treatment$
    Letters


  letters_table <-

    data.frame(

      Treatment =
        names(
          letter_values
        ),

      Letters =
        unname(
          letter_values
        ),

      stringsAsFactors =
        FALSE

    )


  letters_table$Treatment <-

    factor(

      letters_table$Treatment,

      levels =
        treatment_order

    )


  letters_table <-

    letters_table %>%

    arrange(
      Treatment
    )


  list(

    tukey =
      tukey_table,

    letters =
      letters_table

  )

}


# ========================================================================
# 19. STORAGE OBJECTS
# ========================================================================

all_descriptive <- list()

all_anova <- list()

all_tukey <- list()

all_letters <- list()

all_models <- list()

all_plots <- list()


# ========================================================================
# 20. MAIN ANALYSIS LOOP
# ========================================================================

for (

  parameter in parameters

) {


  message(

    "Analysing: ",

    parameter

  )


  # ======================================================================
  # DESCRIPTIVE STATISTICS
  # ======================================================================

  descriptive <-

    calculate_descriptive(

      data,

      parameter

    )


  all_descriptive[[parameter]] <-
    descriptive


  # ======================================================================
  # ANOVA
  # ======================================================================

  anova_result <-

    perform_anova(

      data,

      parameter

    )


  model <-
    anova_result$model


  all_models[[parameter]] <-
    model


  all_anova[[parameter]] <-
    anova_result$table


  # ======================================================================
  # TUKEY
  # ======================================================================

  tukey_result <-

    perform_tukey(

      model,

      parameter

    )


  all_tukey[[parameter]] <-
    tukey_result$tukey


  all_letters[[parameter]] <-
    tukey_result$letters


  # ======================================================================
  # PREPARE PLOT DATA
  # ======================================================================

  plot_data <-

    descriptive %>%

    left_join(

      tukey_result$letters,

      by = "Treatment"

    )


  # ======================================================================
  # CALCULATE Y-AXIS RANGE
  # ======================================================================

  max_value <-

    max(

      plot_data$Mean +
        plot_data$SE,

      na.rm = TRUE

    )


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <-

      max(

        plot_data$Mean,

        na.rm = TRUE

      )

  }


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <- 1

  }


  # ----------------------------------------------------------------------
  # 1.25 × MAXIMUM
  # ----------------------------------------------------------------------

  y_max <-

    1.25 *
    max_value


  # ======================================================================
  # SIGNIFICANCE LETTER POSITION
  # ======================================================================

  letter_offset <-

    0.035 *
    y_max


  plot_data <-

    plot_data %>%

    mutate(

      Letter_y =

        Mean +

        SE +

        letter_offset

    )


  # ======================================================================
  # KEEP LETTERS INSIDE PANEL
  # ======================================================================

  plot_data$Letter_y <-

    pmin(

      plot_data$Letter_y,

      0.95 *
        y_max

    )


  # ======================================================================
  # CREATE INDIVIDUAL BAR CHART
  # ======================================================================

  p <-

    ggplot(

      plot_data,

      aes(

        x = Treatment,

        y = Mean,

        fill = Treatment

      )

    ) +


    # --------------------------------------------------------------------
    # BARS
    # --------------------------------------------------------------------

    geom_col(

      width = 0.66,

      colour = "black",

      linewidth = 0.5

    ) +


    # --------------------------------------------------------------------
    # ERROR BARS
    # --------------------------------------------------------------------

    geom_errorbar(

      aes(

        ymin =
          pmax(
            0,
            Mean - SE
          ),

        ymax =
          Mean + SE

      ),

      width = 0.14,

      linewidth = 0.55,

      colour = "black"

    ) +


    # --------------------------------------------------------------------
    # SIGNIFICANCE LETTERS
    # --------------------------------------------------------------------

    geom_text(

      aes(

        x =
          Treatment,

        y =
          Letter_y,

        label =
          Letters

      ),

      family = "sans",

      size = 3.0,

      fontface = "bold",

      colour = "black",

      vjust = 0

    ) +


    # --------------------------------------------------------------------
    # PANEL LABEL
    # UPPER-LEFT INSIDE PANEL
    # --------------------------------------------------------------------

    annotate(

      "text",

      x = 0.55,

      y =
        0.985 *
        y_max,

      label =
        panel_labels[parameter],

      hjust = 0,

      vjust = 1,

      family = "sans",

      fontface = "bold",

      size = 3.2

    ) +


    # --------------------------------------------------------------------
    # FIVE-COLOUR CVD-FRIENDLY PALETTE
    # --------------------------------------------------------------------

    scale_fill_manual(

      values =
        treatment_colors,

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # Y-AXIS
    # --------------------------------------------------------------------

    scale_y_continuous(

      limits = c(

        0,

        y_max

      ),

      expand =
        expansion(
          mult =
            c(
              0,
              0
            )
        ),

      n.breaks = 5

    ) +


    # --------------------------------------------------------------------
    # SAME X-AXIS ACROSS ALL PANELS
    # --------------------------------------------------------------------

    scale_x_discrete(

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # LABELS
    # --------------------------------------------------------------------

    labs(

      x = "Treatment",

      y =
        y_labels[parameter]

    ) +


    # --------------------------------------------------------------------
    # CLEAN PUBLICATION THEME
    # --------------------------------------------------------------------

    theme_classic(

      base_size = 9,

      base_family = "sans"

    ) +


    theme(

      # ------------------------------------------------------------------
      # FONT
      # ------------------------------------------------------------------

      text =
        element_text(

          family = "sans",

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS
      # ------------------------------------------------------------------

      axis.text.x =
        element_text(

          size = 8,

          colour = "black",

          angle = 0,

          hjust = 0.5,

          margin =
            margin(
              t = 3
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS
      # ------------------------------------------------------------------

      axis.text.y =
        element_text(

          size = 8,

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.x =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              t = 5
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.y =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              r = 5
            )

        ),


      # ------------------------------------------------------------------
      # NO GRIDLINES
      # ------------------------------------------------------------------

      panel.grid.major =
        element_blank(),

      panel.grid.minor =
        element_blank(),


      # ------------------------------------------------------------------
      # NO LEGEND
      # ------------------------------------------------------------------

      legend.position =
        "none",


      # ------------------------------------------------------------------
      # PANEL BREATHING ROOM
      # ------------------------------------------------------------------

      plot.margin =
        margin(

          t = 5,

          r = 5,

          b = 5,

          l = 5,

          unit = "mm"

        ),


      # ------------------------------------------------------------------
      # AXIS LINE WIDTH ≥ 0.5 PT
      # ------------------------------------------------------------------

      axis.line =
        element_line(

          linewidth = 0.5,

          colour = "black"

        ),

      axis.ticks =
        element_line(

          linewidth = 0.5,

          colour = "black"

        )

    )


  # ======================================================================
  # STORE PLOT
  # ======================================================================

  all_plots[[parameter]] <-
    p


  # ======================================================================
  # SAVE INDIVIDUAL PDF
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(

          parameter,

          ".pdf"

        )

      ),

    plot =
      p,

    device =
      cairo_pdf,

    family =
      "sans",

    width =
      105,

    height =
      65,

    units =
      "mm"

  )


  # ======================================================================
  # SAVE INDIVIDUAL TIFF — 600 DPI
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(

          parameter,

          "_600dpi.tiff"

        )

      ),

    plot =
      p,

    device =
      "tiff",

    dpi =
      600,

    width =
      105,

    height =
      65,

    units =
      "mm",

    compression =
      "lzw"

  )

}


# ========================================================================
# 21. COMBINE ALL 12 PANELS
# ========================================================================

combined_figure <-

  wrap_plots(

    all_plots,

    ncol = 3,

    nrow = 4

  )


# ========================================================================
# 22. MULTI-PANEL SPACING
# ========================================================================

combined_figure <-

  combined_figure &

  theme(

    panel.spacing =
      unit(
        10,
        "mm"
      ),

    plot.margin =
      margin(

        t = 5,

        r = 5,

        b = 5,

        l = 5,

        unit = "mm"

      )

  )


# ========================================================================
# 23. SAVE COMBINED PDF
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_l.pdf"

    ),

  plot =
    combined_figure,

  device =
    cairo_pdf,

  family =
    "sans",

  width =
    180,

  height =
    230,

  units =
    "mm"

)


# ========================================================================
# 24. SAVE COMBINED TIFF — 600 DPI
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_l_600dpi.tiff"

    ),

  plot =
    combined_figure,

  device =
    "tiff",

  dpi =
    600,

  width =
    180,

  height =
    230,

  units =
    "mm",

  compression =
    "lzw"

)


# ========================================================================
# 25. COMBINE DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  bind_rows(
    all_descriptive
  )


# ========================================================================
# 26. COMBINE ANOVA RESULTS
# ========================================================================

anova_final <-

  bind_rows(
    all_anova
  )


# ========================================================================
# 27. COMBINE TUKEY RESULTS
# ========================================================================

tukey_final <-

  bind_rows(
    all_tukey
  )


# ========================================================================
# 28. COMBINE SIGNIFICANCE LETTERS
# ========================================================================

letters_final <-

  bind_rows(

    lapply(

      names(all_letters),

      function(parameter_name) {

        all_letters[[parameter_name]] %>%

          mutate(

            Parameter =
              parameter_name

          ) %>%

          select(

            Parameter,

            everything()

          )

      }

    )

  )


# ========================================================================
# 29. ROUND DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  descriptive_final %>%

  mutate(

    Mean =
      round(
        Mean,
        4
      ),

    SD =
      round(
        SD,
        4
      ),

    SE =
      round(
        SE,
        4
      ),

    Minimum =
      round(
        Minimum,
        4
      ),

    Maximum =
      round(
        Maximum,
        4
      )

  )


# ========================================================================
# 30. ROUND ANOVA RESULTS
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Sum_Squares =
      round(
        Sum_Squares,
        4
      ),

    Mean_Square =
      round(
        Mean_Square,
        4
      ),

    F_value =
      round(
        F_value,
        4
      ),

    p_value =
      signif(
        p_value,
        4
      ),

    Eta_squared =
      round(
        Eta_squared,
        4
      )

  )


# ========================================================================
# 31. SIGNIFICANCE AT α = 0.05 ONLY
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Significance = case_when(

      is.na(p_value)
      ~ "",

      p_value < 0.05
      ~ "Significant",

      TRUE
      ~ "Not significant"

    )

  )


# ========================================================================
# 32. PUBLICATION-READY ANOVA TABLE
# ========================================================================

anova_publication <-

  anova_final %>%

  select(

    Parameter,

    Source,

    df,

    Sum_Squares,

    Mean_Square,

    F_value,

    p_value,

    Significance,

    Eta_squared

  )


# ========================================================================
# 33. EXPORT STATISTICAL RESULTS TO EXCEL
# ========================================================================

write.xlsx(

  list(

    Descriptive_Statistics =
      descriptive_final,

    ANOVA_Publication =
      anova_publication,

    Tukey_HSD =
      tukey_final,

    Significance_Letters =
      letters_final,

    Treatment_Coding =
      data.frame(

        Original_Treatment = c(

          "Control 1",

          "Control 2",

          "Zinc Oxide NP",

          "Iron oxide NP",

          "Zinc/Iron Oxide NP"

        ),

        Analysis_Label = c(

          "C1",

          "C2",

          "A",

          "B",

          "C"

        )

      )

  ),

  file =
    file.path(

      output_dir,

      "Statistical_Results.xlsx"

    ),

  overwrite =
    TRUE

)


# ========================================================================
# 34. EXPORT CSV TABLES
# ========================================================================

write.csv(

  descriptive_final,

  file.path(

    output_dir,

    "Statistics",

    "Descriptive_Statistics.csv"

  ),

  row.names = FALSE

)


write.csv(

  anova_publication,

  file.path(

    output_dir,

    "Statistics",

    "ANOVA_Publication_Table.csv"

  ),

  row.names = FALSE

)


write.csv(

  tukey_final,

  file.path(

    output_dir,

    "Statistics",

    "Tukey_HSD.csv"

  ),

  row.names = FALSE

)


write.csv(

  letters_final,

  file.path(

    output_dir,

    "Statistics",

    "Significance_Letters.csv"

  ),

  row.names = FALSE

)


# ========================================================================
# 35. PRINT TREATMENT CODING
# ========================================================================

cat("\n")
cat("============================================================\n")
cat("TREATMENT CODING\n")
cat("============================================================\n\n")

cat("Control 1  = C1\n")
cat("Control 2  = C2\n")
cat("Zinc Oxide NP = A\n")
cat("Iron oxide NP    = B\n")
cat("Zinc/Iron Oxide NP  = C\n\n")


# ========================================================================
# 36. PRINT ANALYSIS INFORMATION
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS INFORMATION\n")
cat("============================================================\n\n")

cat("Parameters analysed: 12\n")
cat("Replications per treatment: 8\n")
cat("Significance level: α = 0.05\n")
cat("Statistical transformation: None\n")
cat("Error bars: Mean ± SE\n\n")


# ========================================================================
# 37. PRINT COMPLETION MESSAGE
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")


cat(

  "Combined PDF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_l.pdf"

  ),

  "\n\n"

)


cat(

  "Combined 600-dpi TIFF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_l_600dpi.tiff"

  ),

  "\n\n"

)


cat(

  "Statistical workbook:\n",

  file.path(

    output_dir,

    "Statistical_Results.xlsx"

  ),

  "\n\n"

)


# ========================================================================
# 38. PRINT ANOVA TABLE
# ========================================================================

cat("============================================================\n")
cat("PUBLICATION-READY ANOVA TABLE\n")
cat("============================================================\n\n")

print(
  anova_publication
)

# ========================================================================
# 39. PRINT DESCRIPTIVE STATISTICS
# ========================================================================

cat("\n")
cat("============================================================\n")
cat("DESCRIPTIVE STATISTICS\n")
cat("============================================================\n\n")

print(
  descriptive_final
)


# ========================================================================
# END OF SCRIPT
# ========================================================================



R script for proline and phenolic assay summary statistics, ANOVA, post hoc and bar charts



# ========================================================================
# 1. INSTALL REQUIRED PACKAGES
# ========================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "multcompView",
  "patchwork",
  "openxlsx",
  "grid"
)

installed_packages <- rownames(
  installed.packages()
)

for (pkg in required_packages) {

  if (!pkg %in% installed_packages) {

    install.packages(pkg)

  }

}


# ========================================================================
# 2. LOAD PACKAGES
# ========================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(multcompView)
library(patchwork)
library(openxlsx)
library(grid)


# ========================================================================
# 3. SELECT EXCEL FILE
# ========================================================================
#
# Select your Excel workbook when the file browser opens.
#
# ========================================================================

file_path <- file.choose()


# ========================================================================
# 4. CREATE OUTPUT FOLDERS
# ========================================================================

output_dir <- "NP_Phenolics_Proline_Results"

dir.create(
  output_dir,
  showWarnings = FALSE
)

dir.create(
  file.path(
    output_dir,
    "Individual_Figures"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    output_dir,
    "Statistics"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)


# ========================================================================
# 5. READ EXCEL FILE
# ========================================================================
#
# The worksheet should be named Sheet1.
#
# ========================================================================

raw_data <- read_excel(

  file_path,

  sheet = "Sheet1"

)


# ========================================================================
# 6. CHECK REQUIRED COLUMNS
# ========================================================================

required_columns <- c(

  "Replicate",

  "Treatment",

  "Phenolics",

  "Proline"

)


missing_columns <- setdiff(

  required_columns,

  names(raw_data)

)


if (
  length(missing_columns) > 0
) {

  stop(

    paste(

      "The following required columns are missing:",

      paste(
        missing_columns,
        collapse = ", "
      )

    )

  )

}


# ========================================================================
# 7. CLEAN AND RECODE TREATMENTS
# ========================================================================
#
# Original Excel names:
#
# Control 1
# Control 2
# Zinc Oxide NP
# Iron oxide NP
# Zinc/Iron Oxide NP
# Final analysis labels:
#
# C1
# C2
# A
# B
# C
#
# ========================================================================

data <- raw_data %>%

  mutate(

    Treatment_original =

      trimws(
        as.character(
          Treatment
        )
      ),

    Treatment = case_when(

      Treatment_original == "Control 1"
        ~ "C1",

      Treatment_original == "control 1"
        ~ "C1",

      Treatment_original == "Control1"
        ~ "C1",

      Treatment_original == "C1"
        ~ "C1",

      Treatment_original == "Control 2"
        ~ "C2",

      Treatment_original == "control 2"
        ~ "C2",

      Treatment_original == "Control2"
        ~ "C2",

      Treatment_original == "C2"
        ~ "C2",

      Treatment_original == "Zinc Oxide NP"
        ~ "A",

      Treatment_original == "Zinc Oxide NP"
        ~ "A",

      Treatment_original == "Zinc Oxide NP"
        ~ "A",

      Treatment_original == "Zinc Oxide NP"
        ~ "A",

      Treatment_original == "A"
        ~ "A",

      Treatment_original == "Iron oxide NP"
        ~ "B",

      Treatment_original == "Iron oxide NP"
        ~ "B",

      Treatment_original == "Iron oxide NP"
        ~ "B",

      Treatment_original == "Iron oxide NP"
        ~ "B",

      Treatment_original == "B"
        ~ "B",

      Treatment_original == "Zinc/Iron Oxide NP"
        ~ "C",

      Treatment_original == "Zinc/Iron Oxide NP"
        ~ "C",

      Treatment_original == "Zinc/Iron Oxide NP"
        ~ "C",

      Treatment_original == "Zinc/Iron Oxide NP"
        ~ "C",

      Treatment_original == "C"
        ~ "C",

      TRUE
        ~ NA_character_

    ),

    Treatment = factor(

      Treatment,

      levels = c(

        "C1",

        "C2",

        "A",

        "B",

        "C"

      )

    )

  ) %>%

  mutate(

    across(

      all_of(
        c(
          "Phenolics",
          "Proline"
        )
      ),

      as.numeric

    )

  )


# ========================================================================
# 8. VERIFY TREATMENT RECODING
# ========================================================================

if (
  any(
    is.na(
      data$Treatment
    )
  )
) {

  invalid_treatments <-

    unique(

      data$Treatment_original[
        is.na(data$Treatment)
      ]

    )


  stop(

    paste(

      "The following treatment names were not recognised:",

      paste(
        invalid_treatments,
        collapse = ", "
      ),

      "\n\nExpected Excel treatment names:",

      "\nControl 1",

      "\nControl 2",

      "\nZinc Oxide NP",

      "\nIron oxide NP",

      "\nZinc/Iron Oxide NP"

    )

  )

}


# ========================================================================
# 9. CHECK REPLICATIONS
# ========================================================================

replication_check <-

  data %>%

  count(
    Treatment
  )


print(
  replication_check
)


# ========================================================================
# 10. TREATMENT ORDER
# ========================================================================

treatment_order <- c(

  "C1",

  "C2",

  "A",

  "B",

  "C"

)


# ========================================================================
# 11. PARAMETERS
# ========================================================================

parameters <- c(

  "Phenolics",

  "Proline"

)


# ========================================================================
# 12. EXACT Y-AXIS LABELS
# ========================================================================

y_labels <- c(

  Phenolics =
    "Total Phenolic Content (mg GAE/g)",

  Proline =
    "Proline Content (μmoles g⁻¹ FW)"

)


# ========================================================================
# 13. PANEL LABELS
# ========================================================================

panel_labels <- c(

  Phenolics =
    "(a)",

  Proline =
    "(b)"

)


# ========================================================================
# 14. FIVE-COLOUR CVD-FRIENDLY PALETTE
# ========================================================================
#
# C1 = black
# C2 = blue
# A  = orange
# B  = green
# C  = purple
#
# ========================================================================

treatment_colors <- c(

  "C1" =
    "#000000",

  "C2" =
    "#0072B2",

  "A" =
    "#E69F00",

  "B" =
    "#009E73",

  "C" =
    "#CC79A7"

)


# ========================================================================
# 15. STANDARD ERROR FUNCTION
# ========================================================================

calculate_se <- function(x) {

  x <- x[
    !is.na(x)
  ]


  if (
    length(x) <= 1
  ) {

    return(
      NA_real_
    )

  }


  sd(x) /
    sqrt(
      length(x)
    )

}


# ========================================================================
# 16. DESCRIPTIVE STATISTICS FUNCTION
# ========================================================================

calculate_descriptive <- function(

  df,

  parameter

) {

  df %>%

    group_by(
      Treatment
    ) %>%

    summarise(

      n =
        sum(
          !is.na(
            .data[[parameter]]
          )
        ),

      Mean =
        mean(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SD =
        sd(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SE =
        calculate_se(
          .data[[parameter]]
        ),

      Minimum =
        min(
          .data[[parameter]],
          na.rm = TRUE
        ),

      Maximum =
        max(
          .data[[parameter]],
          na.rm = TRUE
        ),

      .groups =
        "drop"

    ) %>%

    mutate(

      Parameter =
        parameter

    ) %>%

    select(

      Parameter,

      Treatment,

      n,

      Mean,

      SD,

      SE,

      Minimum,

      Maximum

    )

}


# ========================================================================
# 17. ONE-WAY ANOVA FUNCTION
# ========================================================================

perform_anova <- function(

  df,

  parameter

) {

  analysis_data <-

    df %>%

    select(

      Treatment,

      Value =
        all_of(parameter)

    ) %>%

    filter(

      !is.na(Treatment),

      !is.na(Value)

    )


  # ----------------------------------------------------------------------
  # ONE-WAY ANOVA
  # ----------------------------------------------------------------------

  model <-

    aov(

      Value ~ Treatment,

      data = analysis_data

    )


  # ----------------------------------------------------------------------
  # ANOVA TABLE
  # ----------------------------------------------------------------------

  anova_table <-

    summary(
      model
    )[[1]]


  # ----------------------------------------------------------------------
  # ETA-SQUARED
  # ----------------------------------------------------------------------

  SS_treatment <-

    anova_table[
      "Treatment",
      "Sum Sq"
    ]


  SS_total <-

    sum(
      anova_table[
        ,
        "Sum Sq"
      ]
    )


  eta_squared <-

    SS_treatment /
    SS_total


  # ----------------------------------------------------------------------
  # PUBLICATION-READY TABLE
  # ----------------------------------------------------------------------

  publication_table <-

    data.frame(

      Parameter =
        parameter,

      Source =
        rownames(
          anova_table
        ),

      df =
        anova_table[
          ,
          "Df"
        ],

      Sum_Squares =
        anova_table[
          ,
          "Sum Sq"
        ],

      Mean_Square =
        anova_table[
          ,
          "Mean Sq"
        ],

      F_value =
        anova_table[
          ,
          "F value"
        ],

      p_value =
        anova_table[
          ,
          "Pr(>F)"
        ],

      Eta_squared =
        NA_real_

    )


  publication_table$Eta_squared[

    publication_table$Source ==
      "Treatment"

  ] <-

    eta_squared


  list(

    model =
      model,

    table =
      publication_table

  )

}


# ========================================================================
# 18. TUKEY HSD + SIGNIFICANCE LETTERS
# ========================================================================
#
# Significance level:
#
# α = 0.05
# Confidence level = 0.95
#
# ========================================================================

perform_tukey <- function(

  model,

  parameter

) {


  # ----------------------------------------------------------------------
  # TUKEY HSD
  # ----------------------------------------------------------------------

  tukey_result <-

    TukeyHSD(

      model,

      "Treatment",

      conf.level = 0.95

    )


  # ----------------------------------------------------------------------
  # TUKEY TABLE
  # ----------------------------------------------------------------------

  tukey_table <-

    as.data.frame(

      tukey_result$Treatment

    )


  tukey_table$Comparison <-

    rownames(
      tukey_table
    )


  rownames(
    tukey_table
  ) <- NULL


  tukey_table <-

    tukey_table %>%

    select(

      Comparison,

      everything()

    ) %>%

    mutate(

      Parameter =
        parameter,

      Alpha =
        0.05

    ) %>%

    select(

      Parameter,

      Alpha,

      everything()

    )


  # ----------------------------------------------------------------------
  # SIGNIFICANCE LETTERS
  # ----------------------------------------------------------------------

  letter_object <-

    multcompLetters4(

      model,

      tukey_result

    )


  letter_values <-

    letter_object$
    Treatment$
    Letters


  letters_table <-

    data.frame(

      Treatment =
        names(
          letter_values
        ),

      Letters =
        unname(
          letter_values
        ),

      stringsAsFactors =
        FALSE

    )


  letters_table$Treatment <-

    factor(

      letters_table$Treatment,

      levels =
        treatment_order

    )


  letters_table <-

    letters_table %>%

    arrange(
      Treatment
    )


  list(

    tukey =
      tukey_table,

    letters =
      letters_table

  )

}


# ========================================================================
# 19. STORAGE OBJECTS
# ========================================================================

all_descriptive <- list()

all_anova <- list()

all_tukey <- list()

all_letters <- list()

all_models <- list()

all_plots <- list()


# ========================================================================
# 20. MAIN ANALYSIS LOOP
# ========================================================================

for (

  parameter in parameters

) {


  message(

    "Analysing parameter: ",

    parameter

  )


  # ======================================================================
  # 20A. DESCRIPTIVE STATISTICS
  # ======================================================================

  descriptive <-

    calculate_descriptive(

      data,

      parameter

    )


  all_descriptive[[parameter]] <-
    descriptive


  # ======================================================================
  # 20B. ONE-WAY ANOVA
  # ======================================================================

  anova_result <-

    perform_anova(

      data,

      parameter

    )


  model <-
    anova_result$model


  all_models[[parameter]] <-
    model


  all_anova[[parameter]] <-
    anova_result$table


  # ======================================================================
  # 20C. TUKEY HSD
  # ======================================================================

  tukey_result <-

    perform_tukey(

      model,

      parameter

    )


  all_tukey[[parameter]] <-
    tukey_result$tukey


  all_letters[[parameter]] <-
    tukey_result$letters


  # ======================================================================
  # 20D. PREPARE PLOT DATA
  # ======================================================================

  plot_data <-

    descriptive %>%

    left_join(

      tukey_result$letters,

      by = "Treatment"

    )


  # ======================================================================
  # 20E. CALCULATE Y-AXIS RANGE
  # ======================================================================

  max_value <-

    max(

      plot_data$Mean +
        plot_data$SE,

      na.rm = TRUE

    )


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <-

      max(

        plot_data$Mean,

        na.rm = TRUE

      )

  }


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <- 1

  }


  # ======================================================================
  # 1.25 × MAXIMUM VALUE
  # ======================================================================

  y_max <-

    1.25 *
    max_value


  # ======================================================================
  # 20F. SIGNIFICANCE LETTER POSITION
  # ======================================================================

  letter_offset <-

    0.035 *
    y_max


  plot_data <-

    plot_data %>%

    mutate(

      Letter_y =

        Mean +

        SE +

        letter_offset

    )


  # ======================================================================
  # 20G. KEEP LETTERS INSIDE PANEL
  # ======================================================================

  plot_data$Letter_y <-

    pmin(

      plot_data$Letter_y,

      0.95 *
        y_max

    )


  # ======================================================================
  # 20H. CREATE INDIVIDUAL BAR CHART
  # ======================================================================

  p <-

    ggplot(

      plot_data,

      aes(

        x = Treatment,

        y = Mean,

        fill = Treatment

      )

    ) +


    # --------------------------------------------------------------------
    # BARS
    # --------------------------------------------------------------------

    geom_col(

      width = 0.66,

      colour = "black",

      linewidth = 0.5

    ) +


    # --------------------------------------------------------------------
    # ERROR BARS = MEAN ± SE
    # --------------------------------------------------------------------

    geom_errorbar(

      aes(

        ymin =
          pmax(
            0,
            Mean - SE
          ),

        ymax =
          Mean + SE

      ),

      width = 0.14,

      linewidth = 0.55,

      colour = "black"

    ) +


    # --------------------------------------------------------------------
    # SIGNIFICANCE LETTERS
    # --------------------------------------------------------------------

    geom_text(

      aes(

        x =
          Treatment,

        y =
          Letter_y,

        label =
          Letters

      ),

      family = "sans",

      size = 3.0,

      fontface = "bold",

      colour = "black",

      vjust = 0

    ) +


    # --------------------------------------------------------------------
    # PANEL LABEL
    # TOP-LEFT, INSIDE PANEL
    # --------------------------------------------------------------------

    annotate(

      "text",

      x = 0.55,

      y =
        0.985 *
        y_max,

      label =
        panel_labels[parameter],

      hjust = 0,

      vjust = 1,

      family = "sans",

      fontface = "bold",

      size = 3.2

    ) +


    # --------------------------------------------------------------------
    # FIVE-COLOUR CVD-FRIENDLY PALETTE
    # --------------------------------------------------------------------

    scale_fill_manual(

      values =
        treatment_colors,

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # Y-AXIS
    # --------------------------------------------------------------------

    scale_y_continuous(

      limits = c(

        0,

        y_max

      ),

      expand =
        expansion(

          mult = c(

            0,

            0

          )

        ),

      n.breaks = 5

    ) +


    # --------------------------------------------------------------------
    # SAME X-AXIS ACROSS ALL PANELS
    # --------------------------------------------------------------------

    scale_x_discrete(

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # AXIS LABELS
    # --------------------------------------------------------------------

    labs(

      x = "Treatment",

      y =
        y_labels[parameter]

    ) +


    # --------------------------------------------------------------------
    # CLEAN PUBLICATION THEME
    # --------------------------------------------------------------------

    theme_classic(

      base_size = 9,

      base_family = "sans"

    ) +


    theme(

      # ------------------------------------------------------------------
      # FONT
      # ------------------------------------------------------------------

      text =
        element_text(

          family = "sans",

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.x =
        element_text(

          size = 8,

          colour = "black",

          angle = 0,

          hjust = 0.5,

          margin =
            margin(
              t = 3
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.y =
        element_text(

          size = 8,

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.x =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              t = 5
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.y =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              r = 5
            )

        ),


      # ------------------------------------------------------------------
      # NO GRIDLINES
      # ------------------------------------------------------------------

      panel.grid.major =
        element_blank(),

      panel.grid.minor =
        element_blank(),


      # ------------------------------------------------------------------
      # NO LEGEND
      # ------------------------------------------------------------------

      legend.position =
        "none",


      # ------------------------------------------------------------------
      # PANEL BREATHING ROOM
      # ------------------------------------------------------------------

      plot.margin =
        margin(

          t = 5,

          r = 5,

          b = 5,

          l = 5,

          unit = "mm"

        ),


      # ------------------------------------------------------------------
      # LINE WIDTH
      # ------------------------------------------------------------------

      axis.line =
        element_line(

          linewidth = 0.5,

          colour = "black"

        ),

      axis.ticks =
        element_line(

          linewidth = 0.5,

          colour = "black"

        )

    )


  # ======================================================================
  # 20I. STORE PLOT
  # ======================================================================

  all_plots[[parameter]] <-
    p


  # ======================================================================
  # 20J. SAVE INDIVIDUAL PDF
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(

          parameter,

          ".pdf"

        )

      ),

    plot =
      p,

    device =
      cairo_pdf,

    family =
      "sans",

    width =
      105,

    height =
      65,

    units =
      "mm"

  )


  # ======================================================================
  # 20K. SAVE INDIVIDUAL TIFF — 600 DPI
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(

          parameter,

          "_600dpi.tiff"

        )

      ),

    plot =
      p,

    device =
      "tiff",

    dpi =
      600,

    width =
      105,

    height =
      65,

    units =
      "mm",

    compression =
      "lzw"

  )

}


# ========================================================================
# 21. COMBINE THE TWO PARAMETERS
# ========================================================================

combined_figure <-

  wrap_plots(

    all_plots,

    ncol = 2,

    nrow = 1

  )


# ========================================================================
# 22. MULTI-PANEL SPACING
# ========================================================================

combined_figure <-

  combined_figure &

  theme(

    panel.spacing =
      unit(
        10,
        "mm"
      ),

    plot.margin =
      margin(

        t = 5,

        r = 5,

        b = 5,

        l = 5,

        unit = "mm"

      )

  )


# ========================================================================
# 23. SAVE COMBINED PDF
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_b.pdf"

    ),

  plot =
    combined_figure,

  device =
    cairo_pdf,

  family =
    "sans",

  width =
    180,

  height =
    85,

  units =
    "mm"

)


# ========================================================================
# 24. SAVE COMBINED TIFF — 600 DPI
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_b_600dpi.tiff"

    ),

  plot =
    combined_figure,

  device =
    "tiff",

  dpi =
    600,

  width =
    180,

  height =
    85,

  units =
    "mm",

  compression =
    "lzw"

)


# ========================================================================
# 25. COMBINE DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  bind_rows(
    all_descriptive
  )


# ========================================================================
# 26. COMBINE ANOVA RESULTS
# ========================================================================

anova_final <-

  bind_rows(
    all_anova
  )


# ========================================================================
# 27. COMBINE TUKEY RESULTS
# ========================================================================

tukey_final <-

  bind_rows(
    all_tukey
  )


# ========================================================================
# 28. COMBINE SIGNIFICANCE LETTERS
# ========================================================================

letters_final <-

  bind_rows(

    lapply(

      names(all_letters),

      function(parameter_name) {

        all_letters[[parameter_name]] %>%

          mutate(

            Parameter =
              parameter_name

          ) %>%

          select(

            Parameter,

            everything()

          )

      }

    )

  )


# ========================================================================
# 29. ROUND DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  descriptive_final %>%

  mutate(

    Mean =
      round(
        Mean,
        4
      ),

    SD =
      round(
        SD,
        4
      ),

    SE =
      round(
        SE,
        4
      ),

    Minimum =
      round(
        Minimum,
        4
      ),

    Maximum =
      round(
        Maximum,
        4
      )

  )


# ========================================================================
# 30. ROUND ANOVA RESULTS
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Sum_Squares =
      round(
        Sum_Squares,
        4
      ),

    Mean_Square =
      round(
        Mean_Square,
        4
      ),

    F_value =
      round(
        F_value,
        4
      ),

    p_value =
      signif(
        p_value,
        4
      ),

    Eta_squared =
      round(
        Eta_squared,
        4
      )

  )


# ========================================================================
# 31. SIGNIFICANCE AT α = 0.05 ONLY
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Significance = case_when(

      is.na(p_value)
        ~ "",

      p_value < 0.05
        ~ "Significant",

      TRUE
        ~ "Not significant"

    )

  )


# ========================================================================
# 32. PUBLICATION-READY ANOVA TABLE
# ========================================================================

anova_publication <-

  anova_final %>%

  select(

    Parameter,

    Source,

    df,

    Sum_Squares,

    Mean_Square,

    F_value,

    p_value,

    Significance,

    Eta_squared

  )


# ========================================================================
# 33. EXPORT STATISTICAL RESULTS TO EXCEL
# ========================================================================

write.xlsx(

  list(

    Descriptive_Statistics =
      descriptive_final,

    ANOVA_Publication =
      anova_publication,

    Tukey_HSD =
      tukey_final,

    Significance_Letters =
      letters_final,

    Treatment_Coding =
      data.frame(

        Original_Treatment = c(

         
          "Control 1",

          "Control 2",

          "Zinc Oxide NP",

          "Iron oxide NP",

          "Zinc/Iron Oxide NP"

        ),

        Analysis_Label = c(

          "C1",

          "C2",

          "A",

          "B",

          "C"

        )

      )

  ),

  file =
    file.path(

      output_dir,

      "Statistical_Results.xlsx"

    ),

  overwrite =
    TRUE

)


# ========================================================================
# 34. EXPORT CSV TABLES
# ========================================================================

write.csv(

  descriptive_final,

  file.path(

    output_dir,

    "Statistics",

    "Descriptive_Statistics.csv"

  ),

  row.names = FALSE

)


write.csv(

  anova_publication,

  file.path(

    output_dir,

    "Statistics",

    "ANOVA_Publication_Table.csv"

  ),

  row.names = FALSE

)


write.csv(

  tukey_final,

  file.path(

    output_dir,

    "Statistics",

    "Tukey_HSD.csv"

  ),

  row.names = FALSE

)


write.csv(

  letters_final,

  file.path(

    output_dir,

    "Statistics",

    "Significance_Letters.csv"

  ),

  row.names = FALSE

)


# ========================================================================
# 35. PRINT TREATMENT CODING
# ========================================================================

cat("\n")
cat("============================================================\n")
cat("TREATMENT CODING\n")
cat("============================================================\n\n")

cat("Control 1 = C1\n")
cat("Control 2 = C2\n")
cat("Zinc Oxide NP = A\n")
cat("Iron oxide NP = B\n")
cat("Zinc/Iron Oxide NP = C\n\n")


# ========================================================================
# 36. PRINT ANALYSIS INFORMATION
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS INFORMATION\n")
cat("============================================================\n\n")

cat("Parameters analysed: 2\n")
cat("Replications per treatment: 3\n")
cat("Total observations: 15\n")
cat("Significance level: α = 0.05\n")
cat("Statistical transformation: None\n")
cat("Error bars: Mean ± SE\n\n")


# ========================================================================
# 37. PRINT COMPLETION MESSAGE
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")


cat(

  "Combined PDF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_b.pdf"

  ),

  "\n\n"

)


cat(

  "Combined 600-dpi TIFF:\n",

  file.path(

    output_dir,

    "Figure_1_a_to_b_600dpi.tiff"

  ),

  "\n\n"

)


cat(

  "Statistical workbook:\n",

  file.path(

    output_dir,

    "Statistical_Results.xlsx"

  ),

  "\n\n"

)


# ========================================================================
# 38. PRINT ANOVA TABLE
# ========================================================================

cat("============================================================\n")
cat("PUBLICATION-READY ANOVA TABLE\n")
cat("============================================================\n\n")

print(
  anova_publication
)


# ========================================================================
# 39. PRINT DESCRIPTIVE STATISTICS
# ========================================================================

cat("\n")
cat("============================================================\n")
cat("DESCRIPTIVE STATISTICS\n")
cat("============================================================\n\n")

print(
  descriptive_final
)


# ========================================================================
# END OF SCRIPT
# ========================================================================



R script for chlorophyll, carotenoid, MDA, CAT, POD, assay summary statistics, ANOVA, post hoc and bar charts



========================================================================
# 1. INSTALL REQUIRED PACKAGES
# ========================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "multcompView",
  "patchwork",
  "openxlsx",
  "grid"
)

installed_packages <- rownames(
  installed.packages()
)

for (pkg in required_packages) {

  if (!pkg %in% installed_packages) {

    install.packages(pkg)

  }

}


# ========================================================================
# 2. LOAD PACKAGES
# ========================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(multcompView)
library(patchwork)
library(openxlsx)
library(grid)


# ========================================================================
# 3. SELECT EXCEL FILE
# ========================================================================
#
# Select your Excel workbook when the file browser opens.
#
# ========================================================================

file_path <- file.choose()


# ========================================================================
# 4. CREATE OUTPUT FOLDERS
# ========================================================================

output_dir <- "NP_Biochemical_Results"

dir.create(
  output_dir,
  showWarnings = FALSE
)

dir.create(
  file.path(
    output_dir,
    "Individual_Figures"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    output_dir,
    "Statistics"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)


# ========================================================================
# 5. READ EXCEL DATA
# ========================================================================
#
# Your worksheet should be named Sheet1.
#
# ========================================================================

raw_data <- read_excel(
  file_path,
  sheet = "Sheet1"
)


# ========================================================================
# 6. CHECK REQUIRED COLUMNS
# ========================================================================

required_columns <- c(

  "Replicate",

  "Treatment",

  "Chlorophyll_A",

  "Chlorophyll_B",

  "Carotenoid",

  "Total_Chlorophyll",

  "Total_Chlorophyll_FW",

  "Total_Carotenoid_FW",

  "CAT",

  "POD",

  "MDA"

)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {

  stop(
    paste(
      "The following required columns are missing:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}


# ========================================================================
# # 7. RECODE TREATMENTS
# ========================================================================
#
# Original Excel labels:
#
# Control 1
# Control 2
# Zinc Oxide NP
# Iron oxide NP
# Zinc/Iron Oxide NP

#
# Analysis/figure labels:
#
# C1
# C2
# A
# B
# C
#
# 
# ========================================================================

data <- raw_data %>%

  mutate(

    Treatment_original =

      trimws(
        as.character(
          Treatment
        )
      ),

    Treatment = case_when(

      Treatment_original %in% c(
        "Control 1",
        "control 1",
        "Control1"
      )

      ~ "C1",


      Treatment_original %in% c(
        "Control 2",
        "control 2",
        "Control2"
      )

      ~ "C2",


      Treatment_original %in% c(
        "Zinc Oxide NP",
        "Zinc Oxide NP",
        "Zinc Oxide NP",
        
        "A"
      )

      ~ "A",


      Treatment_original %in% c(
        "Iron oxide NP",
        "Iron oxide NP",
        "Iron oxide NP",
        "B"
      )

      ~ "B",


      Treatment_original %in% c(
        "Zinc/Iron Oxide NP",
        "Zinc/Iron Oxide NP",
        "Zinc/Iron Oxide NP",
                "C"
      )

      ~ "C",


      TRUE ~ NA_character_

    ),


    Treatment = factor(

      Treatment,

      levels = c(

        "C1",

        "C2",

        "A",

        "B",

        "C"

      )

    )

  ) %>%

  mutate(

    across(

      all_of(

        c(


          "Chlorophyll_A",
          "Chlorophyll_B",
          "Carotenoid",
          "Total_Chlorophyll",
          "Total_Chlorophyll_FW",
          "Total_Carotenoid_FW",
          "CAT",
          "POD",
          "MDA"

        )

      ),

      as.numeric

    )

  )


# ========================================================================
# 8. VERIFY TREATMENTS
# ========================================================================

if (
  any(
    is.na(
      data$Treatment
    )
  )
) {

  invalid_treatments <- unique(
    data$Treatment_original[
      is.na(data$Treatment)
    ]
  )

  stop(
    paste(
      "The following treatment names were not recognised:",
      paste(
        invalid_treatments,
        collapse = ", "
      ),
      "\n\nExpected Excel treatment names:",
      "\nControl 1",
      "\nControl 2",
      "\nZinc Oxide NP",
      "\nIron oxide NP",
      "\nZinc/Iron Oxide NP"
    )
  )

}


# ========================================================================
# 9. CHECK REPLICATION NUMBERS
# ========================================================================

replication_check <- data %>%
  count(
    Treatment
  )

print(
  replication_check
)


# ========================================================================
# 10. PARAMETERS TO ANALYSE
# ========================================================================

parameters <- c(

  "Chlorophyll_A",

  "Chlorophyll_B",

  "Carotenoid",

  "Total_Chlorophyll",

  "Total_Chlorophyll_FW",

  "Total_Carotenoid_FW",

  "CAT",

  "POD",

  "MDA"

)


# ========================================================================
# 11. EXACT Y-AXIS LABELS
# ========================================================================

y_labels <- c(

  Chlorophyll_A =
    "Chlorophyll A (μg⋅mL⁻¹)",

  Chlorophyll_B =
    "Chlorophyll B (μg⋅mL⁻¹)",

  Carotenoid =
    "Carotenoid (μg⋅mL⁻¹)",

  Total_Chlorophyll =
    "Total Chlorophyll (μg⋅mL⁻¹)",

  Total_Chlorophyll_FW =
    "Total Chlorophyll Fresh weight (mg.g FW⁻¹)",

  Total_Carotenoid_FW =
    "Total Carotenoid Fresh weight (mg.g FW⁻¹)",

  CAT =
    "CAT (U/g FW)",

  POD =
    "POD (U/g FW)",

  MDA =
    "MDA (nmol·g⁻¹ FW)"

)


# ========================================================================
# 12. PANEL LABELS
# ========================================================================

panel_labels <- c(

  Chlorophyll_A =
    "(a)",

  Chlorophyll_B =
    "(b)",

  Carotenoid =
    "(c)",

  Total_Chlorophyll =
    "(d)",

  Total_Chlorophyll_FW =
    "(e)",

  Total_Carotenoid_FW =
    "(f)",

  CAT =
    "(g)",

  POD =
    "(h)",

  MDA =
    "(i)"

)


# ========================================================================
# 13. TREATMENT ORDER
# ========================================================================

treatment_order <- c(

  "C1",

  "C2",

  "A",

  "B",

  "C"

)


# ========================================================================
# 14. FIVE-COLOUR CVD-FRIENDLY PALETTE
# ========================================================================
#
# Five high-contrast colours suitable for simultaneous comparison.
#
# C1 = black
# C2 = blue
# A  = orange
# B  = green
# C  = purple
#
# ========================================================================

treatment_colors <- c(

  "C1" =
    "#000000",

  "C2" =
    "#0072B2",

  "A" =
    "#E69F00",

  "B" =
    "#009E73",

  "C" =
    "#CC79A7"

)


# ========================================================================
# 15. STANDARD ERROR FUNCTION
# ========================================================================

calculate_se <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (
    length(x) <= 1
  ) {

    return(
      NA_real_
    )

  }

  sd(x) /
    sqrt(
      length(x)
    )

}


# ========================================================================
# 16. DESCRIPTIVE STATISTICS FUNCTION
# ========================================================================

calculate_descriptive <- function(

  df,

  parameter

) {

  df %>%

    group_by(
      Treatment
    ) %>%

    summarise(

      n =
        sum(
          !is.na(
            .data[[parameter]]
          )
        ),

      Mean =
        mean(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SD =
        sd(
          .data[[parameter]],
          na.rm = TRUE
        ),

      SE =
        calculate_se(
          .data[[parameter]]
        ),

      Minimum =
        min(
          .data[[parameter]],
          na.rm = TRUE
        ),

      Maximum =
        max(
          .data[[parameter]],
          na.rm = TRUE
        ),

      .groups =
        "drop"

    ) %>%

    mutate(

      Parameter =
        parameter

    ) %>%

    select(

      Parameter,

      Treatment,

      n,

      Mean,

      SD,

      SE,

      Minimum,

      Maximum

    )

}


# ========================================================================
# 17. ONE-WAY ANOVA FUNCTION
# ========================================================================

perform_anova <- function(

  df,

  parameter

) {

  analysis_data <-

    df %>%

    select(

      Treatment,

      Value =
        all_of(parameter)

    ) %>%

    filter(

      !is.na(Treatment),

      !is.na(Value)

    )


  # ----------------------------------------------------------------------
  # ONE-WAY ANOVA
  # ----------------------------------------------------------------------

  model <-

    aov(

      Value ~ Treatment,

      data = analysis_data

    )


  # ----------------------------------------------------------------------
  # ANOVA TABLE
  # ----------------------------------------------------------------------

  anova_table <-

    summary(
      model
    )[[1]]


  # ----------------------------------------------------------------------
  # ETA-SQUARED
  # ----------------------------------------------------------------------

  SS_treatment <-

    anova_table[
      "Treatment",
      "Sum Sq"
    ]

  SS_total <-

    sum(
      anova_table[
        ,
        "Sum Sq"
      ]
    )

  eta_squared <-

    SS_treatment /
    SS_total


  # ----------------------------------------------------------------------
  # PUBLICATION-READY TABLE
  # ----------------------------------------------------------------------

  publication_table <-

    data.frame(

      Parameter =
        parameter,

      Source =
        rownames(
          anova_table
        ),

      df =
        anova_table[
          ,
          "Df"
        ],

      Sum_Squares =
        anova_table[
          ,
          "Sum Sq"
        ],

      Mean_Square =
        anova_table[
          ,
          "Mean Sq"
        ],

      F_value =
        anova_table[
          ,
          "F value"
        ],

      p_value =
        anova_table[
          ,
          "Pr(>F)"
        ],

      Eta_squared =
        NA_real_

    )


  publication_table$Eta_squared[

    publication_table$Source ==
      "Treatment"

  ] <-

    eta_squared


  list(

    model =
      model,

    table =
      publication_table

  )

}


# ========================================================================
# 18. TUKEY HSD + SIGNIFICANCE LETTERS
# ========================================================================
#
# α = 0.05
# Confidence level = 0.95
#
# ========================================================================

perform_tukey <- function(

  model,

  parameter

) {

  # ----------------------------------------------------------------------
  # TUKEY HSD
  # ----------------------------------------------------------------------

  tukey_result <-

    TukeyHSD(

      model,

      "Treatment",

      conf.level = 0.95

    )


  # ----------------------------------------------------------------------
  # TUKEY TABLE
  # ----------------------------------------------------------------------

  tukey_table <-

    as.data.frame(

      tukey_result$Treatment

    )


  tukey_table$Comparison <-

    rownames(
      tukey_table
    )


  rownames(
    tukey_table
  ) <- NULL


  tukey_table <-

    tukey_table %>%

    select(

      Comparison,

      everything()

    ) %>%

    mutate(

      Parameter =
        parameter,

      Alpha =
        0.05

    ) %>%

    select(

      Parameter,

      Alpha,

      everything()

    )


  # ----------------------------------------------------------------------
  # COMPACT LETTER DISPLAY
  # ----------------------------------------------------------------------

  letter_object <-

    multcompLetters4(

      model,

      tukey_result

    )


  letter_values <-

    letter_object$
    Treatment$
    Letters


  letters_table <-

    data.frame(

      Treatment =
        names(
          letter_values
        ),

      Letters =
        unname(
          letter_values
        ),

      stringsAsFactors =
        FALSE

    )


  letters_table$Treatment <-

    factor(

      letters_table$Treatment,

      levels =
        treatment_order

    )


  letters_table <-

    letters_table %>%

    arrange(
      Treatment
    )


  list(

    tukey =
      tukey_table,

    letters =
      letters_table

  )

}


# ========================================================================
# 19. STORAGE OBJECTS
# ========================================================================

all_descriptive <- list()

all_anova <- list()

all_tukey <- list()

all_letters <- list()

all_models <- list()

all_plots <- list()


# ========================================================================
# 20. MAIN ANALYSIS LOOP
# ========================================================================

for (
  parameter in parameters
) {

  message(
    "Analysing: ",
    parameter
  )


  # ======================================================================
  # 20A. DESCRIPTIVE STATISTICS
  # ======================================================================

  descriptive <-

    calculate_descriptive(

      data,

      parameter

    )


  all_descriptive[[parameter]] <-

    descriptive


  # ======================================================================
  # 20B. ONE-WAY ANOVA
  # ======================================================================

  anova_result <-

    perform_anova(

      data,

      parameter

    )


  model <-
    anova_result$model


  all_models[[parameter]] <-
    model


  all_anova[[parameter]] <-
    anova_result$table


  # ======================================================================
  # 20C. TUKEY HSD
  # ======================================================================

  tukey_result <-

    perform_tukey(

      model,

      parameter

    )


  all_tukey[[parameter]] <-

    tukey_result$tukey


  all_letters[[parameter]] <-

    tukey_result$letters


  # ======================================================================
  # 20D. PREPARE PLOT DATA
  # ======================================================================

  plot_data <-

    descriptive %>%

    left_join(

      tukey_result$letters,

      by = "Treatment"

    )


  # ======================================================================
  # 20E. CALCULATE INDEPENDENT Y-AXIS RANGE
  # ======================================================================

  max_value <-

    max(

      plot_data$Mean +
        plot_data$SE,

      na.rm = TRUE

    )


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <-

      max(

        plot_data$Mean,

        na.rm = TRUE

      )

  }


  if (

    !is.finite(max_value) ||

    max_value <= 0

  ) {

    max_value <- 1

  }


  # ----------------------------------------------------------------------
  # 1.25 × MAXIMUM VALUE
  # ----------------------------------------------------------------------

  y_max <-

    1.25 *
    max_value


  # ======================================================================
  # 20F. SIGNIFICANCE LETTER POSITION
  # ======================================================================

  letter_offset <-

    0.035 *
    y_max


  plot_data <-

    plot_data %>%

    mutate(

      Letter_y =

        Mean +

        SE +

        letter_offset

    )


  # ======================================================================
  # 20G. KEEP LETTERS INSIDE PANEL
  # ======================================================================

  plot_data$Letter_y <-

    pmin(

      plot_data$Letter_y,

      0.95 *
        y_max

    )


  # ======================================================================
  # 20H. CREATE BAR CHART
  # ======================================================================

  p <-

    ggplot(

      plot_data,

      aes(

        x = Treatment,

        y = Mean,

        fill = Treatment

      )

    ) +


    # --------------------------------------------------------------------
    # BARS
    # --------------------------------------------------------------------

    geom_col(

      width = 0.66,

      colour = "black",

      linewidth = 0.5

    ) +


    # --------------------------------------------------------------------
    # ERROR BARS = MEAN ± SE
    # --------------------------------------------------------------------

    geom_errorbar(

      aes(

        ymin =
          pmax(
            0,
            Mean - SE
          ),

        ymax =
          Mean + SE

      ),

      width = 0.14,

      linewidth = 0.55,

      colour = "black"

    ) +


    # --------------------------------------------------------------------
    # TUKEY SIGNIFICANCE LETTERS
    # --------------------------------------------------------------------

    geom_text(

      aes(

        x =
          Treatment,

        y =
          Letter_y,

        label =
          Letters

      ),

      family = "sans",

      size = 3.0,

      fontface = "bold",

      colour = "black",

      vjust = 0

    ) +


    # --------------------------------------------------------------------
    # PANEL LABEL
    # UPPER-LEFT, INSIDE THE PANEL
    # --------------------------------------------------------------------

    annotate(

      "text",

      x = 0.55,

      y =
        0.985 *
        y_max,

      label =
        panel_labels[parameter],

      hjust = 0,

      vjust = 1,

      family = "sans",

      fontface = "bold",

      size = 3.2

    ) +


    # --------------------------------------------------------------------
    # FIVE-COLOUR CVD-FRIENDLY PALETTE
    # --------------------------------------------------------------------

    scale_fill_manual(

      values =
        treatment_colors,

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # Y-AXIS
    # --------------------------------------------------------------------

    scale_y_continuous(

      limits = c(

        0,

        y_max

      ),

      expand =
        expansion(

          mult = c(

            0,

            0

          )

        ),

      n.breaks = 5

    ) +


    # --------------------------------------------------------------------
    # SAME X-AXIS ACROSS ALL PANELS
    # --------------------------------------------------------------------

    scale_x_discrete(

      limits =
        treatment_order,

      drop = FALSE

    ) +


    # --------------------------------------------------------------------
    # AXIS LABELS
    # --------------------------------------------------------------------

    labs(

      x = "Treatment",

      y =
        y_labels[parameter]

    ) +


    # --------------------------------------------------------------------
    # CLEAN PUBLICATION THEME
    # --------------------------------------------------------------------

    theme_classic(

      base_size = 9,

      base_family = "sans"

    ) +


    theme(

      # ------------------------------------------------------------------
      # FONT
      # ------------------------------------------------------------------

      text =
        element_text(

          family = "sans",

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.x =
        element_text(

          size = 8,

          colour = "black",

          angle = 0,

          hjust = 0.5,

          margin =
            margin(
              t = 3
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TEXT
      # ------------------------------------------------------------------

      axis.text.y =
        element_text(

          size = 8,

          colour = "black"

        ),


      # ------------------------------------------------------------------
      # X-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.x =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              t = 5
            )

        ),


      # ------------------------------------------------------------------
      # Y-AXIS TITLE
      # ------------------------------------------------------------------

      axis.title.y =
        element_text(

          size = 9,

          colour = "black",

          margin =
            margin(
              r = 5
            )

        ),


      # ------------------------------------------------------------------
      # NO GRIDLINES
      # ------------------------------------------------------------------

      panel.grid.major =
        element_blank(),

      panel.grid.minor =
        element_blank(),


      # ------------------------------------------------------------------
      # REMOVE LEGEND
      # ------------------------------------------------------------------

      legend.position =
        "none",


      # ------------------------------------------------------------------
      # PANEL BREATHING ROOM
      # ------------------------------------------------------------------

      plot.margin =
        margin(

          t = 5,

          r = 5,

          b = 5,

          l = 5,

          unit = "mm"

        ),


      # ------------------------------------------------------------------
      # AXIS LINE WIDTH >= 0.5 PT
      # ------------------------------------------------------------------

      axis.line =
        element_line(

          linewidth = 0.5,

          colour = "black"

        ),

      axis.ticks =
        element_line(

          linewidth = 0.5,

          colour = "black"

        )

    )


  # ======================================================================
  # 20I. STORE PLOT
  # ======================================================================

  all_plots[[parameter]] <-
    p


  # ======================================================================
  # 20J. SAVE INDIVIDUAL PDF
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(
          parameter,
          ".pdf"
        )

      ),

    plot =
      p,

    device =
      cairo_pdf,

    family =
      "sans",

    width =
      105,

    height =
      65,

    units =
      "mm"

  )


  # ======================================================================
  # 20K. SAVE INDIVIDUAL TIFF — 600 DPI
  # ======================================================================

  ggsave(

    filename =
      file.path(

        output_dir,

        "Individual_Figures",

        paste0(

          parameter,

          "_600dpi.tiff"

        )

      ),

    plot =
      p,

    device =
      "tiff",

    dpi =
      600,

    width =
      105,

    height =
      65,

    units =
      "mm",

    compression =
      "lzw"

  )

}


# ========================================================================
# 21. COMBINE ALL NINE PANELS
# ========================================================================

combined_figure <-

  wrap_plots(

    all_plots,

    ncol = 3,

    nrow = 3

  )


# ========================================================================
# 22. MULTI-PANEL SPACING
# ========================================================================

combined_figure <-

  combined_figure &

  theme(

    panel.spacing =
      unit(
        10,
        "mm"
      ),

    plot.margin =
      margin(

        t = 5,

        r = 5,

        b = 5,

        l = 5,

        unit = "mm"

      )

  )


# ========================================================================
# 23. SAVE COMBINED PDF
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_i.pdf"

    ),

  plot =
    combined_figure,

  device =
    cairo_pdf,

  family =
    "sans",

  width =
    180,

  height =
    190,

  units =
    "mm"

)


# ========================================================================
# 24. SAVE COMBINED TIFF — 600 DPI
# ========================================================================

ggsave(

  filename =
    file.path(

      output_dir,

      "Figure_1_a_to_i_600dpi.tiff"

    ),

  plot =
    combined_figure,

  device =
    "tiff",

  dpi =
    600,

  width =
    180,

  height =
    190,

  units =
    "mm",

  compression =
    "lzw"

)


# ========================================================================
# 25. COMBINE DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  bind_rows(
    all_descriptive
  )


# ========================================================================
# 26. COMBINE ANOVA RESULTS
# ========================================================================

anova_final <-

  bind_rows(
    all_anova
  )


# ========================================================================
# 27. COMBINE TUKEY RESULTS
# ========================================================================

tukey_final <-

  bind_rows(
    all_tukey
  )


# ========================================================================
# 28. COMBINE SIGNIFICANCE LETTERS
# ========================================================================

letters_final <-

  bind_rows(

    lapply(

      names(all_letters),

      function(parameter_name) {

        all_letters[[parameter_name]] %>%

          mutate(

            Parameter =
              parameter_name

          ) %>%

          select(

            Parameter,

            everything()

          )

      }

    )

  )


# ========================================================================
# 29. ROUND DESCRIPTIVE STATISTICS
# ========================================================================

descriptive_final <-

  descriptive_final %>%

  mutate(

    Mean =
      round(
        Mean,
        4
      ),

    SD =
      round(
        SD,
        4
      ),

    SE =
      round(
        SE,
        4
      ),

    Minimum =
      round(
        Minimum,
        4
      ),

    Maximum =
      round(
        Maximum,
        4
      )

  )


# ========================================================================
# 30. ROUND ANOVA RESULTS
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Sum_Squares =
      round(
        Sum_Squares,
        4
      ),

    Mean_Square =
      round(
        Mean_Square,
        4
      ),

    F_value =
      round(
        F_value,
        4
      ),

    p_value =
      signif(
        p_value,
        4
      ),

    Eta_squared =
      round(
        Eta_squared,
        4
      )

  )


# ========================================================================
# 31. SIGNIFICANCE AT α = 0.05 ONLY
# ========================================================================

anova_final <-

  anova_final %>%

  mutate(

    Significance = case_when(

      is.na(p_value)
      ~ "",

      p_value < 0.05
      ~ "Significant",

      TRUE
      ~ "Not significant"

    )

  )


# ========================================================================
# 32. PUBLICATION-READY ANOVA TABLE
# ========================================================================

anova_publication <-

  anova_final %>%

  select(

    Parameter,

    Source,

    df,

    Sum_Squares,

    Mean_Square,

    F_value,

    p_value,

    Significance,

    Eta_squared

  )


# ========================================================================
# 33. EXPORT STATISTICAL RESULTS TO EXCEL
# ========================================================================

write.xlsx(

  list(

    Descriptive_Statistics =
      descriptive_final,

    ANOVA_Publication =
      anova_publication,

    Tukey_HSD =
      tukey_final,

    Significance_Letters =
      letters_final,

    Treatment_Coding =
      data.frame(

        Original_Treatment = c(

          "Control 1",

          "Control 2",

          "Zinc Oxide NP",

          "Iron oxide NP",

          "Zinc/Iron Oxide NP"

        ),

        Analysis_Label = c(

          "C1",

          "C2",

          "A",

          "B",

          "C"

        )

      )

  ),

  file =
    file.path(

      output_dir,

      "Statistical_Results.xlsx"

    ),

  overwrite =
    TRUE

)


# ========================================================================
# 34. EXPORT CSV TABLES
# ========================================================================

write.csv(

  descriptive_final,

  file.path(

    output_dir,

    "Statistics",

    "Descriptive_Statistics.csv"

  ),

  row.names = FALSE

)


write.csv(

  anova_publication,

  file.path(

    output_dir,

    "Statistics",

    "ANOVA_Publication_Table.csv"

  ),

  row.names = FALSE

)


write.csv(

  tukey_final,

  file.path(

    output_dir,

    "Statistics",

    "Tukey_HSD.csv"

  ),

  row.names = FALSE

)


write.csv(

  letters_final,

  file.path(

    output_dir,

    "Statistics",

    "Significance_Letters.csv"

  ),

  row.names = FALSE

)


# ========================================================================
# 35. PRINT TREATMENT CODING
# ========================================================================

cat("\n")
cat("============================================================\n")
cat("TREATMENT CODING\n")
cat("============================================================\n\n")

cat("Control 1 = C1\n")
cat("Control 2 = C2\n")
cat("Zinc Oxide NP = A\n")
cat("Iron oxide NP = B\n")
cat("Zinc/Iron Oxide NP = C\n\n")


# ========================================================================
# 36. PRINT ANALYSIS INFORMATION
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS INFORMATION\n")
cat("============================================================\n\n")

cat("Parameters analysed: 9\n")
cat("Replications per treatment: 3\n")
cat("Total observations: 15\n")
cat("Significance level: α = 0.05\n")
cat("Statistical transformation: None\n")
cat("Error bars: Mean ± SE\n\n")


# ========================================================================
# 37. PRINT COMPLETION MESSAGE
# ========================================================================

cat("============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(

  "Combined PDF:\n",

  file.path(
    output_dir,
    "Figure_1_a_to_i.pdf"
  ),

  "\n\n"

)


cat(

  "Combined 600-dpi TIFF:\n",

  file.path(
    output_dir,
    "Figure_1_a_to_i_600dpi.tiff"
  ),

  "\n\n"

)


cat(

  "Statistical workbook:\n",

  file.path(
    output_dir,
    "Statistical_Results.xlsx"
  ),

  "\n\n"

)


# ========================================================================
# 38. PRINT ANOVA TABLE
# ========================================================================

cat(
  "============================================================\n"
)

cat(
  "PUBLICATION-READY ANOVA TABLE\n"
)

cat(
  "============================================================\n\n"
)

print(
  anova_publication
)


# ========================================================================
# 39. PRINT DESCRIPTIVE STATISTICS
# ========================================================================

cat("\n")
cat(
  "============================================================\n"
)

cat(
  "DESCRIPTIVE STATISTICS\n"
)

cat(
  "============================================================\n\n"
)

print(
  descriptive_final
)


# ========================================================================
# END OF SCRIPT
# ========================================================================

