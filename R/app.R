library(shiny)
library(ggplot2)

# ------------------------------------------------------------
# Bootstrap coefficient table for an lm model (pairs bootstrap)
# ------------------------------------------------------------
boot_lm_coef <- function(formula, data, B = 200, conf = 0.95, seed = 12345) {
  set.seed(seed)
  n <- nrow(data)

  fit0 <- lm(formula, data = data)
  coef0 <- coef(fit0)
  terms0 <- names(coef0)
  p <- length(terms0)

  boot_mat <- matrix(NA_real_, nrow = B, ncol = p)
  colnames(boot_mat) <- terms0

  for (b in seq_len(B)) {
    id <- sample.int(n, size = n, replace = TRUE)
    db <- data[id, , drop = FALSE]
    fb <- lm(formula, data = db)
    cb <- coef(fb)
    boot_mat[b, ] <- cb[terms0]
  }

  se_boot <- apply(boot_mat, 2, sd, na.rm = TRUE)
  alpha <- 1 - conf
  ci_lo <- apply(boot_mat, 2, quantile, probs = alpha/2, na.rm = TRUE)
  ci_hi <- apply(boot_mat, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)

  data.frame(
    term = terms0,
    estimate = unname(coef0),
    se_ols = summary(fit0)$coefficients[, "Std. Error"],
    se_boot = unname(se_boot),
    ci_lo = unname(ci_lo),
    ci_hi = unname(ci_hi),
    row.names = NULL
  )
}

# ------------------------------------------------------------
# Added-variable plot with bootstrap SE + pointwise band
# ------------------------------------------------------------
avplot_boot <- function(mod, var, data, B = 200, grid_n = 120,
                        conf = 0.95, seed = 2025, show_band = TRUE,
                        pch = 16, cex = 0.7) {

  yname <- all.vars(formula(mod))[1]
  terms_mod <- attr(terms(mod), "term.labels")
  if (!(var %in% terms_mod)) stop(sprintf("'%s' not in model.", var))
  others <- setdiff(terms_mod, var)

  f_y <- as.formula(paste(yname, "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))
  f_x <- as.formula(paste(var,   "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))

  y_res <- resid(lm(f_y, data = data))
  x_res <- resid(lm(f_x, data = data))

  av_fit <- lm(y_res ~ x_res)
  slope_hat <- coef(av_fit)[2]

  set.seed(seed)
  n <- nrow(data)
  xgrid <- seq(min(x_res), max(x_res), length.out = grid_n)

  slope_b <- numeric(B)
  pred_mat <- matrix(NA_real_, nrow = grid_n, ncol = B)

  for (b in seq_len(B)) {
    id <- sample.int(n, size = n, replace = TRUE)
    db <- data[id, , drop = FALSE]

    yb_res <- resid(lm(f_y, data = db))
    xb_res <- resid(lm(f_x, data = db))

    fit_b <- lm(yb_res ~ xb_res)
    cb <- coef(fit_b)
    slope_b[b] <- cb[2]
    pred_mat[, b] <- cb[1] + cb[2] * xgrid
  }

  slope_se_boot <- sd(slope_b, na.rm = TRUE)
  alpha <- 1 - conf
  slope_ci <- quantile(slope_b, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  if (show_band) {
    band_lo <- apply(pred_mat, 1, quantile, probs = alpha/2, na.rm = TRUE)
    band_hi <- apply(pred_mat, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  }

  plot(
    x_res, y_res,
    xlab = paste0(var, " | others"),
    ylab = paste0(yname, " | others"),
    pch = pch, cex = cex
  )

  if (show_band) {
    polygon(
      c(xgrid, rev(xgrid)), c(band_lo, rev(band_hi)),
      border = NA, col = adjustcolor("grey70", 0.5)
    )
  }

  abline(av_fit, lwd = 2)

  mtext(
    sprintf("Slope=%.3f; Boot SE=%.3f", slope_hat, slope_se_boot),
    side = 3, line = 0.2, cex = 0.85
  )

  invisible(list(slope_hat = slope_hat, slope_se_boot = slope_se_boot, slope_ci = slope_ci))
}

# ------------------------------------------------------------
# Helper: coefficient comparison aligned across two models
# ------------------------------------------------------------
coef_compare <- function(modA, modB) {
  ca <- coef(modA)
  cb <- coef(modB)
  nm <- sort(unique(c(names(ca), names(cb))))
  data.frame(
    term = nm,
    base = unname(ca[nm]),
    compare = unname(cb[nm]),
    delta = unname(cb[nm] - ca[nm]),
    row.names = NULL
  )
}

# ------------------------------------------------------------
# Helper: check nesting (A nested in B) based on terms
# ------------------------------------------------------------
is_nested_lm <- function(modA, modB) {
  # same response?
  yA <- all.vars(formula(modA))[1]
  yB <- all.vars(formula(modB))[1]
  if (!identical(yA, yB)) return(FALSE)

  # compare term labels (ignoring order)
  tA <- attr(terms(modA), "term.labels")
  tB <- attr(terms(modB), "term.labels")
  all(tA %in% tB)
}

# ------------------------------------------------------------
# Helper: K-fold CV metrics for lm(formula, data)
# ------------------------------------------------------------
cv_lm_metrics <- function(formula, data, K = 5, seed = 1) {
  set.seed(seed)
  n <- nrow(data)
  # random fold assignments
  folds <- sample(rep(seq_len(K), length.out = n))

  yname <- all.vars(formula)[1]
  y <- data[[yname]]

  pred <- rep(NA_real_, n)

  for (k in seq_len(K)) {
    idx_te <- which(folds == k)
    idx_tr <- which(folds != k)
    fit <- lm(formula, data = data[idx_tr, , drop = FALSE])
    pred[idx_te] <- predict(fit, newdata = data[idx_te, , drop = FALSE])
  }

  resid <- y - pred
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  mae  <- mean(abs(resid), na.rm = TRUE)

  # out-of-sample R^2 (defined vs mean(y))
  ss_res <- sum(resid^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - ss_res / ss_tot

  data.frame(
    metric = c("RMSE", "MAE", "CV_R2"),
    value  = c(rmse, mae, r2),
    row.names = NULL
  )
}

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  titlePanel("Compare models (uploaded CSV): coefficient shifts + partial plots + bootstrap SE"),

  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Upload CSV", accept = c(".csv")),
      tags$small("Upload any CSV. Choose outcome and predictors below."),

      tags$hr(),

      selectInput("ycol", "Outcome (y)", choices = character(0)),

      uiOutput("pred_ui"),

      tags$hr(),
      checkboxInput("use_boot", "Use bootstrap (for model SE + plot bands)", TRUE),
      conditionalPanel(
        condition = "input.use_boot == true",
        sliderInput("B", "Bootstrap resamples (B)", min = 200, max = 1000, value = 200, step = 200)
      ),

      tags$hr(),
      checkboxInput("use_cv", "Compute K-fold CV (prediction comparison)", FALSE),
      conditionalPanel(
        condition = "input.use_cv == true",
        sliderInput("K", "Number of folds (K)", min = 3, max = 10, value = 5, step = 1),
        numericInput("cv_seed", "CV seed", value = 1, min = 1, step = 1)
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel(
          "Summaries + bootstrap SE",
          tags$h4("Model A coefficient table (OLS SE + bootstrap SE)"),
          tableOutput("bootA_tbl"),
          tags$h4("Model B coefficient table (OLS SE + bootstrap SE)"),
          tableOutput("bootB_tbl")
        ),

        tabPanel(
          "Coefficient shifts",
          tags$h4("Raw coefficients aligned (B - A shows change)"),
          tableOutput("coef_tbl")
        ),

        tabPanel(
          "Model comparison tests",
          tags$h4("Information criteria (always available)"),
          tableOutput("ic_tbl"),

          tags$hr(),
          tags$h4("Nested-model test (only if Model A is nested in Model B)"),
          verbatimTextOutput("nested_txt"),
          tableOutput("anova_tbl"),

          tags$hr(),
          tags$h4("Cross-validation (optional)"),
          tableOutput("cv_tbl")
        ),

        tabPanel(
          "Partial plots (side-by-side)",
          plotOutput("pp_side", height = "1100px")
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {

  # 1) Read uploaded data
  data_reactive <- reactive({
    req(input$datafile)
    dat <- read.csv(input$datafile$datapath, stringsAsFactors = FALSE, check.names = TRUE)

    # drop incomplete rows (keeps residualization/bootstraps stable)
    dat <- dat[complete.cases(dat), , drop = FALSE]

    validate(
      need(nrow(dat) >= 30, "Need at least 30 complete rows."),
      need(ncol(dat) >= 2, "Need at least 2 columns (1 outcome + >=1 predictor).")
    )
    dat
  })

  # 2) Populate/update outcome choices AFTER upload (and preserve selection if possible)
  observeEvent(data_reactive(), {
    dat <- data_reactive()
    cols <- names(dat)

    default_y <- if ("y" %in% cols) "y" else cols[1]
    current_y <- isolate(input$ycol)

    updateSelectInput(
      session, "ycol",
      choices = cols,
      selected = if (!is.null(current_y) && current_y %in% cols) current_y else default_y
    )
  }, ignoreInit = TRUE)

  # 3) Predictor UI depends on chosen y
  output$pred_ui <- renderUI({
    dat <- data_reactive()
    req(input$ycol)

    cols <- names(dat)
    preds <- setdiff(cols, input$ycol)

    tagList(
      checkboxGroupInput(
        "varsA", "Model A predictors",
        choices = preds,
        selected = intersect(c("x1", "x2", "x3"), preds)
      ),
      checkboxGroupInput(
        "varsB", "Model B predictors",
        choices = preds,
        selected = intersect(c("x1", "x2", "x3", "x4"), preds)
      )
    )
  })

  # 4) Keep predictor choices synced when y changes (and preserve selections)
  observeEvent(input$ycol, {
    dat <- data_reactive()
    preds <- setdiff(names(dat), input$ycol)

    selA <- intersect(isolate(input$varsA), preds)
    selB <- intersect(isolate(input$varsB), preds)

    updateCheckboxGroupInput(session, "varsA", choices = preds, selected = selA)
    updateCheckboxGroupInput(session, "varsB", choices = preds, selected = selB)
  }, ignoreInit = TRUE)

  # 5) Validate numeric columns for lm()
  validate_numeric <- reactive({
    dat <- data_reactive()
    req(input$ycol)

    y_ok  <- is.numeric(dat[[input$ycol]]) || is.integer(dat[[input$ycol]])
    xA_ok <- length(input$varsA) == 0 || all(sapply(input$varsA, function(v) is.numeric(dat[[v]]) || is.integer(dat[[v]])))
    xB_ok <- length(input$varsB) == 0 || all(sapply(input$varsB, function(v) is.numeric(dat[[v]]) || is.integer(dat[[v]])))

    validate(
      need(y_ok, "Outcome (y) must be numeric for lm()."),
      need(xA_ok, "All Model A predictors must be numeric."),
      need(xB_ok, "All Model B predictors must be numeric.")
    )
    TRUE
  })

  # 6) Formulas
  formA <- reactive({
    validate_numeric()
    y <- input$ycol
    vars <- input$varsA
    if (length(vars) == 0) as.formula(paste(y, "~ 1")) else as.formula(paste(y, "~", paste(vars, collapse = " + ")))
  })

  formB <- reactive({
    validate_numeric()
    y <- input$ycol
    vars <- input$varsB
    if (length(vars) == 0) as.formula(paste(y, "~ 1")) else as.formula(paste(y, "~", paste(vars, collapse = " + ")))
  })

  # 7) Models
  modelA <- reactive({
    lm(formA(), data = data_reactive())
  })
  modelB <- reactive({
    lm(formB(), data = data_reactive())
  })

  # 8) Bootstrap tables
  bootA <- reactive({
    dat <- data_reactive()
    if (isTRUE(input$use_boot)) {
      boot_lm_coef(formA(), dat, B = input$B, seed = 111)
    } else {
      fit <- modelA()
      co <- summary(fit)$coefficients
      data.frame(
        term = rownames(co),
        estimate = co[, "Estimate"],
        se_ols = co[, "Std. Error"],
        se_boot = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        row.names = NULL
      )
    }
  })

  bootB <- reactive({
    dat <- data_reactive()
    if (isTRUE(input$use_boot)) {
      boot_lm_coef(formB(), dat, B = input$B, seed = 222)
    } else {
      fit <- modelB()
      co <- summary(fit)$coefficients
      data.frame(
        term = rownames(co),
        estimate = co[, "Estimate"],
        se_ols = co[, "Std. Error"],
        se_boot = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        row.names = NULL
      )
    }
  })

  output$bootA_tbl <- renderTable(bootA(), digits = 4)
  output$bootB_tbl <- renderTable(bootB(), digits = 4)

  output$coef_tbl <- renderTable({
    coef_compare(modelA(), modelB())
  }, rownames = FALSE, digits = 4)

  # -----------------------------
  # Model comparison outputs
  # -----------------------------

  # AIC/BIC always
  output$ic_tbl <- renderTable({
    mA <- modelA()
    mB <- modelB()
    data.frame(
      model = c("Model A", "Model B"),
      df    = c(df.residual(mA), df.residual(mB)),
      AIC   = c(AIC(mA), AIC(mB)),
      BIC   = c(BIC(mA), BIC(mB)),
      row.names = NULL
    )
  }, digits = 4)

  # Nested check + ANOVA (partial F-test) if nested
  output$nested_txt <- renderText({
    mA <- modelA()
    mB <- modelB()
    nested <- is_nested_lm(mA, mB)
    if (nested) {
      "Model A appears nested in Model B (same outcome, all terms in A are contained in B). Partial F-test is shown below."
    } else {
      "Models do NOT appear nested (or outcomes differ). Partial F-test via anova(A, B) is not valid/meaningful here."
    }
  })

  output$anova_tbl <- renderTable({
    mA <- modelA()
    mB <- modelB()
    if (!is_nested_lm(mA, mB)) return(NULL)

    # anova gives the partial F-test for nested lm models
    an <- anova(mA, mB)
    # make it a plain data.frame for tableOutput
    out <- as.data.frame(an)
    out$Model <- rownames(out)
    out <- out[, c("Model", setdiff(names(out), "Model")), drop = FALSE]
    rownames(out) <- NULL
    out
  }, digits = 4)

  # Optional K-fold CV comparison
  output$cv_tbl <- renderTable({
    req(input$use_cv)
    dat <- data_reactive()

    resA <- cv_lm_metrics(formA(), dat, K = input$K, seed = input$cv_seed)
    resB <- cv_lm_metrics(formB(), dat, K = input$K, seed = input$cv_seed)

    merge(
      transform(resA, model = "Model A"),
      transform(resB, model = "Model B"),
      by = c("metric"),
      suffixes = c("_A", "_B"),
      all = TRUE
    ) |> within({
      # Clean layout: metric, Model A, Model B
      Model_A <- value_A
      Model_B <- value_B
      value_A <- NULL
      value_B <- NULL
    }) |> subset(select = c(metric, Model_A, Model_B))
  }, digits = 4)

  # -----------------------------
  # Partial plots (side-by-side)
  # -----------------------------
  output$pp_side <- renderPlot({
    dat <- data_reactive()
    modA <- modelA()
    modB <- modelB()

    varsA <- attr(terms(modA), "term.labels")
    varsB <- attr(terms(modB), "term.labels")

    vars_union <- sort(unique(c(varsA, varsB)))
    if (length(vars_union) == 0) {
      plot.new()
      text(0.5, 0.5, "No predictors selected in either model.")
      return(invisible(NULL))
    }

    par(mfrow = c(length(vars_union), 2), mar = c(4, 4, 3, 1))

    for (v in vars_union) {

      # Model A
      if (v %in% varsA) {
        if (isTRUE(input$use_boot)) {
          avplot_boot(modA, v, dat, B = input$B, seed = 1000 + match(v, vars_union), show_band = TRUE)
        } else {
          others <- setdiff(varsA, v)
          yname <- all.vars(formula(modA))[1]
          f_y <- as.formula(paste(yname, "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))
          f_x <- as.formula(paste(v, "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))
          y_res <- resid(lm(f_y, data = dat))
          x_res <- resid(lm(f_x, data = dat))
          fit <- lm(y_res ~ x_res)
          plot(x_res, y_res, xlab = paste0(v, " | others"), ylab = paste0(yname, " | others"), pch = 16, cex = 0.7)
          abline(fit, lwd = 2)
        }
        title(main = paste0("Model A: ", v))
      } else {
        plot.new()
        text(0.5, 0.5, paste0("Model A:\n", v, " not included"))
        title(main = paste0("Model A: ", v))
      }

      # Model B
      if (v %in% varsB) {
        if (isTRUE(input$use_boot)) {
          avplot_boot(modB, v, dat, B = input$B, seed = 2000 + match(v, vars_union), show_band = TRUE)
        } else {
          others <- setdiff(varsB, v)
          yname <- all.vars(formula(modB))[1]
          f_y <- as.formula(paste(yname, "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))
          f_x <- as.formula(paste(v, "~", if (length(others) == 0) "1" else paste(others, collapse = " + ")))
          y_res <- resid(lm(f_y, data = dat))
          x_res <- resid(lm(f_x, data = dat))
          fit <- lm(y_res ~ x_res)
          plot(x_res, y_res, xlab = paste0(v, " | others"), ylab = paste0(yname, " | others"), pch = 16, cex = 0.7)
          abline(fit, lwd = 2)
        }
        title(main = paste0("Model B: ", v))
      } else {
        plot.new()
        text(0.5, 0.5, paste0("Model B:\n", v, " not included"))
        title(main = paste0("Model B: ", v))
      }
    }
  })
}

shinyApp(ui, server)









