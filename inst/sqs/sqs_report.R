###############################################################################
# sqs_report.R
#
# Utilidades INDEPENDIENTES para leer, consolidar y visualizar los logs CSV
# generados por el filtro SQS en internalFunctions_1.r
# 
# NO depende de internalFunctions_1.r. No altera ningún pipeline.
# Uso típico:
#   source("sqs_report.R")
#   df   <- load_sqs_logs(ref_dir = "results/Ref_myrun",
#                         query_dir = "results/Ref_myrun/Query_ALL")
#   tbl  <- summarize_sqs(df)
#   p1   <- plot_sqs_scatter(df)
#   p2   <- plot_sqs_heatmap(df)
#   save_sqs_report(df, out_dir = "results/SQS_report")
###############################################################################

# -----------------------------------------------------------------------------
# Carga y consolidación
# -----------------------------------------------------------------------------

#' Carga uno o varios CSV de log SQS y devuelve un data.frame consolidado.
#'
#' @param ref_dir   Directorio que contiene *_SQS_log.csv producidos en la ref.
#'                  Puede ser NULL.
#' @param query_dir Directorio que contiene *_SQS_log.csv producidos en la query.
#'                  Puede ser NULL.
#' @param files     Alternativamente, vector de rutas directas a CSVs.
#' @return data.frame con columnas:
#'   stage, file.id, marker, sample, sqs, threshold, apply_sqs, action, passed
load_sqs_logs <- function(ref_dir = NULL, query_dir = NULL, files = NULL){
  # acepta los tres nombres: legado (_SQS_log.csv) y nuevos (_SQS_ref_log.csv, _SQS_query_log.csv)
  pat <- "_SQS(_ref|_query)?_log\\.csv$"
  paths <- character(0)
  if (!is.null(ref_dir) && dir.exists(ref_dir)) {
    paths <- c(paths, list.files(ref_dir, pattern = pat,
                                 full.names = TRUE, recursive = FALSE))
  }
  if (!is.null(query_dir) && dir.exists(query_dir)) {
    paths <- c(paths, list.files(query_dir, pattern = pat,
                                 full.names = TRUE, recursive = FALSE))
  }
  if (!is.null(files)) {
    paths <- c(paths, files[file.exists(files)])
  }
  paths <- unique(paths)
  if (length(paths) == 0) {
    warning("No se encontraron CSV de log SQS en las rutas indicadas.")
    return(data.frame())
  }
  dfs <- lapply(paths, function(p) {
    tryCatch(
      utils::read.csv(p, stringsAsFactors = FALSE, na.strings = c("","NA")),
      error = function(e) { warning("No se pudo leer ", p, ": ", conditionMessage(e)); NULL }
    )
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0) return(data.frame())
  df <- do.call(rbind, dfs)
  # columna booleana derivada
  df$passed <- df$sqs >= df$threshold
  # tipos consistentes
  df$stage      <- as.character(df$stage)
  df$marker     <- as.character(df$marker)
  df$sample     <- as.character(df$sample)
  df$sqs        <- suppressWarnings(as.numeric(df$sqs))
  df$threshold  <- suppressWarnings(as.numeric(df$threshold))
  df$apply_sqs  <- as.character(df$apply_sqs)
  df$action     <- as.character(df$action)
  df
}

# -----------------------------------------------------------------------------
# Resumen tabular
# -----------------------------------------------------------------------------

#' Resumen por marcador (y opcionalmente por stage) del SQS.
#'
#' @param df     data.frame producido por load_sqs_logs().
#' @param by_stage TRUE para agrupar también por stage (ref/query).
#' @return data.frame ordenado por %pass ascendente.
summarize_sqs <- function(df, by_stage = TRUE){
  if (nrow(df) == 0) return(df)
  grp_cols <- if (by_stage) c("stage","marker") else "marker"
  # split por claves
  key <- do.call(paste, c(df[, grp_cols, drop = FALSE], sep = "||"))
  split_df <- split(df, key)
  out <- do.call(rbind, lapply(split_df, function(sub) {
    base <- sub[1, grp_cols, drop = FALSE]
    data.frame(
      base,
      n_samples   = nrow(sub),
      n_pass      = sum(sub$passed, na.rm = TRUE),
      pct_pass    = round(100 * mean(sub$passed, na.rm = TRUE), 2),
      sqs_min     = round(min(sub$sqs,    na.rm = TRUE), 3),
      sqs_median  = round(stats::median(sub$sqs, na.rm = TRUE), 3),
      sqs_mean    = round(mean(sub$sqs,   na.rm = TRUE), 3),
      sqs_max     = round(max(sub$sqs,    na.rm = TRUE), 3),
      threshold   = sub$threshold[1],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out <- out[order(out$pct_pass, out$sqs_median), , drop = FALSE]
  out
}

# -----------------------------------------------------------------------------
# Plots: ggplot2 si está disponible; si no, base R
# -----------------------------------------------------------------------------

.have_ggplot <- function() requireNamespace("ggplot2", quietly = TRUE)

#' Scatter/jitter: SQS por marcador, coloreado por stage, con línea de umbral.
#'
#' @param df       data.frame de load_sqs_logs().
#' @param markers  vector opcional de marcadores a incluir (por defecto todos).
#' @return objeto ggplot o NULL (base R se dibuja directamente por side effect).
plot_sqs_scatter <- function(df, markers = NULL){
  if (nrow(df) == 0) { message("df vacío"); return(invisible(NULL)) }
  if (!is.null(markers)) df <- df[df$marker %in% markers, , drop = FALSE]
  thr <- stats::median(df$threshold, na.rm = TRUE)
  if (.have_ggplot()) {
    ggplot2::ggplot(df, ggplot2::aes(x = marker, y = sqs, colour = stage)) +
      ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 1.8) +
      ggplot2::geom_hline(yintercept = thr, linetype = "dashed") +
      ggplot2::labs(x = NULL, y = "SQS",
                    title = "Stain Quality Score por marcador",
                    subtitle = paste0("Umbral = ", thr,
                                      "; puntos bajo la línea = muestras suspensas")) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    op <- graphics::par(mar = c(7, 4, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    markers_u <- sort(unique(df$marker))
    stages_u  <- sort(unique(df$stage))
    col_pal <- c("ref" = "#185FA5", "query" = "#993C1D")
    col_pal <- col_pal[stages_u]
    graphics::plot(NA, xlim = c(0.5, length(markers_u) + 0.5),
                   ylim = range(df$sqs, na.rm = TRUE),
                   xaxt = "n", xlab = "", ylab = "SQS",
                   main = "Stain Quality Score por marcador")
    graphics::axis(1, at = seq_along(markers_u), labels = markers_u, las = 2)
    for (s in stages_u) {
      sub <- df[df$stage == s, , drop = FALSE]
      x <- match(sub$marker, markers_u) + stats::runif(nrow(sub), -0.2, 0.2)
      graphics::points(x, sub$sqs, col = col_pal[s], pch = 16)
    }
    graphics::abline(h = thr, lty = 2)
    graphics::legend("topright", legend = stages_u, col = col_pal,
                     pch = 16, bty = "n")
    invisible(NULL)
  }
}

#' Heatmap sample × marker con SQS.
#'
#' @param df     data.frame de load_sqs_logs().
#' @param stage  "ref", "query" o NULL (concatena todo con sufijo).
#' @return objeto ggplot (si disponible) o NULL (base R image por side effect).
plot_sqs_heatmap <- function(df, stage = NULL){
  if (nrow(df) == 0) { message("df vacío"); return(invisible(NULL)) }
  if (!is.null(stage)) df <- df[df$stage == stage, , drop = FALSE]
  if (nrow(df) == 0) { message("sin filas para stage=", stage); return(invisible(NULL)) }
  thr <- stats::median(df$threshold, na.rm = TRUE)
  
  # etiqueta única por fila si hay ambos stages
  df$.row_label <- if (length(unique(df$stage)) > 1) paste(df$sample, df$stage, sep = "|") else df$sample
  
  if (.have_ggplot()) {
    ggplot2::ggplot(df, ggplot2::aes(x = marker, y = .row_label, fill = sqs)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::scale_fill_gradient2(low = "#A32D2D", mid = "#F1EFE8",
                                    high = "#0F6E56",
                                    midpoint = thr, name = "SQS") +
      ggplot2::labs(x = NULL, y = NULL,
                    title = "Heatmap SQS (sample × marker)",
                    subtitle = paste0("Umbral = ", thr)) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    # pivot manual
    markers_u <- sort(unique(df$marker))
    rows_u    <- sort(unique(df$.row_label))
    mat <- matrix(NA_real_, nrow = length(rows_u), ncol = length(markers_u),
                  dimnames = list(rows_u, markers_u))
    for (i in seq_len(nrow(df))) {
      mat[df$.row_label[i], df$marker[i]] <- df$sqs[i]
    }
    op <- graphics::par(mar = c(7, 10, 3, 2))
    on.exit(graphics::par(op), add = TRUE)
    graphics::image(t(mat), axes = FALSE,
                    main = "Heatmap SQS (sample × marker)")
    graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)),
                   labels = colnames(mat), las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)),
                   labels = rownames(mat), las = 2, cex.axis = 0.7)
    invisible(NULL)
  }
}

#' Histograma de SQS global, con línea de umbral.
plot_sqs_hist <- function(df, bins = 40){
  if (nrow(df) == 0) return(invisible(NULL))
  thr <- stats::median(df$threshold, na.rm = TRUE)
  if (.have_ggplot()) {
    ggplot2::ggplot(df, ggplot2::aes(x = sqs, fill = stage)) +
      ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = bins) +
      ggplot2::geom_vline(xintercept = thr, linetype = "dashed") +
      ggplot2::labs(x = "SQS", y = "frecuencia",
                    title = "Distribución global del SQS",
                    subtitle = paste0("Umbral = ", thr)) +
      ggplot2::theme_minimal(base_size = 11)
  } else {
    graphics::hist(df$sqs, breaks = bins, main = "Distribución global del SQS",
                   xlab = "SQS")
    graphics::abline(v = thr, lty = 2)
    invisible(NULL)
  }
}

# -----------------------------------------------------------------------------
# Export todo-en-uno
# -----------------------------------------------------------------------------

#' Guarda en un directorio: CSV consolidado + tabla resumen + 3 plots (PDF).
#'
#' @param df       data.frame de load_sqs_logs().
#' @param out_dir  Directorio donde se escribirán los ficheros.
#' @param width    Ancho de los PDFs en pulgadas.
#' @param height   Alto de los PDFs en pulgadas.
#' @return out_dir (invisible).
save_sqs_report <- function(df, out_dir, width = 10, height = 7){
  if (nrow(df) == 0) { message("df vacío; nada que guardar."); return(invisible(NULL)) }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, file.path(out_dir, "sqs_consolidated.csv"), row.names = FALSE)
  utils::write.csv(summarize_sqs(df, by_stage = TRUE),
                   file.path(out_dir, "sqs_summary_by_stage_marker.csv"), row.names = FALSE)
  utils::write.csv(summarize_sqs(df, by_stage = FALSE),
                   file.path(out_dir, "sqs_summary_by_marker.csv"), row.names = FALSE)
  
  .save_one <- function(fname, plt){
    grDevices::pdf(file.path(out_dir, fname), width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (inherits(plt, "ggplot")) print(plt) else if (is.function(plt)) plt()
  }
  if (.have_ggplot()) {
    .save_one("sqs_scatter.pdf", plot_sqs_scatter(df))
    .save_one("sqs_heatmap.pdf", plot_sqs_heatmap(df))
    .save_one("sqs_hist.pdf",    plot_sqs_hist(df))
  } else {
    .save_one("sqs_scatter.pdf", function() plot_sqs_scatter(df))
    .save_one("sqs_heatmap.pdf", function() plot_sqs_heatmap(df))
    .save_one("sqs_hist.pdf",    function() plot_sqs_hist(df))
  }
  invisible(out_dir)
}
