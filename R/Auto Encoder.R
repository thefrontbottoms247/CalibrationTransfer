#' Autoencoder-Based Calibration Transfer (Robust + Weighted)
#'
#' Trains a denoising, L2-regularized autoencoder on parent and child spectra,
#' with optional wavelength weighting for chemically informative peaks,
#' then applies the learned mapping to correct child modeling spectra.
#'
#' @param parent.trans Numeric matrix or data frame of parent calibration spectra.
#' @param child.trans  Numeric matrix or data frame of child calibration spectra.
#' @param child.mod    Numeric matrix or data frame of child modeling spectra.
#' @param latent_dim Integer. Number of latent features to learn (default = 20).
#' @param epochs Integer. Number of training epochs (default = 200).
#' @param batch_size Integer. Batch size for training (default = 32).
#' @param peaks Optional logical or numeric index vector of peak wavelengths.
#'              If provided, the loss will weight those wavelengths higher.
#' @param peak_weight Numeric multiplier for peak weighting (default = 5).
#'
#' @return Invisibly returns training history. Assigns globally:
#' \itemize{
#'   \item \code{auto_latent} — latent encodings of all training data.
#'   \item \code{child_to_parent} — reconstructed child calibration spectra.
#'   \item \code{child.mod_cor} — corrected child modeling spectra.
#' }
#' @import keras3
#' @export
AutoEncode <- function(parent.trans, child.trans, child.mod,
                       latent_dim = 20, epochs = 200, batch_size = 32,
                       peaks = NULL, peak_weight = 5) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Install keras3: install.packages('keras3')")

  parent.trans <- as.matrix(parent.trans[, sapply(parent.trans, is.numeric)])
  child.trans  <- as.matrix(child.trans[,  sapply(child.trans,  is.numeric)])
  child.mod    <- as.matrix(child.mod[,    sapply(child.mod,    is.numeric)])
  if (ncol(parent.trans) != ncol(child.trans))
    stop("Parent and child calibration spectra must have same number of columns.")

  # --- Joint scaling ---
  x_all <- scale(rbind(parent.trans, child.trans))
  mu <- attr(x_all, "scaled:center")
  sd <- attr(x_all, "scaled:scale")

  n <- nrow(parent.trans)
  m <- ncol(parent.trans)

  # --- Optional peak-weighted loss ---
  # --- Optional peak-weighted loss (dimension-safe) ---
  weighted_mse <- NULL
  if (!is.null(peaks)) {

    # Ensure peaks refers to the actual column indices of the numeric matrix
    valid_cols <- sapply(parent.trans, is.numeric)
    m <- sum(valid_cols)

    # Build weight vector exactly matching model input size
    w_vec <- rep(1, m)
    if (is.logical(peaks)) peaks <- which(peaks)
    peaks <- peaks[peaks <= m]                # safety clip
    w_vec[peaks] <- peak_weight

    weighted_mse <- function(y_true, y_pred) {
      tf <- tensorflow::tf
      w <- tf$constant(matrix(w_vec, nrow = 1L), dtype = tf$float32)
      # broadcast weights across batch dimension
      err <- tf$square(y_true - y_pred)
      tf$reduce_mean(tf$reduce_mean(err * w, axis = 2L), axis = 1L)
    }
  }


  # --- Define denoising Autoencoder ---
  input <- keras3::layer_input(shape = m)
  encoded <- input |>
    keras3::layer_gaussian_noise(stddev = 0.02) |>
    keras3::layer_dense(units = 1024, activation = "relu",
                        kernel_regularizer = keras3::regularizer_l2(1e-4)) |>
    keras3::layer_dropout(0.15) |>
    keras3::layer_dense(units = 512, activation = "relu",
                        kernel_regularizer = keras3::regularizer_l2(1e-4)) |>
    keras3::layer_dense(units = latent_dim, activation = "linear",
                        activity_regularizer = keras3::regularizer_l2(1e-5))

  decoded <- encoded |>
    keras3::layer_dense(units = 512, activation = "relu",
                        kernel_regularizer = keras3::regularizer_l2(1e-4)) |>
    keras3::layer_dropout(0.15) |>
    keras3::layer_dense(units = 1024, activation = "relu",
                        kernel_regularizer = keras3::regularizer_l2(1e-4)) |>
    keras3::layer_dense(units = m, activation = "linear")

  autoencoder <- keras3::keras_model(inputs = input, outputs = decoded)

  # --- Compile with adaptive loss selection ---
  keras3::compile(
    autoencoder,
    optimizer = keras3::optimizer_adam(learning_rate = 1e-4),
    loss = if (is.null(weighted_mse)) "mse" else weighted_mse
  )

  # --- Callbacks for robust convergence ---
  callbacks <- list(
    keras3::callback_early_stopping(monitor = "val_loss",
                                    patience = 20,
                                    restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(monitor = "val_loss",
                                          factor = 0.5,
                                          patience = 8,
                                          min_lr = 1e-6)
  )

  # --- Train with validation monitoring ---
  hist <- keras3::fit(
    autoencoder,
    x = x_all,
    y = x_all,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = 0.2,
    callbacks = callbacks,
    verbose = 1
  )

  # --- Encode/Decode training data ---
  encoder <- keras3::keras_model(inputs = input, outputs = encoded)
  encoded_data <- keras3::predict(encoder, x_all)
  reconstructed <- keras3::predict(autoencoder, x_all)
  recon_child <- reconstructed[(n + 1):(2 * n), ]

  # --- Correct child.mod using learned mapping ---
  z <- sweep(child.mod, 2, mu, "-")
  z <- sweep(z, 2, sd, "/")
  z <- keras3::predict(autoencoder, z)
  child.mod_cor <- sweep(z, 2, sd, "*")
  child.mod_cor <- sweep(child.mod_cor, 2, mu, "+")

  # --- Export to Global Environment ---
  assign("auto_latent", encoded_data, envir = .GlobalEnv)
  assign("child_to_parent", recon_child, envir = .GlobalEnv)
  assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)
  assign("ae_autoencoder", autoencoder, envir = .GlobalEnv)
  assign("ae_encoder", encoder, envir = .GlobalEnv)
  assign("ae_center", mu, envir = .GlobalEnv)
  assign("ae_scale", sd, envir = .GlobalEnv)

  message("Autoencoder calibration transfer completed.
           Corrected spectra available as 'child.mod_cor'.")
  invisible(hist)
}
