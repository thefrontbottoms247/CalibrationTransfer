#' Autoencoder-Based Calibration Transfer (Unified)
#'
#' Trains a nonlinear autoencoder on parent and child calibration spectra,
#' then applies the learned mapping to correct the child modeling spectra
#' into the parent instrument domain.
#'
#' @param parent.trans Numeric matrix or data frame of parent calibration spectra.
#' @param child.trans  Numeric matrix or data frame of child calibration spectra.
#' @param child.mod    Numeric matrix or data frame of child modeling spectra.
#' @param latent_dim Integer. Number of latent features to learn (default = 20).
#' @param epochs Integer. Number of training epochs (default = 200).
#' @param batch_size Integer. Batch size for training (default = 32).
#'
#' @return Invisibly returns training history. Assigns globally:
#' \itemize{
#'   \item \code{auto_latent} — latent encodings of all training data.
#'   \item \code{child_to_parent} — reconstructed child calibration spectra.
#'   \item \code{child.mod_cor} — corrected child modeling spectra.
#' }
#'
#' @import keras3
#' @importFrom stats predict
#' @export
AutoEncode <- function(parent.trans, child.trans, child.mod,
                           latent_dim = 20, epochs = 200, batch_size = 32) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Install keras3: install.packages('keras3')")

  # --- Ensure numeric matrices ---
  parent.trans <- as.matrix(parent.trans[, sapply(parent.trans, is.numeric)])
  child.trans  <- as.matrix(child.trans[,  sapply(child.trans,  is.numeric)])
  child.mod    <- as.matrix(child.mod[,    sapply(child.mod,    is.numeric)])
  if (ncol(parent.trans) != ncol(child.trans))
    stop("Parent and child calibration spectra must have same number of columns.")

  # --- Scale all transfer data together ---
  x_all <- scale(rbind(parent.trans, child.trans))
  mu <- attr(x_all, "scaled:center")
  sd <- attr(x_all, "scaled:scale")

  n <- nrow(parent.trans); m <- ncol(parent.trans)

  # --- Define Autoencoder ---
  input <- keras3::layer_input(shape = m)
  encoded <- input |>
    keras3::layer_dense(units = latent_dim, activation = "relu") |>
    keras3::layer_dense(units = latent_dim, activation = "relu")
  decoded <- encoded |>
    keras3::layer_dense(units = m, activation = "linear")

  autoencoder <- keras3::keras_model(inputs = input, outputs = decoded)
  keras3::compile(autoencoder,
                  optimizer = keras3::optimizer_adam(learning_rate = 0.001),
                  loss = "mse"
  )

  # --- Train Autoencoder ---
  hist <- keras3::fit(
    autoencoder,
    x = x_all,
    y = x_all,
    epochs = epochs,
    batch_size = batch_size,
    verbose = 1
  )

  # --- Encode/Decode training data ---
  encoder <- keras3::keras_model(inputs = input, outputs = encoded)
  encoded_data <- stats::predict(encoder, x_all)
  reconstructed <- stats::predict(autoencoder, x_all)
  recon_child <- reconstructed[(n + 1):(2 * n), ]

  # --- Correct child.mod using learned mapping ---
  z <- sweep(child.mod, 2, mu, "-")
  z <- sweep(z, 2, sd, "/")
  z <- as.matrix(z)
  child.mod_cor <- as.array(autoencoder(z))

  # --- Optional: unscale back to parent intensity space ---
  child.mod_cor <- sweep(child.mod_cor, 2, sd, "*")
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
