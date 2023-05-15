#' Function to identify outliers
#' @export
#' @import multtest
#' @param resid  = null
#' @param alpha = 0.05
#' @source METHOD 2: Bonferroni-Holm using studentized residuals (BH-ST)
#' Bernal-Vasquez, AM., Utz, HF. & Piepho, HP. Theor Appl Genet (2016)
#' 129: 787. https://doi.org/10.1007/s00122-016-2666-6
#' @return outliers
#' @author Samuel Fernandes
outlier <- function(resid, alpha = 0.05) {
  # Produce externally studentized residuals
  studresid <- resid / sd(resid, na.rm = TRUE)
  # Calculate adjusted p-values
  rawp.BHStud = 2 * (1 - pnorm(abs(studresid)))
  #Produce a Bonferroni-Holm tests for the adjusted p-values
  #The output is a list
  test.BHStud <-
    multtest::mt.rawp2adjp(rawp.BHStud, proc = c("Holm"), alpha = alpha)
  #Create vectors/matrices out of the list of the BH tests
  adjp = cbind(test.BHStud[[1]][, 1])
  bholm = cbind(test.BHStud[[1]][, 2])
  index = test.BHStud$index
  # Condition to flag outliers according to the BH test
  out_flag = ifelse(bholm < alpha, "OUTLIER ", ".")
  #Create a matrix with all the output of the BH test
  BHStud_test = cbind(adjp, bholm, index, out_flag)
  #Order the file by index
  BHStud_test2 = BHStud_test[order(index), ]
  #Label colums
  names = c("rawp", "bholm", "index", "out_flag")
  colnames(BHStud_test2) <- names
  #Create a final file, with the data and the test and the labels for the outliers
  # Take a look at the outliers
  outliers_BH <-
    as.numeric(BHStud_test2[which(BHStud_test2[, "out_flag"] != "."), "index"])
  if (length(outliers_BH) < 1) {
    cat('No outlier detected \n')
  }
  if (length(outliers_BH) == 1) {
    cat('1 outlier detected! \n')
  }
  if (length(outliers_BH) > 1) {
    cat(length(outliers_BH), 'outliers detected! \n')
  }
  return(outliers_BH)
}