#' Metabolic Plasticity Statistics
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @param MetDiv_df a data frame from MetDiv()
#'
#' @description MetStats calculates sample (HRj) and Kullback–Leibler divergence (Divj)
#'
#' @details
#' Metabolic Plasticity Statistics calculates HRj and Divj considering Hj calculations per species
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame with the all the species-related statistical parameters: Samples, Hj, numb_peaks, HRj, and Divj.
#'
#' @export
#'
#' @examples
#' Stats <- MetStats (Data, MetDiv_df)


MetStats <- function(Data, MetDiv_df) {

  # Mssg
  print("MetStats calculates the divergence associated with Shannon entropy (Hj) ")

  #Warnings
  if (nrow(Data)<ncol(Data)) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Data) == FALSE) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(MetDiv_df) == FALSE) {
    warning('MetPlast functions takes as MetDiv_df argument a data frame')
  }

  #Stop
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }

  if (exists("HRj_fc") == FALSE) {
    stop('Internal function HRj_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  if (exists("Pij_fc") == FALSE) {
    stop('Internal function Pij_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  # Calculate divergence
  Pij <- Pij_fc (Data)
  HRj <- apply (Pij, MARGIN = 2, FUN = HRj_fc)

  # Calculate Kullback–Leibler divergence
  Divj <- HRj - MetDiv_df$Hj

  table <- cbind(MetDiv_df, HRj, Divj)
  MetStat_df <- as.data.frame(table)
  print(MetStat_df)
}

#' Metabolome Specialization index Weight Values
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @param Si_df a data frame with the Si values -column- for all the detected compounds - rows-, calculated from Si_index()
#'
#' @description Dj_index_weight calculates the contribution of METaboLITE SPECialization factor (Pij*Si) to the Metabolome specialization index.
#'
#' @details
#' Metabolite specialization factor (Pij*Si) is defined as product of the a metabolite specilization index (Si) and the frequency of the metabolite in a given sample (Pij).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame with the Pij.Si values per species and compound;
#'
#' @export
#'
#' @examples Dj_weight_df <- Dj_index_weight (Data, Si_df)


Dj_index_weight <- function (Data, Si_df){
  # Mssg
  print("THIS PARAMETER CAN NOT BE EXTRAPOLATED TO OTHER DATA SETS OR SUBSETS")

  #Warnings
  if (nrow(Data)<ncol(Data)) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Data) == FALSE) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Si_df) == FALSE) {
    warning('MetPlast functions takes as Si_df argument a data frame')
  }

  #Stop
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }

  if (exists("Pij_by_Pi_fc") == FALSE) {
    stop('Internal function Pij_by_Pi_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  # Calculate metabolite specialization
  Pij_by_Pi <- Pij_by_Pi_fc (Data)
  Si_value <- Si_df$Si

  # Si weight on Dj
  Pij <- Pij_fc (Data)
  Pij_Si <- sweep(Pij, MARGIN = 1, Si_value, FUN = "*")
  Pij_Si_t <- t(Pij_Si)
  Pij_Si_df <- as.data.frame (Pij_Si_t)
  Compounds <- colnames(Pij_Si_df)
  Pij_Si_df <- setDT(Pij_Si_df, keep.rownames = "Sample")
  Tidy_table <- pivot_longer(Pij_Si_df, all_of(Compounds), names_to = "Compounds", values_to = "Pij_Si")
  Tidy_table$Sample <- sub("\\_?\\.?\\d+", "", Tidy_table$Sample)
  Tidy_table$Sample <- sub("\\..", ".", Tidy_table$Sample)
  print (Tidy_table)
}

#' Metabolome Specialization index Weight plot
#'
#' @param Dj_weight_df
#'
#' @description Dj_index_weight plots the contribution of metabolite specialization factor (Pij*Si) to the Metabolome Specialization (Dj) index.
#'
#' @details
#' Metabolome Specialization index Weight plot generates a graph that shows the contribution of each metabolite to the Metabolome Specialization (Dj) index.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dot plot
#'
#' @export
#'
#' @examples Dj_weight_plot <- Dj_index_weight_plot (Dj_weight_df)


Dj_index_weight_plot <- function(Dj_weight_df){

  #Warnings
  if (is.data.frame(Dj_weight_df) == FALSE) {
    warning('Dj_index_weight_plot functions takes as Dj_weight_df argument a data frame')
  }

  #Stops
  if (hasArg(Dj_weight_df) == FALSE) {
    stop('Dj_weight_df should be provided. Please, run and store Dj_index_weight(). Dj_weight_df <- Dj_index_weight()')
  }

  ggplot2::ggplot(Dj_weight_df, aes(x=.data[["Sample"]], y = .data[["Pij_Si"]], color = .data[["Sample"]])) + geom_point() + geom_boxplot()

}

#' Metabolic Plasticity Parameters
#'
#' @param Stats_df a data frame with statistical values for each samples, calculated from MetStats()
#' @param Dj_df  a data frame with Metabolome specialization index δj calculated from Dj_index()
#'
#' @description MetPar integrates Hj, δj, numb_peaks, HRj, and Divj parameters
#'
#' @details
#' Metabolic Plasticity Parameters extracts and summarizes Hj, δj and statistical information.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame
#'
#' @export
#'
#' @examples
#' MetPar <- MetPar_df (Stat, Dj_df)

MetPar_df <- function(Stats_df, Dj_df){

  # Warnings

  if (is.data.frame(Stats_df) == FALSE) {
    warning('MetPlast functions takes as Stats_df argument a data frame')
  }

  if (is.data.frame(Dj_df) == FALSE) {
    warning('MetPlast functions takes as Dj_df argument a data frame')
  }

  # Stop

  if (hasArg("Stats_df") == FALSE) {
    stop('Please, check Stats_df generated by MetStats(). Stats_df must be a data frame with five variables: Sample, numb_peaks, Hj, HRj, and Divj')
  }

  if (hasArg("Dj_df") == FALSE) {
    stop('Please, check Dj_df generated by Dj_index_weight(). Dj_df must be a data frame with two variables: Sample, and Dj')
  }
  # Extract the Hj and Dj parameters per species
  MetPar_df <- cbind(Stats_df, Dj = Dj_df$Dj)

  print(MetPar_df)
}


#' Metabolic Plasticity Parameters Plot
#'
#' @param MetPar_df a data frame with all the metabolic plasticity paramenters obtained by MetPar_df ()
#' @param x_element column names from MetPar_df: "Divj", "Sample", or "Hj"
#' @param ref_line a logical argument. If TRUE, draws a x=y line. DEFAULT= FALSE
#'
#' @description MetPar allows to graphically evaluate the relationship between Hj, δj, numb_peaks, HRj, and Divj parameters
#'
#' @details
#' Metabolic Plasticity Parameters Plot allows to generate plots to evaluate the relatioship between the Metabolome Specialization (δj) index and other metabolic parameters such as: Divj, Sample, Hj.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dot plot
#'
#' @export
#'
#' @examples
#' MetPar_plot <- MetPar_plot (MetPar_df, Divj, ref_line = "TRUE")


MetPar_plot <- function (MetPar_df, x_element, ref_line = FALSE) {
  #mssg
  if(ref_line != TRUE) {
    print("DEFAULT, ref_line=FALSE. Please, for x=y reference line, select = TRUE")
  }
  #Warnings
  if (is.data.frame(MetPar_df) == FALSE) {
    warning('MetPar_plot functions takes MetPar_df as argument a data frame')
  }
  if (hasArg(x_element) == FALSE) {
    warning('Please, provide y element: Hj, Divj, Sample')
  }
  #Stops
  if (hasArg(MetPar_df) == FALSE) {
    stop('MetPar_df should be provided. Please, run and store MetPar_df(). MetPar_df <- MetPar_df()')
  }
  p <- ggplot2::ggplot(MetPar_df, aes(x={{ x_element }}, y = .data[["Dj"]], color =.data[["Sample"]])) + geom_point()

  if(ref_line == "TRUE") {p = p + expand_limits(x = 0, y = 0) + geom_abline(slope=1, intercept=0)
    }
  print (p)
  }
