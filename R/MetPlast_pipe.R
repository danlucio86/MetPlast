#' Metabolic Diversity Pipeline
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @param Data_raw a data frame with the levels of the detected metabolites without replacing NAs values
#' @param ... any other argument
#'
#' @description MetDiv_pipe calculates and summarises all the calculation and plots related with Metabolic Diversity.
#'
#' @details MetDiv_pipe pipes the functions MetDiv(), MetDiv_plot(), and MetDiv_Supp_plot(). It takes as an input the original data set.
#' All the outputs can be easily unlisted, and extracted, and plotted.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 5 objects:
#' 1. A data frame with the Hj value and number of peaks (compounds) per species;
#' 2. A boxplot despicting the variation of Hj per species;
#' 3. A dot plot despicting the dependency between the number of peaks and Hj
#' 4. A dot plot despicting the dependency beween the number of peaks and Species;
#' 5. A grid with all the plots
#'
#'
#' @export
#'
#' @examples MetDiv_pipe <- MetDiv_pipe (Data, Data_raw)

MetDiv_pipe <- function(Data, Data_raw, ...){
  #Mssg
  print("This function calculates Shannon entropy as a direct measure of Metabolic Diversity")

  print("This function creates plots depicting the relatioship between the different variables")

  #Warnings

  if (nrow(Data)<ncol(Data)) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Data) == FALSE) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  #Stop
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }

  if (exists("Pij_fc") == FALSE) {
    stop('Internal function Pij_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  if (exists("Hj_fc") == FALSE) {
    stop('Internal function numb_peaks_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }


  # Calculate Hj
  Pij <- Pij_fc (Data)
  Hj <- apply(Pij, MARGIN = 2, FUN = Hj_fc )

  # Number of peaks
  numb_peaks <- numb_peaks_fc (Data_raw)

  # Output generation
  Hj <- cbind (Hj, numb_peaks)
  Hj_df <- as.data.frame (Hj)
  Hj_Sample_df <- setDT(Hj_df, keep.rownames = "Sample")
  Hj_Sample_df$Sample <- sub("\\_?\\.?\\d+", "", Hj_Sample_df$Sample)
  Hj_Sample_df$Sample <- sub("\\..", ".", Hj_Sample_df$Sample)

  # Graphical output
  p1 <- ggplot(Hj_Sample_df, aes(Hj, Sample, color = Sample)) + geom_boxplot() + geom_point() + theme(legend.position = "none")
  print(p1)
  p2 <- ggplot(Hj_Sample_df, aes(numb_peaks, Hj, color = Sample)) + geom_point()
  print(p2)
  p3 <- ggplot(Hj_Sample_df, aes(numb_peaks, Sample, color = Hj)) + geom_point() + theme(legend.position = "none")
  print(p3)
  grid_plot <- grid.arrange (p1, p2, p3)

  # List
  list <- list(grid_plot, Hj_Sample_df)
  print (list)
}


#' Metabolite (Si) and Metabolome (δj) Specialization Index Pipeline
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @param ... any other argument
#' @description MetSpec_pipe calculates METabolite SPECialization (Si) and METabolome SPECialization (δj) indexes based on Shannon entropy.
#'
#' @details
#' Metabolic specificity (Si) is defined as the specificity of a particular metabolite (i) among a set of samples (j).
#' Metabolome specialization δj is measured for each jth sample, as the average of the metabolites specificity.
#' All the outputs can be easily unlisted, and extracted, and plotted.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 4 objects:
#' 1. A data frame with the Si per compound;
#' 2. A data frame with the δj per Species;
#' 3. A dot plot depicting the Si value of each compound
#' 4. A dot plot depicting the δj per species
#'
#'
#' @export
#'
#' @examples Dj <- MetSpec_pipe (Data)

MetSpec_pipe <- function (Data, ...){
  # Mssg
  print("THIS PARAMETER CAN NOT BE EXTRAPOLATED TO OTHER DATA SETS")

  print("This function creates plots depicting the relatioship between the different variables")


  #Warnings
  if (nrow(Data)<ncol(Data)) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Data) == FALSE) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  #Stop
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }

  if (exists("Pij_by_Pi_fc") == FALSE) {
    stop('Internal function Pij_by_Pi_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  if (exists("Si_fc") == FALSE) {
    stop('Internal function Si_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  # Calculate metabolite specialization
  Pij_by_Pi <- Pij_by_Pi_fc (Data)
  Si <- apply(Pij_by_Pi, MARGIN = 1, FUN = Si_fc)

  # Ouput metabolite specialization
  Si_df <- data.frame (Si)
  Si_df <- setDT(Si_df, keep.rownames = "Compound")

  # Graphical output
  p1 <- ggplot(Si_df, aes(Compound, Si)) + geom_point()
  print(p1)

  # Calculate metabolome specialization
  Pij <- Pij_fc (Data)
  Pij_Si <- sweep(Pij, MARGIN = 1, Si, FUN = "*")
  Pij_Si_df <- data.frame(Pij_Si)
  Dj <- apply(Pij_Si_df, MARGIN = 2,  FUN = sum)

  #Output metabolome specilization
  Dj_df <- data.frame (Dj)
  Dj_Sample_df <- setDT(Dj_df, keep.rownames = "Sample")
  Dj_Sample_df$Sample <- sub("\\_?\\.?\\d+", "", Dj_Sample_df$Sample)
  Dj_Sample_df$Sample <- sub("\\..", ".", Dj_Sample_df$Sample)
  print(Dj_Sample_df)

  # Graphical output
  p2 <- ggplot(Dj_Sample_df, aes(Sample, Dj)) + geom_boxplot() + geom_point()
  print (p2)
  grid_plot <- grid.arrange (p1, p2)

  #List
  list <- list(grid_plot, Si_df, Dj_Sample_df)
  print(list)
}


#' Metabolite Specialization Analysis Pipeline
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @param ... any other argument
#'
#' @description MetliteSpec_pipe calculates the contribution of Metabolite Specialization factor (Pij*Si) to the Metabolome Specialization (δj) index.
#'
#' @details
#' Metabolite specialization factor (Pij*Si) is defined as product of the a metabolite specilization index (Si) and the frequency of the metabolite in a given sample (Pij).
#' All the outputs can be easily unlisted, and extracted, and plotted.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 4 objects:
#' 1. A data frame with the Pij*Si values per species and compound;
#' 2. A data frame with the highest Pij*Si compound value per species;
#' 3. A dot plot depicting the Pij*Si value of each species
#' 4. A dot plot depicting the highest Pij*Si value per species
#'
#'
#' @export
#'
#' @examples MetliteSpec_pipe <- MetliteSpec_pipe (Data)


MetliteSpec_pipe <- function (Data, ...){

  # Mssg
  print("THIS PARAMETER CAN NOT BE EXTRAPOLATED TO OTHER DATA SETS ")

  print("This function creates plots depicting the relatioship between the different variables")


  #Warnings
  if (nrow(Data)<ncol(Data)) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  if (is.data.frame(Data) == FALSE) {
    warning('MetPlast functions takes as Data argument a data frame with compounds as observations -rows- and individual samples as variables -columns-')
  }

  #Stop
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }


  if (exists("Pij_by_Pi_fc") == FALSE) {
    stop('Internal function Pij_by_Pi_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  if (exists("Pij_fc") == FALSE) {
    stop('Internal function Pij_fc is not uploaded. Please, check MetPlast_int.R functions are properly uploaded.')
  }

  # Calculate metabolite specialization
  Pij_by_Pi <- Pij_by_Pi_fc (Data)
  Si <- apply(Pij_by_Pi, MARGIN = 1, FUN = Si_fc)

  # Si weight on Dj
  Pij <- Pij_fc (Data)
  Pij_Si <- sweep(Pij, MARGIN = 1, Si, FUN = "*")
  Pij_Si_t <- t(Pij_Si)
  Pij_Si_df <- as.data.frame (Pij_Si_t)
  Compounds <- colnames(Pij_Si_df)
  Pij_Si_df <- setDT(Pij_Si_df, keep.rownames = "Sample")
  Dj_index_weight<- pivot_longer(Pij_Si_df, all_of(Compounds), names_to = "Compounds", values_to = "Pij_Si")
  Dj_index_weight$Sample <- sub("\\_?\\.?\\d+", "", Dj_index_weight$Sample)
  Dj_index_weight$Sample <- sub("\\..", ".", Dj_index_weight$Sample)
  print (Dj_index_weight)
  Dj_index_weight_filtered <- my_filter_fc (Dj_index_weight)
  print (Dj_index_weight_filtered)

  # Graphical Output
  p <- ggplot(Dj_index_weight, aes(Sample, Pij_Si, color = Compounds)) + geom_point() + theme(legend.position = "none")
  print(p)
  q <- ggplot (Dj_index_weight_filtered, aes(Sample, Pij_Si, color = Compounds)) + geom_point() + geom_text(aes(label= Compounds), hjust=0,vjust=0, size=3) + theme(legend.position = "none")
  print(q)
  grid_plot <- grid.arrange (p,q)

  # List
  list <- list (grid_plot, Dj_index_weight, Dj_index_weight_filtered)
  print (list)
}

#' Metabolic Plasticity Parameters Pipeline
#'
#' @param Stats_df a data frame with statistical values for each samples, calculated from MetStats()
#' @param Dj_df  a data frame with Metabolome specialization index δj calculated from Dj_index()
#' @param ... any other argument
#'
#' @description MetPar_pipe integrates Hj, δj, numb_peaks, HRj, and Divj parameters, as well as graphically evaluate the relationship between some of them.
#'
#' @details
#' Metabolic Plasticity Parameters Pipe extracts and summarizes Hj, δj and statistical information. It generates dot plots to evaluate the relatioship between the Metabolome Specialization (δj) index and other metabolic parameters such as: Divj, and Hj.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It return a list with 3 objects:
#' 1. A data frame with all the Metabolic Parameters per sample;
#' 2. A dot plot depicting the relationship between the Hj and δj values of each sample
#' 3. A dot plot depicting the relationship between Divj and δj values of each sample
#'
#'
#' @export
#'
#' @examples
#' MetPar <- MetPar_pipe (Stats_df, Dj_df)

MetPar_pipe<- function(Stats_df, Dj_df, ...){

  #Message
  print("Well, it looks you got everything you wanted. Enjoy :)")

  # Warnings

  if (is.data.frame(Stats_df) == FALSE) {
    warning('MetPlast functions takes as Stats_df argument a data frame')
  }

  if (is.data.frame(Dj_df) == FALSE) {
    warning('MetPlast functions takes as Dj_df argument a data frame')
  }

  # Stop

  if (hasArg(Stats_df) == FALSE) {
    stop('Please, check Stats_df generated by MetStats(). Stats_df must be a data frame with five variables: Sample, numb_peaks, Hj, HRj, and Divj')
  }

  if (hasArg(Dj_df) == FALSE) {
    stop('Please, check Dj_df generated by Dj_index_weight(). Dj_df must be a data frame with two variables: Sample, and Dj')
  }

    # Extract the Hj and Dj parameters per species
    MetPar_df <- cbind(Stats_df, Dj = Dj_df$Dj)
    print(MetPar_df)

    # Evaluate the dependence between Hj and Dj
    p1 <- ggplot2::ggplot(MetPar_df, aes(Divj, y = Dj, color ="Sample")) + geom_point() + expand_limits(x = 0, y = 0) + geom_abline(slope=1, intercept=0)
    p2 <- ggplot2::ggplot(MetPar_df, aes(Hj, y = Dj, color ="Sample")) + geom_point()
    grid_plot <- grid.arrange (p1, p2)
    print(grid_plot)

    list <- list (MetPar_df, grid_plot)
    print(list)
  }
