#' Metabolite frequency per sample internal function
#'
#' @param Data
#' @description Internal function to calculate the frequency of metabolite i in species j (Pij)
#'
#' @details
#' The frequency is calculated dividing the intensity of a given metabolite (i) by the sum of intensities in a given species (j).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dataframe with metabolite frequency

Pij_fc <- function (Data) {Total_species <- apply (Data, MARGIN=2, FUN = sum)
Pij <- sweep(Data, MARGIN = 2, Total_species, FUN = '/')}

#' Metabolite expression identity internal function
#'
#' @param Data
#' @description The expression identity of a given MS/MS regarding frequencies among considered samples.
#'
#' @details
#' The metabolite expression identity is calculated dividing the metabolite frequency per sample by the frequency of the given metabolite in the whole dataset
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dataframe with metabolite expression identity

Pij_by_Pi_fc <- function (Data) {Total_Compound <- apply (Data, MARGIN=1, FUN = sum)
Pij <- Pij_fc (Data)
n_samples <- ncol(Data)
Pi <- (1/n_samples)*(apply (Pij, MARGIN=1, FUN=sum))
Pij_divided_Pi <- sweep(Pij, MARGIN = 1, Pi, FUN = "/")
}

#' Metabolic diversity internal function
#'
#' @param x
#' @description It performs internal metabolic diversity (Hj) calculations.
#'
#' @details
#' The metabolic diversity internal function performs an internal calculation necessary to obtained the metabolic diversity per species (Hj).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dataframe

Hj_fc <- function(x){-sum(x * log2(x), na.rm=TRUE)}

#' Metabolic diversity divergence internal function
#'
#' @param x
#' @description It performs internal metabolic diversity divergence (HRj) calculations.
#'
#' @details
#' The metabolic diversity divergence internal function performs an internal calculation necessary to obtained the metabolic diversity divergence (HRj).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dataframe

HRj_function <- function(x){
  Pij <- Pij_fc (Data)
  n_samples <- ncol(Data)
  Pi <- (1/n_samples)*(apply (Pij, MARGIN=1, FUN=sum))
  -sum(x*log2(Pi))
}

#' Metabolite Specialization internal function
#'
#' @param x
#' @description It performs internal metabolite specialization (Si) calculations.
#'
#' @details
#' The metabolite specialization internal function performs an internal calculation necessary to obtained the metabolite specialization (Si).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a dataframe

Si_fc <- function(x){n_samples <- ncol(Data)
(1/n_samples)*(sum(x*log2(x), na.rm =TRUE))}

#' Number of detected compounds function
#'
#' @param x
#' @description It counts the number of detected compounds per species
#'
#' @details
#' Number of detected compounds function counts the number of compound per species that have a value higher to the mininum. (Initially, NA values are replaced by the minimum value of the whole data set divided by 1000000 )
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#'
#' @return It returns a vector with the number of peaks - e.g. detected compounds - per species.

numb_peaks_fc <- function(Data) {min_data <- min(Data)
Data_numb_peaks <- apply(Data, MARGIN = 2, FUN = function(x){x != min_data})
numb_peaks <- apply (Data_numb_peaks, MARGIN = 2, FUN = sum)}

#' Metabolite contribution to metabolome specialization index (δj) function
#'
#' @param Table
#' @description It selects the metabolites that contribute the most to the metabolome specialization index (δj).
#'
#' @details
#' Metabolite contribution to metabolome specialization index (δj) function select the metabolite with the highest Pij*Si value per species.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a table


my_filter_function <- function(table) {
  table %>%
    group_by(Species) %>%
    filter(Pij_Si == max (Pij_Si))
}



#' Metabolic Diversity Index (Hj)
#'
#' @param Data
#' @description MetDiv calculates METabolic Diversity (Hj) index based on Shannon entropy.
#'
#' @details
#' The metabolic profile diversity is defined as the Shannon entropy using MS/MS metabolite frequency distribution in a sample (j) -Hj index- (Martinez et al paper). Hj can take any value between zero when only one metabolite is detected up to log2(m), where all m metabolites are detected and accumulates at the same frequency: 1/m (Martinez and Reyes-Valdes, 2008).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 5 objects:
#' 1. A data frame with the Hj value and number of peaks (compounds) per species;
#' 2. A boxplot despicting the variation of Hj per species;
#' 3. A point plot despicting the dependency between the number of peaks and Hj
#' 4. A point plot despicting the dependency beween the number of peaks and Species;
#' 5. A grid with all the plots
#'
#' @export
#'
#' @examples Hj <- MetDiv (Data)

MetDiv <- function(Data){
  #Mssg
  print("To use MetPar function store this calculation as Hj <- MetDiv (Data)")

  #WARNING
  if (any(is.na(Data))) {
    stop('NA values MUST be replace with Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000')
  }

  # Calculate Hj
  Pij <- Pij_fc (Data)
  Hj <- apply(Pij, MARGIN = 2, FUN = Hj_fc )

  # Number of peaks
  numb_peaks <- numb_peaks_fc (Data)

  # Output generation
  Hj <- cbind (Hj, numb_peaks)
  Hj_df <- as.data.frame (Hj)
  Hj_Species_df <- setDT(Hj_df, keep.rownames = "Species")
  Hj_Species_df$Species <- sub("\\_?\\.?\\d+", "", Hj_Species_df$Species)
  Hj_Species_df$Species <- sub("\\..", ".", Hj_Species_df$Species)


  # Graphical output
  p1 <- ggplot(Hj_Species_df, aes(Hj, Species, color = Species)) + geom_boxplot() + geom_point() + theme(legend.position = "none")
  print(p1)
  p2 <- ggplot(Hj_Species_df, aes(numb_peaks, Hj, color = Species)) + geom_point()
  print(p2)
  p3 <- ggplot(Hj_Species_df, aes(numb_peaks, Species, color = Hj)) + geom_point() + theme(legend.position = "none")
  print(p3)
  grid_plot <- grid.arrange (p1, p2, p3)

  # List
  list <- list(grid_plot, Hj_Species_df)
  print (list)
}

#' Metabolite (Si) and Metabolome (δj) Specialization Index
#'
#' @param Data
#' @description MetSpec calculates METabolite SPECialization (Si) and METabolome SPECialization (δj) indexes based on Shannon entropy.
#'
#' @details
#' Metabolic specificity (Si) is defined as the specificity of a particular MS/MS metabolite (i) among a set of samples (j). Metabolome specialization δj is measured for each jth sample, as the average of the MS/MS specificities.
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 4 objects:
#' 1. A data frame with the Si per compound;
#' 2. A data frame with the δj per Species;
#' 3. A point plot despicting the Si value of each compound
#' 4. A point plot despicting the δj per species
#'
#' @export
#'
#' @examples Dj <- MetSpec (Data)


MetSpec <- function (Data){
  # Mssg
  print("THIS PARAMETER CAN NOT BE EXTRAPOLATED TO OTHER DATA SETS OR SUBSETS")

  #Mssg
  print("To use MetPar function store this calculation as Dj <- MetSpec (Data)")


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
  Dj_Species_df <- setDT(Dj_df, keep.rownames = "Species")
  Dj_Species_df$Species <- sub("\\_?\\.?\\d+", "", Dj_Species_df$Species)
  Dj_Species_df$Species <- sub("\\..", ".", Dj_Species_df$Species)
  print(Dj_Species_df)

  # Graphical output
  p2 <- ggplot(Dj_Species_df, aes(Species, Dj)) + geom_boxplot() + geom_point()
  print (p2)
  grid_plot <- grid.arrange (p1, p2)

  #List
  list <- list(grid_plot, Si_df, Dj_Species_df)
  print(list)
}

#' Metabolite Specialization Analysis
#'
#' @param Data
#' @description MetliteSpec calculates the contribution of METaboLITE SPECialization factor (Pij.Si) to the Metabolome specialization index.
#'
#' @details
#' Metabolite specialization factor (Pij.Si) is defined as product of the a metabolite specilization index (Si) and the frequency of the metabolite in a given sample (Pij).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a list with 4 objects:
#' 1. A data frame with the Pij.Si values per species and compound;
#' 2. A data frame with the highest Pij.Si compound value per species;
#' 3. A point plot despicting the Pij.Si value of each species
#' 4. A point plot despicting the highest Pij.Si value per species
#'
#' @export
#'
#' @examples MetliteSpec <- MetliteSpec (Data)


MetliteSpec <- function (Data){
  # Mssg
  print("THIS PARAMETER CAN NOT BE EXTRAPOLATED TO OTHER DATA SETS OR SUBSETS")

  # Calculate metabolite specialization
  Pij_by_Pi <- Pij_by_Pi_fc (Data)
  Si <- apply(Pij_by_Pi, MARGIN = 1, FUN = Si_fc)

  # Si weight on Dj
  Pij <- Pij_fc (Data)
  Pij_Si <- sweep(Pij, MARGIN = 1, Si, FUN = "*")
  Pij_Si_t <- t(Pij_Si)
  Pij_Si_df <- as.data.frame (Pij_Si_t)
  Compounds <- colnames(Pij_Si_df)
  Pij_Si_df <- setDT(Pij_Si_df, keep.rownames = "Species")
  Tidy_table <- pivot_longer(Pij_Si_df, all_of(Compounds), names_to = "Compounds", values_to = "Pij_Si")
  Tidy_table$Species <- sub("\\_?\\.?\\d+", "", Tidy_table$Species)
  Tidy_table$Species <- sub("\\..", ".", Tidy_table$Species)
  print (Tidy_table)
  Table_filtered <- my_filter_function (Tidy_table)
  print (Table_filtered)

  # Graphical Output
  p <- ggplot(Tidy_table, aes(Species, Pij_Si, color = Compounds)) + geom_point() + theme(legend.position = "none")
  print(p)
  q <- ggplot (Table_filtered, aes(Species, Pij_Si, color = Compounds)) + geom_point() + geom_text(aes(label= Compounds), hjust=0,vjust=0, size=3) + theme(legend.position = "none")
  print(q)
  grid_plot <- grid.arrange (p,q)

  # List
  list <- list (grid_plot, Tidy_table, Table_filtered)
  print (list)
}

#' Metabolic Plasticity Parameters Summary
#'
#' @param Data
#' @description MetPar integrates Hj, and δj indexes calculations.
#'
#' @details
#' Metabolic Plasticity Parameters Summary extracts and summarizes Hj and δj information from MetDiv and MetSpec calculations per species
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame
#'
#' @export
#'
#' @examples
#' Hj <- MetDiv (Data)
#' Dj <- MetDiv (Data)
#' MetPar <- MetPar (Data)

MetPar <- function(Data){
  #Mssg
  print("MetDiv table must be store as Hj, and MetSpec as Dj")
  #Error
  if (exists("Hj") == FALSE) {
    stop ("Hj does not exist. Please, execute MetDiv() function and store it as Hj. Hj <- MetDiv(Data)")}
  if (exists("Dj") == FALSE) {
    stop ("Dj does not exist. Please, execute MetSpec() function and store it as Dj. Dj <- MetSpec(Data)")}
  # Extract the Hj and Dj parameters per species
  MetPar_df <- cbind(Hj[[2]], Dj = Dj[[3]]$Dj)
  # Evaluate the dependence between Hj and Dj
  p1 <- ggplot(MetPar_df, aes(Hj, Dj, color= Species)) + geom_point()

  list <- list (MetPar_df, p1)
  print(list)
}

#' Metabolic Plasticity Statistics
#'
#' @param Data
#' @description MetStats calculates sample (HRj) and Kullback–Leibler divergence (Divj)
#'
#' @details
#' Metabolic Plasticity Statistics calculates HRj and DivJ from MetDiv calculations per species
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame
#'
#' @export
#'
#' @examples
#' Hj <- MetDiv (Data)
#' MetPar <- MetPar (Data)

MetStats <- function(Data) {

  #Mssg
  print("MetDiv table must be store as Hj")
  #Error
  if (exists("Hj") == FALSE) {
    stop ("Hj does not exist. Please, execute MetDiv() function and store it as Hj. Hj <- MetDiv(Data)")}

  # Calculate divergence
  Pij <- Pij_fc (Data)
  HRj <- apply (Pij, MARGIN = 2, FUN = HRj_function)

  # Calculate Kullback–Leibler divergence
  Dj <- HRj - Hj[[2]]$Hj

  #Generating the output
  MetStats <- cbind(Hj[[2]], HRj = HRj, Divj = Dj)
  print(MetStats)
}

