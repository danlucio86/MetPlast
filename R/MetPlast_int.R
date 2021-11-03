#' Metabolite frequency per sample - Internal function
#'
#' @param Data a data frame the levels of different detected  metabolites/compounds -rows-, in different samples -columns-.
#'
#' @description Internal function to calculate the frequency of metabolite i in sample j (Pij)
#'
#' @details
#' Calculates the frequency dividing the intensity of a given metabolite (i) by the sum of intensities of all metabolites detected in a given sample (i.e. species) (j).
#'
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#'
#' @return It returns a dataframe with metabolite frequencies

Pij_fc <- function (Data) {Total_samples <- apply (Data, MARGIN=2, FUN = sum)
Pij <- sweep(Data, MARGIN = 2, Total_samples, FUN = '/')}

#' Metabolite expression identity internal function
#'
#' @param Data a data frame the levels of differente detected  metabolites/compounds -rows-, in different samples -columns-.
#' @description The expression identity of a given metabolite profile regarding frequencies among the set of samples (j).
#'
#' @details
#' The metabolite expression identity is calculated dividing the metabolite frequency per sample by the frequency of the given metabolite in the whole dataset
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a data frame with Pij/Pi values

Pij_by_Pi_fc <- function (Data) {Total_Compound <- apply (Data, MARGIN=1, FUN = sum)
Pij <- Pij_fc (Data)
n_samples <- ncol(Data)
Pi <- (1/n_samples)*(apply (Pij, MARGIN=1, FUN=sum))
Pij_divided_Pi <- sweep(Pij, MARGIN = 1, Pi, FUN = "/")
}

#' Metabolic diversity internal function
#'
#' @param x is a numerical value
#' @description It performs internal metabolic diversity (Hj) calculations.
#'
#' @details
#' Performs an internal calculation necessary to obtained the metabolic diversity per sample (Hj).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a value with the Hj value per sample

Hj_fc <- function(x){-sum(x * log2(x), na.rm=TRUE)}

#' Metabolic diversity divergence internal function
#'
#' @param x is a numerical value
#' @description It performs internal metabolic diversity divergence (HRj) calculations.
#'
#' @details
#' Performs an internal calculation necessary to obtained the metabolic diversity divergence (HRj).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a value with the HRj value per sample

HRj_fc <- function(x){
  Pij <- Pij_fc (Data)
  n_samples <- ncol(Data)
  Pi <- (1/n_samples)*(apply (Pij, MARGIN=1, FUN=sum))
  -sum(x*log2(Pi))
}

#' Metabolite Specialization internal function
#'
#' @param x is a numerical value
#' @description It performs internal metabolite specialization (Si) calculations.
#'
#' @details
#' Performs an internal calculation necessary to obtained the metabolite specialization (Si).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a value with the Si value per metabolite (i)

Si_fc <- function(x){n_samples <- ncol(Data)
(1/n_samples)*(sum(x*log2(x), na.rm =TRUE))}

#' Number of detected compounds function
#'
#' @param Data a data frame the levels of different detected  metabolites/compounds -rows-, in different samples -columns-.
#'
#' @description It counts the number of detected compounds per sample
#'
#' @details
#' It counts the number of compound per sample that have a value higher to the mininum.
#' Initially, NA values are replaced by the minimum value of the whole data set divided by 1000000.
#' Thus, the function considers the minimum value of the corrected data set as a NA in the original data set (i.e. absence of the compound)
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a vector with the number of peaks - e.g. detected compounds - per species.

numb_peaks_fc <- function(Data) {min_data <- min(Data)
Data_numb_peaks <- apply(Data, MARGIN = 2, FUN = function(x){x != min_data})
numb_peaks <- apply (Data_numb_peaks, MARGIN = 2, FUN = sum)}

#' Metabolite contribution to metabolome specialization index (δj) function
#'
#' @param Dj_index_weight a data frame with all the Pij*Si values per compound for each species.
#' @description It selects the metabolites that contribute the most to the metabolome specialization index (δj) - i.e. highest Pij*Si
#'
#' @details
#' Metabolite contribution to metabolome specialization index (δj) function select the metabolite with the highest Pij*Si value per samples (i.e. species).
#'
#' @author Lucio D'Andrea, PhD; Prof Aureliano Bombarely
#'
#' @return It returns a table


my_filter_fc <- function(table) {
   table %>%
    dplyr::group_by(Sample) %>%
    dplyr::filter(Pij_Si == max(Pij_Si))
}
