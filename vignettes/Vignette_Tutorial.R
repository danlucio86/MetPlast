## ----setup, include=FALSE-------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -------------------------------------------------------------------------
# devtools::install_github ("danlucio86/MetPlast")

library("MetPlast")

## -------------------------------------------------------------------------
Data <- MetPlast::Data

## ---- echo=TRUE-----------------------------------------------------------
##Data <- read.csv2(file = "Test.csv", header = TRUE, row.names = "Compounds")

library(rmarkdown)
paged_table(head(Data))


## -------------------------------------------------------------------------
Data[is.na(Data)] <- (min(Data, na.rm=TRUE))/1000000

library(rmarkdown)
paged_table(head(Data))

## ---- results='hide'------------------------------------------------------
MetDiv_df <- MetDiv (Data)


## -------------------------------------------------------------------------
print(head((MetDiv_df)))

## -------------------------------------------------------------------------
MetDiv_plot <- MetDiv_plot(MetDiv_df =  MetDiv_df)

plot(MetDiv_plot)


## -------------------------------------------------------------------------
MetDiv_Supp_plot_1 <- MetDiv_Supp_plot(MetDiv_df = MetDiv_df, y_element = Hj)

plot(MetDiv_Supp_plot_1)

## -------------------------------------------------------------------------
MetDiv_Supp_plot_2 <- MetDiv_Supp_plot(MetDiv_df = MetDiv_df, y_element = Sample)

plot(MetDiv_Supp_plot_2)

## -------------------------------------------------------------------------
categories <- MetPlast::categories

## -------------------------------------------------------------------------
## categories <- read.csv2(file="categories.csv")

## -------------------------------------------------------------------------
Data_categ <- cbind (MetDiv_df, categories)

## -------------------------------------------------------------------------
ggplot(Data_categ, aes(Hj, Sample, color = Categories)) + geom_point() + theme(legend.position = "none") 

## ---- results='hide'------------------------------------------------------
Si_df <- Si_index(Data)

## -------------------------------------------------------------------------
print(head((Si_df)))

## -------------------------------------------------------------------------
Si_plot <- Si_index_plot(Si_df = Si_df)

plot(Si_plot)

## ---- results='hide'------------------------------------------------------
Dj_df <- Dj_index(Data = Data)

## -------------------------------------------------------------------------
print(head((Dj_df)))

## -------------------------------------------------------------------------
Dj_plot <- Dj_index_plot(Dj_df = Dj_df)

plot(Dj_plot)


## -------------------------------------------------------------------------
Dj_categ <- cbind (Dj_df, categories)

## -------------------------------------------------------------------------
ggplot(Dj_categ, aes(Dj, Sample, color = Categories)) + geom_point() + theme(legend.position = "none") 

## ---- results='hide'------------------------------------------------------
Dj_weight_df <- Dj_index_weight(Data = Data, Si_df = Si_df)

## -------------------------------------------------------------------------
print(head((Dj_index_weight)))

## -------------------------------------------------------------------------
Dj_weight_plot <- Dj_index_weight_plot (Dj_weight_df = Dj_weight_df)

plot(Dj_weight_plot)

## ---- results='hide'------------------------------------------------------
Stats <- MetStats(Data = Data, MetDiv_df= MetDiv_df)

## -------------------------------------------------------------------------
print(head((Stats)))

## ---- results='hide'------------------------------------------------------
MetPar <- MetPar_df(Stats_df = Stats, Dj_df = Dj_df)


## -------------------------------------------------------------------------
print(head((MetPar)))

## -------------------------------------------------------------------------
# Statistical Analysis 

## ANOVA analysis

library (ggpubr)

compare_means(Hj ~ Sample, data = MetPar)



## -------------------------------------------------------------------------
compare_means(Dj ~ Sample, data = MetPar)

## -------------------------------------------------------------------------
Dj_vs_Hj_plot <- MetPar_plot(MetPar_df = MetPar, x_element = Hj, ref_line = "FALSE")

## -------------------------------------------------------------------------
Divj_vs_Hj_plot <- MetPar_plot(MetPar_df = MetPar, x_element = Divj, ref_line = TRUE)

## ---- results='hide', fig.show='hide'-------------------------------------
MetDiv <- MetDiv_pipe(Data = Data)

## ---- results='hide', fig.show='hide'-------------------------------------
MetSpec <- MetSpec_pipe(Data = Data)

## ---- results='hide', fig.show='hide'-------------------------------------
MetliteSpec_pipe <- MetliteSpec_pipe(Data = Data)

## ---- results='hide', fig.show='hide'-------------------------------------
MetPar_pipe <- MetPar_pipe(Stats_df = Stats, Dj_df = Dj_df)

