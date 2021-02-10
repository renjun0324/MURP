# MURP

A model-based downsampling algorithm, the “minimal unbiased representative points (MURP)”, which offers an optimal, sufficient and unbiased representation of cell population.

## Installation

```r
install.packages("devtools")
devtools::install_github("renjun0324/MURP")
```

## Usage

```r
data(sdata)
result = MURP(Data = sdata, cores = 1, iter = 20, omega = 1/6, seed = 723)
ggplot2::ggsave("bic_plot.pdf", KBicPlot(murpResult = result))
```
