## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ---- echo=FALSE--------------------------------------------------------------
library(RxODE)
## rxSyntaxFunctions <- read.csv("syntax-functions.csv");usethis::use_data(rxSyntaxFunctions)
library(DT)
datatable(rxSyntaxFunctions, rownames = FALSE, filter="top", options = list(pageLength = 10, scrollX=T) )

## ---- echo=FALSE--------------------------------------------------------------
## rxReservedKeywords <- read.csv("reserved-keywords.csv");names(rxReservedKeywords)[1] <- "Reserved Name";usethis::use_data(rxReservedKeywords)
datatable(rxReservedKeywords, rownames = FALSE, filter="top", options = list(pageLength = 10, scrollX=T) )

