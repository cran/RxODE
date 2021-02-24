## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## -----------------------------------------------------------------------------
syntax <- read.csv("syntax-functions.csv");
library(DT)
datatable(syntax, rownames = FALSE, filter="top", options = list(pageLength = 10, scrollX=T) )

## -----------------------------------------------------------------------------
syntax <- read.csv("reserved-keywords.csv");
library(DT)
datatable(syntax, rownames = FALSE, filter="top", options = list(pageLength = 10, scrollX=T) )

