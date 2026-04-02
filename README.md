# Assessing Sustainable Yield of Lake Trout in Lake Louise and Susitna Lake Using Two Different Models and Documenting Angler Satisfaction of the Fishery, 2026

This study will use 2 different models to estimate sustainable yield of lake trout _Salvelinus namaycush_ from Lake Louise and Susitna Lake.  The lake area (LA) model was developed from data compiled from 43 lakes located in Ontario, Canada and predicts yield potential (YP) from the surface area of a give lake.  The Lester model is a newer, more comprehensive model that considers information from 456 lakes throughout Canada including the Yukon and Northwest Territories.  Unlike the LA model, the Lester model estimates maximum sustained yield (MSY), which is much different from YP, and of more use to fishery managers.  The output of both models is in kg biomass, so the mean weight of lake trout subject to removal from Lake Louise and Susitna Lake needs to be estimated to convert express biomass into in terms of total number of fish.  This study will use lake trout catch data and catch sampling from late winter (April), summer (June), and early winter (November–December) to estimate mean weight.  The data collected will be used to apply the Lester model to information from an Alaskan lake, which may provide additional insights to its suitability for management of lake trout in Alaska and suggest further refinements to its adaptation.  Concurrent with the capture and sampling of lake trout, a survey of anglers will also occur.  The survey will document angler demographics and, seasonal fishery use, and assess angler satisfaction with the current lake trout fishery at Lake Louise and Susitna Lake.  Estimates of yield will be used to evaluate the sustainability of the current sport fishery and along with angler survey results be used to address potential regulatory proposals to liberalize or restrict this fishery.

Key Words: lake trout, Salvelinus namaycush, Lake Louise, lake area model, Lester model, yield potential, maximum sustained yield, mean weight, ice fishing, open water fishing.

## Contents

### OP_2026/R

This folder contains all R code relevant to the Operational Plan published in 2026.

* **1_RP_of_MnWt.R** contains a comprehensive simulation of estimation of mean weight
under several different sampling weights and differences in (fish) weight distribution
between samples, as well as different levels of uncertainty in the true harvest
fractions in each sampling period, and finally reporting the reasonable worst-case
scenario used for estimating relative precision in the operational plan.

* **2_RP_of_LAM** contains a simulation estimating the uncertainty due to estimation
within the Lake Area Model, and the relative precision due to applying the scenario
in **1_RP_of_MnWt.R** to apply estimated mean weight to Yield Potential.
