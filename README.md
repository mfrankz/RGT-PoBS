# Analysis of Rodent Gambling Task (RGT) Data
### This code accompanies a manuscript by Frankot, Young, and Vonder Haar in Perspectives on Behavioral Science. Here, we have compiled behavioral data from 5 preclinical experiments and provide reproducible examples of various techniques used to analyze the data

We will begin by importing two datasets. The first, RGT_data, contains behavioral data collapsed across five experiments. The second, RGT_biome, contains a subset of behavioral data accompanied by biological variables from measurement of the gut microbiome.

```
RGT_data <- read_csv("RGT.Counts.csv") #full behavioral set from 5 experiments
RGT_biome <- read_csv("~/Desktop/School/Lab/WetLab/Microbiome/HFD/Final/Behavior.Analysis.csv") #subset of data with biological variables
```
