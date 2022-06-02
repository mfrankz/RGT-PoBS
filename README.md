# Analysis of Rodent Gambling Task (RGT) Data
### This code accompanies a manuscript by Frankot, Young, and Vonder Haar in Perspectives on Behavioral Science. Here, we have compiled behavioral data from 5 preclinical experiments and provide reproducible examples of various techniques used to analyze the data. Our behavioral outcomes come from the RGT, which concurrently measures optimal, suboptimal, and risky decisions. For more information on the RGT or these datasets, see [our published works]([https://github.com/mfrankz/microbiome/blob/main/ps.rds](https://www.frontiersin.org/articles/10.3389/fnbeh.2022.837654/full)).


# We will begin by importing two datasets. The first, RGT_data, contains behavioral data collapsed across five experiments. The second, RGT_biome, contains a subset of behavioral data accompanied by biological variables from measurement of the gut microbiome. The first section contains an explanation of variables in the dataset. Then, we will perform the following analyses: 
### 1. Correlational analyses
### 2. Mixed-effects modeling
### 3. K-means clustering 

```
RGT_data <- read_csv("RGT_data.csv") #full behavioral set from 5 RGT experiments
RGT_biome <- read_csv("RGT_biome.csv") #subset of RGT data with biological variables (gut microbiome)

```

# Explanation of variables in dataset

### Variables in RGT_data
This dataset contains the count and % level of choice of the 4 options on the RGT
1. Study_ID: study codes for 5 different experiments compiled for this project
2. Subject: unique subject identifiers
3. Injury: TBI (traumatic brain injury) versus Sham (intact/control)
4. Session: Session number. All sessions provided are post injury
5. ChoiceOption: Option on the RGT. Choice 1 is suboptimal, 2 is optimal, and 3 and 4 are risky
6. ChoiceCount: number of times any given option was selected within each subject/session
7. TotChoice: total number of choices made within each subject/session
8. PctChoice: percent choice of any given option (i.e., (ChoiceCount/TotChoice)x100)

### Variables in RGT_biome
This dataset contains behavioral RGT data and a measurement of the gut microbiome
1. Subject: unique subject identifier
2. Injury: TBI (traumatic brain injury) versus Sham (intact/control)
3. Biome_SampleID: unique identifier for microbiome sample. Samples were collected at 4 timepoints
4. Collect_Time: Collection timepoint consisting of pre-injury (Pre) and 3, 30, and 60 days post-injury
5. alpha_diversity: Microbiome measurement of heterogeneity/richness within a sample
6. Week: week of post-injury testing on the RGT. Behavioral outcomes are aggregated across weeks
7. PctOptimal: percent optimal choice (ChoiceOption #2) on the RGT
8. PctPremature: percent premature responses on the RGT. Serves as a measure of motor impulsivity
9. OptimalBase: Baseline (pre-injury) percent optimal choice
10. PrematureBase: Baseline (pre-injury) percent premature



# Correlational analyses

Correlational analyses are useful for analyzing continuous biological variables. Here we will correlate our microbiome measurement (alpha diversity) with behavioral outcomes

```
#prep variables
RGT_biome$Injury<-as.factor(RGT_biome$Injury)
RGT_biome$Collect_Time<-as.factor(RGT_biome$Collect_Time)
RGT_biome$Week <- as.numeric(RGT_biome$Week) 
```
We will first run simple Pearson bivariate correlations between alpha diversity and optimal choice
```
#alpha diversity x optimal choice
cor.test(RGT_biome$alpha_diversity, RGT_biome$PctOptimal)

	Pearson's product-moment correlation

data:  RGT_biome$alpha_diversity and RGT_biome$PctOptimal
t = 4.2171, df = 854, p-value = 2.738e-05
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.07655315 0.20784596
sample estimates:
      cor 
0.1428278 
```
The output shows that the correlation is significant but not particularly strong, r=0.14, p < 0.001
We can get a better idea of the relationship from a plot of the data

```
ggplot(data=RGT_biome, aes(x=alpha_diversity, y=PctOptimal))+
  geom_smooth()+
  theme_classic()
ggplot(data=RGT_biome, aes(x=alpha_diversity, y=PctOptimal))+
  geom_point()+ #add individual-subject data
  geom_smooth()+
  theme_classic()
```
These plots support the conclusion that the correlation across these variables is not particularly robust

