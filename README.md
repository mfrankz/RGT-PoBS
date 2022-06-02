# Analysis of Rodent Gambling Task (RGT) Data
### This code accompanies a manuscript by Frankot, Young, and Vonder Haar in Perspectives on Behavioral Science. Here, we have compiled behavioral data from 5 preclinical experiments and provide reproducible examples of various techniques used to analyze the data. Our behavioral outcomes come from the RGT, which concurrently measures optimal, suboptimal, and risky decisions. For more information on the RGT or these datasets, see [our published works](https://www.frontiersin.org/articles/10.3389/fnbeh.2022.837654/full).


### We will begin by importing two datasets. 

```
RGT_data <- read_csv("RGT_data.csv") #full behavioral set from 5 RGT experiments
RGT_biome <- read_csv("RGT_biome.csv") #subset of RGT data with biological variables (gut microbiome)

```
### The first, RGT_data, contains behavioral data collapsed across five experiments. The second, RGT_biome, contains a subset of behavioral data accompanied by biological variables from measurement of the gut microbiome. Our tutorial for these data consists of the following sections
# 1. Explanation of variables in dataset
# 2. Correlational analyses
# 3. Mixed-effects modeling
# 4. K-means clustering 


# 1. Explanation of variables in dataset

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



# 2. Correlational analyses

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
```
```
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
<img src="https://github.com/mfrankz/RGT-PoBS/blob/main/correlation.png" width="600">

These plots support the conclusion that the correlation across these variables is not particularly robust

Sometimes, more complex regression models are needed to perform correlational analyses. For example, if we wanted to determine whether alpha diversity explains variance in behavior beyond the variance explained by TBI, we would need to perform a linear mixed model. The code for this analysis is provided here, and more information on mixed-effects modeling is provided in the section below. 

```
library(lme4)
#model predicting premature responses using injury as a predictor
m1<-lmer(scale(asin(sqrt(PctPremature/100)))~ 
           Injury*scale(Week)+
           scale(asin(sqrt(PrematureBase/100)))+ 
           (1|Subject), 
           data=subset(RGT_biome, Collect_Time=="D3"))
#model predicting premature responses using injury and alpha diversity as predictors	   
m2<-lmer(scale(asin(sqrt(PctPremature/100))) ~ 
           Injury*scale(alpha_diversity)*scale(Week)+  
           scale(asin(sqrt(PrematureBase/100)))+
           (1|Subject), 
          data=subset(RGT_biome, Collect_Time=="D3"))
anova(m1, m2) #compare models 
```
The output (below) shows that the models are similar, but adding alpha diversity does improve model fit, p = 0.045
```
Data: subset(RGT_biome, Collect_Time == "D3")
Models:
m1: scale(asin(sqrt(PctPremature/100))) ~ Injury * scale(Week) + scale(asin(sqrt(PrematureBase/100))) + (1 | Subject)
m2: scale(asin(sqrt(PctPremature/100))) ~ Injury * scale(alpha_diversity) * scale(Week) + scale(asin(sqrt(PrematureBase/100))) + (1 | Subject)
   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
m1    7 535.87 559.44 -260.94   521.87                       
m2   11 534.16 571.19 -256.08   512.16 9.7128  4    0.04555 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
We will now switch to the RGT_data to perform more complex statistical analyses.


# 2. Mixed-effects modeling

To illustrate why mixed-effects (also called multi-level, random effects) modeling is important for nested data, we will compare a fixed-effects model against 2 mixed models.

```
#subset to stable post-injury data only and prep variables
data=subset(RGT_data, Session>10 & Session<21)
data$ChoiceOption<-as.factor(data$ChoiceOption)
data$ChoiceOption<-relevel(data$ChoiceOption, ref="2")
data$Injury<-as.factor(data$Injury)
```

```
#compare fixed effect versus intercept-only versus random slope+intercept models
library(lme4)
slope=lmer(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session)+
            (1+ChoiceOption|Subject), 
            control = lmerControl(calc.derivs = FALSE),
            data=data)
int=lmer(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session)+
           (0+dummy(ChoiceOption, "2")|Subject), data=data)
fixed=lm(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session), data=data)
anova(slope, int, fixed)
```
The output here shows that the mixed models outperform the fixed effect only model. Among the mixed models, the model that allows both the slope and intercept to vary explains much more variance in behavior
```

      npar     AIC     BIC  logLik deviance  Chisq Df Pr(>Chisq)    
fixed   17 18743.7 18863.6 -9354.8  18709.7                         
int     18 16430.9 16557.9 -8197.4  16394.9 2314.8  1  < 2.2e-16 ***
slope   27  7941.4  8131.9 -3943.7   7887.4 8507.5  9  < 2.2e-16 ***
---
```
Based on these findings, we will look at the regression coefficients for the random slope + intercept model
```
library(lmerTest) #library is necessary for p-value calculation
summary(slope)
ranef(slope)
```
```
#plot random effects
random<-ranef(slope)$Subject
random$Subject<-rownames(random)
random<-tidyr::gather(random, "Variable", "Value", 1:4, factor_key = T)
ggplot(data=random, aes(x=Variable, y=Value))+
  geom_violin()+
  geom_point(size=1, alpha=0.5)+
  ylab("Standardized Value")+
  theme_classic(base_size=20)
ggsave("random_effects.png", width = 25, height = 25, units = "cm")
```
The random effects plot is a quick way of visually the vastly different choice profiles across individual subjects. We will now explore clustering techniques as a means to account for individual subject variability
