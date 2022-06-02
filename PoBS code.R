#script for PoBS GitHub

#import data
setwd("~/Desktop/School/Lab/Stats Projects/PoBS paper/GitHub") #change to your path containing the data files
library(ggplot2)
RGT_data <- readRDS("RGT_data.rds") #full behavioral set from 5 RGT experiments
RGT_biome <- readRDS("RGT_biome.rds") #subset of RGT data with biological variables (gut microbiome)


### Correlational Analyses ###

#prep variables
RGT_biome$Injury<-as.factor(RGT_biome$Injury)
RGT_biome$Collect_Time<-as.factor(RGT_biome$Collect_Time)
RGT_biome$Week <- as.numeric(RGT_biome$Week) 


#simple Pearson correlations

#alpha diversity x optimal choice
cor.test(RGT_biome$alpha_diversity, RGT_biome$PctOptimal)
ggplot(data=RGT_biome, aes(x=alpha_diversity, y=PctOptimal))+
  geom_smooth()+
  theme_classic(base_size = 20)
ggplot(data=RGT_biome, aes(x=alpha_diversity, y=PctOptimal))+
  geom_point()+ #add individual-subject data
  geom_smooth()+
  theme_classic(base_size = 20)
ggsave("correlation.png", width = 25, height = 25, units = "cm")


#alpha diversity x premature responses
cor.test(RGT_biome$alpha_diversity, RGT_biome$PctPremature)
ggplot(data=RGT_biome, aes(x=alpha_diversity, y=PctPremature))+
  geom_smooth()+
  theme_classic()


#more complex regression models
library(lme4)
m1<-lmer(scale(asin(sqrt(PctPremature/100)))~ 
           Injury*scale(Week)+
           scale(asin(sqrt(PrematureBase/100)))+ 
           (1|Subject), 
           data=subset(RGT_biome, Collect_Time=="D3"))
m2<-lmer(scale(asin(sqrt(PctPremature/100))) ~ 
           Injury*scale(alpha_diversity)*scale(Week)+  
           scale(asin(sqrt(PrematureBase/100)))+
           (1|Subject), 
          data=subset(RGT_biome, Collect_Time=="D3"))
anova(m1, m2)





### Mixed-Effects Modeling ###


#subset to stable post-injury data only
data=subset(RGT_data, Session>10 & Session<21)
data$ChoiceOption<-as.factor(data$ChoiceOption)
data$ChoiceOption<-relevel(data$ChoiceOption, ref="2")
data$Injury<-as.factor(data$Injury)

#compare fixed effect versus intercept-only versus random slope+intercept models
slope=lmer(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session)+
            (1+ChoiceOption|Subject), 
            control = lmerControl(calc.derivs = FALSE),
            data=data)
int=lmer(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session)+
           (0+dummy(ChoiceOption, "2")|Subject), data=data)
fixed=lm(scale(asin(sqrt(PctChoice/100)))~ChoiceOption*Injury*scale(Session), data=data)
anova(slope, int, fixed)
library(lmerTest)
summary(slope)
ranef(slope)

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






### k-means clustering ###

#format data
temp<-aggregate(PctChoice~Subject+ChoiceOption, mean, data=data)
wide<-tidyr::spread(temp, ChoiceOption, PctChoice, fill = NA, convert = FALSE)

#plot cluster number
library(factoextra)
fviz_nbclust(wide[,-c(1)], kmeans, method = "wss") +
  ggtitle("Optimal Cluster Number")+
  geom_point(aes(group = 1),size=10)+
  geom_line(aes(group = 1), size = 3) + 
  ylab("Within Sum of Squares")
ggsave("cluster_number.png", width = 25, height = 25, units = "cm")

#create and plot clusters
Clusters<-kmeans(wide[,-c(1)], 5, iter.max=200, nstart=30)
fviz_cluster(Clusters, data=wide[,-c(1)])+
  ggtitle("PCA on Clusters")+
  scale_colour_manual(values = c("#0987D7", "#F46764", "#0DBD51",  "#9A1F09", "#DDB91C"))+
  scale_fill_manual(values = c("#0987D7", "#F46764", "#0DBD51",  "#9A1F09", "#DDB91C"))+
  scale_shape_manual(values=c(19, 15, 17, 18, 3))+
  xlab("Axis 1")+
  ylab("Axis 2")+
  theme_bw(base_size=20)
ggsave("cluster_PCA.png", width = 25, height = 25, units = "cm")


#merge clusters back in with original data for plotting
clust_data<-cbind(temp, cluster=Clusters$cluster)
clust_data$ChoiceOption<-relevel(clust_data$ChoiceOption, ref="1")
ggplot(data=clust_data, aes(x=ChoiceOption, y=PctChoice))+   
  geom_point(size=2, alpha=0.6)+
  stat_summary(aes(group=cluster), size=2, fun=mean, geom="line")+
  facet_wrap(~cluster)+
  ylab("Percent Choice")+
  theme_bw(base_size=20)
ggsave("cluster_profiles.png", width = 25, height = 25, units = "cm")




