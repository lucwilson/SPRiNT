######################################################################################################################################################################
# 
#                             LEMON behav analysis script 
#                 Run regression models and Bayes Factor analysis for
#           Time-resolved parameterization of aperiodic and periodic of brain activity paper
#
#
#
######################################################################################################################################################################
library(dplyr)
library(ggplot2)
library(BayesFactor)
library(lmSupport)
library(sjPlot)
library(Hmisc)
library(effectsize)
library(ggpubr)
colour_scheme= c("#8a9ccb", "#facf93")

###################################################################################
#                                 Import and clean data with R
###################################################################################
EEG_data=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_rmbad.csv', stringsAsFactors = FALSE)
MW= read.csv('~/Desktop/behav_LEMON/MW_Oz_rmbad.csv', stringsAsFactors = FALSE)
demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)

###################################################################################
#                                 Merge and clean data
###################################################################################

# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

MW=MW[MW$subject %in% subj,]

colnames(MW)[2:7]= c("MW_ec_alpha_pow","MW_eo_alpha_pow", "MW_ec_alpha_maxpowfreq_mean", "MW_ec_alpha_maxpowfreq_std",  "MW_eo_alpha_maxpowfreq_mean", "MW_eo_alpha_maxpowfreq_std" )

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data, MW[,-1], tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced

###################################################################################
#                                   Explore/ plot data
###################################################################################

hist(data$eo_ap_slope_mean)
hist(data$eo_ap_slope_std)

hist(data$ec_ap_slope_mean)
hist(data$ec_ap_slope_std)

t.test(data$ec_mean_MAE[data$Age2==1], data$ec_mean_MAE[data$Age2==0])
t.test(data$eo_mean_MAE[data$Age2==1], data$eo_mean_MAE[data$Age2==0])
t.test(data$eo_mean_MAE, data$ec_mean_MAE)

effectsize::cohens_d(data$eo_mean_MAE, data$ec_mean_MAE)

cor.test(data$ec_ap_slope_mean,data$ec_mean_MAE)
cor.test(data$eo_ap_slope_mean,data$eo_mean_MAE)
###################################################################################
#                          Run regression analysis to predict Age
###################################################################################


# run regression eyes-open
lm6=glm(Age2 ~ eo_alpha_cf_mean + eo_alpha_cf_std + eo_alpha_pow_mean + eo_alpha_pow_std+ eo_ap_slope_mean + eo_ap_slope_std , data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table

# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,eo_alpha_cf_mean, eo_alpha_cf_std, eo_alpha_pow_mean, eo_alpha_pow_std, eo_ap_slope_mean, eo_ap_slope_std)
data4BF= data4BF[!is.na(data4BF$eo_alpha_cf_mean),]
colnames(data4BF)[2:7]= c("Eyes-open mean alpha center frequency", "Eyes-open std alpha center frequency", "Eyes-open mean alpha power", "Eyes-open std alpha power","Eyes-open mean aperiodic slope", "Eyes-open std aperiodic slope")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


## plot age effect 
error=confint(lm6)
betas=lm6$coefficients
D_betas=as.data.frame(cbind(error[1:7,],betas))
D_betas$names=c("1.Intercept", "2.Eyes-open mean alpha center frequency", "3.Eyes-open std alpha center frequency", "4.Eyes-open mean alpha power", "5.Eyes-open std alpha power","6.Eyes-open mean aperiodic slope", "7.Eyes-open std aperiodic slope")
D_betas$names=as.factor(D_betas$names)
colnames(D_betas)[1:2]=c("upper", "lower")
D_betas$inter=as.factor(c(1,0,0,0,0,0,0))

ggplot(D_betas, aes(x=betas, y=names)) + geom_rect(data = D_betas, aes(xmin = -Inf, xmax = Inf,ymax=1.5, ymin=0.5, fill = inter, alpha = 0.95)) + scale_fill_manual(values = c("white", "gray53"))+ geom_point(shape=23,size=10, fill= "black")+ geom_line(aes(size=10)) + geom_errorbarh(aes( xmin=lower, xmax=upper), D_betas, height=0, size=2.3) + theme_minimal() +
  geom_vline(xintercept=c(-0,0), linetype="dotted") + scale_y_discrete(position = "left", labels=c("1.Intercept" = "Intercept","2.Eyes-open mean alpha center frequency"= "Eyes-open mean alpha center frequency","3.Eyes-open std alpha center frequency"= "Eyes-open std alpha center frequency","4.Eyes-open mean alpha power"= "Eyes-open mean alpha power", "5.Eyes-open std alpha power" ="Eyes-open std alpha power", "6.Eyes-open mean aperiodic slope" = "Eyes-open mean aperiodic slope", "7.Eyes-open std aperiodic slope" = "Eyes-open std aperiodic slope")) + scale_x_continuous(position = "top")+ labs(y= "", x ="Fixed Effects Coefficients (β)") + theme(legend.position = "none", axis.line.x.top = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(text = element_text(size=20))

colour_scheme= c("#f2b968", "#ffdcab")
data4plot = cbind(data[c(1,52)], stack(data %>% select(eo_alpha_cf_std,  eo_ap_slope_mean)))
data4plot$Age2=as.factor(data4plot$Age2)        
pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values)) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "bar", size= 10, position = pd) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size=2.5, width=0.5, position = pd) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free")

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= Age2)) + geom_jitter(shape=16,size =5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free") + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 



colour_scheme= c("#f2b968", "#ffdcab")
data4plot = cbind(data[c(1,52)], stack(data %>% select(eo_alpha_cf_std)))
data4plot$Age2=as.factor(data4plot$Age2)        

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= Age2)) + geom_jitter(shape=16,size =5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=3, scales = "free") + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=rev(colour_scheme)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

#ggsave('/Users/jasondsc/Desktop/behav_LEMON/figure_new_June2022_eyesopen_effect_std_long2.pdf', device = "pdf", width=8, height=10, dpi=800)

# run same analysis for eyes-closed data
lm6=glm(Age2 ~ ec_alpha_cf_mean + ec_alpha_cf_std + ec_alpha_pow_mean + ec_alpha_pow_std+ ec_ap_slope_mean + ec_ap_slope_std , data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)
tab_model(lm6)


# clean 4 bayes analysis
data4BF= data %>% select(Age2,ec_alpha_cf_mean, ec_alpha_cf_std, ec_alpha_pow_mean, ec_alpha_pow_std, ec_ap_slope_mean, ec_ap_slope_std)
data4BF= data4BF[!is.na(data4BF$ec_alpha_cf_mean),]
colnames(data4BF)[2:7]= c("Eyes-closed mean alpha center frequency", "Eyes-closed std alpha center frequency", "Eyes-closed mean alpha power", "Eyes-closed std alpha power","Eyes-closed mean aperiodic slope", "Eyes-closed std aperiodic slope")

bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


## plot age effect 
error=confint(lm6)
betas=lm6$coefficients
D_betas=as.data.frame(cbind(error[1:7,],betas))
D_betas$names=c("1.Intercept", "2.Eyes-closed mean alpha center frequency", "3.Eyes-closed std alpha center frequency", "4.Eyes-closed mean alpha power", "5.Eyes-closed std alpha power","6.Eyes-closed mean aperiodic slope", "7.Eyes-closed std aperiodic slope")
D_betas$names=as.factor(D_betas$names)
colnames(D_betas)[1:2]=c("upper", "lower")
D_betas$inter=as.factor(c(1,0,0,0,0,0,0))

ggplot(D_betas, aes(x=betas, y=names)) + geom_rect(data = D_betas, aes(xmin = -Inf, xmax = Inf,ymax=1.5, ymin=0.5, fill = inter, alpha = 0.95)) + scale_fill_manual(values = c("white", "gray53"))+ geom_point(shape=23,size=10, fill= "black")+ geom_line(aes(size=10)) + geom_errorbarh(aes( xmin=lower, xmax=upper), D_betas, height=0, size=2.3) + theme_minimal() +
  geom_vline(xintercept=c(-0,0), linetype="dotted") + scale_y_discrete(position = "left", labels=c("1.Intercept" = "Intercept","2.Eyes-closed mean alpha center frequency"= "Eyes-closed mean alpha center frequency","3.Eyes-closed std alpha center frequency"= "Eyes-closed std alpha center frequency","4.Eyes-closed mean alpha power"= "Eyes-closed mean alpha power", "5.Eyes-closed std alpha power" ="Eyes-closed std alpha power", "6.Eyes-closed mean aperiodic slope" = "Eyes-closed mean aperiodic slope", "7.Eyes-closed std aperiodic slope" = "Eyes-closed std aperiodic slope")) + scale_x_continuous(position = "top")+ labs(y= "", x ="Fixed Effects Coefficients (β)") + theme(legend.position = "none", axis.line.x.top = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(text = element_text(size=20))

colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,52)], stack(data %>% select( ec_ap_slope_mean, ec_alpha_cf_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   
pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values)) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "bar", size= 10, position = pd) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size=2.5, width=0.5, position = pd) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free")

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= Age2)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free") + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")



# Let us repeat with STFT data max power peak
# run regression eyes-open
lm6=glm(Age2 ~ eo_alpha_maxpowfreq_mean + eo_alpha_maxpowfreq_std, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table

# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,eo_alpha_maxpowfreq_mean, eo_alpha_maxpowfreq_std)
data4BF= data4BF[!is.na(data4BF$eo_alpha_maxpowfreq_mean),]
colnames(data4BF)[2:7]= c("Eyes-open mean alpha max power", "Eyes-open std alpha max power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# run regression eyes-closed
lm6=glm(Age2 ~ ec_alpha_maxpowfreq_mean + ec_alpha_maxpowfreq_std, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table

# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,ec_alpha_maxpowfreq_mean, ec_alpha_maxpowfreq_std)
data4BF= data4BF[!is.na(data4BF$ec_alpha_maxpowfreq_mean),]
colnames(data4BF)[2:7]= c("Eyes-closed mean alpha max power", "Eyes-closed std alpha max power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)



# let us repeat the above analysis with morlet wavelet 
# run regression eyes-open
lm6=glm(Age2 ~  MW_ec_alpha_pow , data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)



# clean 4 bayes analysis
data4BF= data %>% select(Age2,MW_ec_alpha_pow, )
data4BF= data4BF[!is.na(data4BF$MW_ec_alpha_pow),]
data4BF= na.omit(data4BF)

bf = regressionBF(Age2 ~ MW_ec_alpha_pow  , data4BF, whichModels = "top")
plot(bf)
summary(bf)

# let us repeat the above analysis with morlet wavelet 
# run regression eyes-open
lm6=glm(Age2 ~ MW_eo_alpha_pow, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)


# clean 4 bayes analysis
data4BF= data %>% select(Age2, MW_eo_alpha_pow)
data4BF= data4BF[!is.na(data4BF$MW_eo_alpha_pow),]
data4BF= na.omit(data4BF)

bf = regressionBF(Age2 ~  MW_eo_alpha_pow , data4BF, whichModels = "top")
plot(bf)
summary(bf)


###################################################################################
#                          Run regression analysis to predict Condition
###################################################################################

# stack the data into one matrix for regression
data2= data[!is.na(data$eo_alpha_cf_mean),]
data2 %>% group_by(Age2) %>% summarise(m=mean(ec_ap_slope_mean))
data2 %>% group_by(Age2) %>% summarise(m=mean(eo_ap_slope_mean))
data2 %>% group_by(Age2) %>% summarise(m=mean(ec_alpha_cf_mean))
data2 %>% group_by(Age2) %>% summarise(m=mean(eo_alpha_cf_std))

# prep data to look at condition effect (eyes-open vs eyes-closed)
data4reg = cbind(data[c(1)], stack(data %>% select( ec_ap_slope_mean, eo_ap_slope_mean)))
data4reg = cbind(data4reg, stack(data %>% select( ec_ap_slope_std, eo_ap_slope_std)))
data4reg = cbind(data4reg, stack(data %>% select( ec_alpha_cf_mean, eo_alpha_cf_mean)))
data4reg = cbind(data4reg, stack(data %>% select( ec_alpha_cf_std, eo_alpha_cf_std)))
data4reg = cbind(data4reg, stack(data %>% select( ec_alpha_pow_mean, eo_alpha_pow_mean)))
data4reg = cbind(data4reg, stack(data %>% select( ec_alpha_pow_std, eo_alpha_pow_std)))
data4reg = cbind(data4reg, stack(data %>% select( MW_ec_alpha_pow, MW_eo_alpha_pow)))
data4reg = cbind(data4reg, stack(data %>% select( ec_mean_MAE, eo_mean_MAE)))
data4reg=data4reg[c(1,2, 4, 6, 8, 10, 12, 14, 16)]
data4reg$cond= c(rep(0,178),rep(1,178))

colnames(data4reg)[2:9]= c("ap_slope_mean","ap_slope_std","alpha_cf_mean", "alpha_cf_std", "alpha_pow_mean", "alpha_pow_std", "MW_power", 'MAE')



# run regression condition
lm6=glm(cond ~ alpha_cf_mean + alpha_cf_std + alpha_pow_mean + alpha_pow_std+ ap_slope_mean + ap_slope_std , data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)

# run regression condition
lm6=glm(cond ~ alpha_cf_mean + alpha_cf_std + alpha_pow_mean + alpha_pow_std+ ap_slope_mean + ap_slope_std + MAE , data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)


# remove NA and setup data 4 bayes factor analysis
data4BF= data4reg %>% select(cond, alpha_cf_mean, alpha_cf_std, alpha_pow_mean, alpha_pow_std, ap_slope_mean, ap_slope_std)
data4BF= data4BF[!is.na(data4BF$alpha_cf_mean),]
colnames(data4BF)[2:7]= c("mean alpha center frequency", "sd alpha center frequency", "mean alpha power", "sd alpha power","mean aperiodic slope", "sd aperiodic slope")

bf = regressionBF(cond ~ ., data = data4BF)
plot(bf)
bf = regressionBF(cond ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)

# repeat condition analysis with Morletwavelet power
lm6=glm(cond ~MW_power, data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)

bf = regressionBF(cond ~ MW_power, data = data4reg, whichModels = "top")
plot(bf)
summary(bf)


## plot cond effect 
error=confint(lm6)
betas=lm6$coefficients
D_betas=as.data.frame(cbind(error[1:7,],betas))
D_betas$names=c("1.Intercept", "2.mean alpha center frequency", "3.std alpha center frequency", "4.mean alpha power", "5.std alpha power","6.mean aperiodic slope", "7.std aperiodic slope")
D_betas$names=as.factor(D_betas$names)
colnames(D_betas)[1:2]=c("upper", "lower")
D_betas$inter=as.factor(c(1,0,0,0,0,0,0))

ggplot(D_betas, aes(x=betas, y=names)) + geom_rect(data = D_betas, aes(xmin = -Inf, xmax = Inf,ymax=1.5, ymin=0.5, fill = inter, alpha = 0.95)) + scale_fill_manual(values = c("white", "gray53"))+ geom_point(shape=23,size=10, fill= "black")+ geom_line(aes(size=10)) + geom_errorbarh(aes( xmin=lower, xmax=upper), D_betas, height=0, size=2.3) + theme_minimal() +
  geom_vline(xintercept=c(-0,0), linetype="dotted") + scale_y_discrete(position = "left", labels=c("1.Intercept" = "Intercept","2.mean alpha center frequency"= "mean alpha center frequency","3.std alpha center frequency"= "std alpha center frequency","4.mean alpha power"= "mean alpha power", "5.std alpha power" ="std alpha power", "6.mean aperiodic slope" = "mean aperiodic slope", "7.std aperiodic slope" = "std aperiodic slope")) + scale_x_continuous(position = "top")+ labs(y= "", x ="Fixed Effects Coefficients (β)") + theme(legend.position = "none", axis.line.x.top = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(text = element_text(size=20))


colour_scheme= c("#859fe6", "#f2b968")
data4plot = cbind(data4reg[c(1,9)], stack(data4reg %>% select( ap_slope_mean, alpha_pow_mean, alpha_pow_std)))
data4plot$cond=as.factor(data4plot$cond)   
pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=cond , y = values)) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "bar", size= 10, position = pd) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size=2.5, width=0.5, position = pd) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free")

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=cond , y = values, fill= cond)) + geom_jitter(shape=16, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_minimal() + facet_wrap(~ ind, nrow=3, scales = "free") + scale_x_discrete(labels=c("0" = "eyes-closed","1" = "eyes-open")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")




###################################################################################
#                                 Import and clean data with R
###################################################################################
EEG_data=read.csv('/Users/jasondsc/Desktop/behav_LEMON/LEMON_Oz_FOOOFrevision_rmbad.csv', stringsAsFactors = FALSE)
EEG_data2=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_rmbad.csv', stringsAsFactors = FALSE)
demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)

colnames(EEG_data)[1]= 'subject'

###################################################################################
#                                 Merge and clean data
###################################################################################

# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data,EEG_data2, tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced

colnames(data)[c(2,5)]= c('ap_exp_eo', 'ap_exp_ec')

###################################################################################
#                          Run regression analysis to predict Age
###################################################################################


# run regression eyes-open
lm6=glm(Age2 ~ ap_exp_eo + pk_cfr_eo + pk_amp_eo, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_eo,  pk_cfr_eo, pk_amp_eo)
data4BF= data4BF[!is.na(data4BF$pk_cfr_eo),]
colnames(data4BF)[2:4]= c("Eyes-open aperiodic exponent", "Eyes-openalpha center frequency", "Eyes-open alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# run regression eyes-closed
lm6=glm(Age2 ~ ap_exp_ec + pk_cfr_ec + pk_amp_ec, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_ec,  pk_cfr_ec, pk_amp_ec)
data4BF= data4BF[!is.na(data4BF$pk_cfr_ec),]
colnames(data4BF)[2:4]= c("Eyes-closed aperiodic exponent", "Eyes-closed alpha center frequency", "Eyes-closed alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)



# plot effects 

colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,53)], stack(data %>% select( pk_amp_ec, ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave('/Users/jasondsc/Desktop/behav_LEMON/figure_new_June2022_eyesclosed_effect_fooofVSsprint2.pdf', device = "pdf", width=10, height=8, dpi=800)




colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,53)], stack(data %>% select( pk_amp_ec)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")


colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,53)], stack(data %>% select( ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")



###################################################################################
#                          Run regression analysis to predict Condition
###################################################################################



# prep data to look at condition effect (eyes-open vs eyes-closed)
data4reg = cbind(data[c(1)], stack(data %>% select( pk_cfr_ec, pk_cfr_eo)))
data4reg = cbind(data4reg, stack(data %>% select( pk_amp_ec, pk_amp_eo)))
data4reg = cbind(data4reg, stack(data %>% select( ap_exp_ec, ap_exp_eo)))
data4reg=data4reg[c(1,2,4,6)]
data4reg$cond= c(rep(0,178),rep(1,178))

colnames(data4reg)[2:4]= c("alpha_cf","alpha_amp", "expon")

# run regression condition
lm6=glm(cond ~ alpha_cf + alpha_amp + expon, data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)


# remove NA and setup data 4 bayes factor analysis
data4BF= data4reg %>% select(cond, alpha_cf, alpha_amp, expon)
data4BF= data4BF[!is.na(data4BF$alpha_cf),]
colnames(data4BF)[2:3]= c("alpha cf", "alpha amp", "expon")

bf = regressionBF(cond ~ ., data = data4BF)
plot(bf)
bf = regressionBF(cond ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)







###################################################################################
#                                 Import and clean data with R
###################################################################################
EEG_data=read.csv('/Users/jasondsc/Desktop/behav_LEMON/LEMON_Oz_FOOOFrevision2.csv', stringsAsFactors = FALSE)
EEG_data2=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_final.csv', stringsAsFactors = FALSE)
demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)

colnames(EEG_data)[1]= 'subject'

###################################################################################
#                                 Merge and clean data
###################################################################################

# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data,EEG_data2, tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced

colnames(data)[c(2,5)]= c('ap_exp_eo', 'ap_exp_ec')

###################################################################################
#                          Run regression analysis to predict Age
###################################################################################


# run regression eyes-open
lm6=glm(Age2 ~ ap_exp_eo + pk_cfr_eo + pk_amp_eo, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_eo,  pk_cfr_eo, pk_amp_eo)
data4BF= data4BF[!is.na(data4BF$pk_cfr_eo),]
colnames(data4BF)[2:4]= c("Eyes-open aperiodic exponent", "Eyes-openalpha center frequency", "Eyes-open alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# run regression eyes-closed
lm6=glm(Age2 ~ ap_exp_ec + pk_cfr_ec + pk_amp_ec, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_ec,  pk_cfr_ec, pk_amp_ec)
data4BF= data4BF[!is.na(data4BF$pk_cfr_ec),]
colnames(data4BF)[2:4]= c("Eyes-closed aperiodic exponent", "Eyes-closed alpha center frequency", "Eyes-closed alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)



# corr values of the two methods
cor.test(data$eo_alpha_cf_mean, data$pk_cfr_eo)
cor.test(data$eo_alpha_pow_mean, data$pk_amp_eo)
cor.test(data$eo_ap_slope_mean, data$ap_exp_eo)
cor.test(data$ec_alpha_cf_mean, data$pk_cfr_ec)
cor.test(data$ec_alpha_pow_mean, data$pk_amp_ec)
cor.test(data$ec_ap_slope_mean, data$ap_exp_ec)


cor.test(data$ec_alpha_pow_mean[data$pk_amp_ec >0.8], data$pk_amp_ec[data$pk_amp_ec >0.8])
cor.test(data$ec_alpha_pow_mean[data$pk_amp_ec <0.8], data$pk_amp_ec[data$pk_amp_ec <0.8])

cor.test(data$eo_alpha_pow_mean[data$pk_amp_eo >0.8], data$pk_amp_eo[data$pk_amp_eo >0.8])
cor.test(data$eo_alpha_pow_mean[data$pk_amp_eo <0.8], data$pk_amp_eo[data$pk_amp_eo <0.8])


# plot effects 

colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,45)], stack(data %>% select( pk_amp_ec, ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave('/Users/jasondsc/Desktop/behav_LEMON/figure_new_June2022_eyesclosed_effect_fooofVSsprint.pdf', device = "pdf", width=10, height=8, dpi=800)




colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,55)], stack(data %>% select( pk_amp_ec)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")


colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,55)], stack(data %>% select( ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")



###################################################################################
#                          Run regression analysis to predict Condition
###################################################################################



# prep data to look at condition effect (eyes-open vs eyes-closed)
data4reg = cbind(data[c(1)], stack(data %>% select( pk_cfr_ec, pk_cfr_eo)))
data4reg = cbind(data4reg, stack(data %>% select( pk_amp_ec, pk_amp_eo)))
data4reg = cbind(data4reg, stack(data %>% select( ap_exp_ec, ap_exp_eo)))
data4reg=data4reg[c(1,2,4,6)]
data4reg$cond= c(rep(0,178),rep(1,178))

colnames(data4reg)[2:4]= c("alpha_cf","alpha_amp", "expon")

# run regression condition
lm6=glm(cond ~ alpha_cf + alpha_amp + expon, data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)


# remove NA and setup data 4 bayes factor analysis
data4BF= data4reg %>% select(cond, alpha_cf, alpha_amp, expon)
data4BF= data4BF[!is.na(data4BF$alpha_cf),]
colnames(data4BF)[2:3]= c("alpha cf", "alpha amp", "expon")

bf = regressionBF(cond ~ ., data = data4BF)
plot(bf)
bf = regressionBF(cond ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# parameterized Morlet Wavelet
###################################################################################
#                                 Import and clean data with R
###################################################################################
EEG_data=read.csv('/Users/jasondsc/Desktop/behav_LEMON/MW_Oz.csv', stringsAsFactors = FALSE)
demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)

colnames(EEG_data)[1]= 'subject'

###################################################################################
#                                 Merge and clean data
###################################################################################

# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data, tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced

###################################################################################
#                          Run regression analysis to predict Age
###################################################################################


# run regression eyes-open
lm6=glm(Age2 ~ ap_exp_eo + pk_cfr_eo + pk_amp_eo, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_eo,  pk_cfr_eo, pk_amp_eo)
data4BF= data4BF[!is.na(data4BF$pk_cfr_eo),]
colnames(data4BF)[2:4]= c("Eyes-open aperiodic exponent", "Eyes-openalpha center frequency", "Eyes-open alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# run regression eyes-closed
lm6=glm(Age2 ~ ap_exp_ec + pk_cfr_ec + pk_amp_ec, data, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table


# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_ec,  pk_cfr_ec, pk_amp_ec)
data4BF= data4BF[!is.na(data4BF$pk_cfr_ec),]
colnames(data4BF)[2:4]= c("Eyes-closed aperiodic exponent", "Eyes-closed alpha center frequency", "Eyes-closed alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)



# corr values of the two methods
cor.test(data$eo_alpha_cf_mean, data$pk_cfr_eo)
cor.test(data$eo_alpha_pow_mean, data$pk_amp_eo)
cor.test(data$eo_ap_slope_mean, data$ap_exp_eo)
cor.test(data$ec_alpha_cf_mean, data$pk_cfr_ec)
cor.test(data$ec_alpha_pow_mean, data$pk_amp_ec)
cor.test(data$ec_ap_slope_mean, data$ap_exp_ec)


cor.test(data$ec_alpha_pow_mean[data$pk_amp_ec >0.8], data$pk_amp_ec[data$pk_amp_ec >0.8])
cor.test(data$ec_alpha_pow_mean[data$pk_amp_ec <0.8], data$pk_amp_ec[data$pk_amp_ec <0.8])

cor.test(data$eo_alpha_pow_mean[data$pk_amp_eo >0.8], data$pk_amp_eo[data$pk_amp_eo >0.8])
cor.test(data$eo_alpha_pow_mean[data$pk_amp_eo <0.8], data$pk_amp_eo[data$pk_amp_eo <0.8])


# plot effects 

colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,45)], stack(data %>% select( pk_amp_ec, ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave('/Users/jasondsc/Desktop/behav_LEMON/figure_new_June2022_eyesclosed_effect_fooofVSsprint.pdf', device = "pdf", width=10, height=8, dpi=800)




colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,55)], stack(data %>% select( pk_amp_ec)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")


colour_scheme= c("#859fe6", "#b5c7f7")

data4plot = cbind(data[c(1,55)], stack(data %>% select( ec_alpha_pow_mean)))
data4plot$Age2=as.factor(data4plot$Age2)   

pd <- position_dodge(-0.2)
ggplot(data4plot, aes(x=Age2 , y = values, fill= ind)) + geom_jitter(shape=16, size=5, position=position_jitter(0.2), aes(alpha=0.1)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, width= 0.2) +
  theme_classic2() + facet_wrap(~ ind, nrow=1) + scale_x_discrete(labels=c("0" = "young adults","1" = "older adults")) +
  scale_fill_manual(values=colour_scheme) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")



###################################################################################
#                          Run regression analysis to predict Condition
###################################################################################



# prep data to look at condition effect (eyes-open vs eyes-closed)
data4reg = cbind(data[c(1)], stack(data %>% select( pk_cfr_ec, pk_cfr_eo)))
data4reg = cbind(data4reg, stack(data %>% select( pk_amp_ec, pk_amp_eo)))
data4reg = cbind(data4reg, stack(data %>% select( ap_exp_ec, ap_exp_eo)))
data4reg=data4reg[c(1,2,4,6)]
data4reg$cond= c(rep(0,178),rep(1,178))

colnames(data4reg)[2:4]= c("alpha_cf","alpha_amp", "expon")

# run regression condition
lm6=glm(cond ~ alpha_cf + alpha_amp + expon, data4reg, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL)


# remove NA and setup data 4 bayes factor analysis
data4BF= data4reg %>% select(cond, alpha_cf, alpha_amp, expon)
data4BF= data4BF[!is.na(data4BF$alpha_cf),]
colnames(data4BF)[2:3]= c("alpha cf", "alpha amp", "expon")

bf = regressionBF(cond ~ ., data = data4BF)
plot(bf)
bf = regressionBF(cond ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


###################################################################################
#                                 Import and clean data with R
###################################################################################
EEG_data=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_pkdet_2.csv', stringsAsFactors = FALSE)
EEG_data2=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_pkdet.csv', stringsAsFactors = FALSE)

demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)

colnames(EEG_data)[1]= 'subject'

###################################################################################
#                                 Merge and clean data
###################################################################################

# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data, tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced


###################################################################################
#                          Run regression analysis to predict Age
###################################################################################


# run regression eyes-open
lm6=glm(Age2 ~ eo_alpha_pow_mean, data, family = 'binomial')
lm6=glm(Age2 ~ eo_alpha_cf_mean + eo_alpha_cf_std + eo_alpha_pow_mean + eo_alpha_pow_std+ eo_ap_slope_mean + eo_ap_slope_std + eo_alpha_pkdet , data, family = 'binomial')
lm6=glm(Age2 ~ eo_alpha_pkdet*eo_alpha_pow_mean +eo_alpha_pkdet + eo_alpha_cf_mean + eo_alpha_cf_std + eo_alpha_pow_mean + eo_alpha_pow_std+ eo_ap_slope_mean + eo_ap_slope_std , data, family = 'binomial')



# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2, eo_alpha_pkdet,  eo_alpha_cf_mean,  eo_alpha_cf_std, eo_alpha_pow_mean, eo_alpha_pow_std, eo_ap_slope_mean, eo_ap_slope_std)
data4BF= data4BF[!is.na(data4BF$eo_alpha_cf_mean),]
colnames(data4BF)[2:8]= c("alpha peak detection", "mean alpha center frequency", "sd alpha center frequency", "mean alpha power", "sd alpha power","mean aperiodic slope", "sd aperiodic slope")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)


# run regression eyes-open
lm6=glm(Age2 ~ ec_alpha_pkdet, data, family = 'binomial')
lm6=glm(Age2 ~ ec_alpha_pkdet + ec_alpha_cf_mean + ec_alpha_cf_std + ec_alpha_pow_mean + ec_alpha_pow_std+ ec_ap_slope_mean + ec_ap_slope_std , data, family = 'binomial')

lm6=glm(Age2 ~ ec_alpha_pkdet*ec_alpha_pow_mean +ec_alpha_pkdet + ec_alpha_cf_mean + ec_alpha_cf_std + ec_alpha_pow_mean + ec_alpha_pow_std+ ec_ap_slope_mean + ec_ap_slope_std , data, family = 'binomial')

summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table



# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,, ap_exp_ec,  pk_cfr_ec, pk_amp_ec)
data4BF= data4BF[!is.na(data4BF$pk_cfr_ec),]
colnames(data4BF)[2:4]= c("Eyes-closed aperiodic exponent", "Eyes-closed alpha center frequency", "Eyes-closed alpha power")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)






###################################################################################
#                                 Import and clean block data
###################################################################################
EEG_data=read.csv('~/Desktop/behav_LEMON/tFOOOF_Oz_blockrevision.csv', stringsAsFactors = FALSE)
demo=read.csv('~/Desktop/behav_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv',stringsAsFactors = FALSE)


# get ids of subjects and remove "bad subjects" see manuscript for details 
subj=unique(EEG_data$subject) # get unqiue list of part ids
EEG_data=EEG_data[!subj %in% c( 'sub-010050'),]  #'sub-010044'
subj=subj[!subj %in% c( 'sub-010050')]  #'sub-010044',

tt5=demo[demo$ID %in% subj,]
tt5=tt5[order(tt5$ID),]

data=cbind(EEG_data, tt5[,-1]) # merge datasets 

data$Age2=plyr::mapvalues(data$Age, unique(data$Age), c(1, 0, 0, 1, 0, 1, 1, 1, 0)) # codde age groups according to actual distance: 1 - 20-25
data$Age2=as.numeric(data$Age2)                                                       # and 12- 75-80, this way values are equally spaced

data_eo= cbind(data[,c(1, 227:247)] ,stack(data[,2:11]), stack(data[,12:21]), stack(data[,22:31]), stack(data[,32:41]), stack(data[,42:51]),
               stack(data[,52:61]),stack(data[,62:71]),stack(data[,72:81]) )

data_ec= cbind(data[,c(1, 227:247)] ,stack(data[,82:96]),stack(data[,97:111]),stack(data[,112:126]),stack(data[,127:141]),
               stack(data[,142:156]), stack(data[,157:171]), stack(data[,172:186]), stack(data[,187:201]) )


colnames(data_eo)[c(23, 25, 27, 29, 31, 33, 35, 37)]= c('eo_ap_slope_mean', 'eo_ap_slope_std', 'eo_ap_off_mean', 'eo_ap_off_std', 'eo_alpha_cf_mean', 'eo_alpha_cf_std', 'eo_alpha_pow_mean', 'eo_alpha_pow_std')

colnames(data_ec)[c(23, 25, 27, 29, 31, 33, 35, 37)]= c('ec_ap_slope_mean', 'ec_ap_slope_std', 'ec_ap_off_mean', 'ec_ap_off_std', 'ec_alpha_cf_mean', 'ec_alpha_cf_std', 'ec_alpha_pow_mean', 'ec_alpha_pow_std')

write.csv(data_eo, '~/Desktop/behav_LEMON/stacked_eo_tFOOOF_Oz_blockrevision.csv')
write.csv(data_ec, '~/Desktop/behav_LEMON/stacked_ec_tFOOOF_Oz_blockrevision.csv')



###################################################################################
#                          Run regression analysis to predict Age
###################################################################################
library(lme4)

# run regression eyes-open
lm6=glmer(Age2 ~ (1| subject) + eo_alpha_cf_mean + eo_alpha_cf_std + eo_alpha_pow_mean + eo_alpha_pow_std+ eo_ap_slope_mean + eo_ap_slope_std , data_eo, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table

# remove NA and setup data 4 bayes factor analysis
data4BF= data %>% select(Age2,eo_alpha_cf_mean, eo_alpha_cf_std, eo_alpha_pow_mean, eo_alpha_pow_std, eo_ap_slope_mean, eo_ap_slope_std)
data4BF= data4BF[!is.na(data4BF$eo_alpha_cf_mean),]
colnames(data4BF)[2:7]= c("Eyes-open mean alpha center frequency", "Eyes-open std alpha center frequency", "Eyes-open mean alpha power", "Eyes-open std alpha power","Eyes-open mean aperiodic slope", "Eyes-open std aperiodic slope")
bf = regressionBF(Age2 ~ ., data = data4BF)
plot(bf)
bf = regressionBF(Age2 ~ ., data = data4BF, whichModels = "top")
plot(bf)
summary(bf)



# run regression eyes-open
lm6=glmer(Age2 ~ (1| subject) + ec_alpha_cf_mean + ec_alpha_cf_std + ec_alpha_pow_mean + ec_alpha_pow_std+ ec_ap_slope_mean + ec_ap_slope_std , data_ec, family = 'binomial')
summary(lm6)
confint(lm6)
tab_model(lm6, transform = NULL) # make APA table
