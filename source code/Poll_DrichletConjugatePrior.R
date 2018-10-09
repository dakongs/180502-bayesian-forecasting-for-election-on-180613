#### -- Conjugate Prior : Drichlet Distribution -- ####

library(MCMCpack); library(lubridate); library(reshape2); library(ggplot2);
library(ggmcmc); library(plyr); library(readxl);

rm(list = ls())
setwd('Google Drive/Foresight/20180613/')
list.files()

#### -- 1. Read and process data -- ####
df    <- read.csv('2018_kyungnam.csv')

survey        <- df[,c('조사끝', '표본크기','김경수','김태호', '기타', '무응답')]
survey$조사끝 <- as.Date(as.character(df$조사끝), "%y%m%d")

Candidate1    <- 'KimKyungSoo'
Candidate2    <- 'KimTaeHo'
year <- 2018
diff_real <- 0.099

Candidate1_init <- unlist(strsplit(Candidate1, '[[:lower:]]*'))
Candidate1_init <- paste0(Candidate1_init[Candidate1_init != ''], collapse = '')
Candidate2_init <- unlist(strsplit(Candidate2, '[[:lower:]]*'))
Candidate2_init <- paste0(Candidate2_init[Candidate2_init != ''], collapse= '')

names(survey) <- c('Date', 'N', Candidate1, Candidate2, 'Extra', 'NonResponse')

survey$N1     <- survey[,'N'] * survey[,Candidate1] / 100
survey$N2     <- survey[,'N'] * survey[,Candidate2] / 100 
survey$N3     <- survey[,'N'] * survey[,'Extra'] / 100 



#### -- 2. Set parameter -- ####

## 1) Prior : No information
alpha_noinfo     <- c(1, 1, 1)

## 2) Prior : President supporting ratio
# 2017: N - 1004, Rulling - 49%, Opposition - 13%, Extra - 14.3%
alpha_party_2018 <- c(1004*.49, 1004*.13, 1004*.14)


#### -- 3. See trend -- ####
var          <- melt(data = survey, 
                     id.vars = 'Date', 
                     measure.vars = c(Candidate1, Candidate2)) 
var$variable <- as.factor(var$variable)
names(var)   <- c('Date', 'Poll', 'value')

var %>%
  ggplot(aes(Date , value)) +
  geom_point(aes(colour = Poll)) + 
  stat_smooth(aes(colour = Poll), method = 'loess', alpha=0.5) +
  scale_y_continuous('Polls Ratio', limits = c(10, 80)) + 
  scale_color_discrete('Candidates', labels = c(Candidate1, Candidate2)) + 
  scale_color_manual(values = c('blue', 'red')) + 
  theme_bw() + 
  ggtitle(sprintf('%s vs %s (%s Elction)', Candidate1, Candidate2, year)) 
       

#### -- 4. Bayesian by Conjugate Prior -- ####
survey        <- survey[order(survey$Date, decreasing = FALSE), ]

bayes_diffs   <- freq_diff   <- ci   <- c()

## Step 1) Prior Select
# Prior : No information
alpha         <- alpha_noinfo 
# Prior : President supporting ratio
# alpha         <- alpha_party_2018
election_date <- as.Date('20180613', "%Y%m%d")

## Step 2) Calcuate Posterior
for(i in 1:nrow(survey)){

  poll_date   <- survey[i, ]$Date
  obs         <- unlist(survey[i, c('N1', 'N2', 'N3')])
  
  post        <- MCmultinomdirichlet(obs, alpha, mc = 10000)

  bayes_diffs <- append(bayes_diffs, round(mean(post[ , 1] - post[ , 2]), 3))
  
  ## Step 3) Weight Select
  # No weight
  # alpha       <- alpha + obs # No weight
  # Period Weight
  alpha       <- alpha + obs * 7 / as.numeric(election_date - poll_date)
  
  p1 <- obs[1] / sum(obs)
  p2 <- obs[2] / sum(obs)
  
  conf_interval <- qnorm(0.975) * 1 / sqrt(sum(obs)) * sqrt(p1 * (1 - p1) + p2 * (1 - p2) + 2 * p1 * p2)
  freq_diff <- append(freq_diff, p1 - p2)
  ci <- append(ci, conf_interval)

  }



## Step 4) Draw Histogram
diff_dist <- data.frame(diffs_val = as.numeric(post[ , 1] - post[ , 2]))
mdiff     <- mean(diff_dist$diffs_val)

n1_post        <- as.data.frame(as.numeric(post[,1]))
n2_post        <- as.data.frame(as.numeric(post[,2]))
n1_post$class  <- Candidate1_init
n2_post$class  <- Candidate2_init
names(n1_post) <- c('prob', 'class')
names(n2_post) <- c('prob', 'class')

post_hist      <- rbind(n1_post, n2_post)
post_hist %>%
  ggplot(aes(prob)) + 
  geom_histogram(aes(fill=class), position = 'identity', alpha=0.7) +
  scale_fill_manual(values = c('blue', 'red')) + 
  xlab('expected probability including no response') + 
  theme_bw()

# Plot Histogram
diff_dist %>%
  ggplot(aes(diffs_val)) +
  geom_histogram(color='black', fill='white') + 
  xlab(sprintf('expected probability(%s - %s)', Candidate1_init, Candidate2_init)) + 
  theme_bw() +
  geom_vline(xintercept = quantile(diff_dist$diffs_val, c(0.05, 0.95)), color = 'red')
#  geom_text(x = quantile(diff_dist$diffs_val, 0.25), y = 200, label = paste('25%', 'quantile', sep = '\n'), color = 'red') +
#  geom_text(x = quantile(diff_dist$diffs_val, 0.75), y = 200, label = paste('75%', 'quantile', sep = '\n'), color = 'red')

quantile(diff_dist$diffs_val, c(0.05, 0.95))

#### -- 5. Result -- ####
res <- rbind(data.frame(Date = survey$Date, kind = 'poll',  diff_ratio = freq_diff,   ci = ci),
             data.frame(Date = survey$Date, kind = 'bayes', diff_ratio = bayes_diffs, ci = 0))

res %>%
  ggplot(aes(x = Date, y = diff_ratio, color = kind)) + 
  geom_point() + 
  geom_line() +
  geom_hline(yintercept = 0.0, linetype = 2) +
  ylab(sprintf('%s - %s', Candidate1_init, Candidate2_init)) + 
  scale_color_discrete('Inference Methods') +
  scale_color_manual(values = c('red', 'blue', 'b lack')) + 
  theme_bw() +
  ggtitle(sprintf('%s vs %s (%s Election)', Candidate1, Candidate2, year)) +
  geom_hline(yintercept = diff_real, linetype = 2, color = 'navy', size = 1.5)
  

