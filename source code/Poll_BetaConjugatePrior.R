#### -- Conjugate Prior : Beta Distribution -- ####

library(MCMCpack); library(lubridate); library(reshape2); library(ggplot2);
library(ggmcmc); library(plyr); library(readxl);

rm(list = ls())
setwd('C:/Users/HK_PARK/Desktop/Poll')
list.files()

#### -- 1. Read and process data -- ####
df    <- read.csv('1.Data/2014_pusan.csv')
df    <- df[-nrow(df), ]

survey        <- df[,c('조사끝', '표본크기','서병수','오거돈', '기타')]
survey$조사끝 <- as.Date(as.character(df$조사끝), "%y%m%d")

Candidate1    <- 'SeoByungSoo'
Candidate2    <- 'OhGuDon'
year <- 2014
diff_real <- 0.014 

Candidate1_init <- unlist(strsplit(Candidate1, '[[:lower:]]*'))
Candidate1_init <- paste0(Candidate1_init[Candidate1_init != ''], collapse = '')
Candidate2_init <- unlist(strsplit(Candidate2, '[[:lower:]]*'))
Candidate2_init <- paste0(Candidate2_init[Candidate2_init != ''], collapse= '')

names(survey) <- c('Date', 'N', Candidate1, Candidate2, 'Extra')

survey$N1     <- survey[,'N'] * survey[,Candidate1] / 100
survey$N2     <- survey[,'N'] * survey[,Candidate2] / 100 
survey$Ex_N   <- survey[,'N'] * survey[,'Extra']    / 100 


#### -- 2. Set parameter -- ####

## 1) Prior : No information
alpha_noinfo         <- c(1, 1)                       

## 2) Prior : President supporting ratio
# 2014: N - 808, Good - 46%, Bad - 41%
alpha_president_2014 <- c(808*0.46, 808*0.41) 


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
# alpha         <- alpha_president_2014 
election_date <- as.Date('2014064', "%Y%m%d")

## Step 2) Calcuate Posterior
for(i in 1:nrow(survey)){

  poll_date   <- survey[i, ]$Date
  obs         <- survey[i, c('N1', 'N2')]
  
  posterior   <- MCbinomialbeta(y = obs$N1, n = obs$N1 + obs$N2, 
                              alpha = alpha[1], beta = alpha[2], mc = 10000)

  bayes_diffs <- append(bayes_diffs, mean(posterior - (1-posterior)))
  
  ## Step 3) Weight Select
  # No weight
  alpha       <- alpha + c(obs$N1, obs$N2) # No weight
  # Period Weight
  alpha       <- alpha + c(obs$N1, obs$N2)  * 7 / as.numeric(election_date - poll_date)

  p1          <- obs$N1 / sum(obs)
  p2          <- obs$N2 / sum(obs)
  
  conf_interval <- qnorm(0.975) * 1 / sqrt(sum(obs)) * sqrt(p1 * (1 - p1) + p2 * (1 - p2) + 2 * p1 * p2)
  freq_diff     <- append(freq_diff, p1 - p2)
  ci            <- append(ci, conf_interval)
  
}


## Step 4) Draw Histogram
diff_dist   <- data.frame(diffs_val = as.numeric(posterior - (1-posterior)))

# Plot Histogram
diff_dist %>%
  ggplot(aes(diffs_val)) +
  geom_histogram(color='black', fill='white') + 
  xlab("expected probability(Lee - Nam)") + 
  theme_bw() +
  geom_vline(xintercept = quantile(diff_dist$diffs_val, c(0.25, 0.75)), color = 'red') +
  geom_text(x = quantile(diff_dist$diffs_val, 0.25), y = 200, label = paste('25%', 'quantile', sep = '\n'), color = 'red') +
  geom_text(x = quantile(diff_dist$diffs_val, 0.75), y = 200, label = paste('75%', 'quantile', sep = '\n'), color = 'red')


#### -- 5. Result -- ####
res    <- rbind(data.frame(Date = survey$Date, kind = 'freq',  diff_ratio = freq_diff,   ci = ci),
                data.frame(Date = survey$Date, kind = 'bayes', diff_ratio = bayes_diffs, ci = 0))


res %>%
  ggplot(aes(x = Date, y = diff_ratio, color = kind)) + 
  geom_point() +
  geom_line() +
  # stat_smooth(aes(colour = kind), method = 'loess', alpha=0.5) +
  geom_hline(yintercept = diff_real, linetype = 2, color = 'navy', size= 1.5) + 
  geom_hline(yintercept = 0.0, linetype = 2) +
  ylab(sprintf('%s - %s', Candidate1_init, Candidate2_init)) + 
  scale_color_discrete('Inference Methods') +
  scale_color_manual(values = c('red', 'blue', 'b lack')) + 
  theme_bw() +
  ggtitle(sprintf('%s vs %s (%s Elction)', Candidate1, Candidate2, year)) 


