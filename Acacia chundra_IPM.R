library(dplyr)
library(ggplot2)
library(tidyr)
# if(!requireNamespace('remotes', quietly = TRUE)) {
#   install.packages('remotes')
# }

# remotes::install_github('levisc8/ipmr', build_vignettes = TRUE, force = TRUE)
library(ipmr)


# Read data ---------------------------------------------------------------
chundratrees <- read.csv2("chundratree.csv", sep= ",") %>% na_if("")
str(chundratrees)
chundratrees$size2009 <- as.numeric(as.character(chundratrees$size2009))
chundratrees$size2016 <- as.numeric(as.character(chundratrees$size2016))


#read sprout data
chundrasprouts <- read.csv2("chundrasprouts.csv", sep= ",")

# logsize
chundratrees$log_size_t0<- log(chundratrees$size2009)
chundratrees$log_size_t1<- log(chundratrees$size2016)


#Subset harvest
no_harvest_chundra=subset(chundratrees, chundratrees$harvestall==0)


# 'All' (harvested and not harvested individuals) -------------------------
# chundra <-matrix(NA, ncol=26)  #Dataframe for the parameters you get during bootstrapping +
# 
# colnames(chundra) <- c("g_int", "g_slope", "g_sd", "s_int","s_slope", "m_int","m_slope", "d_mu", "d_sd", "a", "b", "lambda_orgs", "S_g_int", "S_g_slope", "S_g_sd", "S_s_int","S_s_slope", "S_m_int", "S_m_slope", c(1:7)) #Adding names to columns

chundra <- matrix(NA, ncol=31) #Dataframe for the parameters you get during bootstrapping +


colnames(chundra) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda_orgs",  "S_g_int", "S_g_slope", "S_g_sd", "S_s_all", "S_m_int", "S_m_slope",  "S_d_mean", "S_d_sd","S_a", "S_b", c(1:10))   ### Adding names to columns


#trees 
n           <- nrow(chundratrees)
#survival
survivors   <- subset(chundratrees, chundratrees$survival=='1')
chundra[4] <- length(survivors$survival)/n
# S2_all <-   glm(survival ~ log_size_t0 , data = chundratrees, family="binomial")
# summary(S2_all)
# 
# chundra[4]=coef(S2_all)[1]
# chundra[5]=coef(S2_all)[2]

#growth
growers <- subset(chundratrees, chundratrees$survival=='1' & chundratrees$sprout2016=='0')
G1_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G1_all)
chundra[1]=coefficients(G1_all)[1]
chundra[2]=coefficients(G1_all)[2]
chundra[3]=sd(resid(G1_all))
#tree to sprout
M2_mod <- glm(sprout2016 ~ log_size_t0, data = chundratrees,family="binomial") 
summary(M2_mod)
chundra[5] =coef(M2_mod)[1]
chundra[6]=coef(M2_mod)[2]

#sprouts
n2 = nrow(chundrasprouts)
#survival sprouts
survivingsprouts=subset(chundrasprouts, chundrasprouts$survival=='1') 
chundra[7]=length(survivingsprouts$survival)/n2
chundra[8]=0  
chundra[9]=0
chundra[10]=0
#sprout to tree probability
# gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
# chundra[9]=length(gototree$sprout2016)/length(survivingsprouts$survival)
# #size distribution sprout to tree
# gototree2=subset(gototree, gototree$size2016>0)
# if (length(gototree2$size2016)<10) {
#   gumi[8]=1.43193
#   gumi[9]=0.2112621
# }else
# {gumi[8]=mean(gototree2$size2016)
# gumi[9]=sd(gototree2$size2016)
# }
# 
# chundra[8]=0  #surv
# chundra[9]=0  
# chundra[10]=0
# chundra[11]=0

run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,
      family        = 'CD',
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = 0,  #a * (1-b),
      family  = "DD",
      data_list     = data_list, 
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = 0,
      #f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', 'sprout')),
      has_hier_effs = FALSE,
      evict_cor     = FALSE
      # evict_fun     = truncated_distributions('norm','f_d')
    )  %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  # The lower and upper bounds for the continuous state variable and the number of meshpoints for the midpoint rule integration.
  L <- log(0.5)*1.2 #minimum value is negative, hence multiplying by 1.2, if not 0.8
  U <- log(30)  
  n <- 100
  init_pop_vec  <- rpois(100, 2)
  init_sprout   <- 30
  
  ipm <- general_ipm %>%
    define_domains(
      size = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size = init_pop_vec,
        n_sprout  = init_sprout
      )
    ) %>%
    make_ipm(iterations = 100, normalize_pop_size = FALSE)
  
  
  return(ipm)  #return ipm instead of lambda
}

pars_orgs <- list(                           ##Original dataset
  g_int    = chundra[1],
  g_slope  = chundra[2],
  g_sd     = chundra[3],
  s_all    = chundra[4],
  m_int    = chundra[5],
  m_slope  = chundra[6],
  d_mu     = chundra[7],
  d_sd     = chundra[8],
  a        = chundra[9],
  b        = chundra[10]
) 


ipm_temp=run_ipm(pars_orgs) 

chundra[11] = lambda(ipm_temp)


# chundraall[i,6] = lambda(ipm_temp, comp_method = 'eigen')


for(p in 1:10){        # Now for each sample add 0.001 to each parameter, one by one
  
  
  pars_man      <- pars_orgs  ## Copy the original parameters into a manual one
  pars_man[[p]] <- pars_man[[p]] + 0.001 # p is the specific parameter
  
  #with this new list rerun the ipm to get the new lambda
  ipm_man   <- run_ipm(pars_man)
  lambda_man <- lambda(ipm_man)
  sens              <- ((lambda_man - chundra[11]) / 0.001)
  chundra[(11+p)]  <- sens          #all the new sensitivities
  chundra[(21+p)] <- lambda_man   #all the new lambdas
}

chundra

write.csv (chundra, file= "chundraobs14oct.csv")


time_series <- matrix(0, nrow = 100, ncol = 6)
for(i in 2:101) {
  
  time_series[(i-1), ] <- unlist(collapse_pop_state(ipm_temp,
                                                    time_step = i,
                                                    small1 = size <= 1,
                                                    small2 = size > 1 & size <= 1.5,
                                                    medium1= size >1.5 & size <= 2, 
                                                    medium2= size> 2 & size <= 2.5,
                                                    big1 = size > 2.5 & size <= 3,
                                                    big2 = size> 3 & size <= 3.5))
}



#ipm_temp is IPM run with 'All' data

colnames(time_series) <- c("small1","small2", "medium1", "medium2", "big1", "big2")

#get n_size #number of trees
acall1 <- time_series %>% as.data.frame

#get n_sprout #number of sprouts
acall2 <- as.data.frame(t(ipm_temp$pop_state$n_sprout))

# #write csv or cbind
write.csv(acall1,file="all1_12oct.csv")
write.csv(acall2,file="all2_12oct.csv")

#nsize and nsprout in single file
#popprojall <- read.csv2("chundraprojall.csv", sep= ",") %>% na_if("")
#popprojall <- read.csv2("sizeprojchundra.csv", sep= ",") %>% na_if("")
popprojall <- read.csv2("sizeproj13oct.csv", sep= ",") %>% na_if("")

popall <- popprojall  %>% 
  pivot_longer(small1:sprout, names_to = "stages", values_to = "number")
str(popall)
popall$number <-as.numeric(as.character(popall$number))
str(popall)
#Plot
popall %>% 
  subset( time < 140) %>% 
  ggplot( aes(y=log(number),x=time,color=factor(stages)))+
  geom_line(stat="identity",position="identity",
            lwd=1)+
  scale_fill_discrete(name="stages",
                      breaks=c(1, 2, 3, 4, 5, 6, 7),
                      labels=c("small1","small2", "medium1", "medium2", "big1", "big2", "sprout"))+
  xlab("number of years")+ylab("number of individuals") + ggtitle("Chundra All")+
  theme_classic()

# All Bootstrap ---------------------------------------
B = 1000
# chundraall <-matrix(NA, nrow=B, ncol=28)  #Dataframe for the parameters you get during bootstrapping +
# 
# colnames(chundraall) <- c("g_int", "g_slope", "g_sd", "s_int","s_slope", "m_int","m_slope", "d_mu", "d_sd", "a", "b", "lambda_orgs", "S_g_int", "S_g_slope", "S_g_sd", "S_s_int","S_s_slope", "S_m_int", "S_m_slope", c(1:9))  
# ### Adding names to columns  

chundraall <- matrix(NA, nrow=B, ncol=31) #Dataframe for the parameters you get during bootstrapping +


colnames(chundraall) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda_orgs",  "S_g_int", "S_g_slope", "S_g_sd", "S_s_all", "S_m_int", "S_m_slope",  "S_d_mean", "S_d_sd","S_a", "S_b", c(1:10))   ### Adding names to columns

for (i in 1:B){
  
  set.seed(i)
  print(i)
#trees 
n           <- nrow(chundratrees)
x3          <- sample(1:n, n, replace = TRUE)
databoot    <- chundratrees[x3,]
#survival
survivors   <- subset(databoot, databoot$survival=='1')
chundraall[i,4] <- length(survivors$survival)/n
# S2_all <-   glm(survival ~ log_size_t0 , data = chundratrees, family="binomial")
# summary(S2_all)
# 
# chundra[4]=coef(S2_all)[1]
# chundra[5]=coef(S2_all)[2]

#growth
growers <- subset(databoot, databoot$survival=='1' & databoot$sprout2016=='0')
G1_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G1_all)
chundraall[i,1]=coefficients(G1_all)[1]
chundraall[i,2]=coefficients(G1_all)[2]
chundraall[i,3]=sd(resid(G1_all))
#tree to sprout
M2_mod <- glm(sprout2016 ~ log_size_t0, data = databoot,family="binomial") 
summary(M2_mod)
chundraall[i,5] =coef(M2_mod)[1]
chundraall[i,6]=coef(M2_mod)[2]

#sprouts
n2 = nrow(chundrasprouts)
x2 <- sample(1:n2, n2, replace = TRUE)
databoot2 <- chundrasprouts[x2,]
#survival sprouts
survivingsprouts=subset(databoot2, databoot2$survival=='1') 
chundraall[i,7]=length(survivingsprouts$survival)/n2
chundraall[i,8]=0  
chundraall[i,9]=0
chundraall[i,10]=0
#sprout to tree probability
# gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
# chundra[9]=length(gototree$sprout2016)/length(survivingsprouts$survival)
# #size distribution sprout to tree
# gototree2=subset(gototree, gototree$size2016>0)
# if (length(gototree2$size2016)<10) {
#   gumi[8]=1.43193
#   gumi[9]=0.2112621
# }else
# {gumi[8]=mean(gototree2$size2016)
# gumi[9]=sd(gototree2$size2016)
# }
# 
# chundra[8]=0  #surv
# chundra[9]=0  
# chundra[10]=0
# chundra[11]=0

run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,
      family        = 'CD',
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = 0,  #a * (1-b),
      family  = "DD",
      data_list     = data_list, 
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = 0,
      #f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', 'sprout')),
      has_hier_effs = FALSE,
      evict_cor     = FALSE
      # evict_fun     = truncated_distributions('norm','f_d')
    )  %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  # The lower and upper bounds for the continuous state variable and the number of meshpoints for the midpoint rule integration.
  L <- log(0.5)*1.2 #minimum value is negative, hence multiplying by 1.2, if not 0.8
  U <- log(30)  
  n <- 100
  init_pop_vec  <- rpois(100, 2)
  init_sprout   <- 30
  
  ipm <- general_ipm %>%
    define_domains(
      size = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size = init_pop_vec,
        n_sprout  = init_sprout
      )
    ) %>%
    make_ipm(iterations = 100, normalize_pop_size = FALSE)
  
  
  return(ipm)  #return ipm instead of lambda
}

pars_orgs <- list(                           ##Original dataset
  g_int    = chundraall[i,1],
  g_slope  = chundraall[i,2],
  g_sd     = chundraall[i,3],
  s_all    = chundraall[i,4],
  m_int    = chundraall[i,5],
  m_slope  = chundraall[i,6],
  d_mu     = chundraall[i,7],
  d_sd     = chundraall[i,8],
  a        = chundraall[i,9],
  b        = chundraall[i,10]
) 


ipm_temp=run_ipm(pars_orgs) 

chundraall[i,11] = lambda(ipm_temp)


# chundraall[i,6] = lambda(ipm_temp, comp_method = 'eigen')


for(p in 1:10){        # Now for each sample add 0.001 to each parameter, one by one
  
  
  pars_man      <- pars_orgs  ## Copy the original parameters into a manual one
  pars_man[[p]] <- pars_man[[p]] + 0.001 # p is the specific parameter
  
  #with this new list rerun the ipm to get the new lambda
  ipm_man   <- run_ipm(pars_man)
  lambda_man <- lambda(ipm_man)
  sens              <- ((lambda_man - chundraall[i,11]) / 0.001)
  chundraall[i,(11+p)]  <- sens          #all the new sensitivities
  chundraall[i,(21+p)] <- lambda_man   #all the new lambdas
}
  
  
}


y<-  data.frame(mean_est=colMeans(chundraall),
                t(apply(chundraall,2,quantile,c(0.025,0.975))))
write.csv(y,file="chundraboot14oct.csv")



# No harvest --------------------------------------------------------------

chundranoharv <-matrix(NA, ncol=11)  #Dataframe for the parameters you get during bootstrapping +

colnames(chundranoharv) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int","m_slope", "d_mu", "d_sd", "a", "b", "lambda")   ### Adding names to columns  

n           <- nrow(no_harvest_chundra)

survivors   <- subset(no_harvest_chundra, no_harvest_chundra$survival=='1')
chundranoharv[4] <- length(survivors$survival)/n

growers <- subset(no_harvest_chundra, no_harvest_chundra$survival=='1' & no_harvest_chundra$sprout2016=='0')
#growth
G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growers)
chundranoharv[1]=coefficients(G1_noharv)[1]
chundranoharv[2]=coefficients(G1_noharv)[2]
chundranoharv[3]=sd(resid(G1_noharv))
#tree to sprout
M2_noharv <- glm(sprout2016 ~ log_size_t0, data = no_harvest_chundra,family="binomial") 
chundranoharv[5] =coef(M2_noharv)[1]
chundranoharv[6] =coef(M2_noharv)[2]

chundranoharv[7] =0  #surv
chundranoharv[8] =0  
chundranoharv[9]=0
chundranoharv[10]=0

# Run IPM -----------------------------------------------------------------


run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,
      family        = 'CD',
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = 0,  #a * (1-b),
      family  = "DD",
      data_list     = data_list, 
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = 0,
      #f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', 'sprout')),
      has_hier_effs = FALSE,
      evict_cor     = FALSE
      # evict_fun     = truncated_distributions('norm','f_d')
    )  %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  # The lower and upper bounds for the continuous state variable and the number of meshpoints for the midpoint rule integration.
  L <- log(0.5)*1.2
  U <- log(30)
  n <- 100
  init_pop_vec  <- rpois(100, 2)
  init_sprout   <- 30
  
  ipm <- general_ipm %>%
    define_domains(
      size = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size = init_pop_vec,
        n_sprout  = init_sprout
      )
    ) %>%
    make_ipm(iterations = 100, normalize_pop_size = FALSE)
  
  
  return(ipm)  #return ipm instead of lambda
}

pars_orgs <- list(                           ##Original dataset
  g_int     = chundranoharv[1],
  g_slope   = chundranoharv[2],
  g_sd      = chundranoharv[3],
  s_all    = chundranoharv[4],
  m_int    = chundranoharv[5],
  m_slope = chundranoharv[6],
  d_mu    = chundranoharv[7],
  d_sd    = chundranoharv[8],
  # d_mu    = 1.8453,
  # d_sd    = 0.722706,
  a       = chundranoharv[9],
  b       = chundranoharv[10]
) 

ipm_temp3=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

chundranoharv[11] = lambda(ipm_temp3)

chundranoharv
write.csv(chundranoharv, file="chundranoharvobs14oct.csv")


#Bootstrap
B = 1000


chundranoharvest <-matrix(NA,nrow=B, ncol=11)  #Dataframe for the parameters you get during bootstrapping +

colnames(chundranoharvest) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int","m_slope", "d_mu", "d_sd", "a", "b", "lambda")   ### Adding names to columns  

for (i in 1:B){
  
  set.seed( i )
  print(i)
  
  n           <- nrow(no_harvest_chundra)
  x3          <- sample(1:n, n, replace = TRUE)
  databoot    <- no_harvest_chundra[x3,]

survivors   <- subset(databoot, databoot$survival=='1')
chundranoharvest[i,4] <- length(survivors$survival)/n

growers <- subset(databoot, databoot$survival=='1' & databoot$sprout2016=='0')
#growth
G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growers)
chundranoharvest[i,1]=coefficients(G1_noharv)[1]
chundranoharvest[i,2]=coefficients(G1_noharv)[2]
chundranoharvest[i,3]=sd(resid(G1_noharv))
#tree to sprout
M2_noharv <- glm(sprout2016 ~ log_size_t0, data = databoot,family="binomial") 
chundranoharvest[i,5] =coef(M2_noharv)[1]
chundranoharvest[i,6] =coef(M2_noharv)[2]


n2 = nrow(chundrasprouts)
x2 <- sample(1:n2, n2, replace = TRUE)
databoot2 <- chundrasprouts[x2,]
#survival sprouts
survivingsprouts=subset(databoot2, databoot2$survival=='1') 
chundranoharvest[i,7]=length(survivingsprouts$survival)/n2
chundranoharvest[i,8]=0  
chundranoharvest[i,9]=0
chundranoharvest[i,10]=0


# Run IPM -----------------------------------------------------------------


run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,
      family        = 'CD',
      s             = s_all,
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = 0,  #a * (1-b),
      family  = "DD",
      data_list     = data_list, 
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = 0,
      #f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', 'sprout')),
      has_hier_effs = FALSE,
      evict_cor     = FALSE
      # evict_fun     = truncated_distributions('norm','f_d')
    )  %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  # The lower and upper bounds for the continuous state variable and the number of meshpoints for the midpoint rule integration.
  L <- log(0.5)*1.2
  U <- log(30)
  n <- 100
  init_pop_vec  <- rpois(100, 2)
  init_sprout   <- 30
  
  ipm <- general_ipm %>%
    define_domains(
      size = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size = init_pop_vec,
        n_sprout  = init_sprout
      )
    ) %>%
    make_ipm(iterations = 100, normalize_pop_size = FALSE)
  
  
  return(ipm)  #return ipm instead of lambda
}

pars_orgs <- list(                           ##Original dataset
  g_int     = chundranoharvest[i,1],
  g_slope   = chundranoharvest[i,2],
  g_sd      = chundranoharvest[i,3],
  s_all    = chundranoharvest[i,4],
  m_int    = chundranoharvest[i,5],
  m_slope = chundranoharvest[i,6],
  d_mu    = chundranoharvest[i,7],
  d_sd    = chundranoharvest[i,8],
  # d_mu    = 1.8453,
  # d_sd    = 0.722706,
  a       = chundranoharvest[i,9],
  b       = chundranoharvest[i,10]
) 

ipm_temp4=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

chundranoharvest[i, 11] = lambda(ipm_temp4)

}

chundranoharvest


x<-  data.frame(mean_est=colMeans(chundranoharvest),
                t(apply(chundranoharvest,2,quantile,c(0.025,0.975))))
write.csv(x,file="chundrabootnoharvest14oct.csv")
