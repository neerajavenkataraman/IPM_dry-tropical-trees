#IPM implementation for Chloroxylon swietenia

# Ambient harvesting scenario ---------------------------------------------


all <- matrix(NA, ncol=28) #Dataframe for the parameters

colnames(all) <- c("g_int", "g_slope", "g_sd", "s_all", "m_all", "d_mu", "d_sd",
                   "a", "b", "lambda_orgs","S_g_int", "S_g_slope", "S_g_sd", "S_s_all", "S_m_all", "S_d_mu", "S_d_sd", "S_a", "S_b", c(1:9))   
### Adding names to columns
n           <- nrow(trees)
#tree survivorship
#S_cs <- glm(survival ~ log_size_t0, data=trees, family="binomial") 
# summary(S_cs) #not significantly related to size hence modeled as constant
survivors   <- subset(trees, trees$survival=='1') 
all[4] <- length(survivors$survival)/n
gotosprout  <- subset(trees, trees$sprout2016=='1')
#transition from tree to sprout
# M_cs <- glm(sprout2016 ~ log_size_t0, data=trees, family="binomial") 
# summary(M_cs) #not significantly related to size hence modeled as constant
all[5] <- length(gotosprout$sprout2016)/length(survivors$survival) 
#growth of trees
growers     <- subset(trees, trees$survival=='1' & trees$sprout2016=='0')
G1_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G1_all)
all[1]=coefficients(G1_all)[1]
all[2]=coefficients(G1_all)[2]
all[3]=sd(resid(G1_all))

#Sprouts
n2 = nrow(sprouts)
#Sprout survivorship
survivingsprouts=subset(sprouts, sprouts$survival=='1')
all[8]=length(survivingsprouts$survival)/n2
#Sprout to tree transition
gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
all[9]=length(gototree$sprout2016)/length(survivingsprouts$survival)
#Sprout size distribution
gototree2=subset(gototree, gototree$size2016>0)
if(length(gototree2$size2016)<4){
  all[6]=0.5412373
  all[7]=0.1164737
}else{
  all[6]=mean(gototree2$size2016, na.rm=T)
  all[7]=sd(gototree2$size2016, na.rm=T)
}

# Run IPM -----------------------------------------------------------------

run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size, #tree survivorship, growth, (1-treetosprout)
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_all,
      m             = m_all,
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,  #tree survivorship, treetosprout
      family        = 'CD',
      s             = s_all,
      m             = m_all,
      data_list     = data_list,
      states        = list(c('size', "sprout")), 
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b), #Sprout survivorship, 1-treetosprout
      family  = "DD",
      data_list     = data_list,
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d *d_size,
      f_d           = dnorm(size_2, d_mu, d_sd), #size distribution of sprouts
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', 'sprout')), 
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
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
  L <- log(1.1)*0.8
  U <- log(20)
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
  g_int    = all[1],
  g_slope  = all[2],
  g_sd     = all[3],
  s_all    = all[4],
  m_all    = all[5],
  d_mu     = all[6],
  d_sd     = all[7],
  a        = all[7],
  b        = all[8]
) 


ipm_temp=run_ipm(pars_orgs) 

all[10] = lambda(ipm_temp)


for(p in 1:9){        # add 0.001 to each parameter, one by one
  
  
  pars_man      <- pars_orgs  ## Copy the original parameters into a manual one
  pars_man[[p]] <- pars_man[[p]] + 0.001 # p is the specific parameter
  
  #rerun the ipm to get the new lambda
  ipm_man   <- run_ipm(pars_man)
  lambda_man <- lambda(ipm_man)
  
  sens              <- ((lambda_man - all[10]) / 0.001) #Calculate sensitivities
  all[(10+p)]  <- sens          #all the new sensitivities
  all[(19+p)] <- lambda_man   #all the new lambdas
}


# No harvesting scenario --------------------------------------------------

noharv <- matrix(NA, ncol=10) #Dataframe for the parameters


colnames(noharv) <- c("g_int", "g_slope", "g_sd", "s_notH", "m_notH", "d_mu", "d_sd",
                      "a", "b", "lambda"
)   ### Adding names to columns, to keep track


n = nrow(no_harvest)
#tree survivorship
#S_cs <- glm(survival ~ log_size_t0, data=no_harvest, family="binomial") 
# summary(S_cs) #not significantly related to size hence modeled as constant
survivors=subset(no_harvest, no_harvest$survival=='1')
noharv[4]=length(survivors$survival)/n
#transition from tree to sprout
# M_cs <- glm(sprout2016 ~ log_size_t0, data=no_harvest, family="binomial") 
# summary(M_cs) #not significantly related to size hence modeled as constant
gotosprout=subset(no_harvest, no_harvest$sprout2016=='1')
noharv[5]=length(gotosprout$sprout2016)/length(survivors$survival)
#growth of trees
growers=subset(no_harvest, no_harvest$survival=='1' & no_harvest$sprout2016=='0')
G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G1_noharv)
noharv[1]=coefficients(G1_noharv)[1]
noharv[2]=coefficients(G1_noharv)[2]
noharv[3]=sd(resid(G1_noharv))

#Sprouts
n2 = nrow(sprouts)
#Sprout survivorship
survivingsprouts=subset(sprouts, sprouts$survival=='1')
noharv[8]=length(survivingsprouts$survival)/n2
#Sprout survivorship
gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
noharv[9]=length(gototree$sprout2016)/length(survivingsprouts$survival)
#Sprout size distribution
gototree2=subset(gototree, gototree$size2016>0)
if (length(gototree2$size2016)<4){
  noharv[6]=0.5412373
  noharv[7]=0.1164737
}else{noharv[6]=mean(gototree2$size2016)
noharv[7]=sd(gototree2$size2016)
}

# Run IPM -----------------------------------------------------------------


run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_notH,
      m             = m_notH,
      data_list     = data_list,
      states        = list(c('size')),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    
    define_kernel(
      name          = "go_sprout",
      formula       = s * m* d_size,
      family        = 'CD',
      s             = s_notH,
      m             = m_notH,
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b),
      family  = "DD",
      data_list     = data_list, states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d *d_size,
      f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
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
  L <- log(1.1)*0.8
  U <- log(20)
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
  g_int     = noharv[1],
  g_slope   = noharv[2],
  g_sd      = noharv[3],
  s_notH    = noharv[4],
  m_notH    = noharv[5],
  d_mu    = noharv[6],
  d_sd    = noharv[7],
  # d_mu    = 1.8453,
  # d_sd    = 0.722706,
  a       = noharv[8],
  b       = noharv[9]
) 

ipm_temp2=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

noharv[10] = lambda(ipm_temp2)
