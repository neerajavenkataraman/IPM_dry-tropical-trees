#Acacia chundra
# Ambient harvesting scenario ---------------------------------------------

chundra <- matrix(NA, ncol=31) #Dataframe for the parameters 

colnames(chundra) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda_orgs",  "S_g_int", "S_g_slope", "S_g_sd", "S_s_all", "S_m_int", "S_m_slope",  "S_d_mean", "S_d_sd","S_a", "S_b", c(1:10))   ### Adding names to columns


#trees 
n           <- nrow(chundratrees)
#tree survivorship
# S_all <-   glm(survival ~ log_size_t0 , data = chundratrees, family="binomial")
# summary(S_all) #no significant effect of size on survivorship hence modeled as constant
survivors   <- subset(chundratrees, chundratrees$survival=='1')
chundra[4] <- length(survivors$survival)/n
#growth of trees
growers <- subset(chundratrees, chundratrees$survival=='1' & chundratrees$sprout2016=='0')
G_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G1_all)
chundra[1]=coefficients(G_all)[1]
chundra[2]=coefficients(G_all)[2]
chundra[3]=sd(resid(G1_all))
#transition tree to sprout
M_all <- glm(sprout2016 ~ log_size_t0, data = chundratrees,family="binomial") 
summary(M2_mod)
chundra[5] =coef(M_all)[1]
chundra[6]=coef(M_all)[2]

#sprouts
n2 = nrow(chundrasprouts)
#sprout survivorship
#none of the sprouts survive
survivingsprouts=subset(chundrasprouts, chundrasprouts$survival=='1') 
chundra[7]=length(survivingsprouts$survival)/n2
#transition sprout to tree
chundra[8]=0 
#sprout size distribution
chundra[9]=0
chundra[10]=0

# Run IPM -----------------------------------------------------------------
run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%  #general, density independent, deterministic
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


# Ambient harvesting Size structure projection ----------------------------
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

colnames(time_series) <- c("small1","small2", "medium1", "medium2", "big1")
#number of trees
time_series %>% as.data.frame

#number of sprouts: ipm_temp$pop_state$n_sprout
as.data.frame(t(ipm_temp$pop_state$n_sprout))

# No harvesting scenario --------------------------------------------------------------

chundranoharv <-matrix(NA, ncol=11)  #Dataframe for the parameters 

colnames(chundranoharv) <- c("g_int", "g_slope", "g_sd", "s_all", "m_int","m_slope", "d_mu", "d_sd", "a", "b", "lambda")   ### Adding names to columns  

n           <- nrow(no_harvest_chundra)
#tree survivorship
# S_noharv <-   glm(survival ~ log_size_t0 , data = no_harvest_chundra, family="binomial")
# summary(S_noharv) #no significant effect of size on survivorship hence modeled as constant
survivors   <- subset(no_harvest_chundra, no_harvest_chundra$survival=='1')
chundranoharv[4] <- length(survivors$survival)/n
#growth of trees
growers <- subset(no_harvest_chundra, no_harvest_chundra$survival=='1' & no_harvest_chundra$sprout2016=='0')
G_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growers)
chundranoharv[1]=coefficients(G_noharv)[1]
chundranoharv[2]=coefficients(G_noharv)[2]
chundranoharv[3]=sd(resid(G1_noharv))
#transition tree to sprout
M_noharv <- glm(sprout2016 ~ log_size_t0, data = no_harvest_chundra,family="binomial") 
chundranoharv[5] =coef(M_noharv)[1]
chundranoharv[6] =coef(M_noharv)[2]
chundranoharv[7] =0  #sprout survivorship
chundranoharv[8] =0  #transition sprout to tree
chundranoharv[9]=0   #sprout size distribution mean
chundranoharv[10]=0  #sprout size distribution sd

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
  a       = chundranoharv[9],
  b       = chundranoharv[10]
) 

ipm_temp3=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

chundranoharv[11] = lambda(ipm_temp3)
