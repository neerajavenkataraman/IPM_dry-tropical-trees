
# 'Ambient' (harvested and not harvested individuals) ---------------------------------------------------------------------

gumi <- matrix(NA, ncol=34) #Dataframe for the parameters you get during bootstrapping +


colnames(gumi) <- c("g_int", "g_slope", "g_sd", "s_int", "s_slope", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda_orgs",  "S_g_int", "S_g_slope", "S_g_sd", "S_s_int", "S_s_slope", "S_m_int", "S_m_slope",  "S_d_mean", "S_d_sd","S_a", "S_b", c(1:11))   ### Adding names to columns

#trees
n =         nrow(gumitrees)
#survival trees
S_all <-   glm(survival ~ log_size_t0, data = gumitrees, family="binomial")
summary(S_all)
gumi[4]=coef(S_all)[1]
gumi[5]=coef(S_all)[2]
#plot survival
# ggplot(gumitrees, aes(x=log_size_t0, y=survival)) + geom_point() + geom_smooth(method = "glm",  method.args = list(family = "binomial"),  se = FALSE)+theme_classic() + ggtitle("G.gummifera")

#probability of tree to sprout
M_all <- glm(sprout2016 ~ log_size_t0, data= gumitrees, family="binomial")
summary(M_all)
gumi[6] =coef(M_all)[1]
gumi[7]=coef(M_all)[2]
#plot tree to sprout
#ggplot(gumitrees, aes(x=log_size_t0, y=sprout2016)) + geom_point() + geom_smooth(method = "glm",  method.args = list(family = "binomial"),  se = FALSE)+theme_classic()

#growth
growers=subset(gumitrees, gumitrees$survival=='1' & gumitrees$sprout2016=='0')
G_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
summary(G_all)
gumi[1]=coefficients(G_all)[1]
gumi[2]=coefficients(G_all)[2]
gumi[3]=sd(resid(G_all))

#sprouts
n2 = nrow(gumisprouts)
#survival sprouts
survivingsprouts=subset(gumisprouts, gumisprouts$survival=='1') 
gumi[10]=length(survivingsprouts$survival)/n2
#sprout to tree probability
gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
gumi[11]=length(gototree$sprout2016)/length(survivingsprouts$survival)
#size distribution sprout to tree
gototree2=subset(gototree, gototree$size2016>0)
if (length(gototree2$size2016)<10) {
  gumi[8]=1.43193
  gumi[9]=0.2112621
}else
{gumi[8]=mean(gototree2$size2016)
gumi[9]=sd(gototree2$size2016)
}

run_ipm <- function(data_list) {
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = plogis(s_int + s_slope  * size_1),
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
      formula       = s * m * d_size,  #you need d_size because you want to know what sizes go to discrete
      family        = 'CD',
      s             = plogis(s_int + s_slope  * size_1),
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b),
      family  = "DD",
      data_list     = data_list, 
      states  = list(c( "sprout")),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d * d_size,
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
  L <- log(0.6)*1.2 #logsize is negative, if not multiply by 0.8
  U <- log(25)
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
  g_int     = gumi[1],
  g_slope   = gumi[2],
  g_sd      = gumi[3],
  s_int    = gumi[4],
  s_slope  =  gumi[5],
  m_int   =  gumi[6],
  m_slope  = gumi[7],
  d_mu    = gumi[8],
  d_sd    = gumi[9],
  # d_mu    = 1.43193,
  # d_sd    = 0.2112621,
  a       = gumi[10],
  b       = gumi[11]
) 

ipm_temp=run_ipm(pars_orgs) 

gumi[12] = lambda(ipm_temp)

### Run the ipm and calculate the lambda straight away

for(p in 1:11){        # Now for each sample add 0.001 to each parameter, one by one
  
  
  pars_pos <- pars_orgs  ## Copy the original parameters into a manual one
  pars_pos[[p]] <- pars_pos[[p]] + 0.001 # p is the specific parameter
  
  
  #with this "new" list rerun the ipm to get the "new" lambda
  ipm_pos   <- run_ipm(pars_pos)
  lambda_pos <- lambda(ipm_pos)
  
  # save the sensitivity to the right paramter
  sens <- ((lambda_pos - gumi[12]) / 0.001)
  gumi[(12+p)] <- sens          #all the new sensitivities
  gumi[(23+p)] <- lambda_pos   #all the new lambdas
}



# Bootstrap resampling 'Ambient' -----------------------------------------------------------------
B = 1000
gumiall <- matrix(NA, nrow=B, ncol=34) #Dataframe for the parameters you get during bootstrapping +


colnames(gumiall) <- c("g_int", "g_slope", "g_sd", "s_int", "s_slope", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda_orgs",  "S_g_int", "S_g_slope", "S_g_sd", "S_s_int", "S_s_slope", "S_m_int", "S_m_slope",  "S_d_mean", "S_d_sd","S_a", "S_b", c(1:11))   ### Adding names to columns


for (i in 1:B){
  
  set.seed(i)
  #trees
  n =         nrow(gumitrees)
  x3 <-       sample(1:n, n, replace = TRUE)
  databoot <- gumitrees[x3,]
  #survival trees
  S_all <-   glm(survival ~ log_size_t0, data = databoot, family="binomial")
  gumiall[i,4]=coef(S_all)[1]
  gumiall[i,5]=coef(S_all)[2]
  
  #Tree to sprout
  M_all <- glm(sprout2016 ~ log_size_t0, data= databoot, family="binomial")
  gumiall[i,6] =coef(M_all)[1]
  gumiall[i,7]=coef(M_all)[2]
  
  #growth
  growers=subset(databoot, databoot$survival=='1' & databoot$sprout2016=='0')
  G_all  <- lm(log_size_t1 ~ log_size_t0, data = growers)
  gumiall[i,1]=coefficients(G_all)[1]
  gumiall[i,2]=coefficients(G_all)[2]
  gumiall[i,3]=sd(resid(G_all))
  
  #sprouts
  n2 = nrow(gumisprouts)
  x2 <- sample(1:n2, n2, replace = TRUE)
  databoot2 <- gumisprouts[x2,]
  #survival sprouts
  survivingsprouts=subset(databoot2, databoot2$survival=='1') 
  gumiall[i,10]=length(survivingsprouts$survival)/n2
  #sprout to tree
  gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
  gumiall[i,11]=length(gototree$sprout2016)/length(survivingsprouts$survival)
  #size distribution sprout to tree
  gototree2=subset(gototree, gototree$size2016>0)
  if (length(gototree2$size2016)<10) {
    gumiall[i,8]=1.43193
    gumiall[i,9]=0.2112621
  }else
  {gumiall[i,8]=mean(gototree2$size2016)
  gumiall[i,9]=sd(gototree2$size2016)
  }
  
  run_ipm <- function(data_list) {
    general_ipm <- init_ipm("general", "di", "det") %>%
      define_kernel(
        name          = "P",
        formula       = s * g * (1-m) * d_size,
        family        = "CC",
        g             = dnorm(size_2, g_mu, g_sd),
        g_mu          = g_int + g_slope * size_1,
        s             = plogis(s_int + s_slope  * size_1),
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
        s             = plogis(s_int + s_slope  * size_1),
        m             = plogis(m_int + m_slope * size_1),
        data_list     = data_list,
        states        = list(c('size', "sprout")),
        has_hier_effs = FALSE
      ) %>%
      
      
      define_kernel(
        name    = 'stay_sprout',
        formula = a * (1-b),
        family  = "DD",
        data_list     = data_list, 
        states  = list(c( "sprout")),
        evict_cor = FALSE
      ) %>%
      
      define_kernel(
        name          = 'leave_sprout',
        formula       = a * b * f_d * d_size,
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
    L <- log(0.6)*1.2 #logsize is negative, if not multiply by 0.8
    U <- log(25)
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
    g_int     = gumiall[i,1],
    g_slope   = gumiall[i,2],
    g_sd      = gumiall[i,3],
    s_int    = gumiall[i,4],
    s_slope  =  gumiall[i,5],
    m_int   =  gumiall[i,6],
    m_slope  = gumiall[i,7],
    d_mu    = gumiall[i,8],
    d_sd    = gumiall[i,9],
    # d_mu    = 1.43193,
    # d_sd    = 0.2112621,
    a       = gumiall[i,10],
    b       = gumiall[i,11]
  ) 
  
  ipm_temp2=run_ipm(pars_orgs) 
  
  gumiall[i,12] = lambda(ipm_temp2)
  
  ### Run the ipm and calculate the lambda straight away
  
  for(p in 1:11){        # Now for each sample add 0.001 to each parameter, one by one
    
    
    pars_pos <- pars_orgs  ## Copy the original parameters into a manual one
    pars_pos[[p]] <- pars_pos[[p]] + 0.001 # p is the specific parameter
    
    
    #with this "new" list rerun the ipm to get the "new" lambda
    ipm_pos   <- run_ipm(pars_pos)
    lambda_pos <- lambda(ipm_pos)
    #lambda_pos <- run_ipm(pars_pos)
    
    
    # save the sensitivity to the right paramter
    sens <- ((lambda_pos - gumiall[i,12]) / 0.001)
    gumiall[i,(12+p)] <- sens          #all the new sensitivities
    gumiall[i, (23+p)] <- lambda_pos   #all the new lambdas
  }
  
  
}


# No harvesting -----------------------------------------------------------

guminoharvest <- matrix(NA, ncol=12) #Dataframe for the parameters you get during bootstrapping +


colnames(guminoharvest) <- c("g_int", "g_slope", "g_sd", "s_int", "s_slope", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda")   ### Adding names to columns



#trees
n =         nrow(no_harvest)

#Survival trees
S2_noharv<-   glm(survival ~ log_size_t0, data = no_harvest, family="binomial")
summary(S2_noharv)
guminoharvest[4]=coef(S2_noharv)[1]
guminoharvest[5]=coef(S2_noharv)[2]
#probability of tree to sprout
M2_noharv <- glm(sprout2016 ~ log_size_t0, data = no_harvest, family="binomial") 
summary(M2_noharv)
guminoharvest[6] =coef(M2_noharv)[1]
guminoharvest[7]=coef(M2_noharv)[2]
growersnoharv=subset(no_harvest, no_harvest$survival=='1' & no_harvest$sprout2016=='0')
#Growth trees
G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growersnoharv)
summary(G1_noharv)
guminoharvest[1]=coefficients(G1_noharv)[1]
guminoharvest[2]=coefficients(G1_noharv)[2]
guminoharvest[3]=sd(resid(G1_noharv))

#sprouts
#sprouts are not harvested so they remain same for both 'All' and 'Not harvested'
n2 = nrow(gumisprouts)
survivingsprouts=subset(gumisprouts, gumisprouts$survival=='1') 
guminoharvest[10]=length(survivingsprouts$survival)/n2  #sprout survival
gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0') 
#probability of sprout to tree
guminoharvest[11]=length(gototree$sprout2016)/length(survivingsprouts$survival)
#size distribution of sprouts
gototree2=subset(gototree, gototree$size2016>0)
if (length(gototree2$size2016)<10) {
  guminoharvest[8]=1.43193
  guminoharvest[9]=0.2112621
}else
{guminoharvest[8]=mean(gototree2$size2016) 
guminoharvest[9]=sd(gototree2$size2016)
}

### Create a function from the ipm so it is easier to run multiple times for the sensitivities
run_ipm <- function(data_list) {  
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = plogis(s_int + s_slope  * size_1), 
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
      s             = plogis(s_int + s_slope  * size_1),
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b),
      family  = "DD",
      data_list     = data_list,
      states  = list(c('sprout')),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d * d_size,
      f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  
  # The lower and upper bounds for the continuous state variable and the number
  # of meshpoints for the midpoint rule integration. We'll also create the initial
  # population vector from a random uniform distribution
  L <- log(0.6)*1.2 #multiply by 1.2 because logsize is negative
  U <- log(25)
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
  
  return(ipm)
}

pars_orgs <- list(                           ##Original dataset
  g_int     = guminoharvest[1],
  g_slope   = guminoharvest[2],
  g_sd      = guminoharvest[3],
  s_int    = guminoharvest[4],
  s_slope  =  guminoharvest[5],
  m_int   =  guminoharvest[6],
  m_slope  = guminoharvest[7],
  d_mu    = guminoharvest[8],
  d_sd    = guminoharvest[9],
  # d_mu    = 1.43193,
  # d_sd    = 0.2112621,
  a       = guminoharvest[10],
  b       = guminoharvest[11]
) 

ipm_temp3=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

guminoharvest[12] = lambda(ipm_temp3)  
)

# 'Not harvested' --------------------------------------------------------------
guminh <- matrix(NA, ncol=12) #Dataframe for the parameters you get during bootstrapping +


colnames(guminh) <- c("g_int", "g_slope", "g_sd", "s_int", "s_slope", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda")   ### Adding names to columns

#trees
n =         nrow(no_harvest)
#survival trees
S2_noharv<-   glm(survival ~ log_size_t0, data = no_harvest, family="binomial")
summary(S2_noharv)
guminh[4]=coef(S2_noharv)[1]
guminh[5]=coef(S2_noharv)[2]

#gumiall[i,4]=length(survivors$survival)/n
#gotosprout=subset(databoot, databoot$sprout2016=='1')

#Tree to sprout probability
M2_noharv <- glm(sprout2016 ~ log_size_t0, data = no_harvest, family="binomial") 
guminh[6] =coef(M2_noharv)[1]
guminh[7]=coef(M2_noharv)[2]

#Growth trees
growersnoharv=subset(no_harvest, no_harvest$survival=='1' & no_harvest$sprout2016=='0')
G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growersnoharv)
guminh[1]=coefficients(G1_noharv)[1]
guminh[2]=coefficients(G1_noharv)[2]
guminh[3]=sd(resid(G1_noharv))

#sprouts
n2 = nrow(gumisprouts)
survivingsprouts=subset(gumisprouts, gumisprouts$survival=='1')
guminh[10]=length(survivingsprouts$survival)/n2  #sprout survival
gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0')
guminh[11]=length(gototree$sprout2016)/length(survivingsprouts$survival) #sprout to tree
#sprout to tree size distribution
gototree2=subset(gototree, gototree$size2016>0)
if (length(gototree2$size2016)<10) {
  guminh[8]=1.43193
  guminh[9]=0.2112621
}else
{guminh[8]=mean(gototree2$size2016) 
guminh[9]=sd(gototree2$size2016)
}

### Create a function from the ipm so it is easier to run multiple times for the sensitivities
run_ipm <- function(data_list) {  
  general_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = plogis(s_int + s_slope  * size_1), #plogis logistic regression
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
      s             = plogis(s_int + s_slope  * size_1),
      m             = plogis(m_int + m_slope * size_1),
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE
    ) %>%
    
    
    
    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b),
      family  = "DD",
      data_list     = data_list,
      states  = list(c('sprout')),
      evict_cor = FALSE
    ) %>%
    
    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d * d_size,
      f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size', "sprout", "sprout", 'size')
      )
    )
  
  
  # The lower and upper bounds for the continuous state variable and the number
  # of meshpoints for the midpoint rule integration. We'll also create the initial
  # population vector from a random uniform distribution
  L <- log(0.6)*1.2 #because logsize is negative
  U <- log(25)
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
  
  return(ipm)
}

pars_orgs <- list(                           ##Original dataset
  g_int     = guminh[1],
  g_slope   = guminh[2],
  g_sd      = guminh[3],
  s_int    = guminh[4],
  s_slope  =  guminh[5],
  m_int   =  guminh[6],
  m_slope  = guminh[7],
  d_mu    = guminh[8],
  d_sd    = guminh[9],
  # d_mu    = 1.43193,
  # d_sd    = 0.2112621,
  a       = guminh[10],
  b       = guminh[11]
) 

ipm_temp3=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away

guminh[12] = lambda(ipm_temp3)


guminh

# Bootstrap resampling 'No harvesting' -----------------------------------------------------------

B = 1000
guminoharv <- matrix(NA, nrow=B, ncol=12) #Dataframe for the parameters you get during bootstrapping +


colnames(guminoharv) <- c("g_int", "g_slope", "g_sd", "s_int", "s_slope", "m_int", "m_slope", "d_mean", "d_sd",   "a", "b", "lambda")   ### Adding names to columns

for (i in 1:B){
  
  #trees
  n =         nrow(no_harvest)
  x3 <-       sample(1:n, n, replace = TRUE)
  databoot <- no_harvest[x3,]
  #Survival trees
  S2_noharv<-   glm(survival ~ log_size_t0, data = databoot, family="binomial")
  summary(S2_noharv)
  guminoharv[i,4]=coef(S2_noharv)[1]
  guminoharv[i,5]=coef(S2_noharv)[2]
  #probability of tree to sprout
  M2_noharv <- glm(sprout2016 ~ log_size_t0, data = databoot, family="binomial") 
  summary(M2_noharv)
  guminoharv[i,6] =coef(M2_noharv)[1]
  guminoharv[i,7]=coef(M2_noharv)[2]
  growersnoharv=subset(databoot, databoot$survival=='1' & databoot$sprout2016=='0')
  #Growth trees
  G1_noharv  <- lm(log_size_t1 ~ log_size_t0, data = growersnoharv)
  summary(G1_noharv)
  guminoharv[i,1]=coefficients(G1_noharv)[1]
  guminoharv[i,2]=coefficients(G1_noharv)[2]
  guminoharv[i,3]=sd(resid(G1_noharv))
  
  #sprouts
  #sprouts are not harvested so they remain same for both 'All' and 'Not harvested'
  n2 = nrow(gumisprouts)
  x2 <- sample(1:n2, n2, replace = TRUE)
  databoot2 <- gumisprouts[x2,]
  survivingsprouts=subset(databoot2, databoot2$survival=='1') 
  guminoharv[i,10]=length(survivingsprouts$survival)/n2  #sprout survival
  gototree=subset(survivingsprouts, survivingsprouts$sprout2016=='0') 
  #probability of sprout to tree
  guminoharv[i,11]=length(gototree$sprout2016)/length(survivingsprouts$survival)
  #size distribution of sprouts
  gototree2=subset(gototree, gototree$size2016>0)
  if (length(gototree2$size2016)<10) {
    guminoharv[i,8]=1.43193
    guminoharv[i,9]=0.2112621
  }else
  {guminoharv[i,8]=mean(gototree2$size2016) 
  guminoharv[i,9]=sd(gototree2$size2016)
  }
  
  ### Create a function from the ipm so it is easier to run multiple times for the sensitivities
  run_ipm <- function(data_list) {  
    general_ipm <- init_ipm("general", "di", "det") %>%
      define_kernel(
        name          = "P",
        formula       = s * g * (1-m) * d_size,
        family        = "CC",
        g             = dnorm(size_2, g_mu, g_sd),
        g_mu          = g_int + g_slope * size_1,
        s             = plogis(s_int + s_slope  * size_1), 
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
        s             = plogis(s_int + s_slope  * size_1),
        m             = plogis(m_int + m_slope * size_1),
        data_list     = data_list,
        states        = list(c('size', "sprout")),
        has_hier_effs = FALSE
      ) %>%
      
      
      
      define_kernel(
        name    = 'stay_sprout',
        formula = a * (1-b),
        family  = "DD",
        data_list     = data_list,
        states  = list(c('sprout')),
        evict_cor = FALSE
      ) %>%
      
      define_kernel(
        name          = 'leave_sprout',
        formula       = a * b * f_d * d_size,
        f_d           = dnorm(size_2, d_mu, d_sd),
        family        = 'DC',
        data_list     = data_list,
        states        = list(c('size', "sprout")),
        has_hier_effs = FALSE,
        evict_cor     = TRUE,
        evict_fun     = truncated_distributions('norm',
                                                'f_d')
      ) %>%
      
      define_impl(
        make_impl_args_list(
          kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
          int_rule     = c(rep("midpoint", 4)),
          state_start    = c('size', "size", "sprout", "sprout"),
          state_end      = c('size', "sprout", "sprout", 'size')
        )
      )
    
    
    # The lower and upper bounds for the continuous state variable and the number
    # of meshpoints for the midpoint rule integration. We'll also create the initial
    # population vector from a random uniform distribution
    L <- log(0.6)*1.2 #multiply by 1.2 because logsize is negative
    U <- log(25)
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
    
    return(ipm)
  }
  
  pars_orgs <- list(                           ##Original dataset
    g_int     = guminoharv[i,1],
    g_slope   = guminoharv[i,2],
    g_sd      = guminoharv[i,3],
    s_int    = guminoharv[i,4],
    s_slope  =  guminoharv[i,5],
    m_int   =  guminoharv[i,6],
    m_slope  = guminoharv[i,7],
    d_mu    = guminoharv[i,8],
    d_sd    = guminoharv[i,9],
    # d_mu    = 1.43193,
    # d_sd    = 0.2112621,
    a       = guminoharv[i,10],
    b       = guminoharv[i,11]
  ) 
  
  ipm_temp4=run_ipm(pars_orgs) ### Run the ipm and calculate the lambda straight away
  
  guminoharv[i,12] = lambda(ipm_temp4)  
}

