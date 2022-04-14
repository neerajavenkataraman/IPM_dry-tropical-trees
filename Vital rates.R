#Read data
dat <- read.csv("Demogdata.csv")
d <- dat[dat$State == "Tree",] # excluding sprouts in t


# Vital rates -------------------------------------------------------------


# Growth ------------------------------------------------------------------

# Acacia chundra
g.ac <- d[d$Species=="Acacia chundra" & d$survival == 1 & d$sprout2016 != 1,]
table(g.ac$harvest) # 8 were h and 42 noh
mod <- lm(log(size2016) ~ log(size2009)*as.factor(harvest), data = g.ac)
summary(mod)

# Chloroxylon swietenia
g.cs <- d[d$Species=="Chloroxylon swietenia" & d$survival == 1 & d$sprout2016 != 1,]
mod <- lm(log(size2016) ~ log(size2009)*as.factor(harvest), data = g.cs)
summary(mod)

# Gardenia gummifera
g.gg <- d[d$Species=="Gardenia gummifera" & d$survival == 1 & d$sprout2016 != 1,]
mod <- lm(log(size2016) ~log (size2009)*as.factor(harvest), data = g.gg)
summary(mod)


# Survival ----------------------------------------------------------------


# Acacia chundra
s.ac <- d[d$Species == "Acacia chundra", ]
table(s.ac$harvest, s.ac$survival)
mod <- glm(survival ~ log(size2009) * harvest, data = s.ac,
         family = binomial(link = "logit"))
summary(mod)  

# Chloroxylon swietenia
s.cs <- d[d$Species == "Chloroxylon swietenia", ]
table(s.cs$harvest, s.cs$survival)
mod <- glm(survival ~ log(size2009) * harvest, data = s.cs,
         family = binomial(link = "logit"))
summary(mod) 

# Gardenia gummifera
s.gg <- d[d$Species == "Gardenia gummifera", ]
table(s.gg$harvest, s.gg$survival)
mod <- glm(survival ~ log(size2009) * harvest, data = s.gg,
         family = binomial(link = "logit"))
summary(mod) 


# Probability of transition to sprout -------------------------------------

# Acacia chundra
m.ac <- d[d$Species=="Acacia chundra" & d$survival == 1,]
table(m.ac$harvest, m.ac$sprout2016)
mod <- glm(sprout2016 ~ log(size2009) * harvest, data = m.ac,
         family = binomial(link = "logit"))
summary(mod) 

# Chloroxylon swietenia
m.cs <- d[d$Species == "Chloroxylon swietenia" & d$survival == 1, ]
table(m.cs$harvest, m.cs$sprout2016)
mod <- glm(sprout2016 ~ log(size2009) * harvest, data = m.cs,
         family = binomial(link = "logit"))
summary(mod) 

# Gardenia gummifera
m.gg <- d[d$Species == "Gardenia gummifera" & d$survival == 1,  ]
table(m.gg$harvest, m.gg$sprout2016)
mod <- glm(sprout2016 ~ log(size2009) * harvest, data = m.gg,
         family = binomial(link = "logit"))
summary(mod) 




