
#-----------------------------------------------------------------------------------#
# Function for simulating half-synthesized data (strong effect size)

# Bin Yang
#-----------------------------------------------------------------------------------#

#################### functions for data generation ##################################

sim_data = readRDS("sim_data.rds")

generate_strong_data = function(n, seed){
  set.seed(seed)
  num.test = mvrnorm(n, mu = sim_data$mu, sigma = sim_data$sigma)
  anxious2 = rbinom(n, 1, sim_data$anxious2.mean) %>% as.character()
  melancholic = rbinom(n, 1, sim_data$melancholic.mean) %>% as.character()
  anger_attack2 = rbinom(n, 1, sim_data$anger_attack2.mean) %>% as.character()
  axis2_2 = rbinom(n, 1, sim_data$axis2.mean) %>% as.character()
  hypersomnia = rbinom(n, 1, sim_data$hypersomnia.mean) %>% as.character()
  atypical = rbinom(n, 1, sim_data$atypical.mean) %>% as.character()
  chronicity = rbinom(n, 1, sim_data$chronicity.mean) %>% as.character()
  severity1 = rbinom(n, 1, sim_data$severity1.mean) %>% as.character()
  sex = rbinom(n, 1, sim_data$sex.mean) %>% as.character()
  site_p = rmultinom(n, 1, c(1- sim_data$site.mg.prop - sim_data$site.tx.prop - sim_data$site.um.prop)/170, 
                     sim_data$site.mg.prop, 
                     sim_data$site.tx.prop, 
                     sim_data$site.um.prop) %>% t() %>% as_tibble()
  site = site_p %>% 
    mutate(site = case_when(V1 == 1 ~ "cu",
                            V2 == 1 ~ "mg",
                            V3 == 1 ~ "tx",
                            V4 == 1 ~ "um")) %>% 
    dplyr::select(site) %>% pull()
  char.test = bind_cols(site, anxious2, melancholic, anger_attack2, axis2_2, hypersomnia, atypical,
                        chronicity, severity1, sex)
  colnames(char.test) = var.char.names
  w_test = rbinom(n, 1, mean(W))
  y.temp = rbinom(n, 1, 0.5)
  X_test = bind_cols(as_tibble(num.test), char.test, y.temp)
  colnames(X_test)[265] = "Y"
  X_test = model.matrix(Y ~ ., X_test)[,-1] %>% as_tibble()
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  # treatment effect is simulated based on the tree structure learned from real data with effect sizes 
  # estimated with the average of the doubly robust scores
  
  p_rem = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(w_test*0.806 + (1-w_test)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(w_test*0.24 + (1-w_test)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(w_test*0.731 + (1-w_test)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(w_test*0.243 + (1-w_test)*0.685)
  Y = rbinom(n, 1, p_rem)
  final_data = list("X" = X_test, "Y" = Y, "W" = w_test,
                    "num.df" = as_tibble(num.test), "char.df" = char.test)
  return(final_data)
}


generate_weak_data = function(n, seed){
  set.seed(seed)
  num.test = mvrnorm(n, mu = sim_data$mu, sigma = sim_data$sigma)
  anxious2 = rbinom(n, 1, sim_data$anxious2.mean) %>% as.character()
  melancholic = rbinom(n, 1, sim_data$melancholic.mean) %>% as.character()
  anger_attack2 = rbinom(n, 1, sim_data$anger_attack2.mean) %>% as.character()
  axis2_2 = rbinom(n, 1, sim_data$axis2.mean) %>% as.character()
  hypersomnia = rbinom(n, 1, sim_data$hypersomnia.mean) %>% as.character()
  atypical = rbinom(n, 1, sim_data$atypical.mean) %>% as.character()
  chronicity = rbinom(n, 1, sim_data$chronicity.mean) %>% as.character()
  severity1 = rbinom(n, 1, sim_data$severity1.mean) %>% as.character()
  sex = rbinom(n, 1, sim_data$sex.mean) %>% as.character()
  site_p = rmultinom(n, 1, c(1- sim_data$site.mg.prop - sim_data$site.tx.prop - sim_data$site.um.prop)/170, 
                     sim_data$site.mg.prop, 
                     sim_data$site.tx.prop, 
                     sim_data$site.um.prop) %>% t() %>% as_tibble()
  site = site_p %>% 
    mutate(site = case_when(V1 == 1 ~ "cu",
                            V2 == 1 ~ "mg",
                            V3 == 1 ~ "tx",
                            V4 == 1 ~ "um")) %>% 
    dplyr::select(site) %>% pull()
  char.test = bind_cols(site, anxious2, melancholic, anger_attack2, axis2_2, hypersomnia, atypical,
                        chronicity, severity1, sex)
  colnames(char.test) = var.char.names
  w_test = rbinom(n, 1, mean(W))
  y.temp = rbinom(n, 1, 0.5)
  X_test = bind_cols(as_tibble(num.test), char.test, y.temp)
  colnames(X_test)[265] = "Y"
  X_test = model.matrix(Y ~ ., X_test)[,-1] %>% as_tibble()
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  # treatment effect is simulated based on the tree structure learned from real data with effect sizes 
  # estimated with the average of the doubly robust scores
  
  p_rem = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(w_test*0.706 + (1-w_test)*0.172) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(w_test*0.34 + (1-w_test)*0.73) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(w_test*0.631 + (1-w_test)*0.382) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(w_test*0.343 + (1-w_test)*0.585)
  Y = rbinom(n, 1, p_rem)
  final_data = list("X" = X_test, "Y" = Y, "W" = w_test,
                    "num.df" = as_tibble(num.test), "char.df" = char.test)
  return(final_data)
}



###################### setting 1: n = 200, strong effect size ###################

A.owl.1 = vector()
acc.owl.1 = vector()
A.tr.1 = vector()
acc.tr.1 = vector()
A.lasso.1 = vector()
acc.lasso.1 = vector()

for (i in 1:100){
  train_data = generate_strong_data(n = 200, seed = i)
  
  ####### train policy tree
  X_train = train_data$X
  Y_train = train_data$Y
  W_train = train_data$W
  
  owl_df = bind_cols(train_data$num.df, train_data$char.df, ifelse(W_train == 1, 1, -1), Y_train)
  colnames(owl_df)[265] = "W"
  colnames(owl_df)[266] = "Y"
  X_owl = model.matrix(Y ~ W*. , owl_df)[,-1]
  owl.fit = owl(H = as.matrix(scale(X_owl)), AA = ifelse(W_train == 1, 1, -1), RR = Y_train,  
                n = nrow(X_owl), K = 1, loss = "logit.lasso")
  
  
  test_data = generate_strong_data(n = 50000, seed = 123)
  
  X_test = test_data$X
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  ## optimal treatment 
  ########### optimal treatment ##############
  W.0 = 0
  p_rem.0 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.0*0.806 + (1-W.0)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.0*0.24 + (1-W.0)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.0*0.731 + (1-W.0)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.0*0.243 + (1-W.0)*0.685)
  W.1 = 1 
  p_rem.1 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.1*0.806 + (1-W.1)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.1*0.24 + (1-W.1)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.1*0.731 + (1-W.1)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.1*0.243 + (1-W.1)*0.685)
  set.seed(123)
  Y_1 = rbinom(50000, 1, p_rem.1)
  Y_0 = rbinom(50000, 1, p_rem.0)
  
  owl_test_df = bind_cols(test_data$num.df, test_data$char.df, ifelse(test_data$W == 1, 1, -1), test_data$Y)
  colnames(owl_test_df)[265] = "W"
  colnames(owl_test_df)[266] = "Y"
  X_test_owl = model.matrix(Y ~ W*. , owl_test_df)[,-1]
  
  
  pred = predict(owl.fit, H = list(as.matrix(scale(X_test_owl))),n = nrow(X_test_owl), K=1)
  A = pred$treatment[[1]]@x == 1
  acc.owl.1[i] = sum((pi.optimal - A) == 0)/50000
  A.owl.1[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# policy tree ##############
  pi = predict(tr, X_test) - 1
  acc.tr.1[i] = sum((pi.optimal - pi) == 0)/50000
  A = pi == 1
  A.tr.1[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  ############# LASSO ####################
  W.test.1 = 1
  W.test.0 = 0
  y.temp = rbinom(50000, 1, 0.5)
  test.df.1 = bind_cols(test_data$num.df, test_data$char.df, W.test.1, y.temp)
  colnames(test.df.1)[265] = "W"
  colnames(test.df.1)[266] = "Y"
  test.df.1 = model.matrix(Y ~ W*. , test.df.1)[,-1]
  test.df.0 = bind_cols(test_data$num.df, test_data$char.df, W.test.0, y.temp)
  colnames(test.df.0)[265] = "W"
  colnames(test.df.0)[266] = "Y"
  test.df.0 = model.matrix(Y ~ W*. , test.df.0)[,-1]
  
  pred.1 = predict(glm.fit, newdata = test.df.1, type = "prob")
  pred.0 = predict(glm.fit, newdata = test.df.0, type = "prob")
  
  pred.el = ifelse(pred.1$Yes >= pred.0$Yes, 1, 0)
  A = pred.el == 1
  acc.lasso.1[i] = sum((pi.optimal - pred.el) == 0)/50000
  A.lasso.1[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A) ##### value of Lasso ######
  
}



###################### setting 2: n = 200, weak effect size ###################

A.owl.2 = vector()
acc.owl.2 = vector()
A.tr.2 = vector()
acc.tr.2 = vector()
A.lasso.2 = vector()
acc.lasso.2 = vector()

for (i in 1:100){
  train_data = generate_weak_data(n = 200, seed = i)
  
  ####### train policy tree
  X_train = train_data$X
  Y_train = train_data$Y
  W_train = train_data$W
  
  owl_df = bind_cols(train_data$num.df, train_data$char.df, ifelse(W_train == 1, 1, -1), Y_train)
  colnames(owl_df)[265] = "W"
  colnames(owl_df)[266] = "Y"
  X_owl = model.matrix(Y ~ W*. , owl_df)[,-1]
  owl.fit = owl(H = as.matrix(scale(X_owl)), AA = ifelse(W_train == 1, 1, -1), RR = Y_train,  
                n = nrow(X_owl), K = 1, loss = "logit.lasso")
  
  
  test_data = generate_weak_data(n = 50000, seed = 123)
  
  X_test = test_data$X
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  ## optimal treatment 
  ########### optimal treatment ##############
  W.0 = 0
  p_rem.0 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.0*0.806 + (1-W.0)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.0*0.24 + (1-W.0)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.0*0.731 + (1-W.0)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.0*0.243 + (1-W.0)*0.685)
  W.1 = 1 
  p_rem.1 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.1*0.806 + (1-W.1)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.1*0.24 + (1-W.1)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.1*0.731 + (1-W.1)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.1*0.243 + (1-W.1)*0.685)
  set.seed(123)
  Y_1 = rbinom(50000, 1, p_rem.1)
  Y_0 = rbinom(50000, 1, p_rem.0)
  
  owl_test_df = bind_cols(test_data$num.df, test_data$char.df, ifelse(test_data$W == 1, 1, -1), test_data$Y)
  colnames(owl_test_df)[265] = "W"
  colnames(owl_test_df)[266] = "Y"
  X_test_owl = model.matrix(Y ~ W*. , owl_test_df)[,-1]
  
  
  pred = predict(owl.fit, H = list(as.matrix(scale(X_test_owl))),n = nrow(X_test_owl), K=1)
  A = pred$treatment[[1]]@x == 1
  acc.owl.2[i] = sum((pi.optimal - A) == 0)/50000
  A.owl.2[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# policy tree ##############
  pi = predict(tr, X_test) - 1
  acc.tr.1[i] = sum((pi.optimal - pi) == 0)/50000
  A = pi == 1
  A.tr.2[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# LASSO ####################
  W.test.1 = 1
  W.test.0 = 0
  y.temp = rbinom(50000, 1, 0.5)
  test.df.1 = bind_cols(test_data$num.df, test_data$char.df, W.test.1, y.temp)
  colnames(test.df.1)[265] = "W"
  colnames(test.df.1)[266] = "Y"
  test.df.1 = model.matrix(Y ~ W*. , test.df.1)[,-1]
  test.df.0 = bind_cols(test_data$num.df, test_data$char.df, W.test.0, y.temp)
  colnames(test.df.0)[265] = "W"
  colnames(test.df.0)[266] = "Y"
  test.df.0 = model.matrix(Y ~ W*. , test.df.0)[,-1]
  
  pred.1 = predict(glm.fit, newdata = test.df.1, type = "prob")
  pred.0 = predict(glm.fit, newdata = test.df.0, type = "prob")
  
  pred.el = ifelse(pred.1$Yes >= pred.0$Yes, 1, 0)
  A = pred.el == 1
  acc.lasso.2[i] = sum((pi.optimal - pred.el) == 0)/50000
  A.lasso.2[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A) ##### value of Lasso ######
  
}
###################### setting 3: n = 500, strong effect size ###################

A.owl.3 = vector()
acc.owl.3 = vector()
A.tr.3 = vector()
acc.tr.3 = vector()
A.lasso.3 = vector()
acc.lasso.3 = vector()

for (i in 1:100){
  train_data = generate_strong_data(n = 500, seed = i)
  
  ####### train policy tree
  X_train = train_data$X
  Y_train = train_data$Y
  W_train = train_data$W
  
  owl_df = bind_cols(train_data$num.df, train_data$char.df, ifelse(W_train == 1, 1, -1), Y_train)
  colnames(owl_df)[265] = "W"
  colnames(owl_df)[266] = "Y"
  X_owl = model.matrix(Y ~ W*. , owl_df)[,-1]
  owl.fit = owl(H = as.matrix(scale(X_owl)), AA = ifelse(W_train == 1, 1, -1), RR = Y_train,  
                n = nrow(X_owl), K = 1, loss = "logit.lasso")
  
  
  test_data = generate_strong_data(n = 50000, seed = 123)
  
  X_test = test_data$X
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  ## optimal treatment 
  ########### optimal treatment ##############
  W.0 = 0
  p_rem.0 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.0*0.806 + (1-W.0)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.0*0.24 + (1-W.0)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.0*0.731 + (1-W.0)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.0*0.243 + (1-W.0)*0.685)
  W.1 = 1 
  p_rem.1 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.1*0.806 + (1-W.1)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.1*0.24 + (1-W.1)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.1*0.731 + (1-W.1)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.1*0.243 + (1-W.1)*0.685)
  set.seed(123)
  Y_1 = rbinom(50000, 1, p_rem.1)
  Y_0 = rbinom(50000, 1, p_rem.0)
  
  owl_test_df = bind_cols(test_data$num.df, test_data$char.df, ifelse(test_data$W == 1, 1, -1), test_data$Y)
  colnames(owl_test_df)[265] = "W"
  colnames(owl_test_df)[266] = "Y"
  X_test_owl = model.matrix(Y ~ W*. , owl_test_df)[,-1]
  
  
  pred = predict(owl.fit, H = list(as.matrix(scale(X_test_owl))),n = nrow(X_test_owl), K=1)
  A = pred$treatment[[1]]@x == 1
  acc.owl.3[i] = sum((pi.optimal - A) == 0)/50000
  A.owl.3[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# policy tree ##############
  pi = predict(tr, X_test) - 1
  acc.tr.1[i] = sum((pi.optimal - pi) == 0)/50000
  A = pi == 1
  A.tr.3[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# LASSO ####################
  W.test.1 = 1
  W.test.0 = 0
  y.temp = rbinom(50000, 1, 0.5)
  test.df.1 = bind_cols(test_data$num.df, test_data$char.df, W.test.1, y.temp)
  colnames(test.df.1)[265] = "W"
  colnames(test.df.1)[266] = "Y"
  test.df.1 = model.matrix(Y ~ W*. , test.df.1)[,-1]
  test.df.0 = bind_cols(test_data$num.df, test_data$char.df, W.test.0, y.temp)
  colnames(test.df.0)[265] = "W"
  colnames(test.df.0)[266] = "Y"
  test.df.0 = model.matrix(Y ~ W*. , test.df.0)[,-1]
  
  pred.1 = predict(glm.fit, newdata = test.df.1, type = "prob")
  pred.0 = predict(glm.fit, newdata = test.df.0, type = "prob")
  
  pred.el = ifelse(pred.1$Yes >= pred.0$Yes, 1, 0)
  A = pred.el == 1
  acc.lasso.3[i] = sum((pi.optimal - pred.el) == 0)/50000
  A.lasso.3[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A) ##### value of Lasso ######
  
}



###################### setting 4: n = 500, weak effect size ###################

A.owl.4 = vector()
acc.owl.4 = vector()
A.tr.4 = vector()
acc.tr.4 = vector()
A.lasso.4 = vector()
acc.lasso.4 = vector()

for (i in 1:100){
  train_data = generate_weak_data(n = 500, seed = i)
  
  ####### train policy tree
  X_train = train_data$X
  Y_train = train_data$Y
  W_train = train_data$W
  
  owl_df = bind_cols(train_data$num.df, train_data$char.df, ifelse(W_train == 1, 1, -1), Y_train)
  colnames(owl_df)[265] = "W"
  colnames(owl_df)[266] = "Y"
  X_owl = model.matrix(Y ~ W*. , owl_df)[,-1]
  owl.fit = owl(H = as.matrix(scale(X_owl)), AA = ifelse(W_train == 1, 1, -1), RR = Y_train,  
                n = nrow(X_owl), K = 1, loss = "logit.lasso")
  
  
  test_data = generate_weak_data(n = 50000, seed = 123)
  
  X_test = test_data$X
  x1 = X_test$f7.open.alpha
  x2 = X_test$c2.close.theta
  x3 = X_test$qids_eval_total
  
  ## optimal treatment 
  ########### optimal treatment ##############
  W.0 = 0
  p_rem.0 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.0*0.806 + (1-W.0)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.0*0.24 + (1-W.0)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.0*0.731 + (1-W.0)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.0*0.243 + (1-W.0)*0.685)
  W.1 = 1 
  p_rem.1 = ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 <= 0.18, 1, 0)*(W.1*0.806 + (1-W.1)*0.072) + 
    ifelse(x1 <= 0.13, 1, 0)*ifelse(x2 > 0.18, 1, 0)*(W.1*0.24 + (1-W.1)*0.83) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 <= 18, 1, 0)*(W.1*0.731 + (1-W.1)*0.282) +
    ifelse(x1 > 0.13, 1, 0)*ifelse(x3 > 18, 1, 0)*(W.1*0.243 + (1-W.1)*0.685)
  set.seed(123)
  Y_1 = rbinom(50000, 1, p_rem.1)
  Y_0 = rbinom(50000, 1, p_rem.0)
  
  owl_test_df = bind_cols(test_data$num.df, test_data$char.df, ifelse(test_data$W == 1, 1, -1), test_data$Y)
  colnames(owl_test_df)[265] = "W"
  colnames(owl_test_df)[266] = "Y"
  X_test_owl = model.matrix(Y ~ W*. , owl_test_df)[,-1]
  
  
  pred = predict(owl.fit, H = list(as.matrix(scale(X_test_owl))),n = nrow(X_test_owl), K=1)
  A = pred$treatment[[1]]@x == 1
  acc.owl.4[i] = sum((pi.optimal - A) == 0)/50000
  A.owl.4[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# policy tree ##############
  pi = predict(tr, X_test) - 1
  acc.tr.1[i] = sum((pi.optimal - pi) == 0)/50000
  A = pi == 1
  A.tr.4[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A)
  
  ############# LASSO ####################
  W.test.1 = 1
  W.test.0 = 0
  y.temp = rbinom(50000, 1, 0.5)
  test.df.1 = bind_cols(test_data$num.df, test_data$char.df, W.test.1, y.temp)
  colnames(test.df.1)[265] = "W"
  colnames(test.df.1)[266] = "Y"
  test.df.1 = model.matrix(Y ~ W*. , test.df.1)[,-1]
  test.df.0 = bind_cols(test_data$num.df, test_data$char.df, W.test.0, y.temp)
  colnames(test.df.0)[265] = "W"
  colnames(test.df.0)[266] = "Y"
  test.df.0 = model.matrix(Y ~ W*. , test.df.0)[,-1]
  
  pred.1 = predict(glm.fit, newdata = test.df.1, type = "prob")
  pred.0 = predict(glm.fit, newdata = test.df.0, type = "prob")
  
  pred.el = ifelse(pred.1$Yes >= pred.0$Yes, 1, 0)
  A = pred.el == 1
  acc.lasso.4[i] = sum((pi.optimal - pred.el) == 0)/50000
  A.lasso.4[i] = mean(Y_1[A]) * mean(A) + mean(Y_0[!A]) * mean(!A) ##### value of Lasso ######
  
}