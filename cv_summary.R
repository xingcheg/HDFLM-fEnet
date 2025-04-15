## setwd(.....)

make_scad <- function(out, lambda, s){
  out1 <- out[(out$lambda==lambda) & (out$s==s) , c('FPR', 'FNR', 'MNE', 'RER')]
  return(round(apply(out1, 2, quantile, probs = c(0.5, 0.025, 0.975)), 4))
}

make_enet <- function(out, lambda, alpha, df, theta){
  out1 <- out[(out$lambda==lambda) & 
                ( abs(10^(out$alpha)-10^alpha) < 1e-10 ) & 
                (out$df==df) & 
                (abs(10^(out$theta)-10^theta) < 1e-7), c('FPR', 'FNR', 'MNE', 'RER')]
  return(round(apply(out1, 2, quantile, probs = c(0.5, 0.025, 0.975)), 4))
}


#####    Scenario I      #####
#####    rho = 0.75      #####
##### N=100, P=200, Q=10 #####
cv_scad2_auto <- readRDS("cv_scad2_auto.rds")
out_scad2_auto_mean <- round(aggregate(.~lambda+s, data=cv_scad2_auto$output[,-1], mean), 5)
RER_scad2_auto_idx <- order(out_scad2_auto_mean$RER)
out_scad2_auto_mean[RER_scad2_auto_idx[1:10],]
make_scad(cv_scad2_auto$output, 1.0, 3)



cv_enet2_auto <- readRDS("cv_enet2_auto.rds")
out_enet2_auto_mean <- round(aggregate(.~lambda+alpha+theta+df, data=cv_enet2_auto$output[,-1], mean), 5)
RER_enet2_auto_idx <- order(out_enet2_auto_mean$RER)
out_enet2_auto_mean[RER_enet2_auto_idx[1:15],]
make_enet(cv_enet2_auto$output, 0.05, -Inf, 10, -5)


