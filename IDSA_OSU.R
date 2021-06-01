library(plotly)
library(processx)
set.seed(123)

#####Initialization
dir.create("~/IDSA_OSU")
print("This script simulates an epidemic and fits the closed IDSA model to the simulated epidemic.")

osu.test.mat = read.csv("OSU_synthetic.csv")
osu.counts = read.csv("OSU_total_tests.csv")
day.zero = as.Date("2020-08-14")
day.end = as.Date("2020-11-25")
test.mat = matrix(NA, nrow = dim(osu.test.mat)[1], ncol = 2)
test.mat[,2] = osu.test.mat[,3]
test.mat[,1] = osu.test.mat[,2]
n = 11950
n.days = 104


decensor = function(test, p.t, S.0)
{
  flip = 1
  right.end = test[2]
  if (is.na(test[1]))
  {
    left.end = 1
  } else
  {
    left.end = test[1]
    if (left.end >= n.days)
    {
      return(NA)
    }
  }
  if (is.na(test[2]) | right.end > n.days)
  {
    flip.prob = ((1 - sum(p.t[1:left.end])) - (1- sum(p.t)))/(1 - sum(p.t[1:left.end]))
    flip = rbinom(1, 1, flip.prob)
    right.end = n.days
  }
  if (flip == 0)
  {
    return(NA)
  } else
  {
    if (is.na(test[1]))
    {
      prob.0 = (1-S.0)/(1-S.0 + sum(p.t[left.end:right.end]))
      flip.0 = rbinom(1,1, prob.0)
      if (flip.0 == 1)
      {
        return(0)
      }
    }
    true.date = c(left.end:right.end) %*% rmultinom(1, 1, p.t[left.end:right.end]/sum(p.t[left.end:right.end]))
    # if (11 <= true.date[1,1] & true.date[1,1] <= 15)
    # {
    #   true.date[1,1] = trunc(runif(1,11,16))
    # } else if (20 <= true.date[1,1] & true.date[1,1] <= 23)
    # {
    #   true.date[1,1] = trunc(runif(1,20,24))
    # } else if (27 <= true.date[1,1] & true.date[1,1] <= 29)
    # {
    #   true.date[1,1] = trunc(runif(1,27,30))
    # }
    return(true.date[1,1])
  }
}

S.fun = function(inf.dates)
{
  S = rep(0, n.days)
  for (i in 1:n.days)
  {
    S[i] = 1 - length(which(inf.dates <= i))/n
  }
  return(S)
}

p.fun = function(S, S.0)
{
  p.t = rep(0, n.days)
  p.t[1] = max(1/n,S.0 - S[1])
  for (j in 2:n.days)
  {
    p.t[j] = S[j-1] - S[j]
  }
  return(p.t)
}

S.comp = function(beta, t, I, S, I.0, E.0, S.0, S.e.0)
{
  S.e = (1-(1-beta[1]/n)^I.0)*S.0*(1-phi.t[1])^(t-1)
  if (t > 1)
  {
    for (k in 1:(t-1))
    {
      S.e = S.e+(1-(1-beta[k]/n)^I[k])*S[k]*(1-phi.t[1])^(t-k-1)
    }
  }
  S.e = S.e + S[t]
  S.e = S.e
  return(S.e)
}

loglik.eps = function(eps, delta, I, S)
{
  loglik = 0
  loglik = loglik + log((1 - (1-beta[1]/n)^I.0)^delta[1]*((1 - beta[1]/n)^I.0)^(S[1]))
  for (i in 2:n.days)
  {
    loglik = loglik + log((1 - (1-beta[i]/n)^I[i-1])^delta[i]*((1 - beta[i]/n)^I.0)^(S[i]))
  }
  return(loglik)
}

# IDSAplot = function()
# {
#   print("Writing posterior I samples to posterior_I.csv...")
#   write.csv(I.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_I.csv")
#   print("Writing posterior S samples to posterior_S.csv...")
#   write.csv(S.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_S.csv")
#   print("Writing posterior delta samples to posterior_dela.csv...")
#   write.csv(delta.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_delta.csv")
#   print("Writing posterior epsilon samples to posterior_epsilon.csv...")
#   write.csv(eps.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_eps.csv")
#   print("Writing posterior beta samples to posterior_beta.csv...")
#   write.csv(beta.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_beta.csv")
#   print("Writing posterior gamma samples to posterior_gamma.csv...")
#   write.csv(gamma.samp[(.5*nreps+1):nreps,], file = "~/IDSAClosed/posterior_gamma.csv")
#   
#   I.mean = apply(I.samp[(.5*nreps+1):nreps,], 2, mean)
#   I.min = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
#   I.q1 = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .25)
#   I.med = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
#   I.q3 = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .75)
#   I.max = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)
#   
#   beta.mean = apply(beta.samp[(.5*nreps+1):nreps,], 2, mean)
#   beta.min = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
#   beta.q1 = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .25)
#   beta.med = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
#   beta.q3 = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .75)
#   beta.max = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)
#   
#   gamma.med = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
#   gamma.max = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)
#   gamma.min = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
#   
#   print("Generating plots...")
#   
#   t <- list(
#     family = "arial",
#     size = 22,
#     color = "black")
#   
#   t.2 <- list(
#     family = "arial",
#     size = 24,
#     color = "black")
#   
#   
#   fig.on <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
#   fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.min[1:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
#   fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.med[1:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
#   fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.max[1:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
#   fig.on <- fig.on  %>% add_trace(x = c(1:(n.days-4)), y = I.count[1:(n.days-4)], mode = "markers", name = 'True prevalence')
#   fig.on <- fig.on  %>% add_segments(x = 8, xend = 8, y = 0, yend = 450, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   fig.on <- fig.on  %>% add_segments(x = 15, xend = 15, y = 0, yend = 450, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   fig.on <- fig.on  %>% add_segments(x = 29, xend = 29, y = 0, yend = 450, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   
#   
#   fig.on <- fig.on %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
#                               yaxis = list(title = "Infected", range = c(0,400)), font = t,
#                               margin = list(b = 100))
#   fig.on <- fig.on %>% config(mathjax = 'cdn')
#   orca(fig.on, "IDSAClosed/Closed_fit.svg")
#   
#   
#   #####R_t
#   true.Rt = beta.t[2:60]/(1/(I.count/rem.count))*S.count/n
#   true.Rt.smooth = supsmu(c(2:60), true.Rt)
#   R.0.samp = beta.samp[(.5*nreps+1):nreps,]/gamma.samp[(.5*nreps+1):nreps,]*(S.t.samp[(.5*nreps+1):nreps,]/n)
#   R.t.min = apply(R.0.samp, 2, quantile, probs = .025)
#   R.t.q1 = apply(R.0.samp, 2, quantile, probs = .25)
#   R.t.med = apply(R.0.samp, 2, quantile, probs = .5)
#   R.t.q3 = apply(R.0.samp, 2, quantile, probs = .75)
#   R.t.max = apply(R.0.samp, 2, quantile, probs = .975)
#   rem.rate = supsmu(c(2:60),1/(I.count/pos.count))
#   
#   fig.rt <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
#   fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.min[2:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
#   fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.med[2:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
#   fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.max[2:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
#   fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = true.Rt.smooth$y[2:(n.days-4)], mode = "lines", name = 'True R_t')
#   fig.rt <- fig.rt  %>% add_segments(x = 8, xend = 8, y = 0, yend = 450, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   fig.rt <- fig.rt  %>% add_segments(x = 15, xend = 15, y = 0, yend = 450, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   fig.rt <- fig.rt  %>% add_segments(x = 29, xend = 29, y = 0, yend = 450, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
#   
#   fig.rt <- fig.rt %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
#                               yaxis = list(title = "R_t", range = c(0,8)), font = t,
#                               margin = list(b = 100))
#   fig.rt <- fig.rt %>% config(mathjax = 'cdn')
#   orca(fig.rt, "IDSAClosed/Closed_Rt.svg")
# }

########Model fitting Initialization
inf.prop = (1-length(which(is.na(test.mat[,2])))/n)
p.t = rep(inf.prop/n.days, n.days)
p.t.ex = p.t
p.t.unif = rep(1/n.days, n.days)
p.t.null = rep(inf.prop/n.days, n.days)
S.null = .99
I = rep(5, n.days)
I.new = rep(0, n.days)
delta = rep(0, n.days)
eps = rep(1, n.days)
eps.new = rep(0, n.days)
omega = rep(0, n.days)
gamma = rep(.166, n.days)
gamma.new = rep(.166, n.days)
beta = rep(.4, n.days)
phi = rep(.25, n.days)
S = rep(1, n.days)
S.e = rep(1, n.days)
S.t = rep(n, n.days)
S.prod = rep(0, n.days)
e.t = rep(0, n.days)
nreps = 10000
a = c(1,1, rep(10,n.days-2))
b = c(9,9, rep(50, n.days-2))

S.samp = matrix(0, nrow = nreps, ncol = n.days)
S.t.samp = matrix(0, nrow = nreps, ncol = n.days)
beta.samp = matrix(0, nrow = nreps, ncol = n.days)
gamma.samp = matrix(0, nrow = nreps, ncol = n.days)
delta.samp = matrix(0, nrow = nreps, ncol = n.days)
eps.samp = matrix(0, nrow = nreps, ncol = n.days)
I.samp = matrix(0, nrow = nreps, ncol = n.days)
e.samp = matrix(0, nrow = nreps, ncol = n.days)

inf.dates = as.vector(apply(test.mat[1:n,], 1, decensor, p.t=p.t, S.0 = S.null))
ex.dates = inf.dates - rgeom(length(inf.dates), phi.t[1])-1
ex.geom = inf.dates - ex.dates
accept.eps = 0

###Gibbs sampler
print("Running MCMC...")
for (k in 1:nreps)
{
  
  X.unif = runif(1, -.003, .003)
  I.0 = trunc((5/331+X.unif)*n)#max(1, trunc(length(which(inf.dates == 0)) + n*X.unif))#trunc(rbeta(1,1+length(which(inf.dates == 0)/2), 199+n)*n)
  E.0 = trunc((5/331+X.unif)*n)#max(1, trunc(length(which(ex.dates <= 0 & inf.dates > 0)) + n*X.unif))#trunc(rbeta(1,1+length(which(inf.dates == 0)/2), 199+n)*n)
  S.0 = 1 - (I.0+E.0)/n
  S.e.0 = 1 - I.0/n
  delta[1] = length(which(ex.dates == 1))
  omega[1] = length(which(inf.dates == 1))
  gamma[1] = rbeta(1,a[1],b[1])
  gamma.new[1] = rbeta(1,a[1],b[1])
  eps.new[1] = rbinom(1, I.0, gamma[1])
  I[1] = I.0 + omega[1] - eps[1]
  I.new[1] = I.0 + omega[1] - eps.new[1]
  e.t[1] = E.0 + delta[1] - omega[1]
  S.t[1] = n - delta[1]
  beta[1] = I[1]/n*rbeta(1, delta[1] + 1, n - delta[1] + 1)
  S.prod[1] = (1 - beta[1]/n)^I.0
  S[1] = S.0*S.prod[1]
  S.e[1] = S.comp(beta, t=1, I, S, I.0, E.0, S.0, S.e.0)+(S.e.0-S.0)*(1-phi.t[1])
  
  for(h in 2:n.days)
  {
    eps.new[h] = rbinom(1, I.new[h-1], gamma.new[h-1])
    delta[h] = length(which(ex.dates == h))
    omega[h] = length(which(inf.dates == h))
    I[h] = I[h-1] + omega[h] - eps[h]
    I.new[h] = I.new[h-1] + omega[h] - eps.new[h]
    S.t[h] = S.t[h-1] - delta[h]
    e.t[h] = e.t[h-1] + delta[h] - omega[h]
    gamma.new[h] = rbeta(1, eps.new[h] + a[h], I.new[h-1] - eps.new[h] + b[h])
  }
  if (length(which(I <= 0)) > 0)
  {
    flip.eps = 1
  } else {
    ratio.eps = loglik.eps(eps.new, delta, I.new, S.t) - loglik.eps(eps, I, delta, S.t)
    ratio.eps = min(1, exp(ratio.eps), na.rm = T)
    flip.eps = rbinom(1,1,ratio.eps)
  }
  if (flip.eps == 1)
  {
    eps = eps.new
    I = I.new
    accept.eps = accept.eps + 1
  }
  
  for (j in 2:(n.days))
  {
    # delta[j] = length(which(ex.dates == j))
    # omega[j] = length(which(inf.dates == j))
    # S.t[j] = S.t[j-1] - delta[j]
    # e.t[j] = e.t[j-1] + delta[j] - omega[j]
    beta[j] = min(1,n/I[j]*rbeta(1, delta[j] + .5, S.t[j-1] - delta[j] + .5))
    gamma[j] = rbeta(1, eps[j] + a[h], I[j-1] - eps[j] + b[h])
    S.prod[j] = (1 - beta[j-1]/n)^I[j-1]
    S[j] = S[j-1]*S.prod[j]
    S.e[j] = S.comp(beta, t=j, I, S, I.0, E.0, S.0, S.e.0)+(S.e.0-S.0)*(1-phi.t[1])^(j)
  }
  p.t = p.fun(S.e, S.e.0)
  p.t.ex = p.fun(S, S.0)
  p.t[p.t <= 0] = 1/n
  p.t.ex[p.t.ex <= 0] = 1/n
  
  inf.dates = as.vector(apply(test.mat[1:n,], 1, decensor, p.t=p.t, S.0 = S.0))
  ex.dates = inf.dates - ex.geom
  #ex.dates.new = inf.dates - rgeom(length(inf.dates), phi.t[1])-1
  
  for (l in which(inf.dates > 0))
  {
    ex.date.new = inf.dates[l] - rgeom(1, phi.t[1]) - 1
    if (is.na(ex.dates[l]))
    {
      a.ratio = 1
    } else {
      
      if (ex.dates[l] <= 0 & ex.date.new <= 0)
      {
        a.ratio = 1
      }
      if (ex.dates[l] > 0 & ex.date.new <= 0)
      {
        a.ratio = min(1, (S.0*(1-S.0))/(S[ex.dates[l]]*p.t[ex.dates[l]]))
      }
      if (ex.dates[l] <= 0 & ex.date.new > 0)
      {
        a.ratio = min(1, (S[ex.date.new]*p.t.ex[ex.date.new])/(S.0*(1-S.0)))
      }
      if (ex.dates[l] > 0 & ex.date.new > 0)
      {
        a.ratio = min(1, (S[ex.date.new]*p.t.ex[ex.date.new])/(S[ex.dates[l]]*p.t.ex[ex.dates[l]]))
      }
    }
    flip.ex = rbinom(1,1,a.ratio)
    if (flip.ex == 1)
    {
      ex.dates[l] = ex.date.new
    }
    
  }
  ex.geom = inf.dates - ex.dates
  
  
  S.t.samp[k,] = S.t
  S.samp[k,] = S
  I.samp[k,] = I
  eps.samp[k,] = eps
  delta.samp[k,] = delta
  beta.samp[k,] = beta
  gamma.samp[k,] = gamma
}

pos.count = osu.counts$Total_pos/osu.counts$Total
na.rows = which(is.na(pos.count))
badrows = which(osu.counts$Total < 200)
pos.count[badrows] = NA
I.smooth = pos.count
smoothed = supsmu(c(1:n.days), pos.count)

j = 1
for (i in smoothed$x)
{
  I.smooth[i] = smoothed$y[j]
  j = j+1
}


count.se = sqrt(I.smooth*(1-I.smooth)/osu.counts$Total)
I.smooth.max = I.smooth + qnorm(.975)*count.se
I.smooth.min = I.smooth - qnorm(.975)*count.se

I.smooth.max[na.rows] = NA
I.smooth.min[na.rows] = NA
I.smooth.min[I.smooth.min < 0] = 0

I.mean = apply(I.samp[(.5*nreps+1):nreps,], 2, mean)
I.min = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
I.q1 = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .25)
I.med = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
I.q3 = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .75)
I.max = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)

beta.mean = apply(beta.samp[(.5*nreps+1):nreps,], 2, mean)
beta.min = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
beta.q1 = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .25)
beta.med = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
beta.q3 = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .75)
beta.max = apply(beta.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)

gamma.med = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .5)
gamma.max = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)
gamma.min = apply(gamma.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)

t <- list(
  family = "arial",
  size = 22,
  color = "black")

t.2 <- list(
  family = "arial",
  size = 22,
  color = "black")


fig.osu <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
fig.osu <- fig.osu %>% add_trace(x = c(1:(n.days-4)), y = I.min[1:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
fig.osu <- fig.osu %>% add_trace(x = c(1:(n.days-4)), y = I.med[1:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
fig.osu <- fig.osu %>% add_trace(x = c(1:(n.days-4)), y = I.max[1:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
fig.osu <- fig.osu  %>% add_trace(x = c(1:(n.days-4)), y = n*I.smooth[1:(n.days-4)], mode = "markers", name = 'Naive prevalence est')
fig.osu <- fig.osu %>% add_trace(x = c(1:(n.days-4)), y = I.smooth.min[1:(n.days-4)]*n, mode = "none", name = 'Naive prevalence 95% CI', fill = "tonexty", fillcolor = "rgba(255, 116, 134, 0.5)")
fig.osu <- fig.osu %>% add_trace(x = c(1:(n.days-4)), y = I.smooth.max[1:(n.days-4)]*n, mode = "none", name = '2.5%', showlegend = FALSE, fill = "tonexty", fillcolor = "rgba(255, 116, 134, 0.5)")


fig.osu <- fig.osu %>% layout(title = "", legend = list(x=.55, y=1.1, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                              yaxis = list(title = "Infected", range = c(0,750)), font = t,
                              margin = list(b = 100))
fig.osu <- fig.osu %>% config(mathjax = 'cdn')
fig.osu

R.0.samp = beta.samp[(.5*nreps+1):nreps,]/gamma.samp[(.5*nreps+1):nreps,]*(S.t.samp[(.5*nreps+1):nreps,]/n)
R.t.min = apply(R.0.samp, 2, quantile, probs = .025)
R.t.q1 = apply(R.0.samp, 2, quantile, probs = .25)
R.t.med = apply(R.0.samp, 2, quantile, probs = .5)
R.t.q3 = apply(R.0.samp, 2, quantile, probs = .75)
R.t.max = apply(R.0.samp, 2, quantile, probs = .975)

fig.rt <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.min[2:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.med[2:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.max[2:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = true.Rt[2:(n.days-4)], mode = "lines", name = 'True partial R_t')

fig.rt <- fig.rt %>% layout(title = "", legend = list(x=.85, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                            yaxis = list(title = "R_t", range = c(0,8)), font = t,
                            margin = list(b = 100))
fig.rt <- fig.rt %>% config(mathjax = 'cdn')
fig.rt

#####Output
print("Output complete.")