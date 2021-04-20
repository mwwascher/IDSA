library(plotly)
library(processx)

#####Initialization
n = 10000
n.days = 60
set.seed(456)

test.mat = matrix(NA, nrow = n, ncol = 2)
prop.mat = matrix(NA, nrow = n.days, ncol = 2)

beta.t = rep(.12, n.days)
beta.t[1:7] = .35
beta.t[8:14] = .45
beta.t[15:21] = .08
beta.t[22:28] = .08
trace.prob = 1/1.15
gamma.t = rep(.25, n.days)
phi.t = rep(.25, n.days)

S <- vector(mode = "list", length = n.days)
E <- vector(mode = "list", length = n.days)
I <- vector(mode = "list", length = n.days)
R <- vector(mode = "list", length = n.days)
traced <- vector(mode = "list", length = n.days)
test <- vector(mode = "list", length = n.days)
order = sample(c(1:10000), 10000)
schedule = list(order[1:1500], order[1501:3000], order[3001:4500], order[4501:6000], order[6001:7500],
                order[7501:9000], order[9001:10000])


S[[1]] = c(101:n)
E[[1]] = c(51:100)
I[[1]] = c(1:50)
R[[1]] = c()

I.count = c()
S.count = c()
E.count = c()
R.count = c()
trace.count = c()
pos.count = c()
rem.count = c()
test.count = c()
pos.pct = c()
adj.test = rep(1500, n.days-1)
adj.n = rep(n, n.days-1)

#####Simulate the epidemic
for (i in 1:(n.days-1))
{
  
  if ((i %% 7) == 1)
  {
    order = sample(c(1:10000), 10000, replace = FALSE)
    schedule = list(order[1:1500], order[1501:3000], order[3001:4500], order[4501:6000], order[6001:7500],
                    order[7501:9000], order[9001:10000])
  }
  
  test.sched = i %% 7
  if (test.sched == 0)
  {
    test.sched = 7
  }
  test[[i]] = schedule[[test.sched]]
  
  num.S = length(S[[i]])
  num.I = length(I[[i]])
  num.E = length(E[[i]])
  num.R = length(R[[i]])
  num.trace = length(traced[[i]])
  num.pos = length(I[[i]][I[[i]] %in% test[[i]]])
  num.test = length(test[[i]])
  
  I.count = c(I.count, num.I)
  S.count = c(S.count, num.S)
  E.count = c(E.count, num.E)
  R.count = c(R.count, num.R)
  pos.count = c(pos.count, num.pos)
  test.count = c(test.count, num.test)
  
  num.delta = rbinom(1,num.S, (1-(1-beta.t[i]/n)^(num.I)))
  delta.t = sample(S[[i]], num.delta)
  print("I")
  
  num.omega = rbinom(1, num.E, phi.t[i])
  omega.t = sample(E[[i]], num.omega)
  print("E")
  
  adj.test[i] = adj.test[i] - length(R[[i]][R[[i]] %in% test[[i]]])
  pos.pct[i] = length(I[[i]][I[[i]] %in% test[[i]]])/length(test[[i]])
  
  if(i > 2)
  {
    eps.t = I[[i-2]][I[[i-2]] %in% test[[i-2]]]
    num.trace = min(sum(rgeom(length(I[[i]][I[[i]] %in% eps.t]), trace.prob)),length(I[[i]][!I[[i]] %in% eps.t]))
    trace.remove = sample(I[[i]][!I[[i]] %in% eps.t], num.trace)
    traced[[i]] = trace.remove
    eps.t = c(eps.t, trace.remove)
    print("R")
  }
  
  if (i <= 2)
  {
    num.remove = rbinom(1, length(I[[i]]), .1)
    eps.t = sample(I[[i]], num.remove)
    
    S[[i+1]] = S[[i]][!S[[i]] %in% delta.t]
    I[[i+1]] = c(I[[i]][!I[[i]] %in% eps.t], omega.t)
    E[[i+1]] = c(E[[i]][!E[[i]] %in% omega.t], delta.t)  
    R[[i+1]] = c(R[[i]], I[[i]][I[[i]] %in% eps.t])
  } else
  {
    S[[i+1]] = S[[i]][!S[[i]] %in% delta.t]
    I[[i+1]] = c(I[[i]][!I[[i]] %in% eps.t], omega.t)
    E[[i+1]] = c(E[[i]][!E[[i]] %in% omega.t], delta.t)  
    R[[i+1]] = c(R[[i]], I[[i]][I[[i]] %in% eps.t])
  }
  adj.n[i] = adj.n[i] - length(R[[i]])
  rem.count = c(rem.count,length(eps.t))
}

#####Build the test matrix
for (i in 1:(n.days-1))
{
  pop = c(1:n)
  pop.test = pop[pop %in% test[[i]]]
  pop.pos = pop.test[pop.test %in% I[[i]]]
  pop.neg = pop.test[pop.test %in% S[[i]]]
  trace.pos = traced[[i]]
  
  tot.test = length(test[[i]])
  prop.pos = length(c(pop.pos, trace.pos))/tot.test
  prop.neg = length(pop.neg)/tot.test
  
  test.mat[c(pop.pos,traced[[i]]),2] = i
  test.mat[pop.neg,1] = i
  
  prop.mat[i,2] = prop.pos
  prop.mat[i,1] = prop.neg
  
  num.trace = length(traced[[i]])
  trace.count = c(trace.count, num.trace)
}

#####Take a look at the simulated epidemic
# plot(S.count, main = "Susceptibles", xlab = "day", ylab = "S")
# plot(I.count, main = "Infecteds", xlab = "day", ylab = "I", ylim = c(0,400))

#########Model Functions

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

IDSAplot = function()
{
  write.csv(I.samp[(.5*nreps+1):nreps,], file = "posterior_I.csv")
  write.csv(S.samp[(.5*nreps+1):nreps,], file = "posterior_S.csv")
  write.csv(delta.samp[(.5*nreps+1):nreps,], file = "posterior_delta.csv")
  write.csv(eps.samp[(.5*nreps+1):nreps,], file = "posterior_eps.csv")
  write.csv(beta.samp[(.5*nreps+1):nreps,], file = "posterior_beta.csv")
  write.csv(gamma.samp[(.5*nreps+1):nreps,], file = "posterior_gamma.csv")
  
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
  plot(gamma.max, ylim = c(0,.5), type = "l")
  lines(gamma.min)
  points(gamma.med)
  
  t <- list(
    family = "arial",
    size = 22,
    color = "black")
  
  t.2 <- list(
    family = "arial",
    size = 24,
    color = "black")
  
  
  plot(beta.med, type = "l", ylim = c(0,1))
  lines(beta.min, col = "red")
  lines(beta.max, col = "red")
  lines(beta.t, col = "blue")
  
  fig.on <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
  fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.min[1:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
  fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.med[1:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.on <- fig.on %>% add_trace(x = c(1:(n.days-4)), y = I.max[1:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.on <- fig.on  %>% add_trace(x = c(1:(n.days-4)), y = I.count[1:(n.days-4)], mode = "markers", name = 'True prevalence')
  fig.on <- fig.on  %>% add_segments(x = 8, xend = 8, y = 0, yend = 450, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.on <- fig.on  %>% add_segments(x = 15, xend = 15, y = 0, yend = 450, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.on <- fig.on  %>% add_segments(x = 29, xend = 29, y = 0, yend = 450, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  
  
  fig.on <- fig.on %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                              yaxis = list(title = "Infected", range = c(0,400)), font = t,
                              margin = list(b = 100))
  fig.on <- fig.on %>% config(mathjax = 'cdn')
  orca(fig.on, "Closed_fit.svg")
  
  
  #####R_t
  true.Rt = beta.t[2:60]/(1/(I.count/rem.count))*S.count/n
  true.Rt.smooth = supsmu(c(2:60), true.Rt)
  R.0.samp = beta.samp[(.5*nreps+1):nreps,]/gamma.samp[(.5*nreps+1):nreps,]*(S.t.samp[(.5*nreps+1):nreps,]/n)
  R.t.min = apply(R.0.samp, 2, quantile, probs = .025)
  R.t.q1 = apply(R.0.samp, 2, quantile, probs = .25)
  R.t.med = apply(R.0.samp, 2, quantile, probs = .5)
  R.t.q3 = apply(R.0.samp, 2, quantile, probs = .75)
  R.t.max = apply(R.0.samp, 2, quantile, probs = .975)
  rem.rate = supsmu(c(2:60),1/(I.count/pos.count))
  
  fig.rt <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
  fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.min[2:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
  fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.med[2:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = R.t.max[2:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = true.Rt.smooth$y[2:(n.days-4)], mode = "lines", name = 'True R_t')
  fig.rt <- fig.rt  %>% add_segments(x = 8, xend = 8, y = 0, yend = 450, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.rt <- fig.rt  %>% add_segments(x = 15, xend = 15, y = 0, yend = 450, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.rt <- fig.rt  %>% add_segments(x = 29, xend = 29, y = 0, yend = 450, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  
  fig.rt <- fig.rt %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                              yaxis = list(title = "R_t", range = c(0,8)), font = t,
                              margin = list(b = 100))
  fig.rt <- fig.rt %>% config(mathjax = 'cdn')
  orca(fig.rt, "Closed_Rt.svg")
}

########Model fitting Initialization
inf.prop = (1-length(which(is.na(test.mat[,2])))/n)
p.t = rep(inf.prop/n.days, n.days)
p.t.ex = p.t
p.t.unif = rep(1/n.days, n.days)
p.t.null = rep(inf.prop/n.days, n.days)n
S.null = .99
I = rep(5, n.days)
I.new = rep(0, n.days)
delta = rep(0, n.days)
eps = rep(100, n.days)
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
nreps = 2000
a = c(1,1, rep(10,58))
b = c(9,9, rep(50, 58))

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
for (k in 1:nreps)
{
  
  X.unif = runif(1, -.003, .003)
  I.0 = trunc((length(R[[3]])/1500+X.unif)*n)#max(1, trunc(length(which(inf.dates == 0)) + n*X.unif))#trunc(rbeta(1,1+length(which(inf.dates == 0)/2), 199+n)*n)
  E.0 = trunc((length(R[[3]])/1500+X.unif)*n)#max(1, trunc(length(which(ex.dates <= 0 & inf.dates > 0)) + n*X.unif))#trunc(rbeta(1,1+length(which(inf.dates == 0)/2), 199+n)*n)
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
      ratio.eps = min(1, exp(ratio.eps))
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

#####Output
IDSAplot()