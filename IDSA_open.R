library(plotly)
library(processx)

#####Initialize
dir.create("~/IDSAOpen")
print("This script simulates an epidemic and fits the open IDSA model to the simulated epidemic.")

n = 10000
n.off = 30000
n.days = 60
n.tests = 1500
set.seed(789)

test.mat = matrix(NA, nrow = n, ncol = 2)
prop.mat = matrix(NA, nrow = n.days, ncol = 2)

beta.t = rep(.12, n.days)
beta.t[1:7] = .3
beta.t[8:14] = .35
beta.t[14:21] = .08
beta.t[22:28] = .08
trace.prob = .85
gamma.t = rep(.25, n.days)
phi.t = rep(.25, n.days)

beta.off = rep(.12, n.days)
gamma.off = rep(.1, n.days)
beta.cross = rep(.0225, n.days)


S <- vector(mode = "list", length = n.days)
E <- vector(mode = "list", length = n.days)
I <- vector(mode = "list", length = n.days)
R <- vector(mode = "list", length = n.days)

S.off <- vector(mode = "list", length = n.days)
E.off <- vector(mode = "list", length = n.days)
I.off <- vector(mode = "list", length = n.days)
R.off <- vector(mode = "list", length = n.days)

traced <- vector(mode = "list", length = n.days)
test <- vector(mode = "list", length = n.days)
order = sample(c(1:10000), 10000)
schedule = list(order[1:1500], order[1501:3000], order[3001:4500], order[4501:6000], order[6001:7500],
                order[7501:9000], order[9001:10000])

test.off <- vector(mode = "list", length = n.days)
order.off = sample(c(1:30000), 30000)
schedule.off = list(order.off[1:1000], order.off[1001:2000], order.off[2001:3000], order.off[3001:4000], order.off[4001:5000],
                    order.off[5001:6000], order.off[6001:7000], order.off[7001:8000], order.off[8001:9000], order.off[9001:10000],
                    order.off[10001:11000], order.off[11001:12000], order.off[12001:13000], order.off[13001:14000], order.off[14001:15000],
                    order.off[15001:16000], order.off[16001:17000], order.off[17001:18000], order.off[18001:19000], order.off[19001:20000],
                    order.off[20001:21000], order.off[21001:22000], order.off[22001:23000], order.off[23001:24000], order.off[24001:25000],
                    order.off[25001:26000], order.off[26001:27000], order.off[27001:28000], order.off[28001:29000], order.off[29001:30000])


S[[1]] = c(101:n)
E[[1]] = c(51:100)
I[[1]] = c(1:50)
R[[1]] = c()

S.off[[1]] = c(601:30000)
I.off[[1]] = c(1:300)
E.off[[1]] = c(301:600)
R[[1]] = c()

I.count = c()
S.count = c()
E.count = c()
R.count = c()
trace.count = c()

I.count.off = c()
S.count.off = c()
E.count.off = c()
R.count.off = c()
trace.count.off = c()
pos.count = rep(0, n.days-1)

####Simulate the epidemic
print("Simulating epidemic...")
for (i in 1:(n.days-1))
{
  
  if ((i %% 7) == 1)
  {
    order = sample(c(1:10000), 10000, replace = FALSE)
    schedule = list(order[1:1500], order[1501:3000], order[3001:4500], order[4501:6000], order[6001:7500],
                    order[7501:9000], order[9001:10000])
  }
  
  if ((i %% 30) == 1)
    {
    order.off = sample(c(1:30000), 30000)
    schedule.off = list(order.off[1:1000], order.off[1001:2000], order.off[2001:3000], order.off[3001:4000], order.off[4001:5000],
                        order.off[5001:6000], order.off[6001:7000], order.off[7001:8000], order.off[8001:9000], order.off[9001:10000],
                        order.off[10001:11000], order.off[11001:12000], order.off[12001:13000], order.off[13001:14000], order.off[14001:15000],
                        order.off[15001:16000], order.off[16001:17000], order.off[17001:18000], order.off[18001:19000], order.off[19001:20000],
                        order.off[20001:21000], order.off[21001:22000], order.off[22001:23000], order.off[23001:24000], order.off[24001:25000],
                        order.off[25001:26000], order.off[26001:27000], order.off[27001:28000], order.off[28001:29000], order.off[29001:30000])
    }
  
  test.sched = i %% 7
  if (test.sched == 0)
  {
    test.sched = 7
  }
  test[[i]] = schedule[[test.sched]]
  
  test.sched.off = i %% 30
  if (test.sched.off == 0)
  {
    test.sched.off = 30  
  }
  test.off[[i]] = schedule.off[[test.sched.off]]
  
  num.S = length(S[[i]])
  num.I = length(I[[i]])
  num.E = length(E[[i]])
  num.R = length(R[[i]])
  num.trace = length(traced[[i]])
  
  num.S.off = length(S.off[[i]])
  num.I.off = length(I.off[[i]])
  num.E.off = length(E.off[[i]])
  num.R.off = length(R.off[[i]])
  
  I.count = c(I.count, num.I)
  S.count = c(S.count, num.S)
  E.count = c(E.count, num.E)
  R.count = c(R.count, num.R)
  trace.count = c(trace.count, num.trace)
  
  I.count.off = c(I.count.off, num.I.off)
  S.count.off = c(S.count.off, num.S.off)
  E.count.off = c(E.count.off, num.E.off)
  R.count.off = c(R.count.off, num.R.off)
  
  num.delta = rbinom(1,num.S, (1-(1-beta.t[i]/n)^(num.I))) + rbinom(1,num.S, (1-(1-beta.cross[i]/n)^(num.I.off)))
  delta.t = sample(S[[i]], num.delta)
  
  num.delta.off = rbinom(1, num.S.off, (1-(1-beta.off[i]/n.off)^(num.I.off))) + rbinom(1, num.S.off, (1-(1-beta.cross[i]/n.off)^(num.I)))
  delta.t.off = sample(S.off[[i]], num.delta.off)
  
  num.eps.off = rbinom(1, num.I.off, gamma.off[i])
  eps.t.off = sample(I.off[[i]], num.eps.off)
  
  num.omega = rbinom(1, num.E, phi.t[i])
  omega.t = sample(E[[i]], num.omega)
  
  num.omega.off = rbinom(1, num.E.off, phi.t[i])
  omega.t.off = sample(E.off[[i]], num.omega.off)
  
  eps.t = c()
  if(i > 2)
  {
    eps.t = I[[i-2]][I[[i-2]] %in% test[[i-2]]]
    num.trace = min(sum(rgeom(length(I[[i]][I[[i]] %in% eps.t]), trace.prob)),length(I[[i]][!I[[i]] %in% eps.t]))
    trace.remove = sample(I[[i]][!I[[i]] %in% eps.t], num.trace)
    traced[[i]] = trace.remove
    eps.t = c(eps.t, trace.remove)
    pos.count[i] = length(eps.t)
  }
  
  if (i <= 2)
  {
    num.remove = rbinom(1, length(I[[i]]), .1)
    eps.t = sample(I[[i]], num.remove)
    pos.count[i] = length(eps.t)
    
    S[[i+1]] = S[[i]][!S[[i]] %in% delta.t]
    I[[i+1]] = c(I[[i]][!I[[i]] %in% eps.t], omega.t)
    E[[i+1]] = c(E[[i]][!E[[i]] %in% omega.t], delta.t)  
    R[[i+1]] = I[[i]][I[[i]] %in% eps.t]
  } else
  {
    S[[i+1]] = S[[i]][!S[[i]] %in% delta.t]
    I[[i+1]] = c(I[[i]][!I[[i]] %in% eps.t], omega.t)
    E[[i+1]] = c(E[[i]][!E[[i]] %in% omega.t], delta.t)  
    R[[i+1]] = I[[i]][I[[i]] %in% eps.t]
  }
  
  S.off[[i+1]] = S.off[[i]][!S.off[[i]] %in% delta.t.off]
  I.off[[i+1]] = c(I.off[[i]][!I.off[[i]] %in% eps.t.off], omega.t.off)
  E.off[[i+1]] = c(E.off[[i]][!E.off[[i]] %in% omega.t.off], delta.t.off)  
  R.off[[i+1]] = I.off[[i]][I.off[[i]] %in% eps.t.off] 
}

#####Take a look at the epidemic
# plot(S.count, main = "Susceptibles", xlab = "day", ylab = "S")
# plot(I.count, main = "Infecteds", xlab = "day", ylab = "I")
# plot(S.count.off, main = "Susceptibles", xlab = "day", ylab = "S, off campus")
# plot(I.count.off, main = "Infecteds", xlab = "day", ylab = "I, off campus", ylim = c(0,900))
# points(I.count, col = "red")


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
  
  if(i == 1)
    {
    prop.0 = length(c(pop.pos))/tot.test
    }
  
  test.mat[c(pop.pos,traced[[i]]),2] = i
  test.mat[pop.neg,1] = i
  
  prop.mat[i,2] = prop.pos
  prop.mat[i,1] = prop.neg
}

#####Naive background prevlance estimate
prev.off = rep(NA, n.days)
for (i in 1:(n.days-1))
{
  pop.off = c(1:n.off)
  pop.test.off = pop.off[pop.off %in% test.off[[i]]]
  pop.pos.off = pop.test.off[pop.test.off %in% I.off[[i]]]
  
  tot.test.off = length(test.off[[i]])
  prop.pos.off = length(pop.pos.off)/tot.test.off
  prev.off[i] = prop.pos.off
}

#########Functions
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
    if (left.end > n.days)
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

S.comp.group = function(S, S.0, t)
{
  S.E = 0
  if (t > 1)
  {
    for (i in 1:(t))
    {
      if (i == 1)
      {
        S.E = S.E + (S.0-S[1])*S.0*(1-phi.t[1])^(t-i)
      } else {
        S.E = S.E + (S[i-1]-S[i])*S[i-1]*(1-phi.t[1])^(t-i)
      }
    }
  }
  S.E = S.E + S[t]/S.0
  return(S.E)
}

loglik.eps = function(eps, delta, I, S)
{
  loglik = 0
  lambda = S*psi.1/n
  loglik = loglik + log(lambda[1]^delta[1]*exp(-lambda[1])/factorial(delta[1]))
  for (i in 2:n.days)
  {
    loglik = loglik + log(lambda[i]^delta[i]*exp(-lambda[i])/factorial(delta[i]))
  }
  return(loglik)
}

IDSAplot_open = function()
{
  
  print("Writing posterior I_1 samples to posterior_I1.csv...")
  write.csv(I.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_I1.csv")
  print("Writing posterior s_1 samples to posterior_S1.csv...")
  write.csv(S.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_S1.csv")
  print("Writing posterior delta_1 samples to posterior_delta1.csv...")
  write.csv(delta.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_delta1.csv")
  print("Writing posterior epsilon_1 samples to posterior_eps1.csv...")
  write.csv(eps.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_eps1.csv")
  print("Writing posterior phi_1 samples to posterior_phi1.csv...")
  write.csv(psi.1.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_phi1.csv")
  print("Writing posterior gamma_1 samples to posterior_gamma1.csv...")
  write.csv(gamma.samp[(.5*nreps+1):nreps,], file = "~/IDSAOpen/posterior_gamma1.csv")
  
  I.med = apply(I.samp[(.5*nreps+1):nreps,], 2, median)
  I.max = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .975)
  I.min = apply(I.samp[(.5*nreps+1):nreps,], 2, quantile, probs = .025)
  
  S.t.med = apply(S.t.samp[(.5*nreps+1):nreps,], 2, median)
  
  phi.med = apply(psi.1.samp[(.5*nreps+1):nreps,], 2, median)
  print("Generating plots...")
  
  fig.off <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
  fig.off <- fig.off %>% add_trace(x = c(1:(n.days-4)), y = I.min[1:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
  fig.off <- fig.off %>% add_trace(x = c(1:(n.days-4)), y = I.med[1:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.off <- fig.off %>% add_trace(x = c(1:(n.days-4)), y = I.max[1:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.off <- fig.off  %>% add_trace(x = c(1:(n.days-4)), y = I.count[1:(n.days-4)], mode = "markers", name = 'True prevalence')
  fig.off <- fig.off  %>% add_segments(x = 8, xend = 8, y = 0, yend = 1000, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.off <- fig.off  %>% add_segments(x = 15, xend = 15, y = 0, yend = 1000, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.off <- fig.off  %>% add_segments(x = 29, xend = 29, y = 0, yend = 1000, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.off <- fig.off  %>% add_segments(x = 22, xend = 22, y = 0, yend = 1000, name = "chpt 3", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  
  t <- list(
    family = "arial",
    size = 22,
    color = "black")
  
  t.2 <- list(
    family = "arial",
    size = 24,
    color = "black")
  
  fig.off <- fig.off %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                                yaxis = list(title = "Infected", range = c(0,660)), font = t,
                                margin = list(b = 100))
  fig.off <- fig.off %>% config(mathjax = 'cdn')
  orca(fig.off, "IDSAOpen/Open_fit.svg")
  
  
  #####R_t
  gamma.true = pos.count/I.count
  true.Rt = ((beta.t[1:59]*I.count) + (beta.cross[1:59]*I.count.off))/(gamma.true)*(S.count/I.count)/n
  true.Rt[1:2] = ((beta.t[1:2]*I.count[1:2]) + (beta.cross[1:2]*I.count.off[1:2]))/(c(1/6,1/6))*(S.count[1:2]/I.count[1:2])/n
  true.Rt.smooth = supsmu(c(2:60), true.Rt)
  
  #R.0.samp = (psi.1.samp[(.5*nreps+1):nreps,1:59]*I.samp[(.5*nreps+1):nreps,1:59]/n)/gamma.samp[(.5*nreps+1):nreps,1:59]*(S.t.samp[(.5*nreps+1):nreps,1:59]/n)
  #R.0.samp = psi.1.samp[(.5*nreps+1):nreps,1:59]/gamma.samp[(.5*nreps+1):nreps,1:59]*(I.samp[(.5*nreps+1):nreps,1:59]*S.t.samp[(.5*nreps+1):nreps,1:59])/(n*n)
  R.0.samp = (psi.1.samp[(.5*nreps+1):nreps,1:59]/gamma.samp[(.5*nreps+1):nreps,1:59])*(S.t.samp[(.5*nreps+1):nreps,1:59]/I.samp[(.5*nreps+1):nreps,1:59])/n
  R.t.min = apply(R.0.samp, 2, quantile, probs = .025)
  R.t.q1 = apply(R.0.samp, 2, quantile, probs = .25)
  R.t.med = apply(R.0.samp, 2, quantile, probs = .5)
  R.t.q3 = apply(R.0.samp, 2, quantile, probs = .75)
  R.t.max = apply(R.0.samp, 2, quantile, probs = .975)
  
  # R.t.min = supsmu(c(1:59), R.t.min)
  # R.t.min = R.t.min$y
  # R.t.med = supsmu(c(1:59), R.t.med)
  # R.t.med = R.t.med$y
  # R.t.max = supsmu(c(1:59), R.t.max)
  # R.t.max = R.t.max$y
  
  fig.rt <- plot_ly(type = "scatter", mode = "lines", showlegend = TRUE)
  fig.rt <- fig.rt %>% add_trace(x = c(1:(n.days-4)), y = R.t.min[1:(n.days-4)], mode = "lines", name = '2.5%', line = list(color = "rgb(44, 160, 44)"))
  fig.rt <- fig.rt %>% add_trace(x = c(1:(n.days-4)), y = R.t.med[1:(n.days-4)], mode = "lines", name = '50%', line = list(color = "rgb(255, 127, 14)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.rt <- fig.rt %>% add_trace(x = c(1:(n.days-4)), y = R.t.max[1:(n.days-4)], mode = "lines", name = '97.5%', line = list(color = "rgb(214, 39, 40)"), fill = "tonexty", fillcolor = "rgba(168, 216, 234, 0.5)")
  fig.rt <- fig.rt %>% add_trace(x = c(1:(n.days-4)), y = true.Rt.smooth$y[1:(n.days-4)], mode = "lines", name = 'True partial R_t')
  # fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = phi.med[2:(n.days-4)]/5, mode = "lines", name = 'phi/5')
  # fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = I.med[2:(n.days-4)]/30, mode = "lines", name = 'I/30')
  # fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = phi.med[2:(n.days-4)]/5*I.med[2:(n.days-4)]/300, mode = "lines", name = 'phi*I')
  #fig.rt <- fig.rt %>% add_trace(x = c(2:(n.days-4)), y = phi.med[2:(n.days-4)]/5*I.med[2:(n.days-4)]/300*S.t.med[2:(n.days-4)]/n, mode = "lines", name = 'phi*I*S/n')
  fig.rt <- fig.rt  %>% add_segments(x = 8, xend = 8, y = 0, yend = 450, name = TeX("\\beta_t \\text{ changepoint}"), line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.rt <- fig.rt  %>% add_segments(x = 15, xend = 15, y = 0, yend = 450, name = "chpt 2", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.rt <- fig.rt  %>% add_segments(x = 22, xend = 22, y = 0, yend = 450, name = "chpt 3", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  fig.rt <- fig.rt  %>% add_segments(x = 29, xend = 29, y = 0, yend = 450, name = "chpt 1", line = list(color = 'rgb(12, 12, 24)', width = 2, dash = 'dot'), showlegend = FALSE)
  
  fig.rt <- fig.rt %>% layout(title = "", legend = list(x=.75, y=.95, font = t.2), xaxis = list(title = "Day", tickfont = list(size = 12)),
                              yaxis = list(title = "Partial R_t", range = c(0,6)), font = t,
                              margin = list(b = 100))
  fig.rt <- fig.rt %>% config(mathjax = 'cdn')
  fig.rt
  
  
}

#####Smooth background group I
I.off = rep(0, n.days)
for (i in 1:n.days)
  {
  I.off[i] = prev.off[i]*30000
}
I.smoothed = smooth.spline(c(1:59), I.off[1:59], lambda = .0005)
I.smoothed = trunc(I.smoothed$y)

#####Model fit initialization
inf.prop = (1-length(which(is.na(test.mat[,2])))/n)
p.t = rep(inf.prop/n.days, n.days)
S.0 = .99
S.e.0 = .99
p.t.null = rep(1/n.days, n.days)
I = rep(0, n.days)
I.new = rep(0, n.days)
delta = rep(0, n.days)
delta.off = rep(0, n.days)
beta = rep(.4, n.days)
eps = rep(0, n.days)
eps.new = rep(0, n.days)
eps.off = rep(0, n.days)
omega = rep(0, n.days)
omega.off = rep(0, n.days)
gamma = rep(.2, n.days)
gamma.new = rep(.166, n.days)
psi.1 = rep(.4, n.days)
psi.2 = rep(.4, n.days)
phi = rep(.25, n.days)
S = rep(.99, n.days)
S.e = rep(.99, n.days)
S.t = rep(n, n.days)
S.t.off = rep(n.off, n.days)
S.prod = rep(0, n.days)
e.t = rep(0, n.days)
e.t.off = rep(0, n.days)
nreps = 10000
a = c(1,1, rep(10,58))
b = c(9,9, rep(50, 58))

S.samp = matrix(0, nrow = nreps, ncol = n.days)
S.t.samp = matrix(0, nrow = nreps, ncol = n.days)
beta.samp = matrix(0, nrow = nreps, ncol = n.days)
beta.off.samp = matrix(0, nrow = nreps, ncol = n.days)
beta.cross.samp = matrix(0, nrow = nreps, ncol = n.days)
gamma.samp = matrix(0, nrow = nreps, ncol = n.days)
delta.samp = matrix(0, nrow = nreps, ncol = n.days)
eps.samp = matrix(0, nrow = nreps, ncol = n.days)
I.samp = matrix(0, nrow = nreps, ncol = n.days)
e.samp = matrix(0, nrow = nreps, ncol = n.days)
beta.track = matrix(0, nrow = 2, ncol = nreps)
psi.1.samp = matrix(0, nrow = nreps, ncol = n.days)

decensor(test.mat[500,], p.t = p.t, S.0 = .99)
inf.dates = as.vector(apply(test.mat[1:n,], 1, decensor, p.t=p.t, S.0 = .99))  
ex.dates = inf.dates - rgeom(length(inf.dates), phi.t[1])-1
ex.geom = inf.dates - ex.dates
rec.dates = inf.dates + rgeom(length(inf.dates), 1/6)+1
ex.rec = 0
accept.eps = 0
bad.ratio = 0

accept = 0
reject = 0

#####Gibbs sampler
print("Running MCMC...")
for (i in 1:nreps)
  {
  X.unif = runif(1, -.003, .003)
  I.0 = trunc((length(R[[3]])/1500+X.unif)*n)
  E.0 = trunc((length(R[[3]])/1500+X.unif)*n)
  delta[1] = length(which(ex.dates == 1))
  omega[1] = length(which(inf.dates == 1))
  gamma[1] = rbeta(1,a[1],b[1])
  gamma.new[1] = rbeta(1,a[1],b[1])
  eps.new[1] = rbinom(1, I.0, gamma[1])
  eps.off[1] = rbinom(1,trunc(.01*n.off), gamma[1])
  I[1] = I.0 + omega[1] - eps[1]
  I.new[1] = I.0 + omega[1] - eps.new[1]
  delta.off[1] = max(0, I.smoothed[1] + eps.off[1] - trunc(.01*n.off))
  e.t[1] = I.0 + delta[1] - omega[1]
  S.t[1] = n - I.0 - E.0 - delta[1]
  S.t.off[1] = .98*n.off - delta.off[1]
  psi.1[1] = rgamma(1, delta[1] + 1, S.t[1]/n*(1+1))
  psi.2[1] = rgamma(1, delta.off[1] + 1, S.t.off[1]/n.off*(1+1))
  
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
    ratio.eps = loglik.eps(eps.new, delta, I.new, S.t) - loglik.eps(eps, delta, I, S.t)
    ratio.eps = min(1, exp(ratio.eps))
    if (is.na(ratio.eps))
    {
      bad.ratio = bad.ratio + 1
      ratio.eps = 1
    }
    flip.eps = rbinom(1,1,ratio.eps)
  }
  if (flip.eps == 1)
  {
    eps = eps.new
    I = I.new
    accept.eps = accept.eps + 1
  }
  
  
  for (j in 2:(n.days-1))
    {
    eps.off[j] = rbinom(1, I.smoothed[j-1], gamma.off[j])
    delta.off[j] = max(I.smoothed[j] - I.smoothed[j-1] + eps.off[j], 0)
    S.t.off[j] = S.t.off[j-1] - delta[j]
    #eps[j] = rbinom(1, I[j-1], gamma[j-1])
    
    #delta[j] = length(which(ex.dates == j))
    #omega[j] = length(which(inf.dates == j))
    #I[j] = I[j-1] + omega[j] - eps[j]
    #e.t[j] = e.t[j-1] + delta[j] - omega[j]
    #S.t[j] = S.t[j-1] - delta[j]
  }
  
  for (k in 2:n.days)
    {
    psi.1[k] = rgamma(1, delta[k] + .01, S.t[k]/n*(1+.01))
    psi.2[k] = rgamma(1, delta.off[k] + .01, S.t.off[k]/n.off*(1+.01))
    gamma[k] = rbeta(1, eps[j] + a[h], I[j-1] - eps[j] + b[h])
    gamma.off[k] = rbeta(1, eps.off[j] + a[k], I.smoothed[j-1] - eps.off[j] + b[k])
  }

  
  S.0 = (1 - (I.0 + E.0)/n)
  S.e.0 = (1 - I.0/n)
  n.count = 9900
  arrivals = c(0:n.count)
  probs = dpois(c(0:n.count), S.t[k]*psi.1[k]/n)
  prob.vec = probs*(n.count - arrivals)/n.count
  S.prod[1] = sum(prob.vec)
  S[1] = S.0*S.prod[1]
  n.surv = trunc(S[1]*n)
  S.e[1] = S.0+(S.e.0-S.0)*phi.t[1]^(1)
  
  for (k in 2:(n.days-1))
    {
    n.count = S.t[k-1]
    n.surv = trunc(S[k-1]*n)
    arrivals = c(0:n.surv)
    probs = dpois(c(0:n.surv), n.surv*psi.1[k]/n)
    prob.vec = probs*(n.surv - arrivals)/n.surv
    S.prod[k] = sum(prob.vec)
    S[k] = S[k-1]*S.prod[k]
  }
  
  for (h in 2:(n.days - 1))
    {
    S.e[h] = S.0*S.comp.group(S, S.0, h)+(S.e.0-S.0)*phi.t[1]^(h)
  }
  
  
  p.t = p.fun(S.e, S.0)
  p.t.ex = p.fun(S, S.0)
  p.t[p.t <= 0] = 1/n
  p.t.ex[p.t.ex <= 0] = 1/n
  
  inf.dates = as.vector(apply(test.mat[1:n,], 1, decensor, p.t=p.t, S.0 = S.e.0))
  ex.dates = inf.dates - ex.geom
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
      accept = accept + 1
    } else {
      reject = reject + 1
    }
    
  }
  ex.geom = inf.dates - ex.dates
  
  S.t.samp[i,] = S.t
  S.samp[i,] = S
  I.samp[i,] = I
  delta.samp[i,] = delta
  gamma.samp[i,] = gamma
  psi.1.samp[i,] = psi.1
  eps.samp[i,] = eps
}



#####Prevalence
IDSAplot_open()
print("Output complete.")
