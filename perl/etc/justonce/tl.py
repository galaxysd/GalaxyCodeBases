#!/usr/bin/env python3

import numpy as np

x = np.array([ 23,23,25,27,31,31,32,33,35,36,38,38,39,42 ])
y = np.array([ 0.58,0.35,0.26,0.48,0.21,0.19,0.69,0.65,0.09,0.14,0.3,0.22,0.47,0.26 ])
e = np.array([ 0.096794787,0.075415516,0.070237692,0.052086884,0.043921177,0.047927184,0.115623311,0.068138514,0.0379057,0.043373379,0.079772404,0.068101673,0.056511855,0.067682733 ])
sampleid = '11044888_F'

# e = sqrt( p(1-p)/(n-1) )
import matplotlib
matplotlib.use('pdf')

import matplotlib.pyplot as plt

from scipy import optimize

def squared_loss(theta, x=x, y=y, e=e):
    dy = y - theta[0] - theta[1] * x
    return np.sum(0.5 * (dy / e) ** 2)

theta1 = optimize.fmin(squared_loss, [np.mean(x), np.mean(y)], disp=False)

t = np.linspace(-20, 20)

def huber_loss(t, c=3):
    return ((abs(t) < c) * 0.5 * t ** 2
            + (abs(t) >= c) * -c * (0.5 * c - abs(t)))

def total_huber_loss(theta, x=x, y=y, e=e, c=3):
    return huber_loss((y - theta[0] - theta[1] * x) / e, c).sum()

theta2 = optimize.fmin(total_huber_loss, [np.mean(x), np.mean(y)], disp=False)

xfit = np.linspace(10, 70)

# theta will be an array of length 2 + N, where N is the number of points
# theta[0] is the intercept, theta[1] is the slope,
# and theta[2 + i] is the weight g_i

def log_prior(theta):
    #g_i needs to be between 0 and 1
    if (all(theta[2:] > 0) and all(theta[2:] < 1)):
        return 0
    else:
        return -np.inf  # recall log(0) = -inf

def log_likelihood(theta, x, y, e, sigma_B):
    dy = y - theta[0] - theta[1] * x
    g = np.clip(theta[2:], 0, 1)  # g<0 or g>1 leads to NaNs in logarithm
    logL1 = np.log(g) - 0.5 * np.log(2 * np.pi * e ** 2) - 0.5 * (dy / e) ** 2
    logL2 = np.log(1 - g) - 0.5 * np.log(2 * np.pi * sigma_B ** 2) - 0.5 * (dy / sigma_B) ** 2
    return np.sum(np.logaddexp(logL1, logL2))

def log_posterior(theta, x, y, e, sigma_B):
    return log_prior(theta) + log_likelihood(theta, x, y, e, sigma_B)

ndim = 2 + len(x)  # number of parameters in the model
nwalkers = 50  # number of MCMC walkers
nburn  = 100000  # "burn-in" period to let chains stabilize
nsteps = 150000  # number of MCMC steps to take

# set theta near the maximum likelihood, with 
np.random.seed(0)
starting_guesses = np.zeros((nwalkers, ndim))
starting_guesses[:, :2] = np.random.normal(theta1, 1, (nwalkers, 2))
starting_guesses[:, 2:] = np.random.normal(0.5, 0.1, (nwalkers, ndim - 2))

import emcee
import multiprocessing as mp
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[x, y, e, 50], threads=mp.cpu_count() )
sampler.run_mcmc(starting_guesses, nsteps)

sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
sample = sampler.chain[:, nburn:, :].reshape(-1, ndim)

theta3 = np.mean(sample[:, :2], axis=0)
g = np.mean(sample[:, 2:], 0)
outliers = (g < 0.5)
#plt.show()

plt.errorbar(x, y, e, fmt='.k', ecolor='gray')
plt.plot(xfit, theta1[0] + theta1[1] * xfit, color='gray',label="y = %sx + %s"%(theta1[1],theta1[0]) )
plt.plot(xfit, theta2[0] + theta2[1] * xfit, color='green',label="Huber: y = %sx + %s"%(theta2[1],theta2[0]) )
plt.plot(xfit, theta3[0] + theta3[1] * xfit, color='navy',label="MCMC: y = %sx + %s"%(theta3[1],theta3[0]) )
plt.plot(x[outliers], y[outliers], 'ro', ms=20, mfc='none', mec='red')
plt.title("%s"%(sampleid));
plt.legend(loc='best', frameon=False);
#plt.show()
plt.savefig("%s.pdf"%sampleid)
