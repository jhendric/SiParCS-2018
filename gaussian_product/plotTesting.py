import matplotlib.pyplot as plt
import numpy as np
from math import sqrt as sqrt
from scipy.stats import norm

str_prior_mean = "1"
float_prior_mean = float(str_prior_mean)
str_prior_SD = "1"
float_prior_SD = float(str_prior_SD)

str_ob = "1"
float_ob = float(str_ob)
str_ob_SD = "1"
float_ob_SD = float(str_ob_SD)

float_post_SD = sqrt(pow(pow(float_prior_SD, -2) + pow(float_ob_SD, -2), -1))
float_post_mean = pow(float_post_SD, 2)*(pow(float_prior_SD, -2)*float_prior_mean + pow(float_ob_SD, -2)*float_ob)

fig, ax = plt.subplots()
ax.set_title('gaussian_product')

def gaussian(mean, SD):
    y = norm(loc = mean, scale = SD)
    x = np.arange(mean - 5*SD, mean + 5*SD, .1)
    return x, y

def prior():
    prior_x, prior_y = gaussian(float_prior_mean, float_prior_SD)
    prior_y = prior_y.pdf(prior_x)
    plt.plot(prior_x, prior_y, label = 'Prior', color = 'Green')
    return prior_x, prior_y

def obs():
    obs_x, obs_y = gaussian(float_ob, float_ob_SD)
    obs_y = obs_y.pdf(obs_x)
    plt.plot(obs_x, obs_y, label = 'Obs. Likelihood', ls = 'dashed',
             color = 'Red')
    return obs_x, obs_y

def post():
    post_x, post_y = gaussian(float_post_mean, float_post_SD)
    post_y = post_y.pdf(post_x)
    plt.plot(post_x, post_y, label = 'Posterior', color = 'Blue')
    plt.plot(obs_x, np.multiply(prior_y, obs_y), label = 'Weighted Posterior',
             ls = 'dashed', color = 'Blue')
    return post_x, post_y

prior_x, prior_y = prior()
obs_x, obs_y = obs()
post_x, post_y = post()

ax.legend(loc = 'upper left')
ax.set_ylim(ymin = 0, ymax = max(max(prior_y), max(obs_y), max(post_y)))
ax.set_xlim(xmin = min(min(prior_x), min(obs_x), min(post_x)),
            xmax = max(max(prior_x), max(obs_x), max(post_x)))
fig.tight_layout()
plt.show()
