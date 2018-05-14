#!/Users/hendric/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.widgets import Slider, Button, RadioButtons

def lorenz_63(x, y, z, s=10, r=28, b=2.667):
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

def gaussian(x,mu,sigma):
    return 1.0/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

dx   =  0.0001
xmin = -6.0
xmax =  6.0
ymin =  0.0
ymax =  0.6

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
x = np.arange(xmin, xmax, dx)

mu0    = 0.0
sigma0 = 1.0

y = gaussian(x, 0.0, 1.0)
l, = plt.plot(x, y, lw=2, color='red')
plt.axis([xmin,xmax,ymin,ymax])

axcolor = 'lightgoldenrodyellow'
axfreq  = plt.axes([0.35, 0.6, 0.15, 0.03], facecolor=axcolor)
axamp   = plt.axes([0.35, 0.4, 0.15, 0.03], facecolor=axcolor)

smu    = Slider(axfreq, 'mu',    -1.0,  1.0, valinit=mu0)
ssigma =  Slider(axamp, 'sigma',  0.0, 10.0, valinit=sigma0)



def update(val):
    sigma = ssigma.val
    mu    = smu.val
    l.set_ydata(gaussian(x,mu,sigma))
    fig.canvas.draw_idle()
smu.on_changed(update)
ssigma.on_changed(update)

resetax = plt.axes([0.35, 0.2, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    smu.reset()
    ssigma.reset()
button.on_clicked(reset)

rax = plt.axes([0.35, 0.7, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


def colorfunc(label):
    l.set_color(label)
    fig.canvas.draw_idle()
radio.on_clicked(colorfunc)

plt.show()
