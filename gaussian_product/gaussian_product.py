'''
Replication of gaussian_product.m from the DART_LAB using tkinter

Will Downs, May 2018
'''


'''
update structure is screwed up, probably need to change
'''

from tkinter import *
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt as sqrt
from scipy.stats import norm


def gaussian(mean, SD):
    y = norm(loc = mean, scale = SD)
    x = np.arange(mean - 5*SD, mean + 5*SD, .1)
    return x, y

def plot_initial():

    global canvas
    global fig
    global ax
    
    float_prior_mean = float(prior_mean.get())
    float_prior_SD = float(prior_SD.get())
    float_ob = float(ob.get())
    float_ob_SD = float(ob_SD.get())

    ax.set_title('gaussian_product')

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

    prior_x, prior_y = prior()
    obs_x, obs_y = obs()

    ax.legend(loc = 'upper left')
    ax.set_ylim(ymin = 0, ymax = max(max(prior_y), max(obs_y)))
    ax.set_xlim(xmin = min(min(prior_x), min(obs_x)),
                xmax = max(max(prior_x), max(obs_x)))
    fig.tight_layout()

    canvas = FigureCanvasTkAgg(fig, master=main_frame)
    canvas.get_tk_widget().grid(column = 1, row = 1)

    return fig, ax, float_prior_mean, float_prior_SD, float_ob, float_ob_SD, prior_x, prior_y, obs_x, obs_y

    

def plot_posterior():

    global canvas
    global fig
    global ax
    
    _, _, float_prior_mean, float_prior_SD, float_ob, float_ob_SD, prior_x, prior_y, obs_x, obs_y = plot_initial()

    float_post_SD = sqrt(pow(pow(float_prior_SD, -2) + pow(float_ob_SD, -2), -1))
    float_post_mean = pow(float_post_SD, 2)*(pow(float_prior_SD, -2)*float_prior_mean + pow(float_ob_SD, -2)*float_ob)

    def post():
        
        post_x, post_y = gaussian(float_post_mean, float_post_SD)
        post_y = post_y.pdf(post_x)
        plt.plot(post_x, post_y, label = 'Posterior', color = 'Blue')
        plt.plot(obs_x, np.multiply(prior_y, obs_y), label = 'Weighted Posterior',
                 ls = 'dashed', color = 'Blue')
        return post_x, post_y

    post_x, post_y = post()
    
    ax.legend(loc = 'upper left')
    ax.set_ylim(ymin = 0, ymax = max(max(prior_y), max(obs_y), max(post_y)))
    ax.set_xlim(xmin = min(min(prior_x), min(obs_x), min(post_x)),
                xmax = max(max(prior_x), max(obs_x), max(post_x)))
    fig.tight_layout()
    canvas = FigureCanvasTkAgg(fig, master=main_frame)
    canvas.get_tk_widget().grid(column = 1, row = 1)
    
root = Tk()

'''a mainframe'''

main_frame = ttk.Frame(root, padding = "8") #may need to change padding
main_frame.grid(column = 0, row = 0) #may need to add sticky parameter
main_frame.columnconfigure(0, weight = 1) #weights for whole grid
main_frame.rowconfigure(0, weight = 1) #weights for whole grid

prior_mean = StringVar()
prior_mean.set("1")
prior_SD = StringVar()
prior_SD.set("1")

ob = StringVar()
ob.set("1")
ob_SD = StringVar()
ob_SD.set("1")

fig, ax = plt.subplots()
plot_initial()
canvas = FigureCanvasTkAgg(fig, master=main_frame)
canvas.get_tk_widget().grid(column = 1, row = 1)

'''top right frame 2x2: prior entry and label widgets
entries should be entries, tied to prior mean and SD variables
fixed width widgets, with labels ~3x width of entries'''



prior_frame = ttk.Frame(main_frame, padding = "2") #may need to change padding
prior_frame.grid(column = 2, row = 1, sticky = E)

ttk.Label(prior_frame, text = "Prior Mean").grid(column = 1, row = 1, sticky = E)
ttk.Label(prior_frame, text = "Prior SD").grid(column = 1, row = 2, sticky = E)

prior_mean_entry = ttk.Entry(prior_frame, width = 5, textvariable = prior_mean, validate = "focusout", validatecommand = plot_initial) #may need to change padding
prior_mean_entry.grid(column = 2, row = 1, sticky = E)
prior_SD_entry = ttk.Entry(prior_frame, width = 5, textvariable = prior_SD, validate = "focusout", validatecommand = plot_initial) #may need to change padding
prior_SD_entry.grid(column = 2, row = 2, sticky = E)

'''middle right frame 2x2: obs entry and label widgets
entries are entries, tied to obs mean and SD variables
fixed pixel width widgets, with labels ~3x width of entries'''

ob_frame = ttk.Frame(main_frame, padding = "2") #may need to change padding
ob_frame.grid(column = 2, row = 2, sticky = E)

ttk.Label(ob_frame, text = "Observation").grid(column = 1, row = 1, sticky = E)
ttk.Label(ob_frame, text = "Obs. Error SD").grid(column = 1, row = 2, sticky = E)

ob_entry = ttk.Entry(ob_frame, width = 5, textvariable = ob, validate = "focusout", validatecommand = plot_initial) #may need to change padding
ob_entry.grid(column = 2, row = 1, sticky = E)
ob_SD_entry = ttk.Entry(ob_frame, width = 5, textvariable = ob_SD, validate = "focusout", validatecommand = plot_initial) #may need to change padding
ob_SD_entry.grid(column = 2, row = 2, sticky = E)

'''lower middle right button: plot posterior button makes plot'''

ttk.Button(main_frame, text = "Plot Posterior", command = plot_posterior).grid(column = 2, row = 3)

'''bottom right frame 4x2 (except Posterior label): posterior  labels
posterior label should be left weighted and large font
label: mean = label: postMean
label: SD = label: postSD
label: Weight = label: postWeight'''

post_mean = StringVar()
post_SD = StringVar()
post_weight = StringVar()

post_frame = ttk.Frame(main_frame, padding = "2") #may need to change padding
post_frame.grid(column = 2, row = 4, sticky = E)

ttk.Label(post_frame, text = "Posterior").grid(column = 1, row = 1, sticky = W) #Will need to change font size
ttk.Label(post_frame, text = "Mean = ").grid(column = 1, row = 2, sticky = E)
ttk.Label(post_frame, text = "SD = ").grid(column = 1, row = 3, sticky = E)
ttk.Label(post_frame, text = "Weight = ").grid(column = 1, row = 4, sticky = E)

ttk.Label(post_frame, textvariable = post_mean).grid(column = 2, row = 2, sticky = E)
ttk.Label(post_frame, textvariable = post_SD).grid(column = 2, row = 3, sticky = E)
ttk.Label(post_frame, textvariable = post_weight).grid(column = 2, row = 4, sticky = E)


#TODO: Weights for resizing            


root.mainloop()



