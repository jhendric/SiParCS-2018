'''
Replication of gaussian_product.m from the DART_LAB using tkinter

Will Downs, May 2018
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
from math import pi as pi
from math import exp as exp
from scipy.stats import norm
from decimal import Decimal

class gaussian_product:

    def __init__(self, window, grid_col, grid_row):
        
        self.window = window
        self.window.grid_columnconfigure(0, weight = 1)
        self.window.grid_rowconfigure(0, weight = 1)
        
        #a mainframe
        self.main_frame = ttk.Frame(self.window, padding = "8")
        self.main_frame.grid(column = grid_col, row = grid_row, sticky = "N, S, E, W") 
        self.main_frame.grid_columnconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(0, weight = 1) #weights for whole grid

        self.style = ttk.Style()
        
        #prior variables
        self.prior_mean = StringVar()
        self.prior_mean.set("1")
        self.float_prior_mean = float(self.prior_mean.get())
        self.prior_SD = StringVar()
        self.prior_SD.set("1")
        self.float_prior_SD = float(self.prior_SD.get())
        
        #ob variables
        self.ob = StringVar()
        self.ob.set("1")
        self.float_ob = float(self.ob.get())
        self.ob_SD = StringVar()
        self.ob_SD.set("1")
        self.float_ob_SD = float(self.ob_SD.get())


        #post variables
        self.post_mean = StringVar()
        self.post_SD = StringVar()
        self.post_weight = StringVar()
        
        #plot initialization
        self.plot()

        '''top right frame 2x2: prior entry and label widgets
        entries should be entries, tied to prior mean and SD variables
        fixed width widgets, with labels ~3x width of entries'''

        self.style.configure("Prior.TFrame", background = 'green')
        self.prior_frame = ttk.Frame(self.main_frame, padding = "2", style = "Prior.TFrame")
        self.prior_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        
        self.style.configure("Prior.TLabel", background = 'green', foreground = "white")
        ttk.Label(self.prior_frame, text = "Prior Mean", style = "Prior.TLabel").grid(column = 1, row = 1, sticky = "N, S, E")
        ttk.Label(self.prior_frame, text = "Prior SD", style = "Prior.TLabel").grid(column = 1, row = 2, sticky = "N, S, E")
        
        self.style.configure("My.TEntry", background = 'white')
        self.prior_mean_entry = ttk.Entry(self.prior_frame, width = 5, textvariable = self.prior_mean
                                          , validate = "focusout", validatecommand = self.plot,
                                          style = "My.TEntry") 
        self.prior_mean_entry.grid(column = 2, row = 1, sticky = "N, S, E, W")
        self.prior_SD_entry = ttk.Entry(self.prior_frame, width = 5, textvariable = self.prior_SD
                                        , validate = "focusout", validatecommand = self.plot,
                                        style = "My.TEntry")
        self.prior_SD_entry.grid(column = 2, row = 2, sticky = "N, S, E, W")

        
        '''middle right frame 2x2: obs entry and label widgets
        entries are entries, tied to obs mean and SD variables
        fixed pixel width widgets, with labels ~3x width of entries'''

        self.style.configure("Ob.TFrame", background = '#c20144')
        self.ob_frame = ttk.Frame(self.main_frame, padding = "2", style = "Ob.TFrame")
        self.ob_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")

        self.style.configure("Ob.TLabel", background = '#c20144', foreground = "White")
        ttk.Label(self.ob_frame, text = "Observation", style = "Ob.TLabel").grid(column = 1, row = 1, sticky = "N, S, E")
        ttk.Label(self.ob_frame, text = "Obs. Error SD", style = "Ob.TLabel").grid(column = 1, row = 2, sticky = "N, S, E")

        self.ob_entry = ttk.Entry(self.ob_frame, width = 5, textvariable = self.ob
                                  , validate = "focusout", validatecommand = self.plot,
                                  style = "My.TEntry") 
        self.ob_entry.grid(column = 2, row = 1, sticky = "N, S, E, W")
        self.ob_SD_entry = ttk.Entry(self.ob_frame, width = 5, textvariable = self.ob_SD
                                     , validate = "focusout", validatecommand = self.plot,
                                     style = "My.TEntry") 
        self.ob_SD_entry.grid(column = 2, row = 2, sticky = "N, S, E, W")

        
        '''lower middle right button: plot posterior button makes plot'''
        
        self.style.configure("My.TButton", background = 'white')
        self.button = ttk.Button(self.main_frame, text = "Plot Posterior", style = "My.TButton",
                                 command = lambda : self.plot(True))
        self.button.grid(column = 2, row = 3, sticky = "E, W")
        
        '''bottom right frame 4x2 (except Posterior label): posterior  labels
        posterior label should be left weighted and large font
        label: mean = label: postMean
        label: SD = label: postSD
        label: Weight = label: postWeight'''

        
        self.style.configure("Post.TFrame", background = 'Blue', foreground = "White")
        self.post_frame = ttk.Frame(self.main_frame, padding = "2", style = "Post.TFrame")
        self.post_frame.grid(column = 2, row = 4, sticky = "N, S, E, W")


        self.style.configure("Post.TLabel", background = "Blue", foreground = "White")
        ttk.Label(self.post_frame, text = "Posterior", style = "Post.TLabel").grid(column = 1, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.post_frame, text = "Mean = ", style = "Post.TLabel").grid(column = 1, row = 2, sticky = "N, S, E")
        ttk.Label(self.post_frame, text = "SD = ", style = "Post.TLabel").grid(column = 1, row = 3, sticky = "N, S, E")
        ttk.Label(self.post_frame, text = "Weight = ", style = "Post.TLabel").grid(column = 1, row = 4, sticky = "N, S, E")

        ttk.Label(self.post_frame, textvariable = self.post_mean, style = "Post.TLabel").grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.post_frame, textvariable = self.post_SD, style = "Post.TLabel").grid(column = 2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.post_frame, textvariable = self.post_weight, style = "Post.TLabel").grid(column = 2, row = 4, sticky = "N, S, E, W")


        #resizing settings
        for i in range(1, 3):
            print('Here')
            self.main_frame.grid_columnconfigure(i, weight = 1)
            self.main_frame.grid_rowconfigure(i, weight = 1)
            self.prior_frame.grid_columnconfigure(i, weight = 1)
            self.prior_frame.grid_rowconfigure(i, weight = 1)
            self.ob_frame.grid_columnconfigure(i, weight = 1)
            self.ob_frame.grid_rowconfigure(i, weight = 1)
            self.post_frame.grid_columnconfigure(i, weight = 1)
            self.post_frame.grid_rowconfigure(i, weight = 1)
        for i in range(1, 5):
            self.main_frame.grid_rowconfigure(i, weight = 1)
            self.post_frame.grid_rowconfigure(i, weight = 1)
        
        
    def plot(self, posterior = False):

        
        
        def gaussian(mean, SD):
            
            y = norm(loc = mean, scale = SD)
            x = np.arange(mean - 5*SD, mean + 5*SD, .1)
            return x, y

        def plot_initial():
            
            self.float_prior_mean = float(self.prior_mean.get())
            self.float_prior_SD = float(self.prior_SD.get())
            self.float_ob = float(self.ob.get())
            self.float_ob_SD = float(self.ob_SD.get())

            ax.set_title('gaussian_product')

            def prior():
                prior_x, prior_y = gaussian(self.float_prior_mean, self.float_prior_SD)
                prior_y = prior_y.pdf(prior_x)
                ax.plot(prior_x, prior_y, label = 'Prior', color = 'Green')
                return prior_x, prior_y

            def obs():
                obs_x, obs_y = gaussian(self.float_ob, self.float_ob_SD)
                obs_y = obs_y.pdf(obs_x)
                obs = ax.plot(obs_x, obs_y, label = 'Obs. Likelihood', ls = 'dashed',
                         color = 'Red')
                return obs_x, obs_y

            prior_x, prior_y= prior()
            obs_x, obs_y= obs()
            return prior_x, prior_y, obs_x, obs_y



        def plot_posterior(prior_x, prior_y, obs_x, obs_y):
            
            float_post_SD = sqrt(pow(pow(self.float_prior_SD, -2) + pow(self.float_ob_SD, -2), -1))
            float_post_mean = pow(float_post_SD, 2)*(pow(self.float_prior_SD, -2)*
                                                          self.float_prior_mean + pow(self.float_ob_SD, -2)*self.float_ob)
            self.post_SD.set(str(round(Decimal(float_post_SD), 4)))
            self.post_mean.set(str(round(Decimal(float_post_mean), 4)))

            float_weight = ((1 / (sqrt(2 * pi) * sqrt(pow(self.float_prior_SD, 2) + pow(self.float_ob_SD, 2))))
                            * exp(-0.5 * pow((self.float_ob - self.float_prior_mean), 2) / (pow(self.float_prior_SD
                                                                                      , 2) + pow(self.float_ob_SD, 2))))
            self.post_weight.set(str(round(Decimal(float_weight), 4)))
            
            def post():

                post_x, post_y = gaussian(float_post_mean, float_post_SD)
                post_y = post_y.pdf(post_x)
                ax.plot(post_x, post_y, label = 'Posterior', color = 'Blue')
                ax.plot(post_x, np.multiply(float_weight, post_y), label = 'Weighted Posterior',
                         ls = 'dashed', color = 'Blue')
                return post_x, post_y
            
            post_x, post_y = post()
            return post_x, post_y
    
        
        fig = Figure(figsize = (6,6))
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 4, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)
        
        if not posterior:
            prior_x, prior_y, obs_x, obs_y = plot_initial()
            
            ax.legend(loc = 'upper left')
            ax.set_ylim(ymin = 0, ymax = max(max(prior_y), max(obs_y)))
            ax.set_xlim(xmin = min(min(prior_x), min(obs_x)),
                        xmax = max(max(prior_x), max(obs_x)))
            fig.tight_layout()
            self.post_SD.set('')
            self.post_mean.set('')
            self.post_weight.set('')
            
            
        else:
            
            prior_x, prior_y, obs_x, obs_y = plot_initial()
            
            post_x, post_y = plot_posterior(prior_x, prior_y, obs_x, obs_y)

            ax.legend(loc = 'upper left')
            ax.set_ylim(ymin = 0, ymax = max(max(prior_y), max(obs_y), max(post_y)))
            ax.set_xlim(xmin = min(min(prior_x), min(obs_x), min(post_x)),
                        xmax = max(max(prior_x), max(obs_x), max(post_x)))
            fig.tight_layout()
            
        return True


root = Tk()
style = ttk.Style()
style.theme_use('classic')
widg = gaussian_product(root, 0, 0)
root.mainloop()



