#!/Users/hendric/anaconda3/bin/python


import tkinter
from tkinter import ttk, Button, BOTH

class Window(ttk.Frame):
    
    def __init__(self, master=None):
        ttk.Frame.__init__(self, master)
        
        self.master = master
        self.init_window()
    
    # Creation of init_window
    def init_window(self):
        # title of the master widget
        self.master.title("My GUI")
        
        # allowing the widget to take the full space of the root window space
        self.pack(fill=BOTH, expand=1)
        
        # creating a button instance
        quitButton = Button(self, text="Abort")
        quitButton.place(x=0,y=0)

root = tkinter.Tk()
root.geometry("400x300")

app = Window(root)

root.mainloop()