#!/Users/hendric/anaconda3/bin/python

import tkinter
from tkinter import ttk

class Adder(ttk.Frame):
   """The adders gui and functions"""
   def __init__(self, parent, *args, **kwargs):
      ttk.Frame.__init__(self, parent, *args, **kwargs)
      self.root = parent
      self.init_gui()

   # Inside the Adder class but not inside init_gui()
   def on_quit(self):
       """Exits program."""
       quit()
    
   def calculate(self):
       """Calculates the sum of the two inputted numbers."""
       num1 = int(self.num1_entry.get())
       num2 = int(self.num2_entry.get())
       num3 = num1 + num2
       self.answer_label['text'] = num3

   def init_gui(self):
       """Builds GUI."""
       self.root.title('Number Adder')
       self.root.option_add('*tearOff', 'FALSE')
 
       self.grid(column=0, row=0, sticky='nsew')
 
       self.menubar = tkinter.Menu(self.root)
 
       self.menu_file = tkinter.Menu(self.menubar)
       self.menu_file.add_command(label='Exit', command=self.on_quit)
 
       self.menu_edit = tkinter.Menu(self.menubar)
 
       self.menubar.add_cascade(menu=self.menu_file, label='File')
       self.menubar.add_cascade(menu=self.menu_edit, label='Edit')
 
       self.root.config(menu=self.menubar)
 
       self.num1_entry = ttk.Entry(self, width=5)
       self.num1_entry.grid(column=1, row = 2)
 
       self.num2_entry = ttk.Entry(self, width=5)
       self.num2_entry.grid(column=3, row=2)
 
       self.calc_button = ttk.Button(self, text='Calculate',
               command=self.calculate)
       self.calc_button.grid(column=0, row=3, columnspan=4)
 
       self.answer_frame = ttk.LabelFrame(self, text='Answer',
               height=100)
       self.answer_frame.grid(column=0, row=4, columnspan=4, sticky='nesw')
 
       self.answer_label = ttk.Label(self.answer_frame, text='')
       self.answer_label.grid(column=0, row=0)
 
       # Labels that remain constant throughout execution.
       ttk.Label(self, text='Number Adder').grid(column=0, row=0, columnspan=4)
       ttk.Label(self, text='Number one').grid(  column=0, row=2, sticky='w')
       ttk.Label(self, text='Number two').grid(  column=2, row=2, sticky='w')
 
       ttk.Separator(self, orient='horizontal').grid(column=0, row=1, columnspan=4, sticky='ew')
       
       ttk.Frame
 
       for child in self.winfo_children():
          child.grid_configure(padx=5, pady=5)

if __name__ == '__main__':
    root = tkinter.Tk()
    Adder(root)
    root.mainloop()
