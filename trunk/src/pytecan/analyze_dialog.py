from util import CollectData, WriteCSV

################################################################################

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox

class TkFileDialogExample(Tkinter.Frame):
    def __init__(self, root):
        self.MES = None
        self.N_ROWS = 8
        self.N_COLS = 12
        
        Tkinter.Frame.__init__(self, root)
    
        self.button_input = Tkinter.Button(self, text='select result file ...', command=self.AskOpenFile)
        self.button_input.grid(row=0, column=1, columnspan=2, pady=15)

        self.button_plate = Tkinter.Button(self, text='select plate: ', command=self.SelectPlate)
        self.button_plate.grid(row=0, column=3, columnspan=2, sticky=Tkconstants.E)
        self.list_plate_ids = Tkinter.Listbox(self, height=3)
        self.list_plate_ids.grid(row=0, column=5, columnspan=2, sticky=Tkconstants.W)
        
        self.label_label = Tkinter.Label(self, text='select reading label: ')
        self.label_label.grid(row=0, column=7, columnspan=2, sticky=Tkconstants.E)
        self.list_labels = Tkinter.Listbox(self, height=3)
        self.list_labels.grid(row=0, column=9, columnspan=2, sticky=Tkconstants.W)

        
        for c in xrange(self.N_COLS):
            Tkinter.Label(self, text='%d' % (c+1)).grid(row=1, column=c+1)
        for r in xrange(self.N_ROWS):
            Tkinter.Label(self, text=chr(ord('A') + r)).grid(row=r+1, column=0, 
                                                             sticky=Tkconstants.E, padx=10)

        self.text_matrix = {}
        for r in xrange(self.N_ROWS):
            for c in xrange(self.N_COLS):
                self.text_matrix[r,c] = Tkinter.Entry(self, width=10)
                self.text_matrix[r,c].grid(row=r+1, column=c+1, padx=3, pady=3)

        self.button_savepdf = Tkinter.Button(self, text='save to PDF ...',
                                             command=self.asksaveaspdf,
                                             state=Tkconstants.DISABLED)
        self.button_savepdf.grid(row=self.N_ROWS+1, column=1, columnspan=2, pady=5)

        self.button_savecsv = Tkinter.Button(self, text='save to CSV ...',
                                             command=self.asksaveascsv,
                                             state=Tkconstants.DISABLED)
        self.button_savecsv.grid(row=self.N_ROWS+1, column=3, columnspan=2, pady=5)

        self.button_savecsv = Tkinter.Button(self, text='Quit',
                                             command=self.quit)
        self.button_savecsv.grid(row=self.N_ROWS+1, column=5, columnspan=2, pady=5)
        
    def AskOpenFile(self):
        """
            Opens and parses the robot result file (.tar or .tar.gz)
        """
        tar_fname = tkFileDialog.askopenfilename(
            filetypes=[('gzip files', '.gz'), ('tar files', '.tar'), ('all files', '.*')])
        if tar_fname:
            print "Collecting data from %s ..." % (tar_fname),
            self.MES = CollectData(tar_fname)
            print "[DONE]"
            if self.MES:
                self.button_savepdf.configure(state=Tkconstants.ACTIVE)
                self.button_savecsv.configure(state=Tkconstants.ACTIVE)
                
                # populate the plate ID list
                self.list_plate_ids.delete(0, Tkconstants.END)
                for item in sorted(self.MES.keys()):
                    self.list_plate_ids.insert(Tkconstants.END, item)

    def SelectPlate(self):
        selection = int(self.list_plate_ids.curselection()[0])
        self.selected_plate_id = self.list_plate_ids.get(selection)
        
        # populate the reading label list
        self.list_labels.delete(0, Tkconstants.END)
        for item in sorted(self.MES[self.selected_plate_id].keys()):
            self.list_labels.insert(Tkconstants.END, item)
            
    def GetLabel(self):
        selection = int(self.list_labels.curselection()[0])
        return self.list_labels.get(selection)
    
    def asksaveaspdf(self):
        """
            Writes the plots to a PDF file
        """
        f = tkFileDialog.asksaveasfile(mode='w', filetypes=[('pdf file', '.pdf')])
        f.write("<enter plots here>\n")
        f.write("plate ID = %s\n" % self.selected_plate_id)
        f.write("reading label = %s\n" % self.GetLabel())
        for r in xrange(self.N_ROWS):
            for c in xrange(self.N_COLS):
                f.write("r=%d, c=%d, text=%s\n" % (r, c, self.text_matrix[r,c].get()))
        f.close()

    def asksaveascsv(self):
        """
            Writes the raw data to a CSV file
        """
        f = tkFileDialog.asksaveasfile(mode='w')
        WriteCSV(self.MES, f)
        f.close()

################################################################################

if __name__ == "__main__":
    root = Tkinter.Tk()
    TkFileDialogExample(root).pack()
    root.mainloop()
