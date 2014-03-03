import Tkinter

class SOLO_gui(Tkinter.Frame):
    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
    
        #t = Tkinter.Toplevel()
        #t.wm_title('SOLO GUI')
        #t.grid()
        nrow = 0
        self.nodesLabel = Tkinter.Label(self,text='Nodes')
        self.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky='E')
        self.nodesEntry = Tkinter.Entry(self,width=6)
        self.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky='W')
        self.nodesEntry.insert(0,'5')
        self.trainingLabel = Tkinter.Label(self,text='Training')
        self.trainingLabel.grid(row=nrow,column=2,columnspan=1,sticky='E')
        self.trainingEntry = Tkinter.Entry(self,width=6)
        self.trainingEntry.grid(row=nrow,column=3,columnspan=1,sticky='W')
        self.trainingEntry.insert(0,'500')
        self.factorLabel = Tkinter.Label(self,text='Nda factor')
        self.factorLabel.grid(row=nrow,column=4,columnspan=1,sticky='E')
        self.factorEntry = Tkinter.Entry(self,width=6)
        self.factorEntry.grid(row=nrow,column=5,columnspan=1,sticky='W')
        self.factorEntry.insert(0,'5')
        
        nrow = nrow + 1
        self.learningrateLabel = Tkinter.Label(self,text='Learning')
        self.learningrateLabel.grid(row=nrow,column=2,columnspan=1,sticky='E')
        self.learningrateEntry = Tkinter.Entry(self,width=6)
        self.learningrateEntry.grid(row=nrow,column=3,columnspan=1,sticky='W')
        self.learningrateEntry.insert(0,'0.0')
        self.iterationsLabel = Tkinter.Label(self,text='Iterations')
        self.iterationsLabel.grid(row=nrow,column=4,columnspan=1,sticky='E')
        self.iterationsEntry = Tkinter.Entry(self,width=6)
        self.iterationsEntry.grid(row=nrow,column=5,columnspan=1,sticky='W')
        self.iterationsEntry.insert(0,'1')
        
        nrow = nrow + 1
        self.filestartLabel = Tkinter.Label(self,text='File start date')
        self.filestartLabel.grid(row=nrow,column=0,columnspan=3)
        self.fileendLabel = Tkinter.Label(self,text='File end date')
        self.fileendLabel.grid(row=nrow,column=3,columnspan=3)
        
        nrow = nrow + 1
        self.filestartValue = Tkinter.Label(self,text='2011-01-01 00:30')
        self.filestartValue.grid(row=nrow,column=0,columnspan=3)
        self.fileendValue = Tkinter.Label(self,text='2012-01-01 00:00')
        self.fileendValue.grid(row=nrow,column=3,columnspan=3)
        
        nrow = nrow + 1
        self.startLabel = Tkinter.Label(self, text='Start date (YYYY-MM-DD)')
        self.startLabel.grid(row=nrow,column=0,columnspan=3)
        self.startEntry = Tkinter.Entry(self)
        self.startEntry.grid(row=nrow,column=3,columnspan=3)
        
        nrow = nrow + 1
        self.endLabel = Tkinter.Label(self, text='End date   (YYYY-MM-DD)')
        self.endLabel.grid(row=nrow,column=0,columnspan=3)
        self.endEntry = Tkinter.Entry(self)
        self.endEntry.grid(row=nrow,column=3,columnspan=3)
        
        nrow = nrow + 1
        self.runButton = Tkinter.Button (self, text='Run SOLO', command=self.run_SOLO)
        self.runButton.grid(row=nrow,column=1,columnspan=2)
        self.doneButton = Tkinter.Button (self, text='End SOLO session', command=self.finished_solo)
        self.doneButton.grid(row=nrow,column=3,columnspan=2)
        #t.quitButton = Tkinter.Button (t, text='Quit', command=lambda:quit_solo(t,cf) )
        #t.quitButton.grid(row=4,column=3,columnspan=1)
        #self.wait_window(t)

    def finished_solo(self):
        self.quit()

    def run_SOLO(self):
        self.write_sofm_inf()
        self.write_solo_inf()
        self.write_seqsolo_inf()
        self.quit()

    def write_sofm_inf(self):
        f = open('../solo/inf/sofm_test.inf','w')
        f.write(str(self.nodesEntry.get())+'\n')
        f.write(str(self.trainingEntry.get())+'\n')
        f.write(str(20)+'\n')
        f.write(str(0.01)+'\n')
        f.write(str(1234)+'\n')
        f.write('solo/input/sofm_input.csv'+'\n')
        f.write('solo/output/sofm_1.out'+'\n')
        f.write('solo/output/sofm_2.out'+'\n')
        f.write('solo/output/sofm_3.out'+'\n')
        f.write('solo/output/sofm_4.out'+'\n')
        f.write(str(50)+'\n')
        f.close()

    def write_solo_inf(self):
        f = open('../solo/inf/solo_test.inf','w')
        f.write(str(self.nodesEntry.get())+'\n')
        f.write(str(self.factorEntry.get())+'\n')
        f.write('solo/output/sofm_4.out'+'\n')
        f.write('solo/input/solo_input.csv'+'\n')
        f.write('training'+'\n')
        f.write(str(5678)+'\n')
        f.write(str(0)+'\n')
        f.write('solo/output/eigenValue.out'+'\n')
        f.write('solo/output/eigenVector.out'+'\n')
        f.write('solo/output/accumErr.out'+'\n')
        f.write('solo/output/accumRR.out'+'\n')
        f.write('solo/output/trainProcess.out'+'\n')
        f.write('solo/output/freqTable.out'+'\n')
        f.write('solo/output/hidOutputWt.out'+'\n')
        f.write('solo/output/errorMap.out'+'\n')
        f.write('solo/output/finResult.out'+'\n')
        f.write('solo/output/trainWin.out'+'\n')
        f.write('solo/output/trainWout.out'+'\n')
        f.close()

    def write_seqsolo_inf(self):
        f = open('../solo/inf/seqsolo_test.inf','w')
        f.write(str(self.nodesEntry.get())+'\n')
        f.write(str(0)+'\n')
        f.write(str(self.learningrateEntry.get())+'\n')
        f.write(str(self.iterationsEntry.get())+'\n')
        f.write('solo/output/sofm_4.out'+'\n')
        f.write('solo/input/seqsolo_input.csv'+'\n')
        f.write('simulation'+'\n')
        f.write(str(9100)+'\n')
        f.write(str(0)+'\n')
        f.write('solo/output/eigenValue.out'+'\n')
        f.write('solo/output/eigenVector.out'+'\n')
        f.write('solo/output/trainWout.out'+'\n')
        f.write('solo/output/freqTable.out'+'\n')
        f.write('solo/output/errorMap.out'+'\n')
        f.write('solo/output/finResult.out'+'\n')
        f.write('solo/output/trainingRMSE.out'+'\n')
        f.write('solo/output/seqOut0.out'+'\n')
        f.write('solo/output/seqOut1.out'+'\n')
        f.write('solo/output/seqOut2.out'+'\n')
        f.write('solo/output/seqHidOutW.out'+'\n')
        f.write('solo/output/seqFreqMap.out'+'\n')
        f.write(str(-9999.0)+'\n')
        f.close()

if __name__ == "__main__":
    gui = SOLO_gui()
    main_title = 'SOLO GUI'
    gui.master.title(main_title)
    gui.mainloop()
    gui.master.destroy()