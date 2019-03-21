# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 03:54:05 2019

@author: michael Bogumil
"""
import sys, os
sys.path.append('D:/Google Drive/myRepo/MiLAMP Finder/')


import digestCompClass as MiLAMPFinder
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import Menu
from tkinter import filedialog
import pickle 

class MiLAMPFinderGUI():
    def __init__(self):
        
        # Create instants of MiLAMPFinder
        self.finder = MiLAMPFinder.MiLAMPFinder()

        # =============== Eventually remove ===============
        self.finder.setStart(43124000+1083)
        self.finder.setMethSite([43125347, 43125364, 43125372, 43125377, 43125383, 43125409, 43125411, 43125419]) 
        # =============== **** ===============
                
        # Create main window
        self.mainWin = tk.Tk()
        
        # Name main window
        self.mainWin.title("MiLAMP Finder")
        
        # Choose size of main window
        self.mainWin.resizable(False, False)
        
        # Setup menus for the main window
        self.setupMenus()
        
        # Setup widgets in the main widnow
        self.setupWidgets()
        
    def _quit(self):
        self.mainWin.quit()
        self.mainWin.destroy()
        exit()
        
    def setupMenus(self):
        # Create menu bar
        self.menuBar = Menu(self.mainWin)
        self.mainWin.config(menu=self.menuBar)
        
        # Add items to File menu
        self.fileMenu = Menu(self.menuBar, tearoff = 0)
        self.fileMenu.add_command(label='New')
        self.fileMenu.add_command(label='Open', command=self.loadState)
        self.fileMenu.add_command(label='Save', command=self.saveState)
        self.fileMenu.add_command(label='Quit', command=self._quit)
        self.menuBar.add_cascade(label="File", menu=self.fileMenu)
        
        # Add items to Settings menu
        self.settingsMenu = Menu(self.menuBar, tearoff = 0)
        self.settingsMenu.add_command(label='Temperatures')
        self.settingsMenu.add_command(label='Primer Spacing')
        self.menuBar.add_cascade(label="Settings", menu=self.settingsMenu)
        
        # Add items to Help menu
        self.helpMenu = Menu(self.menuBar, tearoff = 0)
        self.helpMenu.add_command(label='About')
        self.helpMenu.add_command(label='Report Bug')
        self.helpMenu.add_command(label='Help Doc')
        self.menuBar.add_cascade(label="Help", menu=self.helpMenu)
        
    def setupWidgets(self):
        # Setup sequence input box
        self.seqTextBox_w = 70
        self.seqTextBox_h = 20
        self.seqTextBox = scrolledtext.ScrolledText(self.mainWin, width=self.seqTextBox_w, height = self.seqTextBox_h, wrap=tk.WORD)
        self.seqTextBox.grid(column = 0, columnspan = 1,row = 0)
        
        # Setup methylation address input box
        self.methSiteTextBox_w = 10
        self.methSiteTextBox_h = 20
        self.methSiteTextBox = scrolledtext.ScrolledText(self.mainWin, width=self.methSiteTextBox_w, height = self.methSiteTextBox_h, wrap=tk.WORD)
        self.methSiteTextBox.grid(column = 2, columnspan = 1,row = 0)
        
        self.enzChosen = tk.StringVar()
        
        ttk.Label(self.mainWin, text = "Enzyme: ").grid(column = 0, row = 1)
        self.enzSel = ttk.Combobox(self.mainWin, width=40, textvariable=self.enzChosen,state='disabled')
        self.enzSel['values'] = ('Choose Enzyme')
        self.enzSel.grid(column = 0, row = 1)
        self.enzSel.current(0)
        
#        ttk.Label(self.mainWin, text = "Enter Foward Loop Sequence: ").grid(column = 0, row = 3)
#        
#        self.custSeqEntryFText = tk.StringVar()
#        self.custSeqEntryF = ttk.Entry(self.mainWin, width = 25, textvariable=self.custSeqEntryFText)
#        self.custSeqEntryF.grid(column = 1, row = 3)
#        
#        ttk.Label(self.mainWin, text = "Enter Backward Loop Sequence: ").grid(column = 0, row = 4)
#        
#        self.custSeqEntryBText = tk.StringVar()
#        self.custSeqEntryB = ttk.Entry(self.mainWin, width = 25, textvariable=self.custSeqEntryBText)
#        self.custSeqEntryB.grid(column = 1, row = 4)   
        
        self.findResEnzBut = ttk.Button(self.mainWin, text="Find Restrcition Enzymes", command=self.findResEnz)
        self.findResEnzBut.grid(column=1, row=0)
        
        self.chooseResEnzBut = ttk.Button(self.mainWin, text="Generate Primers", command=self.genPrim)
        self.chooseResEnzBut.grid(column=1, row=3)
        
        self.outTextBox_w = 70
        self.outTextBox_h = 20
        self.outTextBox = scrolledtext.ScrolledText(self.mainWin, width=self.outTextBox_w, height = self.outTextBox_h, wrap=tk.WORD)
        self.outTextBox.grid(column = 0, columnspan = 1,row = 3)
        
    
    def findResEnz(self):
        self.textSeq = self.seqTextBox.get('1.0', 'end-1c')
        self.textSeq = ''.join(self.textSeq.split())
        self.finder.loadSeqText(self.textSeq)
    
        check = self.finder.methSiteCheck()
        
        checkList = self.finder.checkSites
        
        print(checkList)
        
        self.finder.setResEnz()
        
        self.finder.numOfMethSitesUsed(-2)
        
        self.finder.findResEnz()
        
        self.resEnz = self.finder.getCurrentEnz()
        
        self.enzOpt = {str(self.resEnz[i][0]) + ', surrounding free base pairs: ' + str(self.resEnz[i][4]) : i for i in range(0,len(self.resEnz))}
        
        self.enzSel['values'] = list(self.enzOpt.keys())
        
        self.enzSel.state(['!disabled', 'readonly'])
    
    def saveState(self):
        if len(self.finder.foldername) > 0:
            foldername = self.finder.foldername
        else:
            foldername = "~/"
        
        self.finder.filename =  filedialog.asksaveasfilename(initialdir = 'D:/Google Drive/myRepo/MiLAMP Finder/',title = "Select file",filetypes = (("pickle","*.pickle"),("all files","*.*")))
        self.finder.filename = self.finder.filename.split('.')[0]+'.pickle'
        filehandler = open(self.finder.filename, 'wb') 
        self.finder.foldername = '/'.join(self.finder.filename.split('/')[0:-1]) + '/'
        pickle.dump(self.finder, filehandler)
        print (self.finder.foldername)
        print (self.finder.filename)
        
    def loadState(self):
        if len(self.finder.foldername) > 0:
            foldername = self.finder.foldername
        else:
            foldername = "~/"
        
        filename =  filedialog.askopenfilename(initialdir = foldername,title = "Select file",filetypes = (("pickle","*.pickle"),("all files","*.*")))
        
        filehandler = open(filename, 'rb')
        self.finder = pickle.load(filehandler)
        
        self.updateSettingsFromLoad()
        
    def updateSettingsFromLoad(self):
        # Set sequence in text box and variable
        self.seqTextBox.delete('1.0', tk.END)
        self.seqTextBox.insert(tk.END, self.finder.getSeqText())
        
        # Update sequence in text box
        # Update settings in settingWindow
        
    
    def genPrim(self):
        self.resEnz = self.finder.getCurrentEnz()
        self.enzOpt = {str(self.resEnz[i][0]) + ', surrounding free base pairs: ' + str(self.resEnz[i][4]) : i for i in range(0,len(self.resEnz))}
        self.finder.setEnz(self.enzOpt[self.enzChosen.get()])
        
        self.finder.setDefaultPrimerParam()
        
        self.finder.genPrimers()
        
        self.finder.genSinglePrimer(primerType='B2', startAdd = 136, endAdd = 202)
        
        self.finder.genSinglePrimer(primerType='B1c', startAdd = 222, endAdd = 318)
        
        #finder.B2Best = 2
        self.finder.something()
        print(self.finder.B2Best)
        
        self.finder.setCustomLoopPrimer(self.custSeqEntryFText.get(), self.custSeqEntryBText.get())
        
        
        
        self.finder.createXIP()    
        
        print("============================================================================")
        print("Restriction Enzyme Used: " + str(self.enzChosen.get()))
        #print("============================================================================")
        #finder.printAddress()
        print("============================================================================")
        self.finder.printPrimerComponentStats()
        print("============================================================================")
        self.finder.printPrimerStats()
        print("============================================================================")
        self.finder.printPrimerPotStats()
        print("============================================================================")



#=================
# Start Main Loop
#=================

if __name__ == '__main__':
    finderGUI = MiLAMPFinderGUI()
    finderGUI.mainWin.mainloop()
