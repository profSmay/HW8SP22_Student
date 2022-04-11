import numpy as np
import PyQt5.QtWidgets as qtw
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from pathlib import Path

import sys
import os

# I built the gui as pump.ui and used pyuic5 to make pump.py
from pump import Ui_Form
from Pump_MVC import Pump_Controller

class PumpCurve_GUI_Class(Ui_Form, qtw.QWidget):  #class for PumpCurve_GUI inherits from Ui_Form & qtw.QWidget
    def __init__(self): #constructor for this class
        super().__init__() #runs constructor from parent classes
        self.setupUi(self) #executes the setupUi from Ui_Form inheritance
        self.AssignSignals() #connects 'clicked' signals from buttons to proper slots
        self.FilePath=os.getcwd() #gets the current working directory

        self.FileName=""

        #these are for matplotlib using pyqt
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(5, 3),tight_layout=True, frameon=True))
        self.ax = self.canvas.figure.add_subplot()
        #add the canvas to the form (can't be done directly in qtdesigner)
        self.GL_Output.addWidget(self.canvas,5,0,1,4)

        self.myPump=Pump_Controller() #create a pump controller object
        self.setViewWidgets() #pass along widgets to myPump for diaplay

        #show the widget
        self.show()

    def AssignSignals(self): #connect the signals from buttons to slots
        self.PB_Exit.clicked.connect(self.Exit)
        self.CMD_Open.clicked.connect(self.ReadAndCalculate)

    def setViewWidgets(self):
        """
        Here is where I pass a list of display widgets to the controller that get passed along to the view so that
        the controller can update the view internally.
        :return:
        """
        w=[self.LE_PumpName, self.LE_FlowUnits, self.LE_HeadUnits, self.LE_HeadCoefs, self.LE_EffCoefs, self.ax, self.canvas]
        self.myPump.setViewWidgets(w)

    def ReadAndCalculate(self):
        '''
        Should open a dialog box to search file system for target file.
        Then, open the file and read the data.
        Then, fill out the form with the data.
        :return:
        '''
        if self.OpenFile()==True:
            f1=open(self.FileName,'r')
            data=f1.readlines()
            f1.close()
            self.myPump.ImportFromFile(data)
            return True
        else:
            return False

    def OpenFile(self):
        '''
        This is the slot to open a dialog and search file system.
        return True
        :return: boolean for if the operation was successful.
        '''
        fname=qtw.QFileDialog.getOpenFileName(self, 'Open file', self.FilePath,'all files (*.*)')
        oTF=len(fname[0])>0
        if oTF:
            self.FileName=fname[0]
            self.FilePath=str(Path(fname[0]).parents[0])+'/'
            self.TE_Filename.setText(self.FileName)
        return oTF

    def Exit(self):
        qapp.exit()

if __name__=="__main__":
    qapp = qtw.QApplication(sys.argv)
    PumpCurve_GUI = PumpCurve_GUI_Class() #object of type PumpCurve_GUI_Class
    qapp.exec_()