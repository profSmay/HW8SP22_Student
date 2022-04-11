import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import PyQt5.QtWidgets as qtw

# importing from previous work on least squares fit
from LeastSquares import LeastSquaresFit_Class

class Pump_Model():
    """
    This is the pump model.  It just stores data.
    """
    def __init__(self): #pump class constructor
        #create some class variables for storing information
        self.PumpName = ""
        self.FlowUnits = ""
        self.HeadUnits = ""

        # place to store data from file
        self.FlowData = np.array([])
        self.HeadData = np.array([])
        self.EffData = np.array([])

        # place to store coefficients for cubic fits
        self.HeadCoefficients = np.array([])
        self.EfficiencyCoefficients = np.array([])

        # create two instances (objects) of least squares class
        self.LSFitHead=LeastSquaresFit_Class()
        self.LSFitEff=LeastSquaresFit_Class()

class Pump_Controller():
    def __init__(self):
        self.Model=Pump_Model()
        self.View =Pump_View()
    
    #region functions to modify data of the model
    def ImportFromFile(self, data):
        """
        This processes the list of strings in data to build the pump model
        :param data: 
        :return: 
        """
        self.Model.PumpName = data[0]
        #data[1] is the units line
        L=data[2].split()
        self.Model.FlowUnits = L[0].strip()
        self.Model.HeadUnits = L[1].strip()

        # extracts flow, head and efficiency data and calculates coefficients
        self.SetData(data[3:])
        self.updateView()
    
    def SetData(self,data):
        '''
        Expects three columns of data in an array of strings with space delimiter
        Parse line and build arrays.
        :param data:
        :return:
        '''
        #erase existing data
        self.Model.FlowData = np.array([])
        self.Model.HeadData = np.array([])
        self.Model.EffData = np.array([])

        #parse new data
        for L in data:
            Cells=L.split() #parse the line into an array of strings
            self.Model.FlowData=np.append(self.Model.FlowData, float(Cells[0].strip())) #remove any spaces and convert string to a float
            self.Model.HeadData=np.append(self.Model.HeadData, float(Cells[1].strip())) #remove any spaces and convert string to a float
            self.Model.EffData=np.append(self.Model.EffData, float(Cells[2].strip())) #remove any spaces and convert string to a float

        #call least square fit for head and efficiency
        self.LSFit()
        
    def LSFit(self):
        '''Fit cubic polynomial using Least Squares'''
        self.Model.LSFitHead.x=self.Model.FlowData
        self.Model.LSFitHead.y=self.Model.HeadData
        self.Model.LSFitHead.LeastSquares(3) #calls LeastSquares function of LSFitHead object

        self.Model.LSFitEff.x=self.Model.FlowData
        self.Model.LSFitEff.y=self.Model.EffData
        self.Model.LSFitEff.LeastSquares(3) #calls LeastSquares function of LSFitEff object
        # lsfit=np.poly1d(np.polyfit(self.Model.LSFitHead.x, self.Model.LSFitHead.y,3))
        # ymodel=lsfit(np.linspace(0,500,500))
        pass
    #endregion

    #region functions interacting with view
    def setViewWidgets(self, w):
        self.View.setViewWidgets(w)

    def updateView(self):
        self.View.updateView(self.Model)
    #endregion

class Pump_View():
    def __init__(self):
        """
        In this constructor, I create some QWidgets as placeholders until they get defined later.
        """
        self.LE_PumpName=qtw.QLineEdit()
        self.LE_FlowUnits=qtw.QLineEdit()
        self.LE_HeadUnits=qtw.QLineEdit()
        self.LE_HeadCoefs=qtw.QLineEdit()
        self.LE_EffCoefs=qtw.QLineEdit()
        self.ax=None
        self.canvas=None

    def updateView(self, Model):
        """
        Put model parameters in the widgets.
        :param Model:
        :return:
        """
        self.LE_PumpName.setText(Model.PumpName)
        self.LE_FlowUnits.setText(Model.FlowUnits)
        self.LE_HeadUnits.setText(Model.HeadUnits)
        self.LE_HeadCoefs.setText(Model.LSFitHead.GetCoeffsString())
        self.LE_EffCoefs.setText(Model.LSFitEff.GetCoeffsString())
        self.DoPlot(Model)

    def DoPlot(self, Model):
        """
        Create the plot.
        :param Model:
        :return:
        """
        headx, heady, headRSq = Model.LSFitHead.GetPlotInfo(3, npoints=500)
        effx, effy, effRSq = Model.LSFitEff.GetPlotInfo(3, npoints=500)

        axes = self.ax
        axes2 = axes.twinx()

        axes.clear()
        axes2.clear()

        axes.plot(headx, heady, linestyle='dashed', color='black', linewidth='2',
                  label="Head" + r'($R^2={:0.3f}$)'.format(headRSq))
        axes2.plot(effx, effy, linestyle='dotted', color='black', linewidth='2',
                   label="Efficiency" + r'($R^2={:0.3f}$)'.format(effRSq))
        axes.plot(Model.FlowData, Model.HeadData, linestyle='none', marker='o', markerfacecolor='white',
                  markeredgecolor='black', markersize=10, label="Head")
        axes2.plot(Model.FlowData, Model.EffData, label="Efficiency", linestyle='none', marker='^',
                   markerfacecolor='white',
                   markeredgecolor='black', markersize=10)

        axes.tick_params(axis='both', direction='in', grid_linewidth=1, grid_linestyle='dashed', grid_alpha=0.5)
        axes.set(xlabel='Flow Rate (' + Model.FlowUnits + ')', ylabel='Head (' + Model.HeadUnits + ')',
                 title="Pump Head and Efficiency Curves")
        axes.legend(loc='center left', fontsize='large')

        axes2.tick_params(axis='both', direction='in', grid_linewidth=1, grid_linestyle='dashed', grid_alpha=0.5)
        axes2.set(ylabel='Efficiency (%)')
        axes2.legend(loc='upper right')
        self.canvas.draw()

    def setViewWidgets(self, w):
        self.LE_PumpName, self.LE_FlowUnits, self.LE_HeadUnits, self.LE_HeadCoefs, self.LE_EffCoefs, self.ax, self.canvas = w

