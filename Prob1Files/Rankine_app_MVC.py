import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Rankine_GUI import Ui_Form
from Rankine_Classes import rankineController
from Calc_state import SatPropsIsobar
from Calc_state import UnitConverter as UC
#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        MainWindow constructor
        """
        super().__init__()  #if you inherit, you generally should run the parent constructor first.
        # Main UI code goes here
        self.setupUi(self)
        self.AssignSlots()
        self.MakeCanvas()
        self.setNewPHigh()
        self.setNewPLow()
        #A tuple containing the widgets that get updated in the View
        self.widgets = (self.le_H1, self.le_H2, self.le_H3, self.le_H4, self.le_TurbineWork, self.le_PumpWork, self.le_HeatAdded, self.le_Efficiency, self.lbl_SatPropHigh, self.lbl_SatPropLow, self.ax, self.canvas)
        self.otherwidgets = (self.lbl_TurbineInletCondition, self.rdo_Quality,self.le_PHigh, self.le_PLow, self.le_TurbineInletCondition, self.lbl_PHigh, self.lbl_PLow, self.lbl_H1Units, self.lbl_H2Units, self.lbl_H3Units, self.lbl_H4Units, self.lbl_TurbineWorkUnits, self.lbl_PumpWorkUnits, self.lbl_HeatAddedUnits)
        self.RC=rankineController(self.widgets)  # instantiate a rankineController object
        self.Calculate()  #calculate using initial values

        # a place to store coordinates from last position on graph
        self.oldXData=0.0
        self.oldYData=0.0
        # End main ui code
        self.show()

    def AssignSlots(self):
        """
        Setup signals and slots for my program
        :return:
        """
        self.btn_Calculate.clicked.connect(self.Calculate)
        self.rdo_Quality.clicked.connect(self.SelectQualityOrTHigh)
        self.rdo_THigh.clicked.connect(self.SelectQualityOrTHigh)
        self.le_PHigh.editingFinished.connect(self.setNewPHigh)
        self.le_PLow.editingFinished.connect(self.setNewPLow)
        self.rb_SI.clicked.connect(self.SetUnits)
        self.rb_English.clicked.connect(self.SetUnits)

    def MakeCanvas(self):
        """
        Create a place to make graph of Rankine cycle
        Step 1:  create a Figure object called self.figure
        Step 2:  create a FigureCanvasQTAgg object called self.canvas
        Step 3:  create an axes object for making plot
        Step 4:  add self.canvas to self.widgetsVerticalLayout which is a Vertical box layout
        :return:
        """
        #Step 1.
        self.figure=Figure(figsize=(1,1),tight_layout=True, frameon=True)
        #Step 2.
        self.canvas=FigureCanvasQTAgg(self.figure)
        #Step 3.
        self.ax = self.figure.add_subplot()
        #Step 4.
        self.widgetsVerticalLayout.addWidget(NavigationToolbar2QT(self.canvas,self))
        self.widgetsVerticalLayout.addWidget(self.canvas)
        #Step 5. attach an event handler for mouse movement on graph
        self.canvas.mpl_connect("motion_notify_event", self.mouseMoveEvent_Canvas)

    #since my main window is a widget, I can customize its events by overriding the default event
    def mouseMoveEvent_Canvas(self, event):
        self.oldXData=event.xdata if event.xdata is not None else self.oldXData
        self.oldYData=event.ydata if event.ydata is not None else self.oldYData
        sUnit='kJ/(kg*K)' if self.rb_SI.isChecked() else 'BTU/(lb*R)'
        TUnit='C' if self.rb_SI.isChecked() else 'F'
        self.setWindowTitle('s:{:0.2f} {}, T:{:0.2f} {}'.format(self.oldXData,sUnit, self.oldYData,TUnit))

    def Calculate(self):
        #use rankineController to update the model based on user inputs
        self.RC.updateModel((self.le_PHigh, self.le_PLow, self.rdo_Quality, self.le_TurbineInletCondition, self.le_TurbineEff))

    def SelectQualityOrTHigh(self):
        """
        Action to take when selecting one of the radio buttons for Quality or THigh
        :return:
        """
        #region Code for P1.1
        if self.rdo_Quality.isChecked():
            self.le_TurbineInletCondition.setText("1.0")
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            SI = self.rb_SI.isChecked()
            PCF = 100 if SI else UC.psi_to_kpa
            Tsat=SatPropsIsobar(float(self.le_PHigh.text())*PCF).TSat
            Tsat=Tsat if SI else UC.C_to_F(Tsat)
            self.le_TurbineInletCondition.setText("{:0.2f}".format(Tsat))
            self.le_TurbineInletCondition.setEnabled(True)
        #endregion
        x=self.rdo_Quality.isChecked()
        self.lbl_TurbineInletCondition.setText(("Turbine Inlet: {}{} =".format('x'if x else 'THigh', '' if x else ('(C)' if SI else '(F)'))))

    def SetUnits(self):
        self.RC.updateUnits(self.widgets, self.otherwidgets, SI=self.rb_SI.isChecked())

    def setNewPHigh(self):
        SI=self.rb_SI.isChecked()
        PCF=100 if SI else UC.psi_to_kpa
        self.lbl_SatPropHigh.setText(SatPropsIsobar(float(self.le_PHigh.text())*PCF,SI=SI).txtOut)
        self.SelectQualityOrTHigh()

    def setNewPLow(self):
        SI=self.rb_SI.isChecked()
        PCF=100 if SI else UC.psi_to_kpa
        self.lbl_SatPropHigh.setText(SatPropsIsobar(float(self.le_PLow.text())*PCF,SI=SI).txtOut)


#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine calculator')
    sys.exit(app.exec())