import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QProgressBar, QPlainTextEdit, \
    QMenuBar, QMenu, QLabel, QTextEdit, QLineEdit, QGridLayout
from PyQt5.QtGui import QMovie
from PyQt5 import QtGui, QtCore
from PyQt5.QtCore import Qt, QTimer, QSize, QRect, QPoint,  QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import solver
import time
import numpy as np
import pyqtgraph as pglt


class CalcTask(QThread):
    calc_info = pyqtSignal()

    def run(self):
        solver.calculate(main)
        self.calc_info.emit()


class External(QThread):
    countChanged = pyqtSignal()
    def run(self):
        self.countChanged.emit()


def onclick(event):
    ind = event.ind[0]
    data = event.artist.get_offsets()
    xdata, ydata = data[ind,:]
    print ((xdata, ydata))


class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.setWindowTitle('Simulator 2phase oil&gas')
        self.setGeometry(100, 100, 850, 800)

        # a figure instance to plot on
        self.figure = plt.figure()
        
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `start calculation` method
        self.button = QPushButton('Start Calculations')
        self.button.clicked.connect(self.onclick1)
        #self.button.setFixedSize(100, 25)

        self.button1 = QPushButton('Stop Calculation')
        #self.button1.setFixedSize(100, 25)
        self.button1.clicked.connect(self.stopcalc)

        # Text area
        self.text = QPlainTextEdit()

        # show graphic time step label
        self.label1 = QLabel('Update graphic time step, s:', self)
        self.label1.setFont(QtGui.QFont("Arial", 12, QtGui.QFont.Bold))
        self.label1.setAlignment(QtCore.Qt.AlignRight)

        self.edittext1 = QLineEdit()
        #self.edittext1.setFixedSize(100, 25)
        self.edittext1.setText('20')

        self.stop = False

        # Progressbar
        self.progbar = QProgressBar()
        self.progbar.setValue(0)

        grid = QGridLayout()
        grid.setSpacing(20)

        # set the layout
        #layout = QVBoxLayout()
        grid.addWidget(self.progbar, 1, 0, 1, 2)
        grid.addWidget(self.button, 2, 0)
        grid.addWidget(self.button1, 2, 1)
        grid.addWidget(self.label1, 3, 0)
        grid.addWidget(self.edittext1, 3, 1)
        grid.addWidget(self.toolbar, 4, 0, 1, 2)
        grid.addWidget(self.canvas, 5, 0, 1, 2)
        grid.addWidget(self.text, 6, 0, 1, 2)

        self.setLayout(grid)

        #self.canvas.mpl_connect('pick_event', onclick)
        self.axes1 = self.figure.add_subplot(141)
        self.axes2 = self.figure.add_subplot(142)
        self.axes3 = self.figure.add_subplot(143)
        self.axes4 = self.figure.add_subplot(144)

        self.plot(0,0)
        self.plot1(0,0)
        self.plot2(0,0)
        self.plot_vel(0,0,0)
        self.load_params()

        #self.figure.tight_layout()


    def plot(self, x, y):
        # instead of ax.hold(False)
        
        self.axes1 = self.figure.add_subplot(141)
        # discards the old graph

        x = np.array(x)
        x = x / 101325
        self.axes1.set(title='Pressure, atm')
        self.axes1.plot(x, y, color='red')
        self.axes1.invert_yaxis()
        self.axes1.grid(color='#dbdbdb', linestyle='dashed')
        self.canvas.draw()


    def plot1(self, x, y):

        # discards the old graph
        self.axes2 = self.figure.add_subplot(142)
        self.axes2.set(title = 'Volume fraction oil phase')
        self.axes2.plot(x, y, label="1", color='black')
        self.axes2.invert_yaxis()
        self.axes2.grid(color='#dbdbdb', linestyle='dashed')
        self.canvas.draw()

    def plot2(self, x, y):
        # discards the old graph
        self.axes3 = self.figure.add_subplot(143)
        self.axes3.set(title='Temperature, K')
        self.axes3.plot(x, y, color='blue')
        self.axes3.invert_yaxis()
        self.axes3.grid(color='#dbdbdb', linestyle='dashed')
        self.canvas.draw()


    def plot_vel(self, v1, v2, z):
        # discards the old graph
        self.axes4 = self.figure.add_subplot(144)
        self.axes4.set(title='Velocities, m/s')
        self.axes4.plot(v1, z, color='green', label='oil')
        self.axes4.plot(v2, z, color='red', label='gas')
        self.axes4.legend()
        self.axes4.invert_yaxis()
        self.axes4.grid(color='#dbdbdb', linestyle='dashed')
        self.canvas.draw()


    def plot_well(self, ht, hb):

        # plt.grid(color='gray', linestyle='dashed')
        plast_cout = len(ht)

        for j in range(plast_cout):
            self.axes1.axhspan(float(ht[j]), float(hb[j]), facecolor='0.5', alpha=0.5)
            self.axes2.axhspan(float(ht[j]), float(hb[j]), facecolor='0.5', alpha=0.5)
            self.axes3.axhspan(float(ht[j]), float(hb[j]), facecolor='0.5', alpha=0.5)
            self.axes4.axhspan(float(ht[j]), float(hb[j]), facecolor='0.5', alpha=0.5)
        self.canvas.draw()


    def load_params(self):
        param_file = "param.txt"
        params = np.loadtxt(param_file, skiprows=20, unpack=True)
        #plt.figure(a, figsize=(6, 8), dpi=200)
        self.ht = params[0]
        self.hb = params[1]
        self.plot_well(self.ht, self.hb)


    def onclick1(self):
        self.task = CalcTask()
        self.task.start()


    def stopcalc(self):
        self.stop = True


    def onCountChanged(self, value):
        self.progbar.setValue(value)


if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    #main.tmr.start(1000)
    #main.resize(100, 50)
    main.show()

    sys.exit(app.exec_())
