    #!/usr/bin/env python

from PyQt4 import QtCore, QtGui
from ui import main
from ui.main import MainWindow
import sys

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    out = app.exec_()
    sys.exit(out)

import os
path = os.getcwd()
os.environ['PYSYN_CDBS']= path+'/CDBS'

import ui.pysynphot
