# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""
import os
path = os.getcwd()
if not(os.environ.has_key('PYSYN_CDBS')):
    os.environ['PYSYN_CDBS']= path+'/CDBS'

from PyQt4.QtGui import QMainWindow, QFileDialog
from PyQt4.QtCore import pyqtSignature

from Ui_main import Ui_MainWindow
import star_phot

class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Class documentation goes here.
    """
    def __init__(self, parent = None):
        """
        Constructor
        """
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.Output.setPlainText('Set input and hit calculate.')
        self.statusbar.showMessage('Loading grid.npz ...')
        self.inmags = None
        self.model_mags = None
        self.pars = None
        self.figs = []
        self.mpl = None
        self.load('ui/grid.npz')
        
    def load(self, filename):
        grid,Teff,logZ,logg,Eb_v = star_phot.load_grid(filename)
        Eb_v = Eb_v[0:10] #Until I have made the whole entire cube
        self.TeffBox.clear()
        self.logZBox.clear()
        self.loggBox.clear()
        self.EbvBox.clear()
        self.TeffBox.addItems(['Teff - free'] + ['%f'%T for T in Teff])
        self.logZBox.addItems(['log Z - free'] +['%f'%Z for Z in logZ])
        self.loggBox.addItems(['log g - free']+['%f'%g for g in logg])
        self.EbvBox.addItems(['E(B-V) - free']+['%f'%e for e in Eb_v])
        self.data = (grid,Teff,logZ,logg,Eb_v)
        self.datafile = filename
        self.statusbar.showMessage('Loaded '+filename)
    
    @pyqtSignature("")
    def on_calculateButton_clicked(self):
        """
        Collect magnitude and constrains selections and pass to calculate routine
        """
        self.statusbar.showMessage('Calculating')
        Tsel, Zsel, gsel, Esel = self.TeffBox.currentIndex(), self.logZBox.currentIndex(), self.loggBox.currentIndex(), self.EbvBox.currentIndex()
        U, B, V, R = self.Ubox.value(), self.Bbox.value(), self.Vbox.value(), self.Rbox.value()
        I, J, H, K = self.Ibox.value(), self.Jbox.value(), self.Hbox.value(), self.Kbox.value()
        dU, dB, dV, dR = self.Ubox_2.value(), self.Bbox_2.value(), self.Vbox_2.value(), self.Rbox_2.value()
        dI, dJ, dH, dK = self.Ibox_2.value(), self.Jbox_2.value(), self.Hbox_2.value(), self.Kbox_2.value()
        grid,Teff,logZ,logg,Eb_v = self.data
        inmags = numpy.array([U, B, V, R, I, J, H, K])
        err = numpy.array([dU, dB, dV, dR, dI, dJ, dH, dK])
        T,Z,g,E, offset, chi, dmags, index = calculate(self.data, Tsel, Zsel, gsel, Esel, inmags,  err)
        model_mags = grid[Teff==T, logZ==Z, logg==g, Eb_v==E, :][0]
        text = "Input mags:\n %s\n"%inmags + "Residuals:\n %s\n"%(model_mags+offset-inmags)
        text += "Offset: %f   std.dev: %f   chi2: %f\n"%(offset,((inmags-model_mags-offset))[inmags>0].std() , chi)
        text += "Teff: %f    log Z: %f\n"%(T, Z)
        text += "log g: %f  E(B-V): %f"%(g, E)
        self.Output.setPlainText(text)
        self.inmags = inmags
        self.err = err
        self.model_mags = model_mags
        self.pars = (T, Z, g, E, offset)
        self.statusbar.showMessage('Calculation finished')
    
    @pyqtSignature("")
    def on_plotButton_clicked(self):
        """
        Plot input and output magnitudes alongside model spectrum, if available.
        Uses a pylab figure for minimum effort.
        """
        if self.mpl == 0 :
            self.statusbar.showMessage('Cannot plot')
            return
        try:
            import matplotlib
            if self.mpl is None:
                matplotlib.use('Qt4Agg')
                self.mpl=1
            import pylab
            fig = pylab.figure()
            self.figs.append(fig)
        except:
            self.statusbar.showMessage("Requires matplotlib with Qt4Agg backend.")
            self.mpl=0
            return
        T, Z, g, E, offset = self.pars
        try: 
            import pysynphot as S
            model = S.Icat('k93models', T,Z,g)
            model2= model*S.Extinction(E)
            model2.convert('Hz')
            model2.convert('fnu')
            w, f = model2.wave, model2.flux* 10**(-0.4*offset)
            ind = (w>1e14)*(w<1e15)
            pylab.loglog(w[ind], f[ind], 'b-')
            mags = self.inmags
            err = self.err
            if numpy.alltrue(err==1):
                err = 0
            zp = numpy.array([1810, 4260, 3640, 3080, 2550, 1600, 1080, 670])*1e-23
            l = numpy.array([0.36, 0.44, 0.55, 0.64, 0.79, 1.26, 1.6, 2.22])*1e-6
            nu = 299792458.0/l
            fnu =zp * 10**(-0.4*mags)
            efnu = (2.51**-err - 1)*fnu
            pylab.errorbar(nu, fnu, efnu, fmt='ko')
            pylab.xlabel('nu (Hz)')
            pylab.ylabel("fnu (erg/s/cm2/Hz)")
            pylab.title("Teff=%6.0f, log Z=%4.2f, log g=%4.2f, E(B-V)=%4.2f "%(T, Z, g, E))
        except:
            raise
        self.statusbar.showMessage('Plot generated')
        

    
    @pyqtSignature("")
    def on_clearButton_clicked(self):
        """
        Clear magnitude values, parameter choices and output text.
        """
        self.Output.clear()
        self.Ubox.setValue(25)
        self.Bbox.setValue(25)
        self.Vbox.setValue(25)
        self.Rbox.setValue(25)
        self.Ibox.setValue(25)
        self.Jbox.setValue(25)
        self.Hbox.setValue(25)
        self.Kbox.setValue(25)
        self.Ubox_2.setValue(1)
        self.Bbox_2.setValue(1)
        self.Vbox_2.setValue(1)
        self.Rbox_2.setValue(1)
        self.Ibox_2.setValue(1)
        self.Jbox_2.setValue(1)
        self.Hbox_2.setValue(1)
        self.Kbox_2.setValue(1)
        self.TeffBox.setCurrentIndex(0)
        self.logZBox.setCurrentIndex(0)
        self.loggBox.setCurrentIndex(0)
        self.EbvBox.setCurrentIndex(0)
        self.statusbar.showMessage('Cleared')
    
    @pyqtSignature("")
    def on_loadButton_clicked(self):
        myfile = str(QFileDialog.getOpenFileName(self,  "Open file",  os.getcwd(), "numpy gzip files files (*npz)"))
        if myfile:
            try:
                self.load(myfile)
            except:
                self.statusbar.showMessage('File Error')

import numpy

def calculate(data, Tsel, Zsel, gsel, Esel, mags,  err):
    mags[mags==0] = numpy.nan
    grid,Teff,logZ,logg,Eb_v = data
    shape = []
    if Tsel:
        indT = numpy.array([Tsel]).reshape(1)-1
        shape.append(1)
    else:
        indT = numpy.arange(len(Teff))
        shape.append(len(Teff))
    if Zsel:
        indZ = numpy.array([Zsel]).reshape(1)-1
        shape.append(1)
    else:
        indZ = numpy.arange(len(logZ))
        shape.append(len(logZ))
    if gsel:
        indg = numpy.array([gsel]).reshape(1)-1
        shape.append(1)
    else:
        indg = numpy.arange(len(logg))
        shape.append(len(logg))
    if Esel:
        indE = numpy.array([Esel]).reshape(1)-1
        shape.append(1)
    else:
        indE = numpy.arange(len(Eb_v))
        shape.append(len(Eb_v))
    shape.append(8)
    select = grid[indT][:, indZ][:, :, indg][:, :, :, indE]
    return star_phot.find_star(select,mags,err,Teff[indT],logZ[indZ],logg[indg],Eb_v[indE], all=1)
