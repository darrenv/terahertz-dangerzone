# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 16:54:57 2015

@author: dvalovcin
"""

import numpy as np
import matplotlib.pylab as plt
import pyqtgraph as pg
from matplotlib import rcParams

rcParams['axes.labelsize'] = 32
rcParams['xtick.labelsize'] = 32
rcParams['ytick.labelsize'] = 32
rcParams['legend.fontsize'] = 32
rcParams['axes.titlesize'] = 32
rcParams['font.size'] = 32
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

rcParams['xtick.major.size'] = 14
rcParams['xtick.major.width'] = 3
rcParams['ytick.major.size'] = 14
rcParams['ytick.major.width'] = 3
#rcParams['text.usetex'] = True

#ax.spines["right"].set_visible(False)
#ax.spines["top"].set_visible(False)
#ax.spines["left"].set_linewidth(2)
#ax.spines["bottom"].set_linewidth(2)
#
#ax.yaxis.set_ticks_position('left')
#ax.xaxis.set_ticks_position('bottom')
#xticks = ax.xaxis.get_major_ticks()
#xticks[0].tick1On = False
#xticks[-1].tick1On = False
#fig.tight_layout(pad=0.1)


import os
import processing as proc


os.chdir('Z:\\Hunter Banks\\Data\\2015\\02-10 HSG attemp Loren 10 nm\\NIRfreq')

filelist = os.listdir('.')
filelist = [i for i in filelist if ('NIRlambda' in i and 'bad' not in i and 'CD' not in i)]
filelist = sorted(filelist, key=lambda s: int(s[9:-4]))

pg.mkQApp()

# Create remote process with a plot window
import pyqtgraph.multiprocess as mp
procee = mp.QtProcess()
rpg = procee._import('pyqtgraph')
win = rpg.PlotWidget()
c = rpg.PlotCurveItem()
b = rpg.ScatterPlotItem()
d = rpg.InfiniteLine()
win.addItem(b)
win.addItem(c)
win.addItem(d)
win.show()

wantedSB = 9

sb = {6:np.empty((0,2)), #dictionary for all the sidebands. Each value will be array of [NIRP, intensity]
      7:np.empty((0,2)),
      8:np.empty((0,2)), 
      9:np.empty((0,2)), 
      10:np.empty((0,2)), 
      11:np.empty((0,2)), 
      12:np.empty((0,2)), 
      13:np.empty((0,2)), 
      14:np.empty((0,2)), 
      15:np.empty((0,2)), 
      16:np.empty((0,2))}

      
for (i, name) in enumerate(filelist[0:50]):
    image = proc.EMCCDimage(name)
    print '\n\ni:{}. File: {}. NIR Frequency/WN: {:.2f}/{:.2f}\n'.format(i, name, image.NIR_lambda, 10000000./image.NIR_lambda)
    image.clean_and_subtract()
    image.integrate()
    image.append_wavelengths()
    image.guess_sidebands(mycutoff=4)
    image.fit_sidebands()
#    print image.sb_fit
    c.setData(image.hsg_data[:,0], image.hsg_data[:,1], _callSync='off')
    b.setData(image.sb_guess[:,0], image.sb_guess[:,1], _callSync='off')
    d.setValue(10000000./(10000000./image.NIR_lambda+42*wantedSB))
    win.setLabel('top', str(image.NIR_lambda))
    nirWN = 10000000./image.NIR_lambda
    
    for (idx, peak) in enumerate(image.sb_guess[:,0]):
        peakWN = 10000000./peak
        sbNumber = int(round((peakWN-nirWN)/42)) #21 was the FEL frequency (x2)
        if sbNumber == wantedSB:
            d.setValue(peak)
        print 'Sideband: {}'.format(sbNumber)
        sb[sbNumber] = np.append(sb[sbNumber], [[image.NIR_lambda, image.sb_guess[idx,1]]], axis=0)
    
    

    







