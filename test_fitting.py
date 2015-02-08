# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 16:26:37 2015

@author: dreadnought

"Brevity required, prurience preferred"
"""

import numpy as np
import matplotlib.pyplot as plt
import processing as proc

file_list = ["nir892.txt"]

data_list = []
for elem in file_list:
    data_list.append(proc.EMCCDimage(elem))
    
#plt.close("all")
spectrum_list = []
for image in data_list:
    image.clean_and_subtract()
    image.integrate()
    image.append_wavelengths()
    image.guess_sidebands(mycutoff=3)
    image.fit_sidebands()
    print image.sb_fit
#    plt.plot(image.hsg_data[:,0], image.hsg_data[:,1])
#    plt.plot(image.sb_guess[:,0], image.sb_guess[:,1], 'o')

#plt.show()