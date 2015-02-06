# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 14:01:23 2015

@author: dreadnought
"""

from __future__ import division
import os, errno
import copy
import numpy as np
import matplotlib.pyplot as plt
import cosmics_hsg as cosmics

class EMCCDimage(object):
    
    def __init__(self, file_name):
        '''
        This does the load text file for you, and reads the header to put the
        appropriate values in the correct place
        It will try to open an already-cleaned version of the background image,
        called <BG_image_name>_clean.txt
        
        '''
        self.fname = file_name
        
        f = open(file_name,'rU')
            
        self.description = f.readline()
        self.description = self.description[1:-3]
        self.parameters = f.readline()
        self.parameters = self.parameters[1:-2] # I should probably delete the comma from LabView
        self.parameters = self.parameters.split(',')
        self.info = f.readline()
        self.info = self.info[1:-1] #
        self.info = self.info.split(',')
        f.close()
        
        # The following come from the second line of the file header:
        self.NIR_power = float(self.parameters[0][5:])
        self.NIR_lambda = float(self.parameters[1][10:])
        self.THz_power = float(self.parameters[2][5:])
        self.THz_freq = float(self.parameters[3][5:]) # Use this to guess at 
                                                      # the right frequency  
                                                      # that is calculated 
                                                      # later from the 
                                                      # sidebands?
        self.temp = None #Need to add this, too.              
        
        # The following come from the third line of the file header:
        rep_rate = float(self.info[0][8:]) # This needs an equals sign!
        exposure = float(self.info[1][8:])
        self.num_FEL = exposure * rep_rate # Get number of FEL shots from rep 
                                           # rate and exposure, but be aware 
                                           # that this does not account for 
                                           # missed pulses!
        self.gain = int(self.info[2][5:])
        self.center = float(self.info[3][14:])
        self.grating = int(self.info[4][8])
        self.background_name = self.info[5][9:] # Just make this value the actual string!
        self.series = self.info[6][7:] # This is to label files that should be summed later
        self.y_min = int(self.info[7][5:])
        self.y_max = int(self.info[8][5:]) # Why do I need a final comma to make this read correctly sometimes?
        self.slits = int(self.info[9][7:]) # I think this may be important for some linewidth studies
        self.CCD_temp = int(self.info[10][9:]) # I think this will be useful early during some noise measurements
        
        # Initializing these attributes for later.  Not sure if necessary.
        self.bg_array = None
        self.bg_array_clean = None
        
        # These come from the various files
        try:
            bg_clean_fname = self.background_name[:-4] + '_clean.txt'
            self.bg_array_clean = np.genfromtxt(bg_clean_fname, comments='#', delimiter=',')
            self.bg_fresh = False
            print "Imported clean array"
        except IOError:
            print "Imported dirty array"
            self.bg_array = np.genfromtxt(self.background_name, comments='#', delimiter=',')
            self.bg_fresh = True
            self.clean_up(self.bg_array)
            self.bg_array_clean = np.array(self.bg_array)
            self.save_clean_bg()
        
        self.raw_array = np.genfromtxt(file_name, comments='#', delimiter=',')
         
        
        # This gets created when remove_cosmics is run
        self.clean_array = None
        
        # These will keep track of what operations have been performed
        self.shot_normalized = False
        self.addenda = [0, self.fname]
        self.subtrahenda = []

    def __str__(self):
        '''
        This will print the description of the file.  I'm not sure what else
        would be useful at this time.
        '''
        return self.description
        
    def __len__(self):
        return len(self.hsg_data)
        
#    def __getitem__(self, *args):
#        '''
#        if single integer is passed, want to get the whole row 
#        '''
#        print 'getitem ', args       
#        if type(args[0]) is int:
#            return self.hsg_data[args[0],:]
#        #This happens if we're doing 2D slicing
#        if len(args[0])==2:
#            return self.hsg_data[args[0]]
        
    def __getslice__(self, *args):
        print 'getslice ', args
        #Simply pass the slice along to the data
        return self.hsg_data[args[0]:args[1]]
        
    def __iter__(self):
        for i in range(len(self.hsg_data)):
            yield self.hsg_data[i]
        return
        
    def __add__(self, other):
        '''
        Add together the image data from self.hsg_data, or add a constant to 
        that np.array.  
        '''
        ret = copy.deepcopy(self)
        #Add a constant offset to the data
        if type(other) in (int, float):
            ret.hsg_data[:,1] = self.hsg_data[:,1] + other
            ret.addenda[0] = ret.addenda[0] + other
        #or add the two counts together
        else:
            if np.isclose(ret.hsg_data[0,0], other.hsg_data[0,0]):
                ret.hsg_data[:,1] = self.hsg_data[:,1] + other.hsg_data[:,1]
                ret.addenda[0] = ret.addenda[0] + other.addenda[0]
                ret.addenda.extend(other.addenda[1:])
                ret.subtrahenda.extend(other.subtrahenda)
            else:
                raise Exception('These are not from the same grating settings')
        return ret
    
    def __sub__(self, other):
        '''
        This subtracts constants or other data sets between self.hsg_data.  I 
        think it even keeps track of what data sets are in the file and how 
        they got there
        '''
        ret = copy.deepcopy(self)
        if type(other) in (int, float):
            ret.hsg_data[:,1] = self.hsg_data[:,1] - other
            ret.addenda[0] = ret.addenda[0] - other
        else:
            if np.isclose(ret.hsg_data[0,0], other.hsg_data[0,0]):
                ret.hsg_data[:,1] = self.hsg_data[:,1] - other.hsg_data[:,1]
                ret.subtrahenda.extend(other.addenda[1:])
                ret.addenda.extend(other.subtrahenda)
            else:
                raise Exception('These are not from the same grating settings')
        return ret
    
    def get_series(self):
        '''
        This just returns the self.series attribute
        '''
        return self.series
    
    def set_ylimits(self, y_min, y_max):
        '''
        This changes y_min and y_max
        '''
        self.y_min = y_min
        self.y_max = y_max
    
    def clean_up(self, my_array, myreadnoise=3.0, mysigclip=5.0, mysigfrac=0.5, myobjlim=5.0, myverbose=True):
        '''
        This is a single operation of cosmic ray removal.

        If EMCCD isn't cold enough, the hot pixels wil be removed as well.  I 
        don't know if this is a bad thing?  I think it should just be a thing
        that doesn't happen.
        
        my_array = the array you want to clean
        myreadnoise = Level (imprecise!) of read noise used in the statistical
                      model built to recognize cosmic arrays.  It's a function
                      of gain on our CCD, I believe.  Dark current should be 
                      separate, technically, but it would look like noise 
                      because of shot noise and crap
        mysigclip = I forget
        mysicfrac = I forget
        myobjlim = I forget
        myverbose = Tells hsg_cosmics to print out what it's doing
        
        returns a clean my_array
        '''
        image_removal = cosmics.cosmicsimage(my_array, self.gain, 
                                             readnoise=myreadnoise, 
                                             sigclip=mysigclip, 
                                             sigfrac=mysigfrac, 
                                             objlim=myobjlim)
        image_removal.run(maxiter=4)
        return image_removal.cleanarray
    
    def clean_and_subtract(self):
        '''
        This subtracts the appropriate background.
        '''
        if self.bg_array_clean is not None:
            print "Background already clean"
            self.clean_array = self.clean_up(self.raw_array) - self.bg_array_clean
        else:
            print "I'm confused..."
            raise Exception("Where the shit is the background image?")
    
    def integrate(self):
        '''
        Integrates over vertical axis of the cleaned image
        '''
        self.hsg_signal = self.clean_array[:,self.y_min:self.y_max].sum(axis=1)
        
    
    def norm_shots(self):
        '''
        This method will normalize to the number of FEL pulses during the exposure time
        '''
        if self.shot_normalized == False: #Gotta keep from doing this more than once
            self.hsg_signal = self.hsg_signal / self.num_FEL
            self.shot_normalized = True
    
    def append_wavelengths(self):
        '''
        This appends the appropriate wavelengths for self.hsg_signal to make 
        np.array([wavelengths, self.hsg_signal]).
        
        Wait until after all the manipulation of the data before doing this?
        '''
        wavelengths = gen_wavelengths(self.center, self.grating)
        self.hsg_data = np.concatenate((wavelengths, self.hsg_signal)).reshape(2,1600).T
        
        #http://docs.scipy.org/doc/numpy/reference/arrays.interface.html        
        self.__array_interface__ = dict(shape=self.hsg_data.shape,
                                        typestr='f',
                                        version=3
                                        )
        
    def stitch_spectra(self):
        '''
        I really hope this is important later!
        '''        
        raise NotImplementedError
        
    def save_spectrum(self, name_str, sweep=None, iternum=''):
        '''
        This will do some saving business, I'm not exactly sure how it will 
        work yet.  I want it to be able to iterate and not overwrite files.
        sweep = what kind of parameter sweep this data set is for.
        iternum = a way to keep track of multiple images
        '''
        temp = list(self.addenda[1:])
        temp.append(str(self.addenda[0]))
        additions = ' + '.join(temp)
        subtractions = ' - '.join(self.subtrahenda)
        additions = additions + ' - ' + subtractions
        
        my_header = '''# Source: {1}, {2}\n# Files in data: {2}\nWavelength,Signal
                       \nnm,arb. u.'''.format(self.fname, self.description, additions)

        # The header should probably(?) contain more stuff?  And it needs to be
        # fixed for Origin's import style
        
        if sweep is None:
            file_name = '{0}{1}.txt'.format(name_str, iternum)
        elif sweep == 'NIR_power':
            file_name = '{0} {1}mW{2}.txt'.format(name_str, self.NIR_power, iternum)
            my_header.append('\n,' + self.NIR_power)
        elif sweep == 'NIR_lambda':
            file_name = '{0} {1}nm{2}.txt'.format(name_str, self.NIR_lambda, iternum)
            my_header.append('\n,' + self.NIR_lambda)
        elif sweep == 'THz_power':
            file_name = '{0} {1}kW{2}.txt'.format(name_str, self.THz_power, iternum)
            my_header.append('\n,' + self.THz_power)
        elif sweep == 'THz_freq':
            file_name = '{0} {1}GHz{2}.txt'.format(name_str, self.THz_freq, iternum)
            my_header.append('\n,' + self.THz_freq)
        else:
            print 'You fucked up!'
            return None
        
        np.savetxt(file_name, self.hsg_data, delimiter=',', header=my_header, comments='')
    
    def save_clean_bg(self):
        '''
        Saves a clean background to save time cleaning the same background over
        and over again.
        '''
        if self.bg_fresh:
            print "Saving me some juicy background"
            filename = self.background_name[:-4] + "_clean.txt"
            parameter_str = ','.join(self.parameters)
            save_info = list(self.info)
            info_str = ','.join(save_info)
            my_header = self.description + '\n' + parameter_str + '\n' + info_str
            np.savetxt(filename, self.bg_array_clean, delimiter=',', header=my_header)
            
    
    def clean_and_save_keepers(self, folder_str, name_str, sweep=None, iternum=''):
        '''
        I want this to save the raw, clean and bg array (the bg array may be 
        saved several times) as well as the processed data.  I think it would
        be best if they were saved in a separate folder.  
        This will only save raw_array, bg_array (if it exists), and 
        clean_array.  
        
        folder_str = name of the subdirectory for files to be saved to.
        name_str = base file name for these images.
        sweep = an explicit call to name files for their purpose in a sweep.
        iternum = index to add to files to keep them from overwriting.  
                  Requires user to increment things externally
        '''
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        if sweep is None:
            image_name = '{0}{1}.txt'.format(name_str, iternum)
        elif sweep == 'NIR_power':
            image_name = '{0} {1}mW{2}.txt'.format(name_str, self.NIR_power, iternum)
        elif sweep == 'NIR_lambda':
            image_name = '{0} {1}nm{2}.txt'.format(name_str, self.NIR_lambda, iternum)
        elif sweep == 'THz_power':
            image_name = '{0} {1}kW{2}.txt'.format(name_str, self.THz_power, iternum)
        elif sweep == 'THz_freq':
            image_name = '{0} {1}GHz{2}.txt'.format(name_str, self.THz_freq, iternum)
        else:
            print 'You fucked up!'
            return None
        
        bg_name = 'bg ' + image_name
        clean_name = 'clean ' + image_name
        
        parameter_str = ','.join(self.parameters)
        
        save_info = list(self.info)
        
        
        save_info[5] = 'BG_image=' + bg_name
        info_str = ','.join(save_info)
        
        my_header = self.description + '\n' + parameter_str + '\n' + info_str
        
        
        np.savetxt(os.path.join(folder_str, image_name), self.raw_array, 
                   delimiter=',', header=my_header)
        
        np.savetxt(os.path.join(folder_str, clean_name), self.clean_array, 
                   delimiter=',', header=my_header)
        
        if self.bg_fresh:
            np.savetxt(os.path.join(folder_str, bg_name), self.bg_array,
                       delimiter=',', header=my_header)
        
        bg_warning = open(os.path.join(folder_str, "Background warning.txt"), "w")
        bg_warning.write("The background file name strings do not increment correctly.  I don't know how to fix this easily at the moment.  ")
        bg_warning.write("The i^th raw_array's background file, j, should be the highest-incremented bg file in the folder such that j < i")
        bg_warning.close()

        
#def open_appropriate_directory(something):
#    '''
#    There needs to exist a good way to open files from a different directory
#    I don't know how to do this right now...
#    '''
#    raise NotImplementedError
    
def gen_wavelengths(center_lambda, grating):
    '''
    This returns a 1600 element list of wavelengths for each pixel in the EMCCD based on grating and center wavelength
    
    grating = which grating, 1 or 2
    center = center wavelength in nanometers
    '''
    b = 0.75
    k = -1.0
    r = 16.0e-6

    if grating == 1:
        d = 1./1800000.
        gamma = 0.2243885861015487
        delta = 1.353959952416065
    elif grating == 2:
        d = 1./1200000.
        gamma = 0.2207676478674860
        delta = 1.352055027245235
    else:
        print "What a dick, that's not a valid grating"
        return None
    
    center = center_lambda*10**-9
    wavelength_list = np.arange(-799.0,801.0)
    
    output = d*k**(-1)*((-1)*np.cos(delta+gamma+(-1)*np.arccos((-1/4)*(1/np.cos((1/2)*gamma))**2*(2*(np.cos((1/2)*gamma)**4*(2+(-1)*d**(-2)*k**2*center**2+2*np.cos(gamma)))**(1/2)+d**(-1)*k*center*np.sin(gamma)))+np.arctan(b**(-1)*(r*wavelength_list+b*np.cos(delta+gamma))*(1/np.sin(delta+gamma))))+(1+(-1/16)*(1/np.cos((1/2)*gamma))**4*(2*(np.cos((1/2)*gamma)**4*(2+(-1)*d**(-2)*k**2*center**2+2*np.cos(gamma)))**(1/2)+d**(-1)*k*center*np.sin(gamma))**2)**(1/2))   
    
    output = (output + center)*10**9
    return output


if __name__ == "__main__":
    test = EMCCDimage('740test.txt')

    test.remove_cosmics()
    test.integrate()
    plt.close()
    plt.plot(test.hsg_signal)
    plt.show()
