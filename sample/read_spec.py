# Module: read_spec
""" Functions for retrieving spectral data

Created by Stephanie Douglas, June 2011
"""


import numpy
import read_spec
from scipy.io.idl import readsav
import asciitable
import pyfits
from scipy.ndimage.filters import uniform_filter
import os
"""
hipath = '/Users/Steph/Documents/school/lowmass/summerAMNH/LdwarfSpectra/'
medpath = '/Users/Steph/Documents/school/lowmass/summerAMNH/LdwarfSpectra/'
lowpath = '/Users/Steph/Documents/school/lowmass/summerAMNH/LdwarfSpectra/'
modelpath= '/Users/Steph/Documents/school/lowmass/summerAMNH/modelSpectra/'
"""
hipath =   '/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/'
medpath =  '/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/'
lowpath =  '/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/'
modelpath= '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'

### hi_index ###
# Opens an IDL save file with high-resolution data, then finds the line
# in the file that corresponds to the input object name and spectral order,
# and optionally a date of observation.  
def hi_index (name,order,date):
    """
    Retrieve the index of an object from obspechigh.save

    Parameters
    ----------
    name: string
         Object identifier
    order: string
         Spectral order ("echelle order").  Although the order is
         given by a number, it must be contained in quotes
    date: string, optional
         Date of observation.  Only necessary if the save file contains more
         than one observation of the object at a particular order.

    Returns
    -------
    index: int
    """
    hifile = hipath + 'obspechigh.save'
    high = readsav(hifile)
    name_array = high.obspechigh.NAME
    order_array = high.obspechigh.ORDER
    date_array = high.obspechigh.DATE
    okn = numpy.equal(name_array,name)
    ln = len(okn)
    oko = numpy.equal(order_array,order)
    okd = numpy.equal(date_array,date)
    index = -1
    i = 0
    while (i<ln):
        j = 0
        while (j<ln):
            if (okn[i] and oko[j] and (i==j)):
                if date:
                    k = 0
                    while (k<ln):
                        if (oko[j] and okd[k] and (j==k)):
                            index = i
                        k = k+1
                else:
                    index = i
            j = j+1
        i = i+1
    return index

### med_index ###
# Opens an IDL save file with high-resolution data, then finds the line
# in the file that corresponds to the input object name and spectral filter  
def med_index (name,filter='n3'):
    """
    Retrieve the index of an object from obspecmed.save

    Parameters
    ----------
    name: string
         Object identifier
    filter: string, optional
         Filter object was observed in (default n3)

    Returns
    -------
    index: int
    """
    medfile = medpath + 'obspecmed.save'
    med = readsav(medfile)
    name_array = med.obspecmed.NAME
    filter_array = med.obspecmed.FILTER
    okn = numpy.equal(name_array,name)
    ln = len(okn)
    oko = numpy.equal(filter_array,filter)
    index = -1
    i = 0
    while (i<ln):
        j = 0
        while (j<ln):
            if (okn[i] and oko[j] and (i==j)):
                index = i
            j = j+1
        i = i+1
    return index

### spexd ###
# Generates a dictionary of filenames (with full paths) for all SpeX Prism
# data.  Returns this dictionary, and filenames can be accessed using
# Emily Rice's naming conventions from her .save files (those accessed above)
def spexd ():
    """
    Returns a dictionary whose keys are object identifiers and whose
    values are the filenames (including the full path, if specified in
    read_spec.py) of all SpeX Prism data.

    Parameters
    ----------
    none

    Returns
    -------
    D : dictionary
    """
    D = {'2m0345':lowpath+'2M0345.fits',
         'lp944':lowpath+'U10201_0339-3525.fits',
         '2m0746':lowpath+'spex_prism_0746+2000_U10668.fits',
         '2m0208':lowpath+'2MASSWJ0208+25_spex.dat',
         '2m1300':lowpath+'U11115_1300+1912.fits',
         'si0921':lowpath+'U20336_0921-2104.fits',
         '2m2057':lowpath+'U12054_2057-0252.fits',
         '2m0015':lowpath+'spex_prism_0015+3516_080908.fits',
         'g1963b':lowpath+'G196-3B_R400.fits',
         'kelu-1':lowpath+'spex_prism_kelu-1_060411.fits',
         '2m1615':lowpath+'spex_prism_1615+3559_090630.txt',
         '2m1506':lowpath+'u11291_1506+1321_050323.fits',
         '2m0036':lowpath+'2M0036.fits',
         'gd165b':lowpath+'spex_prism_gd165b_090629.txt',
         '2m2224':lowpath+'U12128_2224-0158.fits',
         '2m1507':lowpath+'2M1507.fits',
         'DE1228':lowpath+'spex_prism_1228-1547_080114.fits',
         'DE0205':lowpath+'DENIS0205.fits'}
    return D

### spexget ###
# Given a filename (i.e. from the dictionary output by read_spec.spexd())
# returns the data from the file in 3 arrays - wavelength, flux(wavelength), and 
# noise(wavelength)
def spexget (file,normalize=True):
    """
    Retrieves wavelenth, flux, and noise data from SpeX Prism data files. 
    Will get data from either .fits or text files.

    Parameters
    ----------
    file: string
         Filename of the data file to be read

    Returns
    -------
    wav: array
         Wavelength bins froms spectral observation

    flu: array
         flux per wavelength bin

    noi: array
         noise per wavelength bin
    """
    wav = []
    flu = []
    noi = []
    l = len(file)
    ext = file[l-5:l]
    if ext=='.fits':
        print 'Opening a FITS file'
        hdus = pyfits.open(file)
        wav = hdus[0].data[0]
        flu = hdus[0].data[1]
        if len(hdus[0].data)==3:
            noi = hdus[0].data[2]
        else:
            noi.append(0)
    else:
        print 'Opening a text file'
        arr = asciitable.read(file)
        wav = arr['col1']
        flu = arr['col2']
        noi = arr['col3']
    if normalize:
        if len(wav)==564:
            mx = numpy.average(flu[193:197])
            flu = flu/mx
            noi = noi/mx
        else:
            j1 = numpy.where(wav<1.3001)
            j2 = numpy.where(wav>1.2999)
            j = numpy.intersect1d(numpy.asarray(j1),numpy.asarray(j2))
            if j:
                print j
            else:
                j = numpy.arange(1)
                j[0] = j2[0][0]
                print j
            if len(j)==1:
                mx = numpy.average(flu[j-2:j+2])
            else:
                mx = numpy.average(flu[j])
            flu = flu/mx
            noi = noi/mx
    return wav,flu,noi

### higet ###
# Given an index number (i.e. the output of hi_index) returns the data for a
# particular object, spectral order, and observation date in 3 arrays - 
# wavelength, flux(wavelength), and SNR(wavelength)
def higet (index) :
    """
    Retrieves wavelenth, flux, and SNR data from high-
    resolution data found in obspechigh.save

    Parameters
    ----------
    index: integer
         index of desired object and spectral order (i.e. 
         output of hi_index)

    Returns
    -------
    wav: array
         Wavelength bins froms spectral observation

    flu: array
         flux per wavelength bin

    snr: array
         signal-to-noise
    """
    hifile = hipath + 'obspechigh.save'
    high = readsav(hifile)
    arr = high.obspechigh[index]
    wav = arr.W
    flu = arr.F
    snr = arr.SNR
    return wav,flu,snr

### medget ###
# Given an index number (i.e. the output of med_index) returns the data for a
# particular object and filter in 3 arrays -  wavelength, flux(wavelength),
# and SNR(wavelength)
def medget (index,normalize=True) :
    """
    Retrieves wavelenth, flux, and SNR data from medium-
    resolution data found in obspechigh.save

    Parameters
    ----------
    file: index
         index of desired object and spectral order (i.e. 
         output of hi_index)

    Returns
    -------
    wav: array
         Wavelength bins froms spectral observation

    flu: array
         flux per wavelength bin

    snr: array
         signal-to-noise
    """
    medfile = medpath + 'obspecmed.save'
    med = readsav(medfile)
    arr = med.obspecmed[index]
    wav = arr.W
    flu = arr.F
    snr = arr.SNR
    if normalize:
        mx=numpy.average(flu[424:428])
        flu = flu/mx
        snr = snr/mx
    noi = flu*0.05
    return wav,flu,noi

def specname(identifier):
    """
    Returns the name of the object that corresponds 
    to the input identifier
    """
    
    D = {'2m0345':'2MASSW J0345+25',
         'lp944' :'LP 944-20',
         '2m0746':'2MASSW J0746+20AB',
         '2m0208':'2MASSW J0208+25',
         '2m1300':'2MASS 1300+1912',
         'si0921':'SIPS 0921-21',
         '2m2057':'2MASSW J2057-02',
         '2m0015':'2MASSW J0015+35',
         'g1963b':'G196-3B',
         'kelu-1':'Kelu-1',
         '2m1615':'2MASSW J1615+35',
         '2m1506':'2MASSW J1506+13',
         '2m0036':'2MASSW J0036+18',
         'gd165b':'GD 165B',
         '2m2224':'2MASSW J2224-01',
         '2m1507':'2MASSW J1507-16',
         'DE1228':'DENIS-P J1228-15AB',
         'DE0205':'DENIS-P J0205-11AB'}

    name = D[identifier]
    return name

def hidate(identifier):
    """
    Returns a list of all hi-res observation dates for the object 
    that corresponds to the input identifier.  The first date
    in the list is the observation which we are using in our
    analysis.
    """
    
    D = {'2m0345':['06dec00','11jan06','04dec00'],
         'lp944' :False,
         '2m0746':['10jan06','01jan02'],
         '2m0208':['06dec08','25jul00','29jul00'],
         '2m1300':False,
         'si0921':False,
         '2m2057':['29may07','20jul03','25jul00','28jul00'],
         '2m0015':['05dec08','25jul00','29jul00'],
         'g1963b':['21mar08','11jan06','20may06','23apr02'],
         'kelu-1':['20mar08','10jan06','12may03','19may06','31may07'],
         '2m1615':False,
         '2m1506':['21mar08','25jul00','28jul00','30may07'],
         '2m0036':['11dec05','04dec00','05dec00','06dec00'],
         'gd165b':['13may03','20may06'],
         '2m2224':['30may07','25jul00'],
         '2m1507':['29may07','25apr00','25jul00'],
         'DE1228':['30may07'],
         'DE0205':['09oct01','10jan06']}

    date = D[identifier]
    return date


def types():
    """
    Returns a dictionary whose keys are spectral types from
    L0-L7, with some half-subtypes.  The associated values
    are lists of one or more elements, listing all objects
    with that spectral type.
    """
    t = {#'M9':['lp944']',
        'L0':['2m0345'],
        'L0.5':['2m0746'],
        'L1':['2m0208'],
        #'L1/L3',['2m1300'],
        'L1.5':['2m2057'],#'si0921'],
        'L2':['2m0015'],#,'kelu-1'],#'g1963b'],
        'L3':['2m1506'],#'2m1615'],
        'L3.5':['2m0036'],
        'L4':['gd165b'],
        'L4.5':['2m2224'],
        'L5':['2m1507'],
        'L6':['DE1228'],
        'L7':['DE0205']}

    return t

def spectype(identifier):
    """
    Returns the spectral type of the object that corresponds 
    to the input identifier
    """
    
    D = {'2m0345':'L0',
         'lp944' :'M9',
         '2m0746':'L0.5',
         '2m0208':'L1',
         '2m1300':'L1/L3',
         'si0921':'L1.5',
         '2m2057':'L1.5',
         '2m0015':'L2',
         'g1963b':'L2',
         'kelu-1':'L2',
         '2m1615':'L3',
         '2m1506':'L3',
         '2m0036':'L3.5',
         'gd165b':'L4',
         '2m2224':'L4.5',
         '2m1507':'L5',
         'DE1228':'L6',
         'DE0205':'L7'}

    name = D[identifier]
    return name



def modget (temp,grav,fsed=None,model='phoenix-cond',resolution='low',order=None, normalize=True, smooth=True) :
    """
    Returns a model spectrum for given parameters, a given model, 
    and a given resolution.

    All model files should be in the directory specified by the
    modelpath parameter specified at the top of the read_spec.py
    file.

    Parameters
    ----------
    temp: integer
         Effective temperature for desired model spectrum.  Models 
         do not compute to fractions of a Kelvin.

    grav: float
         log(g) for desired model spectrum. 

    model: string, optional
         Name of the general model from which to get the specific 
         synthetic spectrum.  
         Allowed model names: 'phoenix-cond'

    resolution: string, optional
         Resolution level for desired model spectrum.
         Not necessary if there is only one resolution for that model.
         Allowed resolutions: 'low'

    normalize: bool, optional
         If set to true (default), the returned flux will be normalized 
         at an appropriate wavelength for the spectral range.  

    smooth: bool, optional
         If set to true (default), the flux will be smoothed before 
         being normalized.

    Returns
    -------
    wav: array
         Wavelengths for model spectrum, in microns.
         Returns an empty array if a model fulfilling the specified
         parameters is not found.
    
    flu: array
         Model flux at each wavelength.
         Returns an empty array if a model fulfilling the specified
         parameters is not found.
    
    """

    wav = []
    flu = []
    grav2 = ''
    

    if model.lower()=='phoenix-cond':
        if resolution.lower()=='low':
            if temp>1950:
                print 'T>1950 not allowed, using T=1950'
                temp=1950
            elif temp<950:
                print 'T<950 not allowed, using T=1950'
            elif temp%50!=0:
                print 'T='+temp+' not allowed'
                temp = int(temp/50)*50
                print 'using T='+temp
            if grav>6.:
                print 'logg>6.0 not allowed, using logg=6.0'
                grav=6.
            elif grav<3.:
                print 'logg<3.0 not allowed, using logg=3.0'
                grav=3.
            modfile = modelpath + 'modelspeccondlowres.save'
            mod = readsav(modfile)
            temp_array = mod.modelspec.TEFF
            grav_array = mod.modelspec.LOGG
            ln = len(temp_array)
            ok = numpy.equal(temp_array,temp)
            i = numpy.min(numpy.where(ok)[0])
            index = -1
            while (i<ln):
                if ok[i]:
                    if ((grav_array[i] < grav+0.05) and (grav_array[i] > grav-0.05)):
                        index = i
                        break
                i = i+1
            if index>0:
                flu = mod.modelspec.fsyn[index]
                wav = mod.wsyn/10000
            else:
                print 'model not found for given parameters'
                print 'temp = ', temp, 'grav = ', grav
            if normalize:
                if smooth:
                    flu = uniform_filter(flu,20)
                mx = numpy.average(flu[1498:1502])
                flu = flu/mx
        else:
            print 'model not found for resolution ', resolution
    elif model.lower()=='phoenix-dusty':
        if resolution.lower()=='low':
            if temp>4500 :
                print 'T>4500 not allowed, using T=4500'
                temp=2400
            if temp<1800:
                print 'T<1800 not allowed, using T=1800'
                temp=1800
            modfile = modelpath + 'modelspeclowresldwarfs.save'
            mod = readsav(modfile)
            temp_array = mod.modelspec.TEFF
            grav_array = mod.modelspec.LOGG
            ln = len(temp_array)
            ok = numpy.equal(temp_array,temp)
            i = numpy.min(numpy.where(ok))
            index = -1
            while (i<ln):
                if ok[i]:
                    if ((grav_array[i] < grav+0.05) and (grav_array[i] > grav-0.05)):
                        index = i
                        break
                i = i+1
            if index>0:
                flu = mod.modelspec.fsyn[index]
                wav = mod.wsyn/10000
            else:
                print 'model not found for given parameters'
                print 'temp = ', temp, 'grav = ', grav
            if normalize:
                if smooth:
                    flu = uniform_filter(flu,20)
                mx = numpy.average(flu[747:751])
                flu = flu/mx
        elif resolution.lower()=='med' or resolution.lower()=='high':
            modfile = modelpath + 'modelspeclowresldwarfs.save'
            mod = readsav(modfile)
            temp_array = mod.modelspecldwarfs.TEFF
            grav_array = mod.modelspecldwarfs.LOGG
            ln = len(temp_array)
            ok = numpy.equal(temp_array,temp)
            i = numpy.min(numpy.where(ok))
            index = -1
            while (i<ln):
                if ok[i]:
                    if ((grav_array[i] < grav+0.05) and (grav_array[i] > grav-0.05)):
                        index = i
                        break
                i = i+1
            if index>0:
                flu = mod.modelspecldwarfs.fsyn[index]
                wav = mod.wsyn/10000
            else:
                print 'model not found for given parameters'
                print 'temp = ', temp, 'grav = ', grav
            #if normalize:
            #    if smooth:
            #        flu = uniform_filter(flu,20)
            #    mx = numpy.average(flu[747:751])
            #    flu = flu/mx
        else:
            print 'model not found for resolution ', resolution
    elif model.lower()=='marley':
        if grav==4.0:
            grav2='100'
        elif grav==4.5:
            grav2='300'
        elif grav==5.0:
            grav2='1000'
        elif grav==5.5:
            grav2='3000'
        else: 
            print 'no model found with that gravity'
            print 'allowed gravities are 4.0, 4.5, 5.0, and 5.5'
        if resolution.lower()=='low':
            filename='sp_t'+str(int(temp))+'g'+grav2+'f'+str(int(fsed))
            modfile=modelpath+'marley_lowres/'+filename
            if os.path.exists(modfile):
                mod = asciitable.read(modfile,data_start=4)
                wav = mod['col1']
                flu1 = mod['col2']
                flu = flu1*(3e7)/wav**2
                if normalize:
                    if smooth:
                        flu = uniform_filter(flu,50)
                    mx = numpy.average(flu[1498:1502])
                    #mx = numpy.average(flu[16344:16348])
                    flu = flu/mx
            else:
                print 'model file not found '+filename
        elif resolution.lower()=='med':
            filename='sp_t'+str(int(temp))+'g'+grav2+'f'+str(int(fsed))
            modfile=modelpath+'marley_highres/'+filename
            if os.path.exists(modfile):
                mod = asciitable.read(modfile,data_start=4)
                wav = mod['col1']
                flu1 = mod['col2']
                flu = flu1*(3e7)/wav**2
                if normalize:
                    mx = numpy.average(flu[7406:7410])
                    flu = flu/mx
        elif resolution.lower()=='high':
            filename='sp_t'+str(int(temp))+'g'+grav2+'f'+str(int(fsed))
            modfile=modelpath+'marley_highres/'+filename
            w = []
            f= []
            if os.path.exists(modfile):
                mod = asciitable.read(modfile,data_start=4)
                w = mod['col1']
                f1 = mod['col2']
                f = f1*(3e7)/w**2
                if str(order)=='58':
                    wav = numpy.copy(w[5651:7101])
                    flu = numpy.copy(f[5651:7101])
                elif str(order)=='59':
                    wav = numpy.copy(w[7332:8802])
                    flu = numpy.copy(f[7332:8802])
                elif str(order)=='61':
                    wav = numpy.copy(w[10645:12110])
                    flu = numpy.copy(f[10645:12110])
                elif str(order)=='62':
                    wav = numpy.copy(w[12263:13727])
                    flu = numpy.copy(f[12263:13727])
                elif str(order)=='63':
                    wav = numpy.copy(w[13558:15313])
                    flu = numpy.copy(f[13558:15313])
                elif str(order)=='64':
                    wav = numpy.copy(w[15421:16923])
                    flu = numpy.copy(f[15421:16923])
                elif str(order)=='65':
                    wav = numpy.copy(w[16958:18416])
                    flu = numpy.copy(f[16958:18416])
                elif str(order)=='66':
                    wav = numpy.copy(w[16891:18347])
                    flu = numpy.copy(f[16891:18347])
                else:
                    print 'model not found at order '+str(order)
                if normalize:
                    m1 = max(flu)
                    i = 0
                    whr=-1
                    while (i<len(flu)):
                        if (flu[i]<(m1+0.000001)) and (flu[i]>(m1-0.000001)):
                            whr = i
                            break
                        i=i+1
                    mx = numpy.average(flu[whr-2:whr+2])
                    flu = flu/mx
            else:
                print 'model file not found '+filename
        else:
            print 'model not found for resolution ', resolution
    else:
        print 'model ',model, ' not found'
        print 'accepted models are "phoenix-cond","phoenix-dusty", and "marley"'

    return wav,flu


def alldates (iden,order='61'):
    """
    Prints all dates for which high-resolution observations
    were performed at a given spectral order on a given object.

    Parameters
    ----------
    iden: string
         Identifier for the object.

    order: string, optional
         Spectral order for desired observations (defaults to 61)
         Must be a string.
    """
    hifile = hipath + 'obspechigh.save'
    high = readsav(hifile)
    name_array = high.obspechigh.NAME
    order_array = high.obspechigh.ORDER
    date_array = high.obspechigh.DATE
    okn = numpy.equal(name_array,iden)
    ln = len(okn)
    oko = numpy.equal(order_array,order)
    i = 0
    while (i<ln):
        j = 0
        while (j<ln):
            if (okn[i] and oko[j] and (i==j)):
                print date_array[i]
            j = j+1
        i = i+1
