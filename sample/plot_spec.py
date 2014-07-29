# Module: plot_spec
""" Functions for plotting spectral data

Created by Stephanie Douglas, June 2011
"""
import numpy
import matplotlib.pyplot as plt
import read_spec as rs
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage.filters import uniform_filter

def fourplot (iden,hi_date=False):
    """
    Plot low, medium, and high (2 orders: 61 and 65) resolution
    spectra for a single object on one page

    Parameters
    ----------
    iden: string
         object identifier

    hi_date: string, optional
         date of observation for high-res data
    """
    
    w61,f61,snr61 = rs.higet(rs.hi_index(iden,'61',hi_date))
    w65,f65,snr65 = rs.higet(rs.hi_index(iden,'65',hi_date))
    wm,fm,snrm = rs.medget(rs.med_index(iden,'n3'))
    wl,fl,nl = rs.spexget(rs.spexd()[iden])
    name = rs.specname(iden)+" : "+rs.spectype(iden)+" : "+hi_date

    plt.figure(figsize=(7.5,10))

    plt.subplot(411)
    if (len(wl)<570):
        plt.plot(wl[50:540],fl[50:540],'k-')
        plt.xlim(xmin=numpy.min(wl[50:540]-0.01))
    else:
        plt.plot(wl,fl,'k-')
        plt.xlim(numpy.min(wl)-0.01,numpy.max(wl)+0.01)
    if nl[0] and (nl[0]<0.9):
        plt.plot(wl[50:540],nl[50:540],color='gray')
    plt.ylabel("Flux - SpeX Prism",horizontalalignment='center')
    plt.title(name)

    # Cropping off large noise at beginning and end of med-res spectrum
    i=100
    mini=0
    maxi=len(fm)-1
    while i>0:
        if fm[i]>1.5:
            mini=i+1
            break
        if fm[i]<-0.75:
            mini=i+1
            break
        i=i-1
    i=800
    while i<len(fm):
        if fm[i]>1.5:
            maxi=i-1
            break
        if fm[i]<-0.75:
            maxi=i-1
            break
        i=i+1 
    plt.subplot(412)
    plt.plot(wm[mini:maxi],fm[mini:maxi],'k-')
    plt.xlim(xmin=wm[mini]-0.001,xmax=wm[maxi]+0.001)
    plt.ylabel("Flux - NIRSPEC N3",horizontalalignment='center')

#cropping off end of order 65
    i=1000
    maxi = len(f65)-1
    last = f65[maxi]
    while i<len(f65):
        if ((f65[i]>last-0.0001) and (f65[i]<last+0.0001)):
            maxi=i
            break
        i=i+1 
    plt.subplot(413)
    plt.plot(w65[0:maxi],f65[0:maxi],'k-')
    plt.xlim(xmin=w65[0]-0.0001,xmax=w65[maxi]+0.0001)
    plt.ylabel("Flux - NIRSPEC 65",horizontalalignment='center')

#cropping off end of order 61
    i=1000
    maxi = len(f61)-1
    last = f61[maxi]
    while i<len(f61):
        if ((f61[i]>last-0.0001) and (f61[i]<last+0.0001)):
            maxi=i
            break
        i=i+1 
    plt.subplot(414)
    plt.plot(w61[0:maxi],f61[0:maxi],'k-')
    plt.xlim(xmin=w61[0]-0.0001,xmax=w61[maxi]+0.0001)
    plt.ylabel("Flux - NIRSPEC 61",horizontalalignment='center')
    plt.xlabel('wavelength (microns)')

def plotall(output_file):
    """
    Plots all spectra for Ldwarf Project to a single file.

    Parameters
    ----------
    output_file: string
         output filename, including .pdf extension
    """

    pp = PdfPages(output_file)

    ty = rs.types()
    types = ty.keys()
    types.sort()

    i = 0
    while i<len(types):
        j = 0
        while j<len(ty[types[i]]):
            iden = ty[types[i]][j]
            print types[i],iden
            date=rs.hidate(iden)[0]
            fourplot(iden,date)
            pp.savefig()
            plt.close()
            j = j+1
        i=i+1
    pp.close()

def plotalldates(identifier,numplots=2):
    """
    Plots hi-res spectra (orders 65 and 61) for a single
    object over all nights that object was observed
    """
    
    output_file = identifier+"_alldates.pdf"
    pp = PdfPages(output_file)

    dates = rs.hidate(identifier)
    if numplots==2:
        i = 0
        while i<len(dates):
            print dates[i]
            w61,f61,snr61 = rs.higet(rs.hi_index(identifier,'61',dates[i]))
            w65,f65,snr65 = rs.higet(rs.hi_index(identifier,'65',dates[i]))
            name = rs.specname(identifier)+" : "+rs.spectype(identifier)+" : "+dates[i]
            plt.figure(figsize=(7.5,8))

#cropping off end of order 65
            j=1000
            maxi = len(f65)-1
            last = f65[maxi]
            while j<len(f65):
                if ((f65[j]>last-0.0001) and (f65[j]<last+0.0001)):
                    maxi=j
                    break
                j=j+1 
            plt.subplot(211)
            plt.plot(w65[0:maxi],f65[0:maxi],'k-')
            plt.ylabel("Flux - NIRSPEC 65",horizontalalignment='center')
            plt.xlim(xmin=w65[0]-0.0001,xmax=w65[maxi]+0.0001)
            plt.title(dates[i])

#cropping off end of order 65
            j=1000
            maxi = len(f65)-1
            last = f65[maxi]
            while j<len(f65):
                if ((f65[j]>last-0.0001) and (f65[j]<last+0.0001)):
                    maxi=j
                    break
                j=j+1 
            plt.subplot(212)
            plt.plot(w61[0:maxi],f61[0:maxi],'k-')
            plt.ylabel("Flux - NIRSPEC 61",horizontalalignment='center')
            plt.xlim(xmin=w61[0]-0.0001,xmax=w61[maxi]+0.0001)
            plt.xlabel('wavelength (microns)')
            pp.savefig()
            plt.close()
            i = i+1
        pp.close()

    else:
        i = 0
        while i<len(dates):
            print dates[i]
            fourplot(identifier,dates[i])
            pp.savefig()
            plt.close()
            i = i+1
        pp.close()

def lowmodel (iden,temp,grav,fsed=None,model='phoenix-cond',resolution='low'):
    """
    """

    info = 'Model Parameters:\n Teff = '+str(temp)+', Logg = '+str(grav)
    wl,fl,nl = rs.spexget(rs.spexd()[iden])
    wmod,fmod1 = rs.modget(temp,grav,model,resolution)

    fmod = uniform_filter(fmod1,8)
    mmax = numpy.average(fmod[1498:1502])
    lmax = numpy.average(fl[192:196])

    plt.plot(wmod,fmod/mmax,'r-')
    plt.plot(wl,fl/lmax,'k-')
    plt.ylabel("Flux")
    plt.title(rs.specname(iden)+" : "+rs.spectype(iden))
    plt.text(2.,0.8,info,color='r')

def plot_teffs (iden,logg,order='low',tmin=1550,tmax=1950,tstep=100,output_file=None):
    """
    """

    if output_file==None:
        pp = PdfPages(iden+"_logg"+str(logg)+".pdf")
    else:
        pp = PdfPages(output_file)

    check = False
    w=[]
    f=[]
    n=[]
    if order=='low':
        w,f,n = rs.spexget(rs.spexd()[iden])
    temp = tmin
    plt.figure(figsize=(7.5,10))
    snum=1

    while temp<=tmax:
        plt.subplot(3,1,snum)

        if order=='low':
            if rs.mod_index(temp,logg) >= 0:
                check=True
                lowmodel(iden,temp,logg)

        if snum==3:
            plt.xlabel("Wavelength (microns)")
            pp.savefig()
            plt.clf()
            snum = 1
        else:
            snum=snum+1

        temp=temp+tstep
        if temp>tmax:
            if snum!=1:
                pp.savefig()
            plt.close()
            

    pp.close()


def plot_loggs (iden,teff,order='low',gmin=4.,gmax=5.,gstep=0.1,output_file=None):
    """
    """

    if output_file==None:
        pp = PdfPages(iden+"_teff"+str(teff)+".pdf")
    else:
        pp = PdfPages(output_file)

    w=[]
    f=[]
    n=[]
    if order=='low':
        w,f,n = rs.spexget(rs.spexd()[iden])
    logg = gmin
    plt.figure(figsize=(7.5,10))
    snum=1

    while logg<=gmax:
        plt.subplot(3,1,snum)

        if order=='low':
            lowmodel(iden,teff,logg)

        if snum==3:
            plt.xlabel("Wavelength (microns)")
            pp.savefig()
            plt.clf()
            snum = 1
        else:
            snum=snum+1

        logg=logg+gstep
        if logg>gmax:
            if snum!=1:
                pp.savefig()
            plt.close()
            

    pp.close()


def mod_orders(iden,temp,grav,fsed=None,date=None,model='marley',output_file=None):
    """
    """

    if output_file==None:
        if fsed==None:
            output_file = iden+"_t"+str(temp)+"g"+str(grav*10)+model+'.pdf'
        else:
            output_file = iden+"_t"+str(temp)+"g"+str(grav*10)+str(fsed)+'.pdf'

    if date==None:
        date=rs.hidate(iden)[0]

    def ordplot(order):
        w,f,n = rs.higet(rs.hi_index(iden,order,date))
        wm,fm = rs.modget(temp,grav,fsed,model,'high',order,normalize=False)
#        ck = 
        plt.plot(wm,fm,'r-')
        plt.plot(w,f,'k-')
        plt.ylim(ymin=0.0,ymax=1.2)
        plt.xlim(xmin=min(w)-0.0001,xmax=max(w)+0.0001)
        plt.text(min(w)+0.001,0.2,order,color='k')

    plt.figure(figsize=(10,7.5))

    plt.subplot(421)
    ordplot('65')
    plt.xticks((1.165,1.170,1.175,1.180))
    plt.subplot(422)
    ordplot('64')
    plt.subplot(423)
    ordplot('63')
    plt.subplot(424)
    ordplot('62')
    plt.subplot(425)
    ordplot('61')
    plt.subplot(426)
    ordplot('59')
    plt.subplot(427)
    ordplot('58')
    plt.subplot(428)    
    info = 'Model Parameters (Cushing 08):\n Teff = '+str(temp)+'\n Logg = '+str(grav)+'\n fsed = '+str(fsed)
    plt.text(0.2,0.2,info,color='k')

    plt.suptitle(rs.specname(iden)+" : "+rs.spectype(iden)+" : "+date)


