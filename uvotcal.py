#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-



'''
  This code is designed for calibrating Swift UVOT Grism images,
  and writing the calibration files.
  
  List of functions/methods:
  
  
  (xplines, yplines, dislines, lines) = wr86zemaxlines(dis_zmx,wav_zmx,xpix,ypix,wave)
    gets the positions for some prominent spectral lines in WR86 to overlay on the 
    contour plot.
  
  ( (dis, spnet, angle, anker, coef_zmx, pix_zmx, wav_zmx), 
            (anker_uvw1,anker_as, anker_field),
            (bg, bg1, bg2, extimg, spimg, spnetimg, offset) ,
            (C_zero,C_1,C_2,C_3,C_min1), 
            (xpix,ypix, zmxdis,zmxwav, wave,theta)
             (img, xplines, yplines, dislines, lines)) = getFirstOrder(ra,dec,filestub,ext)     
     Returns nearly all that is needed for making a spectral calibration.
     calls the functions listed below. 
            
  (x, y) = findAnker(filestub, ext, RA, DEC, method=None, attfile=None) 
    This routine assumes that an image in a lenticular filter 
    was taken right before or after the grism observation. 
    Enter the file stub (the part that is generic to the observation)
    and the position of the source in RA,DEC (decimal degrees). 
    The position of the anchor is returned in DET pixels with centre
    dt image referenced to [1023.5,1023.5]. 
    returns list of anker position.

  (dis, spnet, bg, bg1, bg2, c, spimg, spnetimg) = extractSpecImg(file,ext,anker,angle)
    This functions similar to the uvot tool uvotimsum, but returns a
    slice of the rotated image including the zero order also. 
    The output currently includes:
      dis = array with the dispersion referenced to the 260nm first order point.
      spnet = array with the first order spectrum net counts
      bg = array with background
      bg1, bg2 = arrays with background on each side of the spectrum
      c = image slice of spectrum and background 
      spimg = the image around the spectrum

   [sumimg, BG, EM, DQ] = uvotGrismCoadd(filelist)
      To coadd first order spectra to the first one in the list, this function 
      was made. Here the input is a list of the files to be used.
      The output is a summed image, background, exposure map, and quality.
       
   Once the observed line positions and zemax model are available, the offset (in pixels) of 260nm in the 
   observed data and the relative scale to the zemax model can be calculated. If only a few data points 
   are present, the order should be low. 
   
       pixZeroOffset, wav_scale, pix_scale = uvotgrism.getScale(wav_obs,pix_obs,wav_zmx,pix_zmx,order=4)
      
   wr86zemaxlines returns the predicted positions of spectral lines predicted by zemax 
   the default lines is specifically for WR86. 
      
      (xplines, yplines, dislines, lines) = wr86zemaxlines(zmxdis[1],zmxwav[1],xpix,ypix,wave)
     
2012-12-23 NPMK added support for zemax flux model       
'''
# Developed by N.P.M. Kuin (MSSL/UCL)

__version__ = '20171221-1.0' 

import math
import numpy as N
from astropy.io import fits as pyfits
import re
import warnings
try:
  import imagestats
except:
  from stsci import imagestats  
import scipy
from scipy import interpolate
import scipy.ndimage as ndimage
from scipy.ndimage import convolve
from scipy.signal import boxcar
from scipy.optimize import leastsq
from pylab import polyfit, polyval, legend
from uvotpy.uvotgetspec import extractSpecImg, boresight, interpol, hydrogen,\
  doublegaussian, singlegaussian, Dsinglegaussian, Ddoublegaussian, bilinear,\
  gaussPlusPoly, DgaussPlusPoly, findBackground 
from zemax import rdzemax, correctAnkPos, read_zemax_flux
from uvotpy import uvotmisc, uvotspec
from uvotpy.uvotmisc import interpgrid, uvotrotvec
from uvotpy import uvotplot
import datetime

def getFirstOrder(RA,DEC,filestub, ext, lfilter='uvw1', lfilt1_ext=None, measured_file=None,clocked=False,
      lfilt2=None,lfilt2_ext=None,test=None,testparam=None,chatter=1,lines='WC',splineorder=3,
      getzmxmode='spline',smooth=50,wpixscale=None,spextwidth=13,calfile=None,xlimits=[15,2100]):
   '''
   Make all the necessary calls to reduce the data ...
   INPUT: RA, DEC
          filestub, like  'sw0005700005'
          ext extension number, like 1
          measured_file (optional, like wr86_23002_1.obs_dis.txt)
          chatter  verbosity of program
   OUTPUT:      regular output:
      
      ( (dis, spnet, angle, anker, coef_zmx, pix_zmx, wav_zmx), 
            (anker_uvw1,anker_as, anker_field,ank_c),
            (bg, bg1, bg2, extimg, spimg, spnetimg, offset) ,
            (C_zero,C_1,C_2,C_3,C_min1), 
            (xpix,ypix, zmxdis,zmxwav, wave,theta), 
            (img, xplines, yplines, dislines, lines), hdr)  
   where        
       dis         # dispersion with zero at 260nm[uvgrism] or at 420nm [vgrism]
       spnet       # background-substracted spectrum from 'spnetimg'
       angle       # rotation-angle used to get 'extimg'
             anker,       # anker position used to rotate around in 'img'
             coef_zmx,    # polynomial coefficients to find zemax wavelength from dis
             pix_zmx,     # interpolated zemax model pixel positions at anker orders [0,1,2,3,-1]
             wav_zmx),    # wavelength for the interpolated pixel positions
            (anker_uvw1,  # anker position in uvw1 DET coordinates derived from RA,DEC, and uvw1 header
             anker_as,    # anker position offset from boresight in arcsec from 'anker_uvw1' and known boresight
             anker_field),# derived field positions from 'anker_as'
            (bg,          # mean background, smoothed, with sources removed
             bg1,         # one-sided background, sources removed, smoothed 
             bg2,         # one-seded background, opposite, sources removed, smoothed
             extimg,   # image extracted of source and background, 201 pixels wide, all orders.
             spimg,    # image centered on first order position
             spnetimg # background-subtracted 'spimg'
             offset   # offset of spectrum from expected position based on 'anker' at 260nm, first order
             C_zero    # polynomial coefficient to zemax model of zero order for the anker position (tbd)
             C_1      # same for first order;  C_2  # same for second order;  C_3  # same for third order
             C_min1   # same for the minus-one order
             xpix,ypix,     # interpolated x,y positions of the zemax model for the 'anker' position and 'wave'
                            # no scaling applied, no offsets applied  
             zmxdis,zmxwav, # for all orders dispersion and wavelengths of valid points 
             wave,          # goes with xpix,ypix
             theta         # angles for each order at 260nm
             img              # original image
             xplines, yplines,# position for selected 'lines'
             dislines,        # dispersion for those 'lines'
             lines          # the selected 'lines' (default from wr86zemaxlines() ) 
             xlimits    # limits for valid zemax points (especially for determining the slope of the spectrum)
          
   '''
   if test != 'caltable':
     # test if u or v grism file and set variable 
     import os
     
     if filestub[:2] != 'sw': filestub = 'sw'+filestub
     
     if os.access(filestub+'ugu_dt.img',os.F_OK):
        ugufile =  filestub+'ugu_dt.img'
     elif os.access(filestub+'ugv_dt.img',os.F_OK):     
        ugufile =  filestub+'ugv_dt.img'
     if chatter > 1 : print("ugufile:",ugufile)
     attfile =  filestub+'pat.fits'
     zmxfile =  filestub+'_'+str(ext)+'_zmx.fit'
     hdr = pyfits.getheader(ugufile,ext)
     wheelpos = hdr['WHEELPOS']
     if chatter > 0:
        print('uvotgrism version :',__version__) 
        print(' Position RA,DEC  : ',RA,' ',DEC)
        print(' grism file       : ',ugufile+'['+str(ext)+']')
        print(' attitude file    : ',attfile)
        print(' wheel position   : ',wheelpos)
        print('======================================')
        
   if ((test == None) ^ (test == 'caltable') ^ (test == "flux") ^ (test == "simple")^(test == 'simple_corrected')):
      if (test == None) ^ (test == 'flux') ^ (test == "simple")^(test == 'simple_corrected'): 
         (anker_uvw1, anker_as, anker_field) = findAnker(filestub, ext, RA, DEC, wheelpos=wheelpos,\
                       lfilter=lfilter, method=None, mode=test, attfile=attfile, extf=lfilt1_ext, chatter=chatter)      
      else:
         anker_field, wheelpos = testparam
                                                 
      xfield = anker_field[0]
      yfield = anker_field[1]
      if chatter > 0:
         print(' input ray coordinate xfield,yfield: ',xfield, yfield)
      C_min1 = []        
      (C_zero, C_1, C_2, C_3, zmxflx, xpix, ypix, zmxdis, zmxwav, wave, theta) = getZmx( \
          xfield, yfield, wheelpos, mode=getzmxmode,chatter=chatter, \
          test=test, kk=splineorder,s=smooth,wpixscale=wpixscale,xlimits=xlimits)
          # note includes align_zemax() call
      # xpix and ypix are the zemax result      
      print("getZmx output: ================")
      print("C_zero:",C_zero)
      print("C_1   :",C_1)
      print("C_2   :",C_2)
      print("C_3   :",C_3)
      #print "C_min1:",C_min1
      print("xpix  :",xpix)
      print("ypix  :",ypix)
      print("zmxdis:",zmxdis)
      print("zmxwav:",zmxwav)
      print("zmxflx:",zmxflx)
      print("wave  :",wave)
      print("theta :",theta)
      print("===============================")  
      coef_zmx = C_1
      pix_zmx = zmxdis[1]
      wav_zmx = zmxwav[1] 
      if wpixscale == None:
         wpixscale = get_wpixscale(xfield,xfield,wheelpos=wheelpos,mode='field')
      if wave.mean() < 10: wave = wave*1e4  # units nm -> A
      q260 = N.where(wave == 2600)
      q420 = N.where(wave == 4200)
      if wheelpos < 500: 
         xdet = xpix[1,q260].squeeze()   
         ydet = ypix[1,q260].squeeze() 
         ang_bs_nom = 28.5 #
         ang_bs_clo = 28.5+6.45
      else: 
         xdet = xpix[1,q420].squeeze()   
         ydet = ypix[1,q420].squeeze() 
         ang_bs_nom = 180.0-150.7
         ang_bs_clo = 180.0-140.75
         
      anker = N.array([xdet,ydet])
      if chatter > 0: print('anker = ',anker)
      if test == 'caltable': print('uncorrected getZmx interpolated anker pos. = ',anker)
      uncorr_anker = anker.copy().squeeze()
      #
      #  refine the anchor position with the empirical correction (bilinear scaling factor - Anker Positie) 
      #
      if test == "flux":
         anker = correctAnkPos(anker[0],anker[1],wheelpos=wheelpos,test="flux")
         anker = N.array(anker).squeeze()
      elif test == "simple":
         anker = N.asarray(anker)        
      elif test == "simple_corrected":
         anker = correctAnkPos(anker[0],anker[1],wheelpos=wheelpos,test="simple_corrected")
         anker = N.asarray(anker)        
      else: 
         anker = correctAnkPos(anker[0],anker[1],wheelpos=wheelpos)
         anker = N.array(anker).squeeze()
         
      # positioncorrection 
      xpix = xpix + anker[0] - uncorr_anker[0]
      ypix = ypix + anker[0] - uncorr_anker[0]
      #
      #       
      if chatter > 0:
             print('=======================================================================')
             if wheelpos < 300:
                print('getFirstOrder: check wave at 2600 = ',N.where(wave == 2600)[0])
                print('the wavelength supposed to be 2600 A is      : ', wave[q260]) 
                print('The anker at 2600 first order (corrected) is : ',anker,' [DET-pix]')
                print('The anker at 2600 first order (corrected) is : ',anker - [77+27,78],' [DET-image]')
             else:
                print('getFirstOrder: check wave at 4200 = ',N.where(wave == 4200)[0])
                print('the wavelength supposed to be 4200 A is      : ', wave[q420])              
                print('The anker at 4200 first order (corrected) is : ',anker,' [DET-pix]')
                print('The anker at 4200 first order (corrected) is : ',anker - [77+27,78],' [DET-image]')
             print(' ')          
             print('The uncorrected anker position is [DET-pix]  : ',uncorr_anker[0],uncorr_anker[1])
             print('The uncorrected anker position is [DET-image]: ',uncorr_anker[0]-104,uncorr_anker[1]-78)
             print('=======================================================================')    
      if len(anker) != 2 :
             print('problem : anker is not found ')
             return None   
             
      if not clocked:
             if  (theta[1] !=  None) :
                angle = abs(theta[1])
             else: angle = ang_bs_nom
      else:
             if  (theta[1] !=  None) :
                angle = abs(theta[1])
             else: angle = ang_bs_clo
             
   elif test == 'cal': 
      Xphi, Yphi = uvotspec.findInputAngle( RA, DEC, filestub, ext,
           wheelpos=wheelpos, lfilter=lfilter, lfilt2=None, lfilt2_ext=None, 
           method=None, attfile=None,chatter=2)
      anker, anker2, C_1, C_2, theta1, calibdat = uvotgrism.getCalData(Xphi,Yphi,wheelpos,calfile=calfile)
      angle = theta1
      anker_field = N.array([Xphi,Yphi])
      #
      theta=N.array([theta1,theta1,theta1,theta1,theta1]) # use the angle from first order everywhere.
      C_0 = N.zeros(3)
      C_3 = N.zeros(3)
      Cmin1 = N.zeros(3)
      zmxdis = N.array([-450,-300,-200,-100,0,150,300,450,600,750,900,1150,1300,1500,1700])
      zmxwav = N.polyval(C_1,zmxdis)    
      coef_zmx = C_1
      pix_zmx = zmxdis[1]
      wav_zmx = zmxwav[1]
      
   else: 
      print("uvotgrism.getFirstOrder: error in test parameter: only none and cal are allowed\n")
      return   
                
   # got the anchor position and the dispersion - next : extract spectrum    
         
   if chatter > 0: print('first order angle = ',angle)        
   if test == 'caltable': return angle, anker       
   #
   #   extract spectrum
   #
   ankerimg = anker - [77+27,77+1] 
   print("anker in img:",ankerimg)
   (dis, spnet, bg, bg1, bg2, extimg, spimg, spnetimg, offset, ank_c) = \
                extractSpecImg(ugufile,ext,ankerimg,angle,spwid=spextwidth)
   #   grism det image 
   img = pyfits.getdata(ugufile, ext)
   
   #   zemax model points (scaled for wavelength comparison) on grism
   if wheelpos < 300:
      H_lines=( 6563.35460346,  4861.74415071,  4340.84299171,  4102.09662716, \
        3970.42438975,  3889.39532057,  3835.72671631,  3798.23761775,      \
        3770.96821946,  3750.48834484,  3734.70346123)
      WC_lines = (1723,1908,2297,2405,2530,2595,2600.,2699,2733,2906,3070,4649)
      WN_lines = (1720,1750,2306,2385,2511,2733,2984,3203,4058,4609,4686,4859,5412,5876,6560)
   else:
       H_lines=( 6563.35460346,  4861.74415071,  4340.84299171,  4102.09662716, \
        3970.42438975,  3889.39532057,  3835.72671631,  3798.23761775,      \
        3770.96821946,  3750.48834484,  3734.70346123)  
       WC_lines = (2699,2733,2906,3070,4649,5093,5473,5696,5803,6741)   
       WN_lines = (2733,2984,3203,4058,4609,4686,4859,5412,5876,6560)
       
   # pixel scale factor for use in plot? 
   if wheelpos < 300: wpixscale = 0.960
   else: wpixscale = 1.00
       
   if ((test == None) | (test == "flux") | (test == "simple")) & (wav_zmx != None):
      from uvotpy.uvotmisc import WC_zemaxlines, hydrogenlines
      reswc = WC_zemaxlines(zmxdis[1],wav_zmx,xpix,ypix,wave,c_offset=(-104,-78),lines=WC_lines,wpixscale=wpixscale,wheelpos=wheelpos)
      resh = hydrogenlines(zmxdis[1],zmxwav[1],xpix,ypix,wave,c_offset=(-104,-78),wpixscale=wpixscale,wheelpos=wheelpos)
      (xhlines, yhlines, hdislines, hlines) = resh
      (xplines, yplines,  dislines, lines)  = reswc
   if (wav_zmx != None):  xplines, yplines,  dislines, lines = 0,0,0,0
   if chatter > 0:
      from pylab import figure,clf, imshow, contour, plot,xlim, ylim
      #   make plot of spectrum
      figure(1); clf(); plot(dis,spnet,ls='steps',label='spectrum'); plot(dis,3.*N.sqrt(bg),'k--',label='3-sigma BG noise')
      legend(loc=0)
      #   make plot of model on image
      xa = N.where( (dis < 1400) & (dis > -300) )
      bga = bg.copy()
      figure(3); clf()
      imshow(img,vmin=bga.mean(),vmax=bga.mean()*10)
      levs = N.array([6,15,27,42,56,70,90,120]) * bg.mean()/7.
      contour(img,levels=levs)
      if ((test == None) | (test == "flux") | (test =="simple")) & (wav_zmx != None):
         plot(xplines,yplines,'oy')
         plot(xhlines,yhlines,'oc')
         xlim(0,2050)
         ylim(0,2050)
#   
   if measured_file != None:
      tab = uvotplot.rdTab(measured_file)
      pixZeroOffset, wav_scale, pix_scale = getScale(tab[:,1],tab[:,0],zmxwav[1],zmxdis[1],order=4)
   
   if (test == None) | (test == "flux") | (test == "simple")| (test == "simple_corrected"):
      g1 = N.max(ank_c[1]-370,0); g2 = N.min( [ank_c[1]+1000,dis.shape[0]-3] )
      if chatter > 5: print(g1,g2,ank_c[1]+1000,dis.shape[0]-3,dis.shape)
      aa = N.arange(g1,g2)
      w1 = polyval(C_1,dis[aa])
      RC1 = polyfit(w1, dis[aa],4)
      print('distance 4649-4200 = ', polyval(RC1,4649.) - polyval(RC1,4200.), '   pix ')
      
      if pix_zmx != None: pix_zmx = wpixscale * pix_zmx
      return ( (dis, spnet, angle, anker, coef_zmx, pix_zmx, wav_zmx), 
            (anker_uvw1,anker_as, anker_field,ank_c),
            (bg, bg1, bg2, extimg, spimg, spnetimg, offset) ,
            (C_zero,C_1,C_2,C_3,), 
            (xpix,ypix,zmxflx, zmxdis,zmxwav, wave,theta), 
            (img, xplines, yplines, dislines, lines), hdr) 
   else: return ( (dis, spnet, angle, anker, anker2, anker_field,ank_c),
            (bg, bg1, bg2, extimg, spimg, spnetimg, offset) ,
            (C_1,C_2, img, H_lines, WC_lines), hdr )        

def align_zemax(xpix,ypix,wheelpos,test=None,chatter=0):
      ''' 
      Align the zemax model using xpix, ypix, to the anchor position
      NOT using fieldpos to align the boresight with fieldpos=[xfield,yfield] = dphi.
      Input filter wheel position
            model xpix, ypix 
      Output xpix, ypix (corrected)
      
      '''
      import numpy as np
            
      #  The zemax model has not yet been corrected for the actual boresight offset 
      if test == "flux":
      #  Rather than lining up the boresight, the flux drop-off at the edges is lined up 2012-12-23
         if wheelpos == 200:
            dxpix = 0. 
            dypix = 0. 
         elif wheelpos == 160: 
            dxpix =    -60.0 # simple model   
            dypix =    0.
         elif wheelpos == 955:
            dxpix =  0.
            dypix =  0.
         elif wheelpos == 1000:  
            dxpix =  0. 
            dypix =  0.
            
      elif test == "simple" or test =="simple_corrected":           
         if wheelpos == 200:
            dxpix =  -60.  # simple model  (-140,0 ?)
            dypix =    0. 
         elif wheelpos == 160: 
            dxpix =    -60.0 # simple model   
            dypix =    0.
         elif wheelpos == 955:
            dxpix =  0.
            dypix =  0.
         elif wheelpos == 1000:  
            dxpix =  0. 
            dypix =  0.
      
      else:
      #  The zemax model has not yet been corrected for the actual boresight offset 
         if wheelpos == 200:
            dxpix = -27 
            dypix = -1 
         elif wheelpos == 160: 
            dxpix =  53.36 
            dypix = 136.20
         elif wheelpos == 955:
            dxpix =  0.
            dypix =  0.
         elif wheelpos == 1000:  
            dxpix =  0. 
            dypix =  0.

      if chatter > 1: print("align_zemax X,Y pixel shift = ",dxpix,",",dypix)        
      xpix = xpix - dxpix
      ypix = ypix - dypix      
      return xpix, ypix      


def align_zemax_field(xfield,yfield,wheelpos,test=None,chatter=0):
      ''' 
      Align the zemax model using fieldpos to align the boresight with fieldpos=[xfield,yfield] = dphi.
      Input filter wheel position
      
      '''
      import numpy as np
      
      # correct the field positions with xfield, yfield numpy arrays
      if ((xfield != None) | (yfield != None)) & (wheelpos == 160) & (test == 'flux'): 
         dphi = np.array([-0.02267,-0.0020953]) # simple model anchor point phi 
         if xfield != None: xfield = xfield + dphi[0]
         if yfield != None: yfield = yfield + dphi[1]
         if chatter > 1: print("field coordinate angular shift (deg) = ",dphi)
      elif ((xfield != None) | (yfield != None)) & (wheelpos == 160) & (test == 'simple_corrected'): 
         dphi = np.array([-0.00199,+0.00080]) # simple model anchor point phi 
         if xfield != None: xfield = xfield + dphi[0]
         if yfield != None: yfield = yfield + dphi[1]
         if chatter > 1: print("field coordinate angular shift (deg) = ",dphi)
      return xfield, yfield


def wrWaveCalFile(wheelpos,xoff=0.0,yoff=0.0,wpixscale=0.960, test=None,tstart = 230531544,wlshift=0.0,\
    zmxmodel="orig",chatter=0,extrapolate=True,resample=1,clobber=False,debug=0):
   ''' Write a fits file with the array of input angles (Xang,Yang) from 
       the boresight in UVW1, and the positions (Xdet,Ydet)_ank of the anker
       at 2600A first order on the UV grism nominal (200) detector image. 
       The positions are from the unscaled zemax model with the following
       correction from correctAnkPos()
       Writes the angle of the first order at the anchor positions.
       Writes the coefficients of a polynomial for the pixel to wavelength 
       conversion. Based on the pixels along the dispersion direction after 
       rotating the image, without distortion or curvature corrections.
       
       if extrapolate = True then the illegal anchor positions and failed 
       dispersion relations are extrapolated from the good array
       
       offsets xoff,yoff to correct boresight position 
          
       the 160  anchor position recommended correction  wlshift=6.1         20090625
       for 1000 anchor position recommended correction  xoff=2.2  yoff=-1.1 20100128
       for 955  anchor position recommended correction wlshift=-8.0 xoff=+3.5 yoff=-2.7 20100419
          but --955  anchor position recommended correction wlshift=-8.0  2011-06-30            
       
 101215  correct fits header keywords (use od DIMn requires that TFORMn has format like '784E', not '1E')
         => need to change second order anchor and dispersion for wheelpos = 160 
 110124  hacks to fix the position of the second order anchor point by using the calibrated distance to the 
         first order anchor point (for wheelpos=160). 
         introduce wpixscale2 for scaling the second order pixels for dispersion        
 110627  the previous version produced wrong results. General cleanup an testing. Fixed error of scaling 
         first order to second order function :(
 121223  change 160(UV clocked) model with rotation (tilt, not image plane) "flux" model. Option zmxmodel="flux".        
 130515  first order anchor and dispersion from "simple_corrected" zemax "flux" model.   
 130604  added array size as variable AZ and resample parameter (is used as power of 2 to grow grid). 
         correctAnkPos(test="simple_corrected") added.
          
   '''    
   version = '130604'
   calversion = '2.6.1'
   ch = chatter
   a = now = datetime.date.today()
   atime = datetime.datetime.now() 
   filtername = '' 
   filtcom = ''
   datestring = a.isoformat()[0:4]+a.isoformat()[5:7]+a.isoformat()[8:10]
   log = open('wrCalFile.log','w')
   log.write("wrCalFile(%i, xoff=%f, yoff=%f, wpixscale=%f, test=%s, tstart=%10i,\n wlshift=%f, zmxmodel=%s, chatter=%i, extrapolate=%s, clobber=%s, debug=%i)\n\n" 
       %(wheelpos,xoff,yoff,wpixscale,test,tstart,wlshift,\
       zmxmodel,chatter,extrapolate,clobber,debug))
   AZ = 28    
   
   if wheelpos == 160:  
      filtername = 'uc160'
      filterheader = 'UGRISM'
      filtcom = 'UV grism Clocked'
      #xoff = -53.36     # aligns zemax model anchor to observed anchor 
      #yoff = -136.20
      if zmxmodel == "flux" :
         xoff, yoff = align_zemax(0,0,160,test="flux")
      elif zmxmodel == "simple_corrected":
         xoff, yoff = 0,0        
      else:
         xoff, yoff = align_zemax(0,0,160,)              
      extname1 = 'WAVECALUGRISM160'
      wave_ank = 0.26
      calfilename = 'swugu0160wcal20041120v002.fits'
   elif wheelpos == 200:  
      filtername = 'ug200'
      filterheader = 'UGRISM'
      filtcom = 'UV grism Nominal'
      extname1 = 'WAVECALUGRISM200'
      wave_ank = 0.26
      calfilename = 'swugu0200wcal20041120v001.fits'
   elif wheelpos == 1000: 
      filtername = 'vg1000'
      filterheader = 'VGRISM'
      filtcom = 'V grism Nominal'
      extname1 = 'WAVECALVGRISM1000'
      wave_ank = 0.42
      calfilename = 'swugv1000wcal20041120v001.fits'
   elif wheelpos == 955:  
      filtername = 'vc955'
      filterheader = 'VGRISM'
      filtcom = 'V grism Clocked'
      extname1 = 'WAVECALVGRISM955'
      wave_ank = 0.42
      calfilename = 'swugv0955wcal20041120v001.fits'
   else:
      print('wrCalFile(wheelpos, xoff, yoff, wpixscale, test) - wheelpos parameter is required input')
      return None   
   if wlshift != 0.:
      wlsh = '_wlshift'+str(0.1*int(wlshift*10)) 
   else: 
      wlsh =''      
      
   filename = 'swwavcal'+datestring+'_'+calversion+'_mssl_'+filtername+wlsh+'.fits'

   #  get zemax model 
   
   if zmxmodel == "orig":
      # this is the first version of the calibration (by august 2011)
      (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, ypix, 
         valid) = rdzemax(wheelpos)
   elif zmxmodel == "flux":
      # using the flux model calculation (never completed)
      (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, ypix, 
         flux, valid) = read_zemax_flux(wheelpos)
   elif zmxmodel == "simple_corrected":
      # for the new wavecal of the clocked uv grism the starting point 
      # is to first scale ALL pixels by a fixed factor and shift the 
      # model so it aligns better with the observed drop in flux 
      # 
      (wave, order, xf, yf, xpix, ypix, flux, valid) = zemax_simply_corrected(wheelpos)
      valid[xpix == -450] = False
      
   # define the size of the axis
   norder,nwave,npos = xpix.shape
                 
   #  find index for the reference wavelength; put data in arrays
   
   q_ank = (N.where(wave == wave_ank) )[0][0]
   xp  = (xpix[1,q_ank,:]).squeeze() + xoff    # anchor 1st order X
   yp  = (ypix[1,q_ank,:]).squeeze() + yoff    #                  Y
   xp2 = (xpix[2,q_ank,:]).squeeze() + xoff    # anchor 2nd order X
   yp2 = (ypix[2,q_ank,:]).squeeze() + yoff    #                  Y 
   if zmxmodel != "simple_corrected":
       xf  = (xfield[1,q_ank,:]).squeeze()         # field coordinate X
       yf  = (yfield[1,q_ank,:]).squeeze()         #                  Y
   wave = wave * 1.0e4                         # wave um -> Angstroem
   xpix2 = xpix.copy()  # used for scaling disp second order pixels
   ypix2 = ypix.copy()

   if chatter > 0: 
      log.write('wrCalfile: index reference wavelength in zemax data = %i\n'%(q_ank))  
   
   #  compute the input angles from the field coordinates 
   #  ****Base the boresight on the old boresight and then 
   #  provide a correction at the end.****
   
   if wheelpos == 200:
      rotangle = -64.6
      scale = 6550.4
      bore = boresight(filter=filtername,r2d=0.,date=000000)  # adopted boresight calibration in raw coordinates
      boreimg = N.array(bore) - [27,1]                        # boresight in image pixel coordinates
      rxf ,ryf = uvotmisc.uvotrotvec(xf,yf,rotangle)          # rotate field coordinates
      xshift, yshift = (1023.5-8-boreimg[0])/scale,(1023.5+4-boreimg[1])/scale # find shift angles
      rx, ry = rxf + xshift, ryf + yshift                   # this is the input angle
   elif (wheelpos == 160) & (zmxmodel == "simple_corrected") :
      rotangle = 0.
      scale = 6550.4
      bore = boresight(filter=filtername, date=000000)      # boresight in DET coordinates
      boreimg = N.array(bore) - N.array([104,78])           # boresight in image pixel coordinates cal. observations
      rxf, ryf = xf, yf                                   # no rotation in correct zemax ("flux") model 
      log.write( "applying align_zemax_field correction - simple_corrected\n")
      rxf, ryf = xfzmx,yfzmx = align_zemax_field(
                 xf,yf,wheelpos,test="simple_corrected",chatter=chatter)
      xshift, yshift = 0.0, 0.0                             # find shift angles  offset on pix matches up things
      rx, ry = rxf + xshift, ryf + yshift                   # this is the input angle        
   elif ((wheelpos == 160) & (zmxmodel != 'simple_corrected'))^ (wheelpos == 955) ^ (wheelpos == 1000) :
      rotangle = -64.6
      scale = 6550.4
      bore = boresight(filter=filtername, date=000000)      # boresight in DET coordinates
      boreimg = N.array(bore) - N.array([104,78])           # boresight in image pixel coordinates cal. observations
      if (zmxmodel == "orig") :
        rxf ,ryf = uvotmisc.uvotrotvec(xf,yf,rotangle)      # rotate field coordinates
        log.write("applying rotation of zemax model \n")
      else:
        rxf, ryf = xf, yf                                   # no rotation in correct zemax ("flux") model       
        log.write( "applying align_zemax_field correction\n")
      xshift, yshift = 0.0, 0.0                             # find shift angles  offset on pix matches up things
      rx, ry = rxf + xshift, ryf + yshift                   # this is the input angle   
   else:
      raise IOError( "error: wrong wheelpos entry")
      
             
   #  correct anchor position for anchor distortion  
   if wheelpos == 200:
      xp1z, yp1z = correctAnkPos(xp+27,yp+1,mode=zmxmodel)
      log.write( 'correctAnkPos done 200\n')
   elif (wheelpos == 160) & (zmxmodel != "simple_corrected"): 
      xp1z, yp1z = correctAnkPos(xp,yp,wheelpos=160,test=zmxmodel)
      log.write( 'correctAnkPos done 160\n')
   elif (wheelpos == 160) & (zmxmodel == "simple_corrected"): 
      xp1z, yp1z = correctAnkPos(xp,yp,wheelpos=160,test=zmxmodel)
      dx1 = xp1z - xp
      dy1 = yp1z - yp
      xp1z[xp == -450.] = N.NaN
      yp1z[yp == -450.] = N.NaN
      xp1z = extrapolateArray(xf,yf,xp1z)
      yp1z = extrapolateArray(xf,yf,yp1z)
      log.write( 'correctAnkPos  simple_corrected 160\n')
   elif wheelpos == 1000:
      xp1z, yp1z = correctAnkPos(xp,yp,wheelpos=1000,mode=zmxmodel) 
      log.write('correctAnkPos done 1000\n')
   elif wheelpos == 955:
      qgd = N.where(xp < 2200.0)
      log.write("******** indices of good points qgd: ",qgd)
      xp1z = xp
      yp1z = yp
      xpi = xp[qgd].copy()
      ypi = yp[qgd].copy()
      for ii in range(len(qgd[0])):
         xpi_in = xpi[ii] ; ypi_in = ypi[ii]
         xpi[ii],ypi[ii] = correctAnkPos(xpi[ii],ypi[ii],wheelpos=955,mode=zmxmodel) 
         if chatter > 0: 
            log.write( "correcting anker position with (%8.2f , %8.2f ) "+\
            "for (%8.2f,%8.2f) \n" % (xpi_in, ypi_in, xpi[ii], ypi[ii]) )
      xp1z[qgd] = xpi
      yp1z[qgd] = ypi
      log.write( 'correctAnkPos done 955 \n' )
   # now load the corrected anchor values back into xp, yp 
   xp, yp = xp1z.reshape(AZ*AZ), yp1z.reshape(AZ*AZ)              
   #
   #  find the slope of the first order spectrum near wave_ank (i.e., 2600 A or 4200 A)
   #
   sp_ang = (ypix[1,q_ank-1,:] - ypix[1,q_ank+1,:]) / ( 
             xpix[1,q_ank-1,:] - xpix[1,q_ank+1,:])
   sp_ang[:] = N.abs(N.arctan(sp_ang.squeeze())/N.pi*180.0)
   sp_ang[sp_ang < 27.0] = N.NaN 
   sp_ang[sp_ang > 37.0] = N.NaN 
   sp_ang = extrapolateArray(xf,yf,sp_ang)
   sp_ang = sp_ang.reshape(AZ,AZ)        
   
   if wheelpos == 200:   
      # fix bad automated angles (use extrapolation later on? ) 
      sp_ang[N.where(N.isnan(sp_ang))] = -28.4
      sp_ang[26,27] = -30.12 ; sp_ang[27,27]=-30.12
      sp_ang[3,27] = -28.6 
      sp_ang[2,26] = sp_ang[2,27] = -28.52
      sp_ang[1,26] = sp_ang[1,27] = -28.48
      sp_ang[0,25] = sp_ang[0,26] = sp_ang[0,27] = -28.38 
      
   sp_ang = sp_ang.reshape(AZ*AZ)

      
   #   The input coordinate angles rx, ry
      
   #  in version 1 (test='orig'), we rotated the field positions by -64.6 deg to get a rectangu-
   #  lar lookup grid [or no rotation for zmxmodel="flux"]
   #  then applied the shift based on DET->field map and observed boresight.
   #  114.0 = 1023.5-8-901.5   25.8=1023.5+4-1001.7  
   rx = rx.reshape(AZ,AZ)
   ry = ry.reshape(AZ,AZ)
   
   #  correct for the angle between the model boresight and the actual boresight
   #  (the first term is the grism boresight, the second the model center of distortion
   #  THE BORESIGHT IS IN PIXEL DETECTOR COORDINATES  correct the zemax pix coords
   #  ... and possibly a similar shift in the "flux" models ??
   borex, borey = boresight(filter=filtername,date=0000)
   
   # select the mean column values to get regular array in the field positions rxs, rys    
   rxm = N.zeros(AZ) ; rym = N.zeros(AZ)
   for i in range(AZ):
      rx1 = rx[:,i]
      rxm[i] = rx1[N.isfinite(rx1)].mean()
      ry1 = ry[i,:]
      rym[i] = ry1[N.isfinite(ry1)].mean()   
   if (zmxmodel == 'orig') & (wheelpos == 160): rxm[27] = 0.164   
   rxs = (N.outer(N.ones(AZ),rxm)).reshape(AZ*AZ)
   rys = (N.outer(rym,N.ones(AZ))).reshape(AZ*AZ)
   if chatter > -1:
      log.write( "making the lookup field position arrays regular - verify the mean error abs(new-old) \n")
      log.write( " x-position : %e\n"%((abs(rx.reshape(AZ*AZ)-rxs)).mean()))
      log.write( " y-position : %e\n"%((abs(ry.reshape(AZ*AZ)-rys)).mean()))
      log.write( " if larger than 1e-6 there is cause for concern \n" )

   # flag anchor coordinates that are clearly wrong to 'None'
   xp [ N.where((xp  > 2350.) & (xp  < -100)) ] = N.NaN 
   yp [ N.where((yp  > 2350.) & (yp  < -100)) ] = N.NaN
   xp2[ N.where((xp2 > 2350.) & (xp2 < -100)) ] = N.NaN
   yp2[ N.where((yp2 > 2350.) & (yp2 < -100)) ] = N.NaN

   #   find new position second order by using the scaled distance to the first order anchor and sp_ang
   #   scaled distance is: 
   dis12 = _orderdist12(xp,yp, None, 1, None, xpix, ypix, mode='cal3',chatter = chatter,wheelpos=wheelpos)
   angle = abs(sp_ang.reshape(AZ*AZ)) * (math.pi / 360.0)   
   # -- neglecting the offset dy2 from second order to angle. 
   #    which would be extra angle atan(dy2/dis12) 
   xp2 = xp - dis12 * N.cos( angle )
   yp2 = yp + dis12 * N.sin( angle )
   xp2[N.isnan(xp)] = N.NaN
   yp2[N.isnan(yp)] = N.NaN
   qv = N.where( (xp2 < -100.) ^ (xp2 > 2300) )
   xp2[qv] = N.NaN
   qv = N.where( (yp2 < -100.) ^ (yp2 > 2300) )
   yp2[qv] = N.NaN
   if (chatter > 3):
      log.write( 'xp2 = %s\n'%(str(xp2)) )
      log.write( 'yp2 = %s\n'%(str(yp2)) )
      log.write( '= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n' )
   
   
   #  determine the wavelength dispersion polynomial coefficients and order 
   ch = 1  # chatter set to one now
   #  define arrays for the dispersion coeeficients
   disp10 = N.zeros(AZ*AZ) ; disp11 = N.zeros(AZ*AZ) 
   disp12 = N.zeros(AZ*AZ) ; disp13 = N.zeros(AZ*AZ)
   disp14 = N.zeros(AZ*AZ) ; disp1n = N.zeros(AZ*AZ) 
   disp20 = N.zeros(AZ*AZ) ; disp21 = N.zeros(AZ*AZ)
   disp22 = N.zeros(AZ*AZ) ; disp2n = N.zeros(AZ*AZ)
   
   #    scale the first order dispersion scale factor WPIXSCL
   #    use now: get_wpixscale(xf[i],yf[i], wheelpos=160, mode='field') xf,yf rotated field coordinates
   #    calculate second order wpixscale2 
   #    use same definition for valid pixels ? qv  
   #    wpixscl2 is applied to zemax and is different to 
   #    first order scaled (now xpix, ypix) 
   #    2011-09-20 wpixscale2 now returns scale at 2000A
   #    getDisp returns the correct C_2 at anchor 2600
   
   # grab index all valid zemax points
   qv = N.where(valid == 0) 
   
   #   initialise wpixscale   
   wpixscl  = N.ones(norder*nwave*AZ*AZ).reshape(norder,nwave,AZ*AZ)
   wpixscl2 = N.ones(norder*nwave*AZ*AZ).reshape(norder,nwave,AZ*AZ)
   if (wheelpos == 160) & (zmxmodel == "simple_corrected"): 
       # make no further scaling to all orders except 1        
       if chatter > 0: print("wpixscale correction with test2=fix2")
       ia = 1          
       xpix_ = xpix.copy()
       ypix_ = ypix.copy()
       for ic in range(AZ*AZ):     # grid on detector
           qv_ = N.where(xpix[1,:,ic] != -450.)
           if len(qv_[0]) == 0: continue
           #
           log.write( "wpixscale call with wave[qv_]="+str(wave[qv_[0]])+'\n')
           if len(qv_[0]) > 0:
               try:
                   log.write("wavw, xpix, ypix get_wpixscale input\n"+str(wave[qv_[0]])
                      +"\n"+str(xpix[1,qv_[0],ic])+"\n"+str(ypix[1,qv_[0],ic])+"\n" ) 
                   # this call applies wpixscale1 :      
                   xpix[1,qv_[0],ic], ypix[1,qv_[0],ic], _w, wpixscale1 = get_wpixscale(
                       xpix[1,qv_[0],ic], ypix[1,qv_[0],ic], w=wave[qv_[0]],wheelpos=wheelpos,
                       mode='fix2',chatter=chatter)
                   log.write("wpixscale: pixscale1= "+str(wpixscale1))    
                   log.write("at anchor=(%7.1f,%7.1f)  \n-------------\n"%(xp[ic],yp[ic]) ) 
               except:
                   # main reason for exception should be no anchor point 
                   log.write( 
                   "pixscale1 not applied for point at anchor=(%7.1f,%7.1f)\n------------\n"
                   %(xp[ic],yp[ic]) ) 
                   pass                
           xpix[1,qv_[0],ic] = xpix[1,qv_[0],ic]+ xpix_[1,q_ank,ic]    
           ypix[1,qv_[0],ic] = ypix[1,qv_[0],ic]+ ypix_[1,q_ank,ic] 
       xpix_ = ''
       ypix_ = ''
       for ib in range(nwave):
           qv_ = N.where(xpix[1,ib,:] != -450.)
           xpix2[2,ib,qv_] = xpix2[2,ib,qv_] * 0.95 
           ypix2[2,ib,qv_] = ypix2[2,ib,qv_] * 0.95
   else:   
       wpixscl  = wpixscl * 0.96
       wpixscl2 = wpixscl2 *0.92
       for ia in range(norder):          
           for ib in range(nwave):         
               for ic in range(AZ*AZ):     # grid on detector
                   wpixscl[ia,ib,ic] = (get_wpixscale(rx.reshape(AZ*AZ)[ic],
                       ry.reshape(AZ*AZ)[ic], wheelpos=wheelpos,mode='rotfield') )  
                   wpixscl2[ia,ib,ic] = get_wpixscale2(xf.reshape(AZ*AZ)[ic],
                       yf.reshape(AZ*AZ)[ic], wheelpos=wheelpos,mode='field')  
       xpix [qv] = xpix[qv] * wpixscl[qv] 
       ypix [qv] = ypix[qv] * wpixscl[qv]
       xpix2[qv] = xpix2[qv] * wpixscl2[qv] 
       ypix2[qv] = ypix2[qv] * wpixscl2[qv]
   if debug == 4:
      return rx,ry,wpixscl
      
   
   #  compute the dispersion coefficients
   
   if wheelpos < 500:
      plimit = [-100,2100]
      dsp0 = 2600.0
      dsp1 = 3.2
      dsp2 = 1.6e-3
   else:   
      plimit = [20,2100]
      dsp0 = 4200.0
      dsp1 = 5.8
      dsp2 = 0.0
   mord = 4 # maximum/optimum order of dispersion coefficient solution
   if wheelpos == 955: mord=3
      
   for ix in range(AZ*AZ):
      test1 = (xpix[1,:,ix] < 2350 ) & (ypix[1,:,ix] < 2200 ) & (xpix[1,:,ix] > -50 ) & (ypix[1,:,ix] > -50 ) 
      test2 = (xpix[2,:,ix] < 2350 ) & (ypix[2,:,ix] < 2200 ) & (xpix[2,:,ix] > -50 ) & (ypix[2,:,ix] > -50 ) 
      test3 = (xpix[1,q_ank,ix] < 2350 ) & (ypix[1,q_ank,ix] < 2200 ) & (xpix[1,q_ank,ix] > -50 ) & (ypix[1,q_ank,ix] > -50 ) 
      test4 = (xpix[2,q_ank,ix] < 2350 ) & (ypix[2,q_ank,ix] < 2200 ) & (xpix[2,q_ank,ix] > -50 ) & (ypix[2,q_ank,ix] > -50 ) 
      q1 = N.where( test1 )    # for wavelengths the points on and close to detector
      q2 = N.where( test2 )
      if chatter > 2: 
         log.write( "\n===========   ============    ===========     ============     ============    ==============\n")
         log.write("PROCESSING 1st order DISPERSION FOR POINT NO="+str(ix)+
           "  ank:["+str(xpix[1,q_ank,ix])+" ,"+str(ypix[1,q_ank,ix])+"  wv:"+str(wave[q1])+'\n')
         log.write("xpix[1,q1]:"+str(xpix[1,q1,ix])+'\n')
         log.write("ypix[1,q1]:"+str(ypix[1,q1,ix])+'\n')
         log.write("xpix[2,q2]:"+str(xpix[2,q2,ix])+'\n')
         log.write("ypix[2,q2]:"+str(ypix[2,q2,ix])+'\n\n')
      # make highest order for given number of wavelength points
      # th is angle
      C = N.zeros(5) 
       
      # if 1st order anker buiten normal range, put default answer for dispersion
      if (not test3):  
         disp10[ix] = dsp0 ; disp11[ix] = dsp1 ; disp12[ix] = 0
         disp13[ix] = 0.   ; disp14[ix] = 0.   ; disp1n[ix] = 0
         xpix[1,q_ank,ix] = N.NAN
         ypix[1,q_ank,ix] = N.NAN          
         if extrapolate:
            disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
      else:      
        if ((len(q1[0]) > 5) & (mord == 4)):
          th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q1,ix],ypix[1,q1,ix],wave[q1], wheelpos=wheelpos, fitorder=4, lim = plimit,chatter=ch)
          if di_1 == None:
             disp14[ix] = 0 ; disp13[ix] = 0 ; disp12[ix] = 0
             disp10[ix] = C_1[1] ; disp11[ix] = C_1[0] ; disp1n[ix] = 0     
             if extrapolate:
               disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
          else:
             if (len(C_1) == 5):
                C[5-len(C_1):] = C_1 
                disp10[ix] = C[4] ; disp11[ix] = C[3] ; disp12[ix] = C[2]
                disp13[ix] = C[1] ; disp14[ix] = C[0] ; disp1n[ix] = len(C_1)-1 ; sp_ang[ix] = th_1
             else:  # this is the rare case that the number of valid points were less than the zemax model 
                # due to points outside the detector not following the smooth spectral dispersion 
                C[5-len(C_1):] = C_1 ; disp1n[ix] = len(C_1)-1 ; sp_ang[ix] = th_1
                disp10[ix] = C[4] ; disp11[ix] = C[3] ; disp12[ix] = C[2]
                disp13[ix] = C[1] ; disp14[ix] = C[0] 
        elif len(q1[0]) > 4:
          th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q1,ix],ypix[1,q1,ix],wave[q1], wheelpos=wheelpos, fitorder=3, lim = plimit,chatter=ch)
          if di_1 == None:
             disp13[ix] = 0 ; disp14[ix] = 0 ; disp12[ix] = 0
             disp10[ix] = C_1[1] ; disp11[ix] = C_1[0] ; disp1n[ix] = 0     
             if extrapolate:
               disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
          else:
             C[4-len(C_1):4] = C_1
             disp10[ix] = C[3] ; disp11[ix] = C[2] ; disp12[ix] = C[1]
             disp13[ix] = C[0] ; disp14[ix] = 0 ; disp1n[ix] = len(C_1)-1 ; sp_ang[ix] = th_1
        elif len(q1[0]) > 3:
          th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q1,ix],ypix[1,q1,ix],wave[q1], fitorder=2,wheelpos=wheelpos,  lim = plimit,chatter=ch)
          if di_1 == None:
             disp13[ix] = 0 ; disp14[ix] = 0 ; disp12[ix] = 0
             disp10[ix] = C_1[1] ; disp11[ix] = C_1[0] ; disp1n[ix] = 0     
             if extrapolate:
               disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
          else:
             C[3-len(C_1):3] = C_1
             disp10[ix] = C[2] ; disp11[ix] = C[1] ; disp12[ix] = C[0]
             disp13[ix] = 0 ; disp14[ix] = 0 ; disp1n[ix] = len(C_1)-1 ; sp_ang[ix] = th_1
        elif len(q1[0]) > 2:
          th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q1,ix],ypix[1,q1,ix],wave[q1], fitorder=1, wheelpos=wheelpos, lim = plimit,chatter=ch)
          if di_1 == None:
             disp13[ix] = 0 ; disp14[ix] = 0 ; disp12[ix] = 0
             disp10[ix] = C_1[1] ; disp11[ix] = C_1[0] ; disp1n[ix] = 0     
             if extrapolate:
               disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
          else:
             C[2-len(C_1):2] = C_1
             disp10[ix] = 0 ; disp11[ix] = 0 ; disp12[ix] = 0
             disp13[ix] = C[1] ; disp14[ix] = C[0] ; disp1n[ix] = 1
        else: 
          disp10[ix] = dsp0 ; disp11[ix] = 2.*dsp1 ; disp12[ix] = 0
          disp13[ix] = 0. ; disp14[ix] = 0. ; disp1n[ix] = 0
          if extrapolate:
             disp10[ix] = disp11[ix] = disp12[ix] = disp13[ix] = disp14[ix] = N.NAN
        first_order_disp = C
                  
        # SECOND ORDER
                
        if wheelpos < 500:
           plimit = [-100,2100]
           dsp0 = 2600.0
           dsp1 = 1.9
           dsp2 = 6.4e-4
        else:   
           plimit = [20,2100]
           dsp0 = 4200.0
           dsp1 = 3.0    # guess
           dsp2 = 0.0
        sp_order=2 
        C = N.zeros(3) 
        ian = 3
        if wheelpos > 300: ian = 5 
         
        # test if second order anchor is inside normal range; if not put default answer
        if (not test4) :  
           disp20[ix] = dsp0 ; disp21[ix] = dsp1 ; disp22[ix] = dsp2 ; disp2n[ix] = 0   
           xpix[2,q_ank,ix] = N.NAN
           ypix[2,q_ank,ix] = N.NAN
           xpix2[2,q_ank,ix] = N.NAN
           ypix2[2,q_ank,ix] = N.NAN
           if extrapolate:
              disp20[ix] = disp21[ix] = disp22[ix] = N.NAN
        else:     # second order anchor 'on' detector
           if len(q2[0]) > ian:     # plenty of wavelengths points in zemax 
             th_2, C_2, di_2, wa_2 = getDisp(xpix2[2,q2,ix],ypix2[2,q2,ix],wave[q2], wheelpos=wheelpos, 
                 spectral_order=sp_order, lim = plimit,chatter=ch)
             log.write("getDisp returns 2nd order: "+str(th_2)+str(C_2)+str(di_2)+str(wa_2)+"\n")
             if di_2 == None:       # bad solution
                disp20[ix] = dsp0 ; disp21[ix] = dsp1 ; disp22[ix] = 0 ; disp2n[ix] = 0     
                if extrapolate:
                   disp20[ix] = disp21[ix] = disp22[ix] = N.NAN
             else:                  # good solution
                C[3-len(C_2):3] = C_2
                disp20[ix] = C[2] ; disp21[ix] = C[1] ; disp22[ix] = C[0] ; disp2n[ix] = len(C_2)-1
           elif len(q2[0]) > ian-1:  # second order one point less on detector
             th_2, C_2, di_2, wa_2 = getDisp(xpix2[2,q2,ix],ypix2[2,q2,ix],wave[q2], fitorder=1,wheelpos=wheelpos,  
                    spectral_order=sp_order,  lim = plimit,chatter=ch)
             if di_2 == None:        # no solution
                disp20[ix] = dsp0 ; disp21[ix] = dsp1 ; disp22[ix] = dsp2 ; disp2n[ix] = 0          
                if extrapolate:
                   disp20[ix] = disp21[ix] = disp22[ix] = N.NAN
             else:                   # store solution
                C[2-len(C_2):2] = C_2
                disp20[ix] = C_2[1] ; disp21[ix] = C_2[0] ; disp22[ix] =0 ; disp2n[ix] = 1
           else:                     # too few second order points found
             disp20[ix] = dsp0 ; disp21[ix] = dsp1 ; disp22[ix] = dsp2 ; disp2n[ix] = 0   
             if extrapolate:
                disp20[ix] = disp21[ix] = disp22[ix] = N.NAN
                
   disp20[disp20 < 2000]=N.NaN          
   disp20[disp20 > 3000]=N.NaN

   if (wheelpos == 160) & (zmxmodel == "simple_corrected"):     
        disp22 = 2.2*disp22
        disp22[disp22 < 3.1e-5] = 3.1e-5
        disp22[disp22 > 6.4e-4] = 6.4e-4
                        
   #  shift the wavelength zero point to optimize this calibration
   #     
   disp10 = disp10 + wlshift    
   
   #  extrapolate dispersion for large rx
   
   if debug == 1:
     return rxs,rys,sp_ang,disp10,disp11,disp12,disp13,disp14,disp20,disp21,disp22,wpixscl
   
   if ((wheelpos == 160) ^ (wheelpos == 200)):
     qv = N.where( (disp10 > 2630.) ^ (disp10 < 2570) )
     disp10[qv] = N.NaN
     disp11[qv] = N.NaN
     disp12[qv] = N.NaN
     disp13[qv] = N.NaN
     disp14[qv] = N.NaN
     qv = N.where( (disp20 > 2800.) ^ (disp20 < 2400) ^ (disp21 < 1.0) ^ (disp21 > 4))
     disp20[qv] = N.NaN
     disp21[qv] = N.NaN
     disp22[qv] = N.NaN
     
   if ((wheelpos == 160) ^ (wheelpos == 200) ^ (wheelpos == 1000) ^ (wheelpos == 955)):
      rx = rxs.reshape(AZ,AZ)
      ry = rys.reshape(AZ,AZ)
      angx = sp_ang.reshape(AZ,AZ)
      disp10x = disp10.reshape(AZ,AZ).copy()
      disp11x = disp11.reshape(AZ,AZ).copy()
      disp12x = disp12.reshape(AZ,AZ).copy()
      disp13x = disp13.reshape(AZ,AZ).copy()
      disp14x = disp14.reshape(AZ,AZ).copy()
      disp20x = disp20.reshape(AZ,AZ).copy()
      disp21x = disp21.reshape(AZ,AZ).copy()
      disp22x = disp22.reshape(AZ,AZ).copy()
      if wheelpos == 160:
      #  fix the bad disp for large rx by extrapolating the inner region
         k = (24,25,26,27)
         disp10x[:,k] = None
         disp11x[:,k] = None
         disp12x[:,k] = None
         disp13x[:,k] = None
         disp14x[:,k] = None
         disp20x[:,k] = None
         disp21x[:,k] = None
         disp22x[:,k] = None    
      #
      xpx = extrapolateArray(rx,ry,xp.copy().reshape(AZ,AZ))
      ypx = extrapolateArray(rx,ry,yp.copy().reshape(AZ,AZ))
      angx = -1*( extrapolateArray(rx,ry,angx) )
      disp10x = extrapolateArray(rx,ry,disp10x)
      disp11x = extrapolateArray(rx,ry,disp11x)
      disp12x = extrapolateArray(rx,ry,disp12x)
      disp13x = extrapolateArray(rx,ry,disp13x)
      disp14x = extrapolateArray(rx,ry,disp14x)
      xp2x = extrapolateArray(rx,ry,xp2.copy().reshape(AZ,AZ))
      yp2x = extrapolateArray(rx,ry,yp2.copy().reshape(AZ,28))
      disp20x = extrapolateArray(rx,ry,disp20x)
      disp21x = extrapolateArray(rx,ry,disp21x)
      disp22x = extrapolateArray(rx,ry,disp22x)   
      
   # define angle as positive 
   
   angx = abs(angx)   
   
   # create primary header 
   
   hdu0 = pyfits.PrimaryHDU()
   hdu0.header.update('TELESCOP','SWIFT   ',comment='Telescope (mission) name')
   hdu0.header.update('INSTRUME','UVOTA   ',comment='Instrument Name')
   hdu0.header.update('DETNAM',wheelpos,comment='Filter Wheel Position')
   hdu0.header.update('FILTER',filterheader,comment='Filter value')
   hdu0.header.update('CREATED','written by uvotcal.wrCalFile version:'+version)
   hdu0.header.update('AUTHOR','NPM Kuin (MSSL)')
   #hdu0.header.update('CALVERS',calversion)
   #hdu0.header.update('ORDERS','12',comment='list of orders included')
   #hdu0.header.update('ROTANGLE',rotangle,comment='input phi_x,phi_y by this angle (deg)')
      
   hdulist=pyfits.HDUList([hdu0])
   #
   # write the results to bintable FITS file
   #
   col1 = pyfits.Column(name='PHI_X',   format='1E',unit='deg',array=rxs.reshape(AZ*AZ))#,comment='X coordinate input angle')
   col2 = pyfits.Column(name='PHI_Y',   format='1E',unit='deg',array=rys.reshape(AZ*AZ))#,comment='Y coordinate input angle')
   col5 = pyfits.Column(name='SP1SLOPE',format='1E',unit='deg',array=angx.reshape(AZ*AZ))#,comment='reference angle spectrum') 
   col3 = pyfits.Column(name='DETX1ANK',format='1E',unit='pixel',array=xpx.reshape(AZ*AZ))#,comment='X coordinate anchor point 1st order')
   col4 = pyfits.Column(name='DETY1ANK',format='1E',unit='pixel',array=ypx.reshape(AZ*AZ))#,comment='Y coordinate anchor point 1st order')
   col6 = pyfits.Column(name='DISP1_0', format='1E',array=disp10x.reshape(AZ*AZ))#,comment='polynomial coefficient 1st order dispersion')
   col7 = pyfits.Column(name='DISP1_1', format='1E',array=disp11x.reshape(AZ*AZ))#,comment='polynomial coefficient 1st order dispersion')
   col8 = pyfits.Column(name='DISP1_2', format='1E',array=disp12x.reshape(AZ*AZ))#,comment='polynomial coefficient 1st order dispersion')
   col9 = pyfits.Column(name='DISP1_3', format='1E',array=disp13x.reshape(AZ*AZ))#,comment='polynomial coefficient 1st order dispersion')
   # enter option to remove empty column for vgrism ?
   if max(disp1n) == 4: 
      col10 = pyfits.Column(name='DISP1_4',format='1E',array=disp14x.reshape(AZ*AZ))#,comment='polynomial coefficient 1st order dispersion')
   col11 = pyfits.Column(name='DISP1_N', format='1I',array=disp1n.reshape(AZ*AZ))#,comment='order polynomial or 0 for extrapolated polynomial')
   col12 = pyfits.Column(name='DETX2ANK',format='1E',unit='pixel',array=xp2x.reshape(AZ*AZ))#,comment='X coordinate anchor point 2nd order') 
   col13 = pyfits.Column(name='DETY2ANK',format='1E',unit='pixel',array=yp2x.reshape(AZ*AZ))#,comment='Y coordinate anchor point 2nd order')
   col14 = pyfits.Column(name='DISP2_0', format='1E',array=disp20x.reshape(AZ*AZ))#,comment='polynomial coefficient 2nd order dispersion')
   col15 = pyfits.Column(name='DISP2_1', format='1E',array=disp21x.reshape(AZ*AZ))#,comment='polynomial coefficient 2nd order dispersion')
   col16 = pyfits.Column(name='DISP2_2', format='1E',array=disp22x.reshape(AZ*AZ))#,comment='polynomial coefficient 2nd order dispersion')
   col17 = pyfits.Column(name='DISP2_N', format='1I',array=disp2n.reshape(AZ*AZ))#,comment='order polynomial or 0 for extrapolated polynomial')
   if max(disp1n) == 4:
      cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,
       col11,col12,col13,col14,col15,col16,col17])
   else:    
      cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,
       col11,col12,col13,col14,col15,col16,col17])
   #
   thdu = pyfits.new_table(cols)
   hdulist.append(thdu)
   thdu.header.update('EXTNAME',extname1,comment='name of this binary extension')
   thdu.header.update('TELESCOP','SWIFT   ',comment='Telescope (mission) name')
   thdu.header.update('INSTRUME','UVOTA   ',comment='Instrument Name')
   thdu.header.update('DETNAM',wheelpos,comment='Filter Wheel Position')
   thdu.header.update('FILTER',filterheader,comment='Filter Value')
   thdu.header.update('ORIGIN','MSSL',comment='Source of FITS file')
   thdu.header.update('CREATOR','N.P.M. Kuin (UCL/MSSL)',comment='Creator')
   thdu.header.update('VERSION',2,comment='Extension version number')
   thdu.header.update('FILENAME',calfilename,comment='File name')
   thdu.header.update('CONTENT','WAVELENGTH_CALIBRATION',comment='File content')
   thdu.header.update('CCLS0001','BCF',comment='Dataset is a Basic Calibration File')
   thdu.header.update('CCNM0001','WAVECAL',comment='Type of Calibration data')
   thdu.header.update('CDES0001',filterheader+' '+str(wheelpos)+' WAVELENGTH CALIBRATION',comment='Description')
   thdu.header.update('CDTP0001','DATA',comment='Calibration file contains data')
   thdu.header.update('CVSD0001','2004-11-20',comment='UTC date when calibration should first be used')
   thdu.header.update('CVST0001','00:00:00',comment='UTC time when calibration should first be used')
   thdu.header.update('CBD10001','FILTER('+filterheader+')',comment='Parameter boundary')
   thdu.header.update('CBD20001','WHEELPOS('+str(wheelpos)+')',comment='Parameter boundary')
   thdu.header.update('ORDERS',12,comment='List of orders included')
   if wheelpos > 500:
      thdu.header.update('HISTORY','second order from Zemax model without further calibration')
      thdu.header.update('BOXG_420',borex,comment='grism boresight X at 420nm first order')
      thdu.header.update('BOYG_420',borey,comment='grism boresight Y at 420nm first order') 
   else:
      thdu.header.update('HISTORY','Preliminary second order calibration')
      thdu.header.update('BOXG_260',borex,comment='Reference grism boresight X at 260nm 1st order')
      thdu.header.update('BOYG_260',borey,comment='Reference grism boresight Y at 260nm 1st order') 

   thdu.header.update('CALVERS',calversion)
   thdu.header.update('CREATED','N.P.M. Kuin (UCL/MSSL) '+atime.isoformat()[0:19])
   thdu.header.update('HISTORY','2013-05 revised zemax model; 1st order recalibrated ')
   thdu.header.add_comment('')   
   thdu.header.add_comment('Based on scaled Zemax OM optical model data.')
   thdu.header.add_comment('Anchor point can be offset by as much as 10 pixels from spectrum, ')
   thdu.header.add_comment('and needs to be projected onto the spectrum. Resulting wavelength ')
   thdu.header.add_comment('accuracy is typically 15A, but can be as large as 50A in extreme ')
   thdu.header.add_comment('cases. ')
   thdu.header.add_comment('')   
   if extrapolate:
      thdu.header.add_comment('Missing anchor points near edges were extrapolated. If the dispersion')
      thdu.header.add_comment('was extrapolated DISP?_N has been set to zero.')
   else:
      thdu.header.add_comment('Missing anchor points near edges were set to -100,-100')
      thdu.header.add_comment('When no dispersion solution, DISP?_N set to zero')
   thdu.header.add_comment('')   
   thdu.header.add_comment('Coordinate system for anchor position:')
   thdu.header.add_comment('The coordinate system is in modified DET coordinates derived from ')
   thdu.header.add_comment('DET coordinates (mm) as defined in the TELDEF file by shifting the')
   thdu.header.add_comment('origen by adding [1100.5,1100.5] pix and converting to pixels using a')
   thdu.header.add_comment('pixel scale of 0.5496"/pixel.')
   thdu.header.add_comment('')  
   thdu.header.update('DATE',atime.isoformat(),comment='Creation date')
   
   thdu.header.add_comment('',after='EXTNAME')
   thdu.header.add_comment('',after='ORIGIN')
   thdu.header.add_comment('',after='CONTENT')   
   thdu.header.add_comment('',before='ORDERS')
   
   hdulist.writeto(filename,clobber=clobber)
   #
   print('wrote file: '+filename)
   log.write( 'wrote file: '+filename+'\n')
   log.close()
   return (rxs,rys,xpx,ypx,angx),(disp10x,disp11x,disp12x,disp13x,disp14x,
           disp1n),(disp20x,disp21x,disp22x,disp2n),(xf,yf,xpix,ypix,valid), wpixscl


   
def _getObservedScale(wav_obs,pix_obs,wav_zmx,pix_zmx,order=4,chatter=1,wheelpos=200):
   '''  Determine the scale factor from comparing the observed data 
   to the zemax model
   
   Use:
   pixZeroOffset, wav_scale, pix_scale = uvotgrism.getScale(wav_obs,pix_obs,wav_zmx,pix_zmx,order=4)
   '''  
   if wheelpos < 500:
      ref_wave = 2600.
   else:
      ref_wave = 4200.
   w_1 = ref_wave-700.
   w_2 = ref_wave+900      
   RevCoefObs = polyfit(wav_obs,pix_obs,order)
   pixZeroOffset = polyval(RevCoefObs,ref_wave)
   if chatter > 0: print('zero point offset along dispersion direction of ',ref_wave,'A Dlam(pix) = ',pixZeroOffset) 
   RevCoefZmx = polyfit(wav_zmx,pix_zmx,order)
   ww = w_1 + N.arange(0.,1.,0.02)*(w_2 - w_1)
   newpixobs = polyval(RevCoefObs,ww)
   newpixzmx = polyval(RevCoefZmx,ww)
   pix_sca_coef = polyfit(newpixobs,newpixzmx,1)
   if chatter > 0: 
     print('the wavelength range used is ',w_1,' to ',w_2)
     print('pixel scale fit pixZmx(pixObs) linear fit coefficients: ',pix_sca_coef)
   pix_scale = pix_sca_coef[0]
   CoefObs = polyfit(pix_obs,wav_obs,order)
   CoefZmx = polyfit(pix_zmx,wav_zmx,order)
   dd = N.array(list(range(50))) * 12 - 300
   newwavobs = polyval(CoefObs,dd)
   newwavzmx = polyval(CoefZmx,dd)
   wav_scale_coef = polyfit(newwavobs,newwavzmx,1)
   if chatter > 0: 
     print('wavelength scale fit wavZmx(wavObs) linear coefficients: ',wav_scale_coef)
   wav_scale = wav_scale_coef[0]
   return pixZeroOffset, wav_scale, pix_scale


def getZmx(Xfield,Yfield,wheelpos,mode='spline',chatter=3,wpixscale=None, \
    wpixscale2=None, kk=3,s=None,test='simple',test2='fix2',xlimits=[10,2100]):
   '''
   This routine, using inputs Xfield and Yfield positions, 
   interpolates the zemax model grid, using either 
   bilinear (mode='bilinear') or 'bisplines' interpolation.
   
   returns the dispersions of all orders and more...  
      C_zero, C_1, C_2, C_3, flux, xpix_, ypix_ , dd, ww, zmxwav, theta,
   after applying the wavelength scaling factor wpixscale to xpix_, ypix_
   and returns the original interpolated model grid 
   points: Xpix[order,wave], Ypix[order,wave] consistent with the anker at 2600, first order.
   kk is the spline fit order
   
   calls align_zemax() to make the correction to line up model (with boresight | overal flux shape)
   
   2013-01-13 changed the domain of the interpolation for Xfield,Yfield since small Xfield was quenched
              due to change anchor offset in field position boresight ("flux" model, 160)
   2013-01-17 updated interpolation; cleaner code; added 'simple' model,C_min1 replaced with flux in output.          
   '''   
   from .zemax import rdzemax, correctAnkPos
   
   if test == "flux":   
      if chatter > 0: print("zemax_flux used")       
      Xfield, Yfield = align_zemax_field(Xfield,Yfield,wheelpos,test="flux",chatter=chatter-1)
      if chatter > 1: print("new Xfield, Yfield : ", Xfield, Yfield)
      zmxwav, zmxorder, xfzmx, yfzmx, xcen, ycen, xpix, \
      ypix, zmxflux, zmxvalid = read_zemax_flux(wheelpos=wheelpos,file=None, 
      rotate=False,chatter=chatter-2, )
      xfzmx = xfzmx[1,3,:]
      yfzmx = yfzmx[1,3,:]
      xpix, ypix = align_zemax(xpix,ypix,wheelpos, test="flux",chatter=chatter)
      zmxflux = N.zeros((5,16,784),dtype=float)
   if test == "simple" or test == "simple_corrected":
      if chatter > 0: print("zemax simple model used: "+test)        
      zmxwav,zmxorder, xfzmx, yfzmx, xpix, ypix, zmxflux, zmxvalid = zemax_simply_corrected(wheelpos)
      wpixscale = 1.0      
      if test == "simple_corrected":
         xfzmx,yfzmx = align_zemax_field(xfzmx,yfzmx,wheelpos,test="simple_corrected",chatter=chatter-1)
         wpixscale = None      
   else:
      if chatter > 0: print("full zemax model used")
      (zmxwav, zmxorder, xfzmx, yfzmx, xcen, ycen, xpix, 
            ypix, zmxvalid) = rdzemax(wheelpos, chatter=chatter-2)
      zmxflux = N.zeros((5,16,784),dtype=float)
      xpix, ypix = align_zemax(xpix,ypix,wheelpos, chatter=chatter-1)
      xfzmx = xfzmx[1,3,:]
      yfzmx = yfzmx[1,3,:]
      
   # align to (the anchor position; overall flux) before proceeding
   #  moved up
   #xpix, ypix = align_zemax(xpix,ypix,wheelpos, flux= test == "flux",chatter=chatter-1)

   print('getZmx: zemax model read in completed')
   # get default anchor wavelength
   if wheelpos < 300:
      k = (N.where(zmxwav == 0.260))[0]
   else: 
      k = (N.where(zmxwav == 0.420))[0]  
   if chatter > 0: print("anchor at k=",k," wave: ",zmxwav[k])
   
   # anchor positions in first order
   xp = (xpix[1,k,:]).squeeze()
   yp = (ypix[1,k,:]).squeeze()
   
   # make final correction to the anchor positions (only in version 2 calibration uv clocked) done in getFirstOrder()
   #if (wheelpos == 160) & (test == "simple_corrected"): 
   #   xp_,yp_ = zemax.correctAnkPos(xp,yp,test=test,chatter=chatter)
   #   dx_ = xp_-xp  # shift
   #   dy_ = yp_-yp
   #   xp = xp+dx_   # corrected anchor
   #   yp = yp+dy_
   #   xpix = xpix+dx_  # corrected pixel positions
   #   ypix = ypix+dy_
    
   n1 = len(zmxorder)
   n2 = len(zmxwav)
   n3 = len(xp)
      
   xpix_ = N.zeros((n1,n2),dtype=N.float64) -999  # values corresponding to Xfield, Yfield 
   ypix_ = N.zeros((n1,n2),dtype=N.float64) -999  #    "        "             "       "
   flux  = N.zeros((n1,n2),dtype=float) 
   
   dmin = 0
   dminx = 0
   dminy = -100
   dmax = 2100
   chat = 2   
   if mode == 'bilinear':
      print('getZmx. using bilinear interpolation on 28,28 grid')
      for j1 in range(n1):         # loop over order
         for j2 in range(n2):      # loop over wavelengths
            xf = (xfzmx[:]).reshape(28,28) 
            yf = (yfzmx[:]).reshape(28,28)
            fx = ( xpix[j1,j2,:]).squeeze().reshape(28,28)   # the correction with align_zemax() is included
            fy = ( ypix[j1,j2,:]).squeeze().reshape(28,28)
            xpix_[j1,j2] = bilinear(Xfield, Yfield, xf[0,:].squeeze(), yf[:,0].squeeze(), fx,chatter=chat)
            ypix_[j1,j2] = bilinear(Xfield, Yfield, xf[0,:].squeeze(), yf[:,0].squeeze(), fy,chatter=chat)
            
   elif mode == 'spline':
      print('getzmx.  using bispline interpolation of order (3=cubic)', kk)
      for j1 in range(n1):         # loop over order
         for j2 in range(n2):      # loop over wavelengths
            # filter 
            q = (xpix[j1,j2,:] > (dminx)) & (xpix[j1,j2,:] < dmax) & (ypix[j1,j2,:] > (dminy)) & (ypix[j1,j2,:] < dmax) 
            if len(N.where(q)[0]) < 17:
               print("getZmx interpolation problem: filter q gives not enough results: ", N.where(q)[0])
               continue
               
            fx = xpix[j1,j2,q]
            fy = ypix[j1,j2,q]     
            xf = xfzmx[q]
            yf = yfzmx[q]
            ff = zmxflux[j1,j2,q]
            
            try:
               tck_x = interpolate.bisplrep(xf,yf,fx,xb=-0.20,xe=0.20,yb=-0.20,ye=0.20,kx=kk,ky=kk,)
               tck_y = interpolate.bisplrep(xf,yf,fy,xb=-0.20,xe=0.20,yb=-0.20,ye=0.20,kx=kk,ky=kk,) 
               tck_f = interpolate.bisplrep(xf,yf,ff,xb=-0.20,xe=0.20,yb=-0.20,ye=0.20,kx=kk,ky=kk,) 
               xpix_[j1,j2] = interpolate.bisplev(Xfield, Yfield, tck_x)
               ypix_[j1,j2] = interpolate.bisplev(Xfield, Yfield, tck_y)
               flux[j1,j2]  = interpolate.bisplev(Xfield, Yfield, tck_f)
            except:
               if chatter > 1: 
                  print("getZmx ERROR in loop (",j1,j2,") valid points len(q) = ",len(N.where(q)[0])) 
                  print(" fx: ",fx.shape," bad :",N.where(fx == -450)) 
                  print(" fy: ",fy.shape," bad :",N.where(fy == -450)) 
                  print(" yf: ",yf.shape ," xf: ",xf.shape) 
            
   else:
      print(' mode not valid') 
      raise RunTimeError
      return None
   
   # now I have the zemax points for the observed det coordinates, but unscaled so far.
   # as (xpix_,ypix_)
   table = N.array([xpix_,ypix_])
   if chatter > 1: 
      print('table interpolated from zemax: wave, DET(X,Y) positions first order')
      for iz in range(n2): print("%7.1f -- (%7.1f,%7.1f)" % (1e4*zmxwav[iz] , xpix_[1,iz], ypix_[1,iz])) 
   
   C_min1 = N.zeros(2)
   C_zero = N.zeros(3)
   C_1 = N.zeros(4)
   C_2 = N.zeros(3)
   C_3 = N.zeros(2)
   di_0,di_1,di_2,di_3,dim1 = 0,0,0,0,0
   wa_0,wa_1,wa_2,wa_3,wam1 = 0,0,0,0,0
   if test =='caltable': 
      theta = N.zeros(5)
      dd =    N.zeros(5)
      ww =    N.zeros(5)
      return C_zero, C_1, C_2, C_3, C_min1, xpix_, ypix_ , dd, ww, zmxwav, theta
   #   
   #     calculate dispersion coefficients  
   #
   #     first order scale pixel size
   #
   if wpixscale == None:  # first order 
       if (test == 'simple_corrected') & (wheelpos == 160):
           if chatter > 0: print("wpixscale correction with test2="+test2) 
           if k != 3: print('getZmx: newpixscale1 call, k not 3')
           xpix = xpix_.copy()
           ypix = ypix_.copy()
           xpix[1,:], ypix[1,:], _w, wpixscale1 = get_wpixscale(xpix_[1,:], ypix_[1,:], 
               w=1e4*zmxwav,wheelpos=wheelpos,mode=test2,chatter=chatter)
           xpix[1,:] = xpix[1,:]+ xpix_[1,k]    
           ypix[1,:] = ypix[1,:]+ ypix_[1,k]    
       else:
          if chatter > 2 : print("using 2009 wpixscale") 
          wpixscale1 = get_wpixscale(Xfield,Yfield,wheelpos=wheelpos,mode='field')
          if chatter > 1:
             print('getZmx: overriding wpixscale argument ')
             print('new wavelength pixel scale = ', wpixscale1)
          xpix = xpix_ * wpixscale1
          ypix = ypix_ * wpixscale1
   else:
      if chatter > 2: print("prescribed wpixscale ",wpixscale)
      wpixscale1 = wpixscale
      xpix = xpix_ * wpixscale1
      ypix = ypix_ * wpixscale1
   xpix1 = xpix.copy()
   ypix1 = ypix.copy()   

   if chatter > 1: 
      print('table pixscale corrected: wave, DET(X,Y) positions first order')
      for iz in range(n2): print("%7.1f -- (%7.1f,%7.1f)" % (1e4*zmxwav[iz] , xpix[1,iz], ypix[1,iz])) 
   #
   #   second order scale pixel size
   #
   if (wheelpos == 160) & (test == 'simple_corrected'): 
      wpixscale2=0.95   
      # finilize by requiring second coefficient of second order to be limited to 
      # the range 3.1e-5 < C_2[0] < 6.4e-4
   if wpixscale2 == None: # second order   
      wpixscale2 = get_wpixscale2(Xfield,Yfield,wheelpos=wheelpos,mode='field')    
      if chatter > 1:
         print('getZmx: overriding wpixscale2 argument ')
         print('new second order wavelength pixel scale [2000A]= ', wpixscale2)
            
        
   if chatter > 0 : 
      print('getZmx: order 1 using wpixscale = ',wpixscale1)
   print('===================================================================================')
   print('\n (getZmx calling getDisp)          order minus one dispersion ---')
   q = ( (xpix_[4,:] != -999) & (ypix_[4,:] != -999) & N.isfinite(xpix_[4,:]) & N.isfinite(ypix_[4,:])) 
   if q.sum() > 3:
      thm1, C_min1, dim1, wam1  = getDisp(xpix[4,q],ypix[4,q],zmxwav[q], fitorder=2,\
       lim = xlimits,chatter=chatter,wheelpos=wheelpos)  
   else:   thm1, C_min1 = (0, [0,0,0]) 
   
   print('===================================================================================')
   print('\n (getZmx calling getDisp)          third  order dispersion ---') 
   q = ( (xpix_[3,:] != -999) & (ypix_[3,:] != -999) & N.isfinite(xpix_[3,:]) & N.isfinite(ypix_[3,:]))
   if q.sum() > 2:
      th_3, C_3, di_3, wa_3    = getDisp(xpix[3,q],ypix[3,q],zmxwav[q], fitorder=1, \
      spectral_order=3, lim = xlimits,chatter=chatter,wheelpos=wheelpos)  
   else:  th_3, C_3    = (0, [0,0])  
   
   xpix2 = xpix_.copy() * wpixscale2
   ypix2 = ypix_.copy() * wpixscale2     

   print('===================================================================================')
   print('\n (getZmx calling getDisp)          second order dispersion ---')
   q = ( (xpix_[2,:] != -999) & (ypix_[2,:] != -999) & N.isfinite(xpix_[2,:]) & N.isfinite(ypix_[2,:])) 
   if q.sum() > 3:
      th_2, C_2, di_2, wa_2    = getDisp(xpix2[2,q],ypix2[2,q],zmxwav[q], fitorder=2, \
        spectral_order=2, lim = xlimits,chatter=chatter,wheelpos=wheelpos)  
   else: th_2, C_2    = (0, [0,0,0])   

   #xpix = xpix_.copy() * wpixscale1
   #ypix = ypix_.copy() * wpixscale1 
   xpix = xpix1
   ypix = ypix1  
        
   print('===================================================================================')
   print('\n                   first order dispersion ---')
   q = ( (xpix_[1,:] != -999) & (ypix_[1,:] != -999)& N.isfinite(xpix_[1,:]) & N.isfinite(ypix_[1,:]) ) 
   if q.sum() > 5:
      th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q],ypix[1,q],zmxwav[q], fitorder=4, \
      spectral_order=1, lim = xlimits,chatter=chatter,wheelpos=wheelpos)
   #   th_1, C_1, di_1, wa_1 = getDisp(xpix[1,q],ypix[1,q],zmxwav[q], order=4, lim = [-100,2100],chatter=chatter,wheelpos=wheelpos)
   else: th_1, C_1    = (0, [0,0,0,0,0])  
   
   print('===================================================================================')
   print('\n                    zeroth order dispersion >>> ')
   q = ( (xpix_[0,:] != -999) & (ypix_[0,:] != -999) & N.isfinite(xpix_[0,:]) & N.isfinite(ypix_[0,:]))
   if q.sum() > 4:
      Z  = getDisp(xpix[0,q],ypix[0,q],zmxwav[q], fitorder=-1, \
      spectral_order=0, lim = xlimits,chatter=chatter,wheelpos=wheelpos)
      if chatter > 5: print(Z)
      th_0, C_zero, di_0, wa_0  = Z
   else:  th_0, C_zero = (0, (N.array([0,0,0]),N.array([0,1])) )  
   print('=========================================done getDisp, done getZmx==========================================')
   
   theta = [th_0,th_1,th_2,th_3,thm1]
   dd = [di_0,di_1,di_2,di_3,dim1]
   ww = [wa_0,wa_1,wa_2,wa_3,wam1]
   # returns dispersion constants of zeroth, first, second, third orders, the unscaled zemax model X(wave),Y(wave),
   #         [position(pix;wave)] for all orders, and corresponding [wave] all orders for all valid points. 
   #         zemax wavelengths (in um), 180-slope of spectrum (in degrees) 
   if (wheelpos == 160) & (test == 'simple_corrected'): 
      # adjust the second coefficient of the second order by * 2.2
      # finalize by requiring second coefficient of the second order 
      # to be limited to the range 3.1e-5 < C_2[0] < 6.4e-4  
      C_2[0] = 2.2*C_2[0]
      if C_2[0] < 3.1e-5: C_2[0] = 3.1e-5
      if C_2[0] > 6.4e-4: C_2[0] = 6.4e-4
   return C_zero, C_1, C_2, C_3, flux, xpix_, ypix_ , dd, ww, zmxwav, theta,   
   
def getDisp(xpix,ypix, wave, thetain=None, wheelpos=None, fitorder=2, \
    spectral_order = 1, lim = [0,2100], chatter=0):
   '''  find the dispersion

   check there are at least 2 points with given wavelengths
     if not return with a default answer.  
   find the position of the anchor point from a fit
   find the angle of the points around the anchor point
   rotate the positions around the anchor point (2600 for UV, 4200 for V)
   fit a polynomial to the dispersion
   
   2010-10-28 New version NPMK (MSSL) 
   2011-09-20 reference wavelength 2nd order (UV) to 2000A
   '''   
   from numpy import array, where, polyfit, polyval, abs, pi, zeros
   from math import atan
   from uvotpy.uvotmisc import uvotrotvec
   from rationalfit import ratfit
   from scipy.interpolate import spline
   
   # wavelengths in Angstrom if in um
   if wave.mean() < 10: wave = wave*1e4
   
   xpix = xpix.flatten()
   ypix = ypix.flatten()
   if wheelpos == 200:
      dxpix = -27 
      dypix = -1 
      ang = 180-151.0
      ref_wave = 2000.
   elif wheelpos == 160:
      # new calibration spring 2013 
      dxpix = 0 # 53.36 +104
      dypix = 0 # 136.20 +78
      ang = 180-144.3
      ref_wave = 2000.
   elif wheelpos == 955: 
      dxpix = 0.
      dypix = 0.
      ang=180-140.0
      ref_wave = 4200.
   elif wheelpos == 1000: 
      dxpix = 0. 
      dypix = 0.0
      ang=180-148.0
      ref_wave = 4200.
      
   xpix -= dxpix
   ypix -= dypix
   
   if chatter > 2:
      print("getDisp input parameters: Xpix, Ypix")
      for i in range(len(xpix)): print("%3i  %7.1f  %5.1f   %5.1f  "%(i,wave[i],xpix[i],ypix[i]))
      print("thetain = ",thetain,"  fitorder = ",fitorder,"  spectral_order = ",spectral_order)
      print("lim = ",lim,"  wheelpos = ", wheelpos)
   
   if wheelpos < 300: 
     w_anch = 2600.
     if spectral_order == 1: 
        defaultdisp = (None, array([3.2,2600.]), None, None)  
     elif spectral_order == 2: 
        w_anch = 2000.
        defaultdisp = (None, array([2.0,2600.]), None, None)
     elif spectral_order == 3: 
        defaultdisp = (None, array([1.0,2600.]), None, None) 
     elif spectral_order == -1: 
        defaultdisp = (None, array([3.2,2600.]), None, None)
     elif spectral_order == 0 :
        defaultdisp = (None,(-array([6.,166.,0.035]),-array([3.2,91.1])), None, None )  # random numbers
     else:
        defaultdisp = (None, None, None, None)  

   if wheelpos > 300: 
     w_anch = 4200.
     if spectral_order == 1:
        defaultdisp = (None, array([5.8,4200.]), None, None) 
     elif spectral_order == 2: 
        defaultdisp = (None, array([2.9,4200.]), None, None)
     elif spectral_order == 3: 
        defaultdisp = (None, array([1.9,4200.]), None, None)
     elif spectral_order == -1: 
        defaultdisp = (None, array([5.8,4200.]), None, None)
     elif spectral_order == 0 :
        defaultdisp = (None,(-array([6.,166.,0.035]),-array([3.2,91.1])), None, None)   # random numbers
     else:
        defaultdisp = (None, None, None, None)
        
   theta, disp, D, wave_ = defaultdisp  # preset output 
   
   N = len(xpix)
   valid = zeros(N,dtype='bool')        
   for i in range(N):
      valid[i] = (xpix[i] > -100.) & (ypix[i] > -100) 
      #valid[i] = (xpix[i] > lim[0]) & & (ypix[i] > lim[0]) & (xpix[i] < lim[1]) & & (ypix[i] < lim[1])

     # iteratively determine which points are good or bad 
     # start with lowest wavelength when on top of detector
     # ??? what to do at other places on the detector???
   val = where(valid)[0]
   M = len(val)
   if M < 2:
      print("getDisp: less than 2 points available: set to default solution")
      return defaultdisp
   elif M == 2:
      # special case: linear extrapolation/interpolation used:
      if chatter > 2: print("getDisp    NOTE: linear fit")
      NP = 2
      xpix_ = xpix[val[0:NP]]
      ypix_ = ypix[val[0:NP]]
      wave_ = wave[val[0:NP]]
      x     = polyval(polyfit(wave_,xpix_,1),w_anch)
      y     = polyval(polyfit(wave_,ypix_,1),w_anch)
      yfit = polyfit(xpix_ - x,ypix_ - y,1)     
      theta = abs(180.0/pi*atan(yfit[-2]))  # slope at anchor
                
   else:
        # now we hunt for a sudden change in the slope of the zemax model spectrum
        # which indicates when Zemax starts getting non-linear (internal reflections)
      slope = polyfit(xpix[val[0:2]],ypix[val[0:2]],1)[0]
      NP = 2
      good = True
      while good:
         for i in range(2,M):
            slope9 = polyfit(xpix[val[i-1:i+1]],ypix[val[i-1:i+1]],1)[0] 
            if abs(slope9 - slope) < 0.04 * abs(slope): 
               slope = slope9
               NP = i
               if i == M-1: good = False 
            else: 
               good = False
               break
            if chatter > 2:   
               print("i=%3i wave=%8.1f slope=%10.4f  slope9=%10.4f" % (i,wave[i],slope,slope9))
               print(good)  
                
      NP += 1          
               
      # so we now know that points xpix[val[0:NP]] are good to go. 
      xpix_ = xpix[val[0:NP]]
      ypix_ = ypix[val[0:NP]]
      wave_ = wave[val[0:NP]]
      # determine the slope of the spectrum on the detector around w_anch.           
      x     = spline(wave_, xpix_, w_anch, order=1)
      #polyval(polyfit(wave_,xpix_,2),w_anch)
      y     = spline(wave_, ypix_, w_anch, order=1)
      #polyval(polyfit(wave_,ypix_,2),w_anch)
      yfit = polyfit(xpix_ - x,ypix_ - y,2)     
      theta = abs(180.0/math.pi*math.atan(yfit[-2]) ) # slope at anchor in degrees
      
   #if ((spectral_order == -1) ^ (spectral_order == 0)): # dispersion points the opposite way! 
   if ((spectral_order == -1) ): # dispersion points the opposite way! 
         theta = 180. - theta
         ang = 180. - ang
         
   if chatter > 2:
      print("getDisp valid points ")
      for i in range(len(xpix_)): print("%3i  %7.1f  %7.1f  %7.1f " % (i,wave_[i],xpix_[i],ypix_[i]))   
      print("x=%7.1f y=%7.1f  w_anch=%7.1f " % (x,y,w_anch))
      print("yfit  = ",yfit)
      print("theta = ",theta)
        
   if ((abs(theta-ang) > 10.)):
      print("theta = %8.2f  ang = %8.2f" % (theta,ang))
      print('getDisp.  message:something looks wrong with determination theta in uvotgrism.getDisp')
      return defaultdisp
      
   print(xpix_ - x)
   print(ypix_ - y)
   print(180.-theta)
   try: 
      D,ry = uvotrotvec( xpix_-x, ypix_-y, 180.-theta )  
   except:
      print("getDisp: fatal problem rotating the data ")  
      raise  
   if chatter > 1:
      print('getDisp. rotated points    #    wave      Dx      ry  ')
      for i1 in range(len(D)): print('                   ',i1,'  ',wave_[i1],'  ',D[i1],'  ',ry[i1])
        
     # now ry should be small : if not, there is a problem with the previous rejection of 
     # bad Zemax model points  
   q = where(abs(ry) > 30)
   if q[0].size > 0: 
      print("\n\n\n getDisp: PROBLEM in codes rejection of bad Zemax model data\n\n\n")  
      return defaultdisp
      
   if ((spectral_order == 1) ^ (spectral_order == 3) ^ (spectral_order == -1)): 
      if fitorder >= len(D) : 
         print("adjusting order of fitting polynomial ")
         fitorder = len(D)-1     
      disp = polyfit(D, wave_, fitorder)      
   # if fitorder == 0 then give default
   # if D[-1] more than 0.4% different from defaultdisp[1][1] then give default 
     
   elif spectral_order == 2:
      if fitorder >= len(D) : 
         print("adjusting order of fitting polynomial ")
         fitorder = len(D)-1     
      disp = polyfit(D, wave_, fitorder)
      if wheelpos < 300:
         # convert dispersion to that at 2600.0 anchor
         dist, disp = change2ndOrderReferencePoint(0.0,disp,oldReferenceWave=2000.,
             newReferenceWave=2600.0)
      
   elif spectral_order == 0: 
      # returns the polynomial coefficients p,q for the nominator and denominator
      # of y = ( poly(p)/poly(q) )
      #return defaultdisp
     try:
      dis1 = ratfit(D,wave_,2,1,) 
      numerator = dis1[2:2+dis1[0]]
      denumerator = dis1[2+dis1[0]:]
      disp = (numerator, denumerator)
     except:
      print('failure to find solution for zeroth order dispersion')
      return defaultdisp 
      pass
        # note:
        # alternatively, use leastsq from scipy.optimize: 
        # with p the list of adjustable variables, i.e.,  p=(a0,y0,sig)
        # params_fit, ier = leastsq(RatFun, (a0,y0,sig,a1,y1,sig) , args=(x,y))
        # where the arguments are the functionname(p,x,y,) to be zero for a perpect fit, 
        # the parameter list and the arguments of the function, i.e., RatFun = y-f(p,x,).
   else:        
      # this should not happen
      print("getDisp: input error in parameter spectral_order = ", spectral_order)
      return defaultdisp
        
   return theta, disp, D, wave_ 
            

def findAnker(filestub, ext, RA, DEC, wheelpos=200, \
     lfilter='uvw1', extf=None, lfilt2=None, lfilt2_ext=None, \
     method=None, mode=None, attfile=None, chatter=2,):
   '''
   provided both a uvw1 and a ugu observation were made 
   in the same observation, use these to determine the 
   anker of 2600(1) in the UV grism image. Assumed is that
   the grism and lenticular filter image have the same 
   extension. Before running this, run uvotgrapcorr on
   the grism image to get a better aspect solution.
   
   Input: 
      filestub = full path + part of filename including obsid.
      ext      = number of extension to be used for grism 
      extf     = number of extension lenticular filter 
      lfilter  = lenticular filter id
      lfilt2   = if two lenticualr filter observations, this is id of the second
      lfilt2_ext  = extension number of second lent. filter observation
      RA       = Right ascension in decimal degrees
      DEC      = Declination in decimal degrees
      method   = None : use the uvw1 image (default)
               = fake : use the UV grism aspect only to find a solution 
                 requires attitude file and makes a fake uvw1 sky file
      attfile  = full path+filename of attitude file
      
   Output: (list of 3 numpy arrays)
      anker    = anker position of 260nm in first order computed 
                 using fixed pixel scale in image coordinates
      anker_as = offset (DX,DY) in arcsec in DET coordinate system of the source from the 
                 needs to be converted to input rays by applying transform. 
      anker_field = offset(theta,phi) in degrees from the axis for 
                 the input field coordinates for the zemax model lookup                  
   '''
   __version__ = '1.1 NPMK 20090914 (MSSL)'
   # npkuin@gmail.com
   if extf == None: extf = ext
   
   anker = (-1,-1)
   if wheelpos < 300:
     gfile = filestub+'ugu_dt.img'
   else:
     gfile = filestub+'ugv_dt.img'
       
   if lfilter == 'wh': ffile = filestub+'uwh_sk.img'
   if lfilter == 'u' : ffile = filestub+'uuu_sk.img'
   if lfilter == 'v' : ffile = filestub+'uvv_sk.img'
   if lfilter == 'b' : ffile = filestub+'ubb_sk.img'
   if lfilter == 'wh' : ffile = filestub+'uwh_sk.img'
   if lfilter == 'uvw1' : ffile = filestub+'uw1_sk.img'
   if lfilter == 'uvw2' : ffile = filestub+'uw2_sk.img'
   if lfilter == 'uvm2' : ffile = filestub+'um2_sk.img'
   
   if lfilter == None: ffile = filestub+'uw1_sk.img'  # make fake lenticular file
   
   if chatter > 2:
      print('grism file = ',gfile)
      print('lenticular filter = ',lfilter,'  filter file = ',ffile)
   
   if ((method == 'fake') ^ (lfilter == None)): 
      from uvotpy.uvotwcs import makewcshdr 
      ffile = makewcshdr(RA,DEC,filestub,ext,attfile,wheelpos=wheelpos,chatter=chatter)
      lfilter = 'uvw1'
      
   hf = pyfits.getheader(ffile,extf)   
   hg = pyfits.getheader(gfile,ext)
   
   if chatter > 1: 
      print('findAnker. wheelpos               = ',wheelpos)
      print('findAnker. grismfile              = ',gfile,'[',ext,']')
      print('findAnker. lenticular filter file = ',ffile,'[',extf,']')
      print('findAnker. lfilter                = ',lfilter)
      
   if lfilt2 != None:
      if lfilt2_ext == None: 
         lfilt2_ext = extf   # if no second extension given, then use the same as the first one
      if lfilt2 == 'wh': ffil2 = filestub+'uwh_sk.img'
      if lfilt2 == 'u' : ffil2 = filestub+'uuu_sk.img'
      if lfilt2 == 'v' : ffil2 = filestub+'uvv_sk.img'
      if lfilt2 == 'b' : ffil2 = filestub+'ubb_sk.img'
      if lfilt2 == 'wh' : ffil2 = filestub+'uwh_sk.img'
      if lfilt2 == 'uvw1' : ffil2 = filestub+'uw1_sk.img'
      if lfilt2 == 'uvw2' : ffil2 = filestub+'uw2_sk.img'
      if lfilt2 == 'uvm2' : ffil2 = filestub+'um2_sk.img'
      if chatter > 2: 
        print('findAnker. lenticular filter # 2  = ',ffil2,'[',lfilt2_ext,']')
        print('findAnker. lfilter2               = ',lfilt2)     
        
   # check for losses in grism image 
   if float(hg['BLOCLOSS']) != 0: print('#### BLOCLOSS = '+repr(hg['BLOCLOSS']))
   if float(hg['STALLOSS']) != 0: print('#### STALLOSS = '+repr(hg['STALLOSS']))
   if float(hg['TOSSLOSS']) != 0: print('#### TOSSLOSS = '+repr(hg['TOSSLOSS']))
   tstart = hg['TSTART']

   print('grism exposure time = ',hg['EXPOSURE'],'  seconds')
   
   RA_PNT  = hg['RA_PNT']
   DEC_PNT = hg['DEC_PNT']
   PA_PNT  = hg['PA_PNT']   # roll angle
   time    = hg['TSTART']   # time observation
   ra_diff  = RA - RA_PNT
   dec_diff = DEC - DEC_PNT
   RAs = repr(RA)
   DECs= repr(DEC)
   exts = repr(ext)
   extfs = repr(extf)
   extfs2 = repr(lfilt2_ext)
   
   from os import getenv,system
   system('/bin/rm radec.txt skyfits.out skyfits.in detmm.txt ')
   system('echo '+RAs+'  '+DECs+' > radec.txt ' )

   CALDB = getenv('CALDB')
   if CALDB == '': 
      print('findAnker. the CALDB environment variable has not been set')
      return None
   HEADAS = getenv('HEADAS')
   if HEADAS == '': 
      print('findAnker. The HEADAS environment variable has not been set')
      print('           which is needed for the uvot Ftools \n')
      return None 
      
   command = HEADAS+'/bin/uvotapplywcs infile=radec.txt outfile=skyfits.out wcsfile=\"'+ffile+'['+extfs+']\" operation=WORLD_TO_PIX'
   if chatter > 1: print(command)
   system( command )
   
   f = open('skyfits.out', "r")
   line = f.read()
   x1, y1 = (line.split())[2:4]
   f.close  
   system( 'echo '+repr(x1)+'  '+repr(y1)+'  > skyfits.in' )
   
   command = HEADAS+'/bin/uvotapplywcs infile=skyfits.in outfile=detmm.txt wcsfile=\"'+ffile+'['+extfs+']\" operation=PIX_TO_WORLD to=D'
   if chatter > 1: print(command)
   system( command )
   
   f = open('detmm.txt', "r")
   line = f.read()
   if chatter > 2: print('detmm: '+line)
   x1, y1 = (line.split())[2:4]
   f.close
   if chatter > 2: print('detector coordinates [mm]: ',x1, y1)
   
   if lfilt2 != None:      
      command = HEADAS+'/bin/uvotapplywcs infile=radec.txt outfile=skyfits.out wcsfile=\"'+ffil2+'['+extfs2+']\" operation=WORLD_TO_PIX'
      if chatter > 1: print(command)
      system( command )
      f = open('skyfits.out', "r")
      line = f.read()
      x2, y2 = (line.split())[2:4]
      f.close  
      system( 'echo '+repr(x2)+'  '+repr(y2)+'  > skyfits.in' )
      command = HEADAS+'/bin/uvotapplywcs infile=skyfits.in outfile=detmm.txt wcsfile=\"'+ffil2+'['+extfs2+']\" operation=PIX_TO_WORLD to=D'
      if chatter > 1: print(command)
      system( command )
      f = open('detmm.txt', "r")
      line = f.read()
      if chatter > 2: print('detmm: '+line)
      x2, y2 = (line.split())[2:4]
      f.close
      if chatter > 2: print('detector coordinates [mm]: ',x2, y2)
            
   # anchor in pixel-DET pixels in lenticular filter and in arcsec from  in lenticular filter

   anker_uvw1det = N.array([float(x1), float(y1)])/0.009075+[1100.5,1100.5]
   if chatter > 0: print('derived coord (pix) with ref.[1100.5,1100.5] in filter = ',lfilter,' = ',anker_uvw1det)
   #anker_uvw1det = anker_uvw1det - [954.60,1044.66]
   borexy = boresight(filter=lfilter,date=000000)
   anker_uvw1det = anker_uvw1det - borexy
   anker_as = anker_uvw1det*0.502 
   if lfilt2 != None:
      anker_uvw1det2 = N.array([float(x2), float(y2)])/0.009075+[1100.5,1100.5]
      if chatter > 0: print('derived coord (pix) with ref.[1100.5,1100.5] in filter2 = ',lfilt2,' = ',anker_uvw1det2)
      #
      borexy2 = boresight(filter=lfilt2,date=00000)
      anker_uvw1det2 = anker_uvw1det2 - borexy
      anker_as2 = anker_uvw1det2*0.502
      anker_as  = 0.5*(anker_as + anker_as2) 
   if lfilt2 == None:   
      if chatter > 0: print('adopted boresight (DET-pix) ref.[1100.5,1100.5] in filter = ',lfilter,' = ',borexy)
      if chatter > 0: print('derived boresight offset lenticular filter ',lfilter,' (DET pix) = ',anker_uvw1det)
      if chatter > 0: print('derived   boresight offset (from lenticular filter)  in (\")  = ', anker_as)
   else:
      if chatter > 0: print('adopted boresight (DET-pix) ref.[1100.5,1100.5] in filter = ',lfilter,' = ',borexy)
      if chatter > 0: print('derived boresight offset lenticular filter ',lfilter,' (DET pix) = ',anker_uvw1det)
      if chatter > 0: print('derived coord (pix) with ref.[1100.5,1100.5] in filter2 = ',lfilt2,' = ',anker_uvw1det2)
      if chatter > 0: print('adopted boresight (DET-pix) ref.[1100.5,1100.5] in filter = ',lfilter,' = ',borexy2)
      if chatter > 0: print('derived boresight offset lenticular filter ',lfilter,' (DET pix) = ',anker_uvw1det2)
      if chatter > 0: print('derived mean boresight offset (from lenticular filter)  in (\")  = ', anker_as)
           
   #
   # 2600(1) anker offset in pixels in grism det image;
   # (*_dt.img) for fixed scale factor (no grism distortion).
   # 
   anker = anker_as/0.5496 
   if chatter > 0: print('derived boresight offset UV grism (0.5496"/pix) in fixed scale DET pixels = ', anker)
   #
   #  define in pixels on _dt image (physical array from extension)
   #
   #anker = anker + [928.53, 1002.69]  -[27,1]
   ankerref = N.array([0,0])
   if   wheelpos == 200: borep =  [901.53, 1001.69]  # IMAGE COORDINATES
   #elif wheelpos == 160: borep =  N.array(boresight(filter='uc160')) - [104,78]   # IMAGE COORDINATES  with uncertainty [8,4]
   elif wheelpos == 160: borep =  N.array(boresight(filter='uc160',date=0000)) - [104,78]   #coordinate
   elif wheelpos == 955: borep =  N.array(boresight(filter='vc955',date=0000)) - [0,0]
   elif wheelpos == 1000: borep =  N.array(boresight(filter='vg1000',date=0000)) - [0,0]
   else: 
      print("no anchor for this wheelpos\n")
      return   
   print('*** before shift and rotation: angles - ', anker_as[0]/3600.,anker_as[1]/3600.)
   anker = anker + borep   
   if wheelpos == 200:
      ankeroff  =  anker_as/0.5496 + [1092.5,1104.5] # [det pixels] : use anchor offset 
   elif wheelpos == 160:
      ankeroff =    anker_as/0.5496 + [1100.5,1100.5]
   elif wheelpos == 955:
      ankeroff =    anker_as/0.5496 + [1100.5,1100.5]
   elif wheelpos == 1000:
      ankeroff =    anker_as/0.5496 + [1100.5,1100.5]
  #       to make the boresight map to phi=(0,0) in the zemax model.  
   if chatter > 0: 
      print('[IMAGE] bore point position of grism used                = ', borep)
      print('derived linear anker position in grism image coordinates = ', anker)   
      #
      # the idea behind this is that the linear position maps 1-1 to the angular offset  
      # use the approximate mapping from the model to get the linear angles and then use those 
      # field coordinates to calculate the model anchor position (in uvotgrism.getZmx).
      # that in turn then has to be offset by the shift of the boresight from the model 
      # boresight at [1092.5,1104.5]
      #
      # print '                                  and original zemax coordinates = ', ankeroff     
      print('RA, DEC used: ' , RA, DEC)
   # now find the field coordinate by applying a shift followed by a rotation of the input angles
   # the result are the field angles as they should be in the absence of distortion. The shift 
   # comes about from the difference in anchor point of the zemax model and that observed.
   # same as zemax.matchDetToField 
   if wheelpos == 200:
      ankx = (anker[0] + 8.0 - 1023.5)/6550.4 
      anky = (anker[1] - 4.0 - 1023.5)/6550.4   
   elif wheelpos == 160:
      ankx =  anker_as[0]/3600. #(anker[0] - 1023.5)/6550.4 
      anky =  anker_as[1]/3600. #(anker[1] - 1023.5)/6550.4   
   elif wheelpos > 300:
      ankx =  anker_as[0]/3600. #(anker[0] - 1023.5)/6550.4 
      anky =  anker_as[1]/3600. #(anker[1] - 1023.5)/6550.4  
   if mode == None:     
      xsk, ysk =   uvotmisc.uvotrotvec(ankx, anky, 64.6)   
      anker_field = N.array([xsk,ysk])  
      rotang = 64.6 
   else:
      anker_field = N.array([ankx,anky])  
      rotang = 0.0 ; xsk, ysk = ankx, anky
   if chatter > 0: 
      print('*** after applying shift and ',rotang,' deg rotation, we find the ')
      print('    corresponding zemax field coordinates : ',xsk, ysk)    
   return anker, anker_as, anker_field
           
           
def _find_wpixscale_xy(wheelpos, sourcelist='sources.160.1.txt'):
   ''' this routine to find the wavelength pix scale across the 
      detector from the first - zeroth order scale ''' 
   f = open(sourcelist)
   fa = f.readlines()
   f.close()
   # run over all observations
   n = len(fa)
   wpixscales = N.zeros(n)
   ankers = N.zeros(2*n).reshape(2,n)
   for k in range(n):
      ras,decs,filestub,exts,filter1,ext1s,filter2,ext2s,xgrasps,ygrasps = fa[k].split(',')
      ra, dec, ext, ext1, ext2 = float(ras),float(decs),float(exts),float(ext1s),float(ext2s)
      xgrasp, ygrasp = float(xgrasps) , float(ygrasps)
      if ext1 == 0: ext1 = ext
      if ext2 == 0: ext2 = None
      if filter2 == 'None': filter2 = None
      if filter1 == 'None': filter1 = 'uvw1'
      anker, anker_as, anker_field = findAnker(filestub, ext, ra, dec, wheelpos=wheelpos, 
        lfilter=filter1, extf=ext1, lfilt2=filter2, lfilt2_ext=ext2,
        method=None, attfile=None, chatter=2, HACK=0)   
      xfield = anker_field[0]
      yfield = anker_field[1]
      (C_zero, C_1, C_2, C_3, C_min1, xpix, ypix, zmxdis, zmxwav, wave, theta) = getZmx( \
          xfield, yfield, wheelpos, mode=None,chatter=chatter, \
          test=None, kk=3,s=50,wpixscale=1.0)
      i260 = N.where( wave == 2600 )
      i500 = N.where( wave == 5000 )  
      dist260_500pix= N.array([xpix[0,i500]-xpix[0,i260],ypix[0,i500]-ypix[0,i260]])     # coordinate distance
      # nu is the observed scale (not corrected):
      # x[260,1] - x[260,0] = xank - xgrasp + xpix[0,i500] - xpix[0,i260] = xank - xgrasp + dist260_500pix[0] < 0
      # x[260,1] - x[260,0] = yank - ygrasp + ypix[0,i500] - ypix[0,i260] = yank - ygrasp + dist260_500pix[1] > 0 
      # en de zemax unscaled distance is:
      # xpix[1,i260] - xpix[0,i260]
      # ypix[1,i260] - xpix[0,i260]
      # but I really want to make a projection on the overall dispersion line theta[1] and then calc the ratio.
      # for now, just use:
      dxobs =  -(anker[0] - xgrasp + dist260_500pix[0] )   
      dyobs =    anker[1] - ygrasp + dist260_500pix[1]
      dxmod =  (xpix[1,i260] - xpix[0,i260])
      dymod = -(ypix[1,i260] - xpix[0,i260])
      disobs = N.sqrt(dxobs**2+dyobs**2)
      dismod = N.sqrt(dxmod**2+dymod**2)
      wscale = disobs/dismod
      wpixscales[k] = wscale
      ankers[0,k] = anker[0]
      ankers[1,k] = anker[1]
   return wpixscales, ankers   
   
def get_wpixscale(x,y,w=None, wheelpos=160, mode='other',chatter = 0):
      ''' return the scale factor for the dispersion 
      correction to zemax over the face of the detector
      as function of position x,y in DET-pix coordinates.
      
      Parameters
      ----------
      x,y : floats, array
         position on detector in det coordinates
      w : array (optional)
         wavelengths in A
         if None, only a scale factor is returned (version 1 calibration)
         only active in uv grism, clocked mode (160)     
      wheelpos : int
         filter wheel position 
      mode : str
         options implemented by mode
         for wheelpos 160, if w given, 4 options, 2013-05-14
         -  mode=fit1: constant map from zemax to observed over detector
         -  mode=fit2: constant quadratic factor in map from zemax to observed
         -  mode=linear: only the linear term in the map from zemax to observed
         -  mode=_anything_else_: bispline fits to linear and quadratic term 
      chatter : int 
         verbosity
                         
      '''
      if chatter > 2: print('uvotcal.wpixscale wheelpos=',wheelpos)
      if wheelpos == 200: 
         return 0.960
      elif wheelpos == 160:
         if w == None:   
             if mode == 'DET-pix':
                tck = [N.array([0.,0.,2048.,2048.]),N.array([0.,0.,2048.,2048.]), 
                N.array([0.92086849,0.93637569,0.94497985,0.97726383]),1,1]
                xx = x - 104
                yy = y - 78  
                return interpolate.bisplev(xx, yy, tck)  
             elif mode == 'DET-pix_interorder0_1':
                # see fitpack.py: knots and coefficients
                tck = [N.array([    0.,     0.,  2048.,  2048.]),
                   N.array([    0.,     0.,  2048.,  2048.]),
                   N.array([ 0.69644746,  0.86675867,  1.20252484,  0.9907974 ]),1, 1]
                # interpolation is valid in the square spanned by corners [308,145] to [1855,1734]
                # and extrapolated beyond    
                xx = x - 104
                yy = y - 78    
                return interpolate.bisplev(xx, yy, tck)
             elif mode == 'field':
                #  major problem with interpolation for xf > 0.1 
                tck = [N.array([-0.175,-0.175,0.175,0.175]),
                   N.array([-0.175,-0.175,0.175,0.175]), 
                   N.array([ 0.95715799,  0.89390723,  0.95062182,  0.97579035]),1,1]
                return interpolate.bisplev(x,y,tck) 
             elif mode == 'rotfield':
                # here the input field coordinates have been rotated by 64.6 degrees
                tck =  [N.array([-0.18, -0.18,  0.18,  0.18]), 
                N.array([-0.18, -0.18,  0.18,  0.18]), 
                N.array([ 0.9190176 ,  0.93401162,  0.94126273,  0.98094215]), 1,1]
                return interpolate.bisplev(x,y,tck) 
             elif mode == 'field_interorder0_1':
                # the input are the zemax system field coordinates (not rotated to align detector axis) 
                tck = [N.array([-0.172, -0.172,  0.172,  0.172]), 
                   N.array([-0.172, -0.172,  0.172,  0.172]), 
                   N.array([ 1.00215967,  0.75247786,  1.21600755,  0.80937101]),1,1]
                return interpolate.bisplev(x,y,tck)   
             else:
                return 0.940 
         else:
             # new calibration with polynomial depending on wavelength w for points x,y (float arrays)
             tck_linear= [N.array([    0.,     0.,  2200.,  2200.]),
                N.array([  -50.,   -50.,  2100.,  2100.]),
                N.array([ 1.04734146,  0.96739897,  0.98335362,  1.02276879]),1,1]
             tck_quadratic = [N.array([    0.,     0.,     0.,  2200.,  2200.,  2200.]),
                N.array([  -50.,   -50.,   -50.,  2100.,  2100.,  2100.]),
                N.array([ -1.69348636,   5.49555444, -19.58392502,  -4.13090029,
                -0.12815537,  16.25640561,   8.43919587,  -0.90921651, -10.81594911]),2,2]
             tck_linonly = [N.array([    0.,     0.,  2200.,  2200.]),
                N.array([  -50.,   -50.,  2100.,  2100.]),
                N.array([ 1.01077821,  0.99491015,  1.02217832,  1.02999323]),1,1]              
             n2 = len(w)
             # apply pixscale to x,y 
             if len(x) == len(y) == n2:
                k = (N.where(w == 2600.))[0]
                if len(k) == 0:
                    raise RuntimeError('get_wpixscale1: anchor index not found')
                else:
                      # anchor image coordinates
                      ankx = x[k] - 140
                      anky = y[k] - 78
                disfactor = N.zeros(3)
                disfactor[1] = interpolate.bisplev(ankx,anky,tck_linear)
                disfactor[0] = interpolate.bisplev(ankx,anky,tck_quadratic)*1.0e-5 
                if mode == 'linear' : 
                   disfactor[1] = interpolate.bisplev(ankx,anky,tck_linonly)
                   disfactor = disfactor[1:]
                   if chatter > 0: print("get_wpixscale: mode='linear'")
                if mode == 'fix2' : 
                   if chatter > 0: print("get_wpixscale: mode='fix2'")
                   disfactor[0]=4.5e-6
                if mode == 'fix1' : 
                   if chatter > 0: print("get_wpixscale: mode='fix1'")
                   disfactor[0:2] = [4.5e-6,1.00]
                wnew = N.polyval(disfactor,w-2600.)+2600.
                if chatter > 0: print(mode, disfactor, wnew)             
                # given x(wnew),y(wnew) compute corrected xnew(w), ynew(w)
                xnew=N.polyval(N.polyfit(wnew-2600.,x-x[k],9),w-2600.) 
                ynew=N.polyval(N.polyfit(wnew-2600.,y-y[k],9),w-2600.) 
                return xnew, ynew, wnew, disfactor   
             else:
                raise RuntimeError('get_wpixscale: array length mismatch') 
                                
      elif wheelpos == 1000: 
         if ((mode == 'other') ^ (mode == 'DET-pix')):
            tck = [N.array([0.,0.,2048.,2048.]),N.array([0.,0.,2048.,2048.]),\
                N.array([0.97337398,0.96848245,0.95029983,0.98082003]),1,1] 
            return interpolate.bisplev(x,y,tck) 
         elif mode == 'rotfield':  # note 20100119: not sure if this is the rotated field coordinate or just field coordinate!!
            tck = [N.array([-0.18,-0.18,0.18,0.18]),N.array([-0.18,-0.18,0.18,0.18]),\
                N.array([0.97613126,0.96686695,0.94780679,0.98290208]),1,1]
            return interpolate.bisplev(x,y,tck) 
         else:   
            return 0.967
      elif wheelpos == 955: 
         if ((mode == 'other') ^ (mode == 'DET-pix')):
            tck = [N.array([  -50.,   -50.,  2200.,  2200.]),N.array([  -50.,   -50.,  2200.,  2200.]), \
                   N.array([ 0.97337996,  0.88165988,  0.98563709,  0.95443548]), 1, 1]
            xx = x - 104
            yy = y - 78  
            return interpolate.bisplev(xx,yy,tck)
         if mode == 'rotfield':  
            # scale from mode DET-pix:
            #anker = array(boresight('vc955',r2d=0.))
            #tck = [(N.array([  -50.,   -50.,  2200.,  2200.])-anker[0])/6550.4, \
            #       (N.array([  -50.,   -50.,  2200.,  2200.])-anker[1])/6550.4, \
            tck = [N.array([-0.17002015, -0.17002015,  0.17347032,  0.17347032]),\
                   N.array([-0.15305936, -0.15305936,  0.19043112,  0.19043112]),\
                   N.array([ 0.97337996,  0.88165988,  0.98563709,  0.95443548]), 1, 1]
         #   # here the field coordinates of the measured ratio were rotated by 70 deg. 
         #   tck = [N.array([-0.17, -0.17,  0.18,  0.18]), \
         #          N.array([-0.17 , -0.17 ,  0.175,  0.175]), \
         #          N.array([ 0.98633387,  0.9603005 ,  0.97220619,  0.89015667]), 1, 1]
            return interpolate.bisplev(x,y,tck)  
         else:
            return 0.953  # mean value                        
      else:
         return 1.00

def zmx2ank(wave,xpix,ypix,wheelpos):
      ''' 
      Calculate the anker position given the zemax model
      
      Snippet of code from getFirstOrder . . . 
      '''
      if wave.mean() < 10: wave = wave*1e4  # units nm -> A
      q260 = N.where(wave == 2600)
      #  The zemax model has not yet been corrected for the actual boresight offset 
      if wheelpos == 200:
        dxpix = -27 
        dypix = -1 
      elif wheelpos == 160: 
         dxpix =  53.36 
         dypix = 136.20
      elif wheelpos == 955:
         dxpix =  0.
         dypix =  0.
      elif wheelpos == 1000:     
         dxpix =  0. 
         dypix =  0.
      xpix = xpix - dxpix
      ypix = ypix - dypix      
      xdet = xpix[1,q260].squeeze()   
      ydet = ypix[1,q260].squeeze() 
      anker = N.array([xdet,ydet])
      uncorr_anker = anker.copy().squeeze()
      #
      #  refine the anchor position with the empirical correction
      #
      if wheelpos == 200:
         anker = correctAnkPos(anker[0],anker[1])
         anker = N.array(anker).squeeze()
      elif wheelpos == 160: 
         anker = correctAnkPos(anker[0],anker[1],wheelpos=160)
         anker = N.array(anker).squeeze()
      else:
         anker = correctAnkPos(anker[0],anker[1],wheelpos=wheelpos)
         anker = N.array(anker).squeeze()
      print('xmx2ank: Anker position correction done with correctAnkPos to zemax model') 
      # positioncorrection 
      xpix = xpix + anker[0] - uncorr_anker[0]
      ypix = ypix + anker[0] - uncorr_anker[0]
      #
      #       
      if chatter > 0:
             print('=======================================================================')
             print('getFirstOrder: check wave at 2600 = ',N.where(wave == 2600)[0])
             print('the wavelength supposed to be 2600 A is      : ', wave[q260]) 
             print('The anker at 2600 first order (corrected) is : ',anker,' [DET-pix]')
             print('The anker at 2600 first order (corrected) is : ',anker - [77+27,78],' [DET-image]')
             print('The uncorrected anker position is [DET-pix]  : ',uncorr_anker[0],uncorr_anker[1])
             print('The uncorrected anker position is [DET-image]: ',uncorr_anker[0]-104,uncorr_anker[1]-78)
             print('=======================================================================')    
      return anker, uncorr_anker, xpix, ypix         

def extrapolateArray(x,y,f):
   '''
   extrapolate f values to points that are NaN for x, y and f, for wrCalFile 
   x, y, f in form shape(28,28)   
   keep the good values
   ''' 
   import numpy
   from scipy import interpolate
   x1, x2, y1, y2 = x.min()-0.002, x.max()+0.002, y.min()-0.002, y.max()+0.002
   xx = x.copy().flatten()
   yy = y.copy().flatten()
   ff = f.copy().flatten()
   q = numpy.isfinite(ff)
   tck = interpolate.bisplrep(xx[q].squeeze(),yy[q].squeeze(),ff[q].squeeze(), xb=x1,xe=x2,yb=y1,ye=y2)  
   qn = numpy.isnan(ff) | numpy.isinf(ff)
   if len(qn) == 0: return ff.flatten()   
   xa = xx[qn].squeeze()
   ya = yy[qn].squeeze()
   fa = xa*0.
   for i in range( len(xa)): 
      xa1 = xa[i]
      ya1 = ya[i]
      fa1 = interpolate.bisplev(xa1, ya1, tck)
      fa[i] = fa1
   ff[qn] = fa      
   return ff.flatten()

def curved_computation(ra,dec,filestub,ext, anchor):
   ''' commands to prepare data for determining curvature 
   anchor=array([ 917.1, 455.3]) observed position in image
   '''
   import uvotgrism
   # cd /Volumes/data/grism/WR52/00057916001/uvot/image
   # get the dispersion etc.
   Y = uvotgrism.getSpec(ra,dec,filestub,ext)
   ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Y
   # now use the observed anchor position rather than the one from the cal file
   IMG = uvotgrism.extractSpecImg(filestub+'ugu_dt.img',1,anchor,angle)
   dis, spnet, bg, bg1, bg2, extimg, spimg, spnetimg, offset, ank_c = IMG
   # get the full background image
   BKG = uvotgrism.findBackground(extimg,-1.0123456789 )
   # the defined positions of the orders may have to be edited in the program file iteratively.
   CURVE = uvotgrism.curved_extraction(extimg,ank_c,anker,160,trackfull=True)
   fitorder, gfit, extras = CURVE
   r1,r2,r3,r4 = curved_sub1(gfit,ank_c)
   # where:
   #(q0,i0,X0,Y0,sig0,a0),(q1,i1,X1,Y1,sig1,a1),(q2,i2,X2,Y2,sig2,a2),(q3,i3,X3,Y3,sig3,a3) = r1,r2,r3,r4
   return r1,r2,r3,r4,(angle,extimg,C_1,C_2),(Y,IMG,BKG,CURVE)   

def curved_sub1(gfit,ank_c):
   from numpy import where
   q0 = where(gfit[0,1,:] == 0 & (gfit[0,3,:] < 160) & (gfit[0,3,:] > 40) )   
   q1 = where((gfit[1,1,:] == 1) & (gfit[1,4,:] < 10.0) & (gfit[1,4,:] > 2.7) & (gfit[1,3,:] > 72.) & (gfit[1,3,:] < 200))
   q2 = where((gfit[2,1,:] == 2) & (gfit[2,4,:] < 6.0) & (gfit[2,4,:] > 2.7) & (gfit[2,3,:] < 200))
   q3 = where((gfit[3,1,:] == 3) & (gfit[3,4,:] < 10) & (gfit[3,4,:] > 3) & (gfit[3,3,:] > 50) & (gfit[3,3,:] < 190) )
   #
   i0 = gfit[0,0,q0][0]
   i1 = gfit[1,0,q1][0]
   i2 = gfit[2,0,q2][0]
   i3 = gfit[3,0,q3][0]
   #
   X0 = i0-ank_c[1]
   X1 = i1-ank_c[1]
   X2 = i2-ank_c[1]
   X3 = i3-ank_c[1]
   Y0 = gfit[0,3,q0][0]
   Y1 = gfit[1,3,q1][0]
   Y2 = gfit[2,3,q2][0]
   Y3 = gfit[3,3,q3][0]
   sig0 = gfit[0,4,q0][0]
   sig1 = gfit[1,4,q1][0]
   sig2 = gfit[2,4,q2][0]
   sig3 = gfit[3,4,q3][0]
   a0 = gfit[0,2,q0][0]
   a1 = gfit[1,2,q1][0]
   a2 = gfit[2,2,q2][0]
   a3 = gfit[3,2,q3][0]
   r0=q0,i0,X0,Y0,sig0,a0
   r1=q1,i1,X1,Y1,sig1,a1
   r2=q2,i2,X2,Y2,sig2,a2
   r3=q3,i3,X3,Y3,sig3,a3
   # (q0,i0,X0,Y0,sig0,a0),(q1,i1,X1,Y1,sig1,a1),(q2,i2,X2,Y2,sig2,a2),(q3,i3,X3,Y3,sig3,a3) = r0,r1,r2,r3
   return r0,r1,r2,r3 
   
 
def _orderdist12(x,y, filestub, ext, hdr, zmx_xpix, zmx_ypix, mode='cal3',chatter = 0,wheelpos=160):
   ''' compute the second order pixel scale factor for the distance between 
       the first order anchor at position x,y (DET pix coordinates) and the 
       second order anchor, as measured relative to the unscaled zemax model. 
   
       when mode = 'cal3' the observed ratios from zemax will be used.
       otherwise, the a calibrated interpolation 2Dpoly will be used.
       
       2011-01-23 NPMK changed call to make call possible without header 
       2013-05-21 NPMK change 160 to directly return dist12 from fit
   '''
   from scipy.interpolate import bisplev, bisplrep
   from numpy import array, where, isfinite, sqrt, asarray, NaN,zeros
   from .cal3 import get_curvedata_asarray

   check = (len(x) == len(y))
   if not check: print("_orderdist12 input arrays x,y different length")
   if hdr != None: 
      wheelpos = hdr['wheelpos']
      crpix = array( [hdr['crpix1'], hdr['crpix2']] )  # centre of image
      cent_ref_2img = array([1100.5,1100.5])-crpix  
      xx = x - cent_ref_2img[0]
      yy = y - cent_ref_2img[1]
   else:
   # convert coordinates to physical image DET coordinates (xx,yy)
      xx = x - 104 # image coordinates anchors
      yy = y - 78  # ditto

   xx = asarray(xx).flatten()
   yy = asarray(yy).flatten() 
   ratio = (xx + yy).copy()   # initialize (especially the NaN values)
   
   if wheelpos == 160:
      outDist12 = zeros(len(xx))   # initialize (especially the NaN values)
      if mode == 'cal3':
         obsid = get_curvedata_asarray('obsid')
         anchor = get_curvedata_asarray('anchor')        
         dist12 = get_curvedata_asarray('dist12').flatten()
         q = isfinite(dist12) 
         X = anchor[0][q]
         Y = anchor[1][q]
         R = dist12[q]
         qq = (xx > 0) & (xx < 2200) & (yy > -100) & (yy < 2100) 
         tck_dist12 = interpolate.bisplrep(X,Y,R,xb=0.,xe=2200.,yb=-100.,ye=2100.,kx=1,ky=1)
      else:
         tck_dist12 = [array([0.,0.,1079.1683,2200.,2200.]),
             array([-100.,-100.,884.96665,2100.,2100.]),
             array([667.46543,501.466869,515.72518,596.66853,
             646.92233886,535.29734785,672.324228,663.72882342,
             683.239]),1,1]
         qq = (xx > 0) & (xx < 2200) & (yy > -100) & (yy < 2100) 
      for i in range(xx.size):          
         if qq[i]: outDist12[i] = bisplev(xx[i],yy[i],tck_dist12)
      qq = qq == False
      for i in range(xx.size):
         if qq[i]: outDist12[i] = NaN    
      return outDist12
      
   elif wheelpos == 200:
      # initial value (0.962_/- 0.0244)
      # average ratio = 0.962 bilinear fit excludes ratios < 0.935 and > 0.985 
      tck = [array([  200.,   200.,  2100.,  2100.]),\
             array([  400.,   400.,  2100.,  2100.]),\
             array([ 0.91687282,  0.9682413 ,  1.00789688,  0.9860058 ]), 1,1]
      # 2011-06-30 replace:          
      #for i in range(xx.size):      
      #  if ( (xx > 200) & (xx < 2100) & (yy > 400) & (yy < 2100)): 
      #      ratio[i] = bisplev(xx[i],yy[i],tck)
      #  else:
      #   ratio[i] = 0.962 
      # with:
      ratio = 0.962*ratio/ratio
      qq = where( (xx > 200) & (xx < 2100) & (yy > 400) & (yy < 2100))
      xx1 = xx[qq]
      yy1 = yy[qq]
      for i in range(xx1.size):
         ratio[qq[0][i]] = bisplev(xx1[i],yy1[i],tck)
      # end replacement 2011-06-30
        
   elif wheelpos == 955:
      #                   needs refinement!
      ratio = 0.92 
   elif wheelpos == 1000:
      #                   needs refinement!
      ratio = 0.96   
   else:
      print("error in wheelpos found in routine uvotcal.orderdist12; assuming no correction") 
           
   orderdistance = sqrt((zmx_ypix[:,3]-zmx_ypix[1,3])**2+( zmx_xpix[:,3]-zmx_xpix[1,3])**2 ) 
   dist12 = orderdistance[2] * ratio 
   if chatter > 0:
      print('orderdistance ratio applied = ',ratio)
      print('orderdistances : ',orderdistance)
      print('distance 1-2   : ',dist12) 
   return dist12     

def _find_2ndOrderDispersion_scale(wheelpos=160,outfile="find_2ndOrderDispersion.out",chatter=0):
   '''Use the measured dispersion and 
   compare to zemax model to find the 
   scale pixel factor for the second order 
   wavelengths.  
   
   2011-07-05 npmk: added output to file to make a check
   2011-09-15 npmk: changed second order ref point to 2000 
   '''
   from .cal3 import get_curvedata_aslist, get_curvedata_asarray, nominaluv, clockeduv
   from numpy import where, array, zeros, isfinite, polyval,arange,polyfit,sqrt,linspace
   from scipy.interpolate import bisplev, bisplrep
   import os
   from pylab import plot, title, imshow, xlabel, legend, ion, figure, clf,subplot,text
   from uvotgrism import getSpec
   from uvotpy.uvotio import fileinfo
   from uvotpy.uvotmisc import uvotrotvec
  
   if wheelpos == 160:
      curves = uv0160
   elif wheelpos == 200:
      curves = nominaluv
   else:
      curves = clockeduv
      print("WARNING this program only works for the UV grism.") 

   obsdir  = get_curvedata_aslist ('obsdir'    ,curvedata=curves)
   obsid   = get_curvedata_aslist ('obsid'     ,curvedata=curves)
   anchor  = get_curvedata_asarray('anchor'    ,curvedata=curves)
   disp2   = get_curvedata_asarray('disp2nd_2000',curvedata=curves)
   dist12  = get_curvedata_aslist('dist12_2000',curvedata=curves)
   ratio12 = get_curvedata_aslist('ratio12',curvedata=curves).flatten()
   
   dist12 = array(dist12)
   q_ = dist12 > 0       
   q = where(q_)[0]
   X = anchor[0][q]      # in image coordinates
   Y = anchor[1][q]
   D = disp2[:,q]
   m = len(X)
   
   fout = open(outfile,"w")
   fout.write("wheelpos %4i,  number of valid observations = %3i\n" % (wheelpos,m) )
   
   ID = list()
   for i in q[0]: ID.append(obsid[i])
   ID = array(ID)
   
   pixscale = zeros(m)
   C_2zmx = list()
   Xdet = zeros(m)      # in image coordinates
   Ydet = zeros(m)     
   xfldin = zeros(m)      # field x-coordinate
   yfldin = zeros(m)      # field y-coordinate
   ratio12_2000 = zeros(m)    # distance anchor of 2000(2nd order)
   
   ext = 1
   msg =""
   
   # loop over all observations
   ipl = -1
   for i in q[0]:
      ipl += 1
      ax = subplot(m/4+1,4,ipl)
      
      filestub = 'sw'+obsid[i]
      os.chdir('/Volumes/data/grism/'+obsdir[i]+'/') 
      print("current dir: ", os.getcwd())           
      radec = get_radec() 
      if radec == None: 
         radec = get_radec(objectid=obsdir[i]) 
         if radec == None:
            print("problem finding ra,dec for obsdir: ", obsdir[i])    
      ra, dec = radec 
      print("new ra,dec = ", ra,dec)
      os.chdir(obsid[i]+'/uvot/image')  
      print("current dir: ", os.getcwd())          
      Y1 = getSpec(ra,dec,filestub, 1, ) 
      ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
         (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Y1
      anker -= array([(1100.5-hdr['crpix1']),(1100.5-hdr['crpix2'])])
      anker2 -= array([(1100.5-hdr['crpix1']),(1100.5-hdr['crpix2'])])
   
      # uvotcal input to get the orderdistances:
      wheelpos = hdr['wheelpos']   
      crpix = array( [hdr['crpix1'], hdr['crpix2']] )  # centre of image
      cent_ref_2img = array([1100.5,1100.5])-crpix 
       
      Xdet[ipl] = X[ipl] + cent_ref_2img[0]    # anchor measurements in DET coordinates
      Ydet[ipl] = Y[ipl] + cent_ref_2img[1]
      
      specfile, lfilt1_, lfilt1_ext_, lfilt2_, lfilt2_ext_, attfile = \
         fileinfo(filestub,ext,)    
      (anker_uvw1, anker_as, anker_field) =  findAnker(filestub, \
         ext, ra, dec, wheelpos=wheelpos,lfilter=lfilt1_,extf=lfilt1_ext_,chatter=chatter)
      C_zero, C_1, C_2, C_3, C_min1, xpix_, ypix_ , dd, ww, zmxwav, theta = \
         getZmx(anker_field[0],anker_field[1],wheelpos,wpixscale2=1.0,chatter=chatter)

      xfldin[ipl] = anker_field[0]      # field coord measurements 
      yfldin[ipl] = anker_field[1]
      
      msg += ("file number %3i: %s anker2-img [%7.1f,%7.1f],  anker_field [%8.4f,%8.4f]\n" \
           % (i,filestub,anker2[0],anker2[1],xfldin[ipl],yfldin[ipl]))
           
      # find position 2000 A reference point in zemax unscaled model data       
      try:   
         zmx_2k_x = polyval( polyfit(zmxwav,xpix_[2,:],3), 2000.0)  # x coordinate
         zmx_2k_y = polyval( polyfit(zmxwav,ypix_[2,:],3), 2000.0)  # y coordinate
      except:    
         zmx_2k_x = polyval( polyfit(zmxwav,xpix_[2,:],2), 2000.0)  # x coordinate
         zmx_2k_y = polyval( polyfit(zmxwav,ypix_[2,:],2), 2000.0)  # y coordinate
      zmx_dist = sqrt( (zmx_2k_x - xpix_[1,3])**2 + (zmx_2k_y-ypix_[1,3])**2)
      
      ratio12_2000[ipl] = dist12[0][i] / zmx_dist
      
      # compare C_2 and D[i] to get pixel scale at measured points
      # BUT Make sure that the zemax values are valid ones.
      C2_meas = D[:,ipl]  
      
      C_2zmx.append(C_2)    
      #pz = arange(-500,260,25)
      p = arange(-250,550,25)
      pz_p, C_2_ = change2ndOrderReferencePoint(0.,C_2)
      
      wzmx = polyval( C_2_   , p )
      wm   = polyval( C2_meas, p )
      C_2inv     = polyfit(wzmx, p, 2)
      C2_measinv = polyfit(wm,   p, 2)
      
      w = arange(1675.,2425.,25.)   # consider range centered around 2000A
      # pixel distances on w scale: 
      p_2 = polyval(C_2inv,w)
      p2m = polyval(C2_measinv,w)
      # now fit a linear relation between the pixel distances to get the scaling factor at 2000A
      pscl = polyfit(p_2,p2m,1)
      
      msg += "pixscale = "+str(pscl)+" ; C_2 zemax = "+str(C_2)+" ; C_2 meas = "+str(C2_meas)+"\n"      
      
      pixscale[ipl] = pscl[-2]
      
      msg = 'obsid=%s  anker=[%5.0f ,%5.0f] - pixscale=%7.3f \n' % (obsid[ipl],X[ipl],Y[ipl],pixscale[ipl],)
      
      ax.plot(p_2,p2m,'--.')
      text(0,0,"%5.2f"%(pscl[-2]))
      
   #C_2zmx = array(C_2zmx)
      
   if wheelpos == 160: 
      good = array([0,1,3,4,5,6,7,8,9,11,12,13,15,16,17,18,])
   else:
      good = arange(len(Xdet))

   # bispline fitting coefficients
          
   tckD2 = bisplrep(xfldin[good],yfldin[good], D[1,good],xb=-0.165,xe=0.165,\
           yb=-0.165,ye=0.165,kx=1,ky=1,s=len(good) )
           
   tckD1 = bisplrep(Xdet[good],Ydet[good], D[1,good],xb=cent_ref_2img[0],xe=2048.+cent_ref_2img[0],\
           yb=cent_ref_2img[1],ye=2048.+cent_ref_2img[1],kx=1,ky=1)
           
   tckS1 = bisplrep(xfldin[good],yfldin[good], pixscale[good],xb=-0.165,xe=0.165,\
           yb=-0.165,ye=0.165,kx=2,ky=2,s=len(good) )              
           
   tckS2 = bisplrep(Xdet[good],Ydet[good], pixscale[good],xb=cent_ref_2img[0],xe=2048.+cent_ref_2img[0],\
           yb=cent_ref_2img[1],ye=2048.+cent_ref_2img[1],kx=2,ky=2)

   msg += "bispline fit of pixscale in fld coord\n"+str(tckS1)+"\n"        
   msg += "bispline fit of pixscale in det coord\n"+str(tckS2)+"\n"        

   fout.write(msg)
   return Xdet, Ydet, xfldin, yfldin, good, tckD1, tckD2, tckS1, tckS2,D, C_2zmx, pixscale, obsid 
   
#
#    The rest is wrong
#
   for i in range(len(good)):
      Dfit = bisplev(Xdet[good[i]],Ydet[good[i]],tckD1)
      C2_meas= D[:,good[i]]
      msg = "[%6.0f,%6.0f]  %9.2e - %9.2e \n" % (Xdet[good[i]],Ydet[good[i]],Dfit,C2_meas[1],)
      fout.write(msg)
      
   # define a regular grid in x and y to calculate the pixel scale factor based on the 
   # D[1,:] fit/C_2       
   msg = ''
   n1 = 8 # 17
   n2 = 8 # 21 
   xdet1 = linspace(300,1600,n1) #110,2150,n1)  # on the grism detector image
   ydet1 = linspace(200,1700,n2) #80,2120,n2)
   pixscl = zeros(n1*n2).reshape(n1,n2)
   pixsclfld = zeros(n1*n2).reshape(n1,n2)
   xfld = zeros(n1*n2).reshape(n1,n2)
   yfld = zeros(n1*n2).reshape(n1,n2)
   xde = zeros(n1*n2).reshape(n1,n2)
   yde = zeros(n1*n2).reshape(n1,n2)
   c2   = zeros(n1*n2*3).reshape(n1,n2,3)
   pixscale = 0.502 # lenticular filter pixel scale (pix -> arcsec)
   gpixscale = 0.5496 # grism filter pixel scale (pix -> arcsec)
   if wheelpos == 160:
      bs = array(boresight(filter='uc160',r2d=0,date=0) )
   elif wheelpos == 200:   
      bs = array(boresight(filter='ug200',r2d=0,date=0) )
   for i1 in range(n1):
      for i2 in range(n2):
      
         # xd,yd distance to boresight in degrees (convert pixels -> arcsec -> degrees)
         # approximate solution det -> field coordinates
         xde[i1,i2] = xdet1[i1]
         xd = (xdet1[i1]-(bs[0])) * gpixscale / 3600.  
         yde[i1,i2] = ydet1[i2]
         yd = (ydet1[i2]-(bs[1])) * gpixscale / 3600.  
         # rotate to line up to zemax model input 
         xr, yr = uvotrotvec(xd, yd, 64.6) 
         xfld[i1,i2] = xr
         yfld[i1,i2] = yr
         if chatter > 0: print('field point coordinates (%4i,%i4):  (%12.5f, %12.5f ) ' % (i1,i2,xr,yr))
         
         C_zero, C_1, C_2, C_3, C_min1, xpix_, ypix_ , dd, ww, zmxwav, theta = \
         getZmx(xfld[i1,i2],yfld[i1,i2],wheelpos,wpixscale2=1.0,chatter=0)
         
         d_new, C_2new = change2ndOrderReferencePoint(0.0,C_2)
         
         # update position (matching field pos from estimate)
         print("initial xde,yde = ",xdet1[i1],ydet1[i1])
         print("actual  xde,yde = ",xpix_[1,3],ypix_[1,3])
         xde[i1,i2] = xpix_[1,3]
         yde[i1,i2] = ypix_[1,3]
         
         if C_2[1] > 3.0:
            C_2 = [4.4e-4,1.53, 2000.] 
         try:
            c2[i1,i2,0] = C_2new[0]
            c2[i1,i2,1] = C_2new[1]
            c2[i1,i2,2] = C_2new[2]
         except:
            print(' problem computing C_2 for i1=%4i and i2=%4i ' % (i1,i2))
            print(' C_2 = ', C_2)
            pass   
         D1 =  bisplev(xpix_[1,3],ypix_[1,3], tckD1 ) 
         D2 =  bisplev(xr,yr,tckD2 )
         print('interpolated observed dispersion scale = ',D1)
         print('                zemax dispersion scale = ',C_2new[1]) 
         msg += '***** coordinate [det] = '+str(xdet1[i1])+','+str(xdet1[i1])+'  [field] '+str(xr)+','+str(yr)+'\n'
         msg += 'interpolated observed dispersion scale [det] = '+str(D1)+'\n'
         msg += 'interpolated observed dispersion scale [field] = '+str(D2)+'\n'
         msg += '                zemax dispersion scale = '+str(C_2new[1]) +'\n'
         pixscl[i1,i2] = D1 /C_2new[1]
         pixsclfld[i1,i2] = D2 /C_2new[1]
   fout.close()  
   print(msg)
   return ID, X,Y, xde, yde, xfld, yfld, pixscl, pixsclfld, c2, tckD1,tckD2      
 

def get_wpixscale2(Xfield,Yfield,wheelpos=160,mode='field'):
   ''' The ratio of the linear coefficient of the dispersion observed 
   to that from the Zemax model  
   
   
   '''
   from numpy import array
   from scipy.interpolate import bisplev
   if wheelpos == 160 :
     if mode == 'physical':
       # interpolating array for wavelength pixels scale second order. 
       # Reference 2000A
       cpk = [array([  104.,   104.,   104.,  2152.,  2152.,  2152.]),\
              array([   79.,    79.,    79.,  2127.,  2127.,  2127.]),\
         array([ 1.96891673, -0.90334084,  2.89798666, -0.24970502,  2.49445747,\
        -0.15633537,  2.07844283,  0.32677755,  1.49149232]), 2, 2]
       #cpk = [array([  500.,   500.,   500.,   500.,  1700.,  1700.,  1700.,  1700.]),\
       # array([  100.,   100.,   100.,   100.,  1800.,  1800.,  1800.,  1800.]),\
       # array([ 1.16117711,  1.15965983,  1.14729946,  1.12701152,  1.22559361,\
       #       1.22259229,  1.2907682 ,  1.10487131,  1.2141921 ,  1.1778114 ,\
       #       1.17837647,  1.12535922,  1.17473496,  1.13558531,  1.24458354,\
       #       0.9985929 ]),3,3]
       scale = bisplev(Xfield,Yfield,cpk)
       scale = 1./scale
     elif mode == 'field':
       # Reference 2000A
       cpk = [array([-0.165, -0.165, -0.165,  0.165,  0.165,  0.165]),\
              array([-0.165, -0.165, -0.165,  0.165,  0.165,  0.165]),\
         array([ 1.41034609,  0.03028612,  1.93963586,  1.08901025,  1.2697106 ,\
         0.83069402,  1.04361485,  0.97337149,  1.01353148]), 2, 2]
       
       #cpk = [array([-0.17, -0.17, -0.17, -0.17,  0.17,  0.17,  0.17,  0.17]),\
       #       array([-0.17, -0.17, -0.17, -0.17,  0.17,  0.17,  0.17,  0.17]),\
       # array([ 0.80475165,  1.93253556,  0.23513829,  2.16020255,  1.69842352,\
       # 0.33280181,  2.38798937, -0.44123452,  0.4135748 ,  2.18689382,\
       # 0.28017754,  2.41073557,  2.00027706, -0.08131805,  2.14549715,\
       #-0.11550307]),3, 3]
       scale = bisplev(Xfield,Yfield,cpk)
       scale = 1./scale
     else:
      scale = 1.02   
   elif wheelpos == 200: 
      scale = (get_wpixscale(Xfield,Yfield,wheelpos=160,mode='field') ) **2
      if mode == 'physical': 
        # Reference 2000A
        cpk = [array([  104.,   104.,   104.,  2152.,  2152.,  2152.]),\
               array([   78.,    78.,    78.,  2126.,  2126.,  2126.]),\
          array([ 2.05106472, -0.68779771,  1.82393057,  0.83775862,  1.40358864,\
          0.6688332 ,  1.06526787,  1.13388682,  1.29621486]), 2, 2]
        
        #cpk = [array([  100.,   100.,  2100.,  2100.]),\
        #       array([  100.,   100.,  2100.,  2100.]),\
        #       array([ 0.11731603,  1.3685309 ,  1.90546364,  0.60702564]),1,1]
        scale = bisplev(Xfield,Yfield,cpk)
      elif mode == 'field':
        # Reference 2000A       
        cpk = [array([-0.165, -0.165, -0.165,  0.165,  0.165,  0.165]),\
               array([-0.165, -0.165, -0.165,  0.165,  0.165,  0.165]),\
           array([-4.39307994,  6.05070218, -3.93411081,  4.78492019, -1.70098785,\
           2.6548429 , -0.17211748,  1.62823572,  0.90859938]), 2, 2]

        #cpk = [array([-0.17, -0.17,  0.17,  0.17]),\
        #       array([-0.17, -0.17,  0.17,  0.17]),\
        #       array([ 1.72874983, -0.43058469,  0.70831056,  1.71966469]),1,1]        
        scale = bisplev(Xfield,Yfield,cpk)
      else:
        scale = 1.018  # mean value near the centre of the field        
   #elif wheelpos == 1000:
   #elif wheelpos == 955:    
   else:
     scale = 1.0    
   return scale     
   
def get_radec(file='radec.usno',objectid=None,chatter=0):
   '''read the decimal ra,dec from a file 
      or look it up using the objectid name
   '''   
   if objectid == None:
     try:
        f = open(file)
        line = f.readline()
        f.close()
        ra,dec = line.split(',')
        ra  = float( ra)
        dec = float(dec)
        if chatter > 0: 
           print("reading from ",file," : ", ra,dec) 
        return ra,dec
     except:
        return None     
   else:
      import os
      # see http://cdsarc.u-strasbg.fr/doc/sesame.htx 
      # using 'sesame' script from cdsclient package 
      #command = "/sciencesoft/bin/sesame -o2 "+objectid+" > radec.sesame"
      command = "sesame -o2 "+objectid+" > radec.sesame"
      if chatter > 0: 
         print(command)
      if not os.system(command):
         os.system('cat radec.sesame')
         f = open('radec.sesame')
         lines = f.readlines() 
         things = lines[1].split()
         f.close()
         command = "scat -c ub1 -ad "+things[0]+" "+things[1]+" > radec.usnofull"
         if chatter > 0: print(command)
         if not os.system(command):
            f = open('radec.usnofull')
            line = f.readline()
            f.close()
            ra,dec = line.split()[1:3]
            f = open('radec.usno','w')
            f.write("%s,%s" % (ra,dec) )
            f.close()
            ra  = float( ra)
            dec = float(dec)
            return ra,dec
         else:
            if chatter > 0: print('get_radec() error call sesame ')
      else:
         if chatter > 0: print("get_radec() error main call ")
         return None, None          
         
def run_fit(ra,dec,obsid,ext,factr=0.94,alt_anchor=None,alt_angle=None,chatter=0,bright=False,calfile=None):
   ''' not sure yet what this will do yet
     test+develop script
   '''
   from pylab import plot,title,imshow,xlabel,legend,ion,figure,contour,clf,sqrt,array,where,asarray,asscalar,xlim
   import optext2
   from uvotgrism import getSpec, extractSpecImg, findBackground, curved_extraction, find_zeroth_orders
   from uvotpy.uvotio import fileinfo
   from uvotpy.uvotplot import plot_ellipsoid_regions
   from convolve import boxcar
   # from uvotcal import findAnker,getZmx,get_wpixscale,orderdist12

   filestub = 'sw'+obsid
   print('input run_fit:',ra,dec,obsid,ext,factr,alt_anchor,alt_angle,chatter)
   print(' filestub = ',filestub)
   Y1 = getSpec(ra,dec,'sw'+obsid, ext, calfile=calfile, chatter=chatter ) 
   ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
           (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Y1
   # change anker and anker2 to image coordinates
   ankerdet = anker.copy()
   anker2det = anker2.copy()       
   anker -= array([(1100.5-hdr['crpix1']),(1100.5-hdr['crpix2'])])
   anker2 -= array([(1100.5-hdr['crpix1']),(1100.5-hdr['crpix2'])])
   #
   maxcounts = hdr['exposure'] * 8.0 / hdr['framtime']
   #
   # uvotcal input to get the orderdistances:
   wheelpos = hdr['wheelpos']   
   specfile, lfilt1_, lfilt1_ext_, lfilt2_, lfilt2_ext_, attfile = \
       fileinfo(filestub,ext,)    
   print("fileinfo: ", specfile, lfilt1_, lfilt1_ext_, lfilt2_, lfilt2_ext_, attfile)
   (anker_uvw1, anker_as, anker_field) =  findAnker(filestub, \
     ext, ra, dec, wheelpos=wheelpos,lfilter=lfilt1_,extf=lfilt1_ext_,chatter=3,)
   Y2 = getZmx(anker_field[0],anker_field[1],wheelpos,)
   C_zero, C_1_, C_2_, C_3, C_min1, xpix_, ypix_ , dd, ww, zmxwav, theta = Y2
   print("from uvotspec.getSpec: angle = ",angle,"   from uvotcal.getZmx",theta)  
   orderdistance = sqrt((ypix_[:,3]-ypix_[1,3])**2+( xpix_[:,3]-xpix_[1,3])**2 ) 
   wpixscale = get_wpixscale(anker_field[0],anker_field[1],wheelpos=wheelpos,mode='field')
   #
   if alt_anchor != None: 
      anker = alt_anchor
   if alt_angle  != None: 
      angle = alt_angle
   print("using angle = ",angle)   
   #
   if wheelpos < 500: gstub =  'ugu_dt.img' 
   else:  gstub =  'ugv_dt.img'  
   Y = extractSpecImg('sw'+obsid+gstub,ext,anker,angle,)
   dis, spnet, bg, bg1, bg2, extimg, spimg, spnetimg, offset, ank_c = Y
   offset = asscalar(offset)
   cval = -1.0123456789
   Z = findBackground(extimg,cval)
   Zsig = Z[-1][-1].std()
   #
   net = extimg-Z[-1][-1]
   smo = Z[-1][-2]-Z[-1][-1]
   var = extimg
   #
   ZOpos = find_zeroth_orders(filestub, ext, wheelpos,clobber="yes", )
   Xim,Yim,Xa,Yb,Thet,b2mag,matched,ondetector = ZOpos
   dims = asarray( img.shape )
   dims = array([dims[1],dims[0]])
   pivot_ori=array([(anker)[0],(anker)[1]])
   if chatter > 1: 
      print('Fig 3: dims = ',dims,'   pivot_ori = ',pivot_ori)
   figure(3)   
   plot_ellipsoid_regions(Xim,Yim,Xa,Yb,Thet,b2mag,matched,ondetector,pivot_ori,pivot_ori,dims,17.,)

   figure(4)
   clf()
   imshow(extimg-Z[-1][-1],vmin=0.01*Zsig,vmax=4*Zsig)
   contour(extimg-Z[-1][-1],levels=[10,20,40,60,80,140,200,300,450])
   xlim(ank_c[2],ank_c[3])
   # add 2 sigma borders
   title(obsid)
   
   # print out some details
   print("run_fit input to curved_extraction: ")
   print("ank_c = ",ank_c)
   print("ankerdet = ",ankerdet, "   anker used for getting extimg = ",anker)
   print("offset = ",offset)
   print("angle = ", angle)
   print("wheelpos = ",wheelpos)
   print(" ^^^^^^^^^^^^^^^^^^^^^ ")
   # if faint use this :
   if not bright:
     fitorder = curved_extraction(extimg,ank_c,ankerdet,wheelpos,trackonly=True,\
         ZOpos=ZOpos,angle=angle,offset=offset, chatter=1)
            # if bright use a simple slit extraction around the curves :
   else:
     fitorder =  curved_extraction(extimg,ank_c,ankerdet, wheelpos,ZOpos=ZOpos,
         angle=angle,offset=offset, chatter=1) 
      
   (present0,present1,present2,present3),(q0,q1,q2,q3), \
              (y0,dlim0L,dlim0U,sig0coef,sp_zeroth),(y1,dlim1L,dlim1U,sig1coef,sp_first),\
              (y2,dlim2L,dlim2U,sig2coef,sp_second),(y3,dlim3L,dlim3U,sig3coef,sp_third),\
              (x,xstart,xend,sp_all,quality)  = fitorder

   print("orders present:",present0,present1,present2,present3)   
   print("sig0coef : ",sig0coef)     
   print("sig1coef : ",sig1coef)     
   print("sig2coef : ",sig2coef)     
   print("sig3coef : ",sig3coef) 
       
   fig4 = figure(4)  
   plot([ank_c[1],ank_c[1]],[0,200],'k',linewidth=2)   

   if present0: plot(q0[0],y0[q0[0]],'w--')
   plot(q1[0],y1[q1[0]],'w--')
   if present2: plot(q2[0],y2[q2[0]],'w--')
   if present3: plot(q3[0],y3[q3[0]],'w--')
   print("lines of predicted order positions drawn in figure 4 ")

   d12 = _orderdist12(anker[0]+104,anker[1]+78, filestub, ext, hdr, xpix_, ypix_, mode='cal3',chatter = chatter)

   fitorder2, fval, fvalerr = optext2.updateFitorder(extimg, fitorder, wheelpos, full=True,
       predict2nd=True, C_1=C_1, C_2=C_2, d12=d12, chatter=chatter)           
   print("updated fitorder")
   
   (present0,present1,present2,present3),(q0,q1,q2,q3), \
              (y0,dlim0L,dlim0U,sig0coef,sp_zeroth),(y1,dlim1L,dlim1U,sig1coef,sp_first),\
              (y2,dlim2L,dlim2U,sig2coef,sp_second),(y3,dlim3L,dlim3U,sig3coef,sp_third),\
              (x,xstart,xend,sp_all,quality)  = fitorder2
   
   print("orders present:",present0,present1,present2,present3)   
   print("sig0coef : ",sig0coef)     
   print("sig1coef : ",sig1coef)     
   print("sig2coef : ",sig2coef)     
   print("sig3coef : ",sig3coef) 
   
   dims2 = asarray(extimg.shape)
   dims2 = array([dims2[1],dims2[0]])
   pivot= array([ank_c[1],ank_c[0]-offset])
   pivot_ori=anker
   if chatter > 1: 
      print('Fig 4: dims = ',dims2,'   pivot_ori = ',pivot_ori,'   pivot = ',pivot,' angle ',angle-180.0)
   plot_ellipsoid_regions(Xim,Yim,Xa,Yb,Thet,b2mag,matched,ondetector,pivot,pivot_ori,dims2,17.,img_angle=angle-180.0,)
   
   print('order distances vanilla zemax : ',orderdistance)
   print('wavelength pixel scale factor : ',wpixscale) 
   print('order distances scaled zemax : ',orderdistance*wpixscale)
   #
   #  enter a new value for the factor (negative = stop) 
   #
   refine = True
   #d12 = orderdist12(anker[0]+104,anker[1]+78, filestub, ext, hdr, xpix_, ypix_, mode='cal3',chatter = 0)
   
   
   if not bright:
      Y3 = optext2.get_initspectrum(net,var,fitorder2,160,anker,C_1=C_1,C_2=C_2,predict2nd=True,chatter=1)
      counts, variance, borderup, borderdown, (fractions,cnts,vars,newsigmas) = Y3 
      
      plot(q0[0],borderup  [0,q0[0]],'k-')
      plot(q0[0],borderdown[0,q0[0]],'k-')
      plot(q1[0],borderup  [1,q1[0]],'w-',lw=1.2)
      plot(q1[0],borderdown[1,q1[0]],'w-',lw=1.2)
      plot(q2[0],borderup  [2,q2[0]],'y-',lw=1.2)
      plot(q2[0],borderdown[2,q2[0]],'y-',lw=1.2)
      plot(q3[0],borderup  [3,q3[0]],'k-',lw=1.2)
      plot(q3[0],borderdown[3,q3[0]],'k-',lw=1.2)
   
   
   while refine:
     print('starting plot 5; the order distance 1-2 is now: ',d12,' pixels')
     
     fig5 = figure(5)
     clf()   
     
     wav1 = polyval(C_1,x[q1[0]])
     wav2 = polyval(C_2,x[q2[0]]-d12)
     
     if not bright:
        qq1 = where( (counts[1,q1[0]] > -100) & (counts[1,q1[0]] < maxcounts) )
        qq2 = where( (counts[2,q2[0]] > -100) & (counts[2,q2[0]] < maxcounts) )
        plot(wav1[qq1[0]],counts[1,q1[0][qq1[0]]],'b',ls='steps')  
        plot(wav2[qq2[0]],counts[2,q2[0][qq2[0]]],'r',ls='steps')
     title(obsid)
     xlabel('$\lambda(\AA)$')
     if not bright:
        # plot variance
        plot(wav2[qq2[0]],sqrt(variance[2,q2[0][qq2[0]]]),'y--') 
        plot(wav1[qq1[0]],sqrt(variance[1,q1[0][qq1[0]]]),'c--')
        legend(['first order counts','second order counts','second order error','first order error'])
     else:
        plot(wav1, sp_first[q1[0]],'b',ls='steps')      
        plot(wav2,sp_second[q2[0]],'r',ls='steps')
        legend(['first order counts','second order counts',])
     #answer = raw_input('give a new value for the order 1-2 distance [negative to stop] : ') 
     #ans = float(answer)
     ans = -1
     if ans < 0: 
        refine = False
     else:
        d12 = ans       
   if bright:
      return Y,Y1,fitorder,fitorder2,fval,fvalerr,d12, ZOpos, Z, cval
   else:              
      return (wav1,C_1,x,q1,qq1,anker),(wav2,C_2,x,q2,qq2,anker2),\
          (counts, variance, orderdistance, d12, wpixscale,fitorder,fitorder2,fval,fvalerr,\
          fractions,cnts,vars), (img,extimg,Z), Y,Y1,Y3,ZOpos,ank_c

def adjust4BSdrift(date,wheelpos):
   ''' 
   provide the correction dx and dy (in pixels) to be applied to the 
   grism anchor position output from the calibration file. The wavelength 
   calibration depends on fixed positions on the detector and therefore 
   is not dependent on the boresight drift. This routine can be useful 
   for comparison to positions derived using uvotgraspcorr.   
   
   This routine is only useful to understand uvotgraspcorr. For the grism 
   calibration, a fixed boresight point is adopted, since in reality, the 
   calibration depends only on the detector position.    
   '''
   date_cal160 = 256000000.
   date_cal200 = 230500000.
   pixscale = 0.54
   dx1 = -7.79 + 10.8*N.exp( -(date - 116e6)/30.8e6 )
   dy1 =  2.58 - 5.07*N.exp( -(date - 125e6)/50.1e6 )
   if wheelpos == 160:
      dx2 = -7.79 + 10.8*N.exp( -(date_cal160 - 116e6)/30.8e6 )
      dy2 =  2.58 - 5.07*N.exp( -(date_cal160 - 125e6)/50.1e6 )
   elif wheelpos == 200: 
      dx2 = -7.79 + 10.8*N.exp( -(date_cal200 - 116e6)/30.8e6 )
      dy2 =  2.58 - 5.07*N.exp( -(date_cal200 - 125e6)/50.1e6 )
   else:
      return N.array( (0,0) ) 
   return N.array( ((dx1-dx2)/pixscale, (dy1-dy2)/pixscale ))

def change2ndOrderReferencePoint(d12,C_2,oldReferenceWave=2600.,
             newReferenceWave=2000.0):
   '''In order to scale the calibration results in a sensible way 
   to the Zemax model values the original anchor point at 2600A in 
   second order is too far out for an accurate scaling. Too many 
   measurements were extrapolated. The lowest we can use Zemax reliably 
   is about 1900A, so I am selecting 2000A as a new reference point 
   for the second order. That comes to about +300 pixels from the 
   first order anchor point. 
   
   The initial measurements were used to determine d12 and C_2 using the 
   2600A in second order as a reference point. This lies at about 
   +630 pixels from the first order anchor. 
   
   This routine will take d12 and C_2 and conver them to the new reference
   value.'''
   
   import numpy as np
   from uvotpy.uvotgetspec import polyinverse
   
   x = np.arange(1150)
   w = polyval(C_2,x-d12)
   C_2inv = polyinverse(C_2,x-d12)
   d_old = polyval(C_2inv,oldReferenceWave) + d12
   d_new = polyval(C_2inv,newReferenceWave) + d_old
   C_2new = polyfit(x-d_new,w,2)    
   return d_new, C_2new


##################
#  the following routine was initially used to compute the distortion 
#
def findLambdaCorr(X, Y, DL, amin=0.4, amax=0.6, tol=1e-3):
   ''' a,b,c = findLambdaCorr(X,Y,DL)
       X,Y position on detector
       DL = delta lambda
       a = slope to be optimised with solution
           DL = b + c*z
       where
           z = -X + aY     
   '''
   #  minimise DLcorr(X,Y,DL,a) with a in range (amin..amax)
   from scipy.optimize.optimize import brent
   from pylab import polyfit
   a = brent(DLcorr,args=(X,Y,DL),brack=(amin,amax),tol=tol)
   coef = polyfit(zfunc(X,Y,a),DL, 1)
   b = coef[1]
   c = coef[0]
   return a,b,c

def DLcorr(a,X,Y,DL):
   from pylab import polyfit
   coef, R,rank,sing,rcond = polyfit(zfunc(X,Y,a),DL, 1, full='true')
   return R/len(X)
def zfunc(X,Y,a):
   return -X + a * Y   
##################


def EffAreaCal(calobsfile,wave,countrate, disp, spectral_order=1,nboxcar=5, 
               nsmooth=3,expo=None,toterr=None,quality=None,write_file=False ):
   '''create effective area  
   
   input countsrate(wave) (in c/s) - from uvot calibration spectrum
   input calobsfile - filename with ascii table (wave,flux) pairs of CALSPEC (or equivalent) 
                      calibration spectrum flux in ergs cm-2 s-1 A-1
   input disp - Numpy polynomial for dispersion
   input keyword spectral order 1 or 2
   input keyword - nboxcar number of data point in CALSPEC file to smooth 
                   prior to interpolation
   input expo - array length wave, with exposure time for each point. Used same boxcar 
        as the effArea data 
   input quality - array length wave, with quality indices      
   input/output wave - wavelength (in A) of the uvot calibration spectrum                  
   output EffArea(wave) (in cm^2) - effective area 
   output wout,eff - nsmooth-smoothed resampled effective area at 2A steps.
   output qual - for each wout a quality value set similar to a neighbor when set 
   
   '''
   import uvotgrism
   from uvotpy import uvotmisc
   import numpy as np
   try:
      from convolve import boxcar
   except:
      from stsci.convolve import boxcar
   nsm = nsmooth
   version = '111117'
   h_planck = 6.626e-27  # erg - s
   lightspeed = 2.9979e10  # cm/sec
   h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom

   hnu = h_c_ang / wave # ergs per photon
   
   dAdpix = uvotgrism.dAngstrom_dpix_wave(wave,disp,sp_order=spectral_order)  # width pixel in angstroms
   
   tab = uvotmisc.rdTab(calobsfile)
   wcal = tab[:,0]
   fcal = tab[:,1] # erg cm-2 sec-1 A-1
   
   fcalbin = uvotgrism.rebin(wcal,fcal,wave, mode='interpolate',N=nboxcar) 
        
   EffArea = (hnu * countrate) / (fcalbin*dAdpix)               

   q = np.isfinite(EffArea)
   wout1 = int(wave[q].min()+0.5)
   #if wout1/2*2 != wout1: wout1 += 1 # start at even number - only needed when writing every 2 A
   wout = np.arange( wout1, int(wave[q].max()+0.5), 1 )  # value every 1 angstroms
   qual = np.zeros( len(wout), dtype=int )
   effa = uvotgrism.interpol(wout, wave[q], boxcar( EffArea[q], (nsmooth,) ))
   expo2 = uvotgrism.interpol(wout, wave[q], boxcar( expo[q], (nsmooth,) ))
   if toterr != None:
      q = np.isfinite(toterr)
      toterror = uvotgrism.interpol(wout, wave[q], boxcar( toterr[q], (nsmooth,) )) 
      # error should be done with RMS averaging
      # overestimate by sqrt(nsmooth) 
   
   # update quality in qual 
   if quality != None:
      qflags = uvotgrism.quality_flags()
      qw1 = wave.searchsorted(wout,side='left')
      qw2 = wave.searchsorted(wout,side='right')
      for i in range(len(wout)-1):
         if (quality[qw1[i]] == qflags['bad']) | (quality[qw2[i]] == qflags['bad']):
            qual[i] = qflags['bad'] 
         elif (quality[qw1[i]] != qflags['good']):
            qual[i] = quality[qw1[i]]
         elif (quality[qw2[i]] != qflags['good']):
            qual[i] = quality[qw2[i]]     

   if write_file:
      f = open('effective_area_'+str(spectral_order)+'.txt','w')
      f.write('# effective area in cm^2 for order='+str(spectral_order)+'  smoothed='+str(nsm)+'\n')
      f.write('# CALSPEC file = '+calobsfile+'\n')
      for i in range(len(wave)):
         f.write( '%8.2f  %13.4e\n'%(wave[i],EffArea[i]) )
      f.close()
      
      f = open('effective_area_'+str(spectral_order)+'.dat','w')
      f.write('# effective area in cm^2 for order='+str(spectral_order)+'  smoothed='+str(nsm)+'\n')
      f.write('# CALSPEC file = '+calobsfile+'\n')
      for i in range(wout.size):
         f.write( '%8.2f  %13.4e\n'%(wout[i],effa[i]) )  
      f.close()
   return wave, EffArea, wout, effa, expo2, toterror, qual
   
def EffAreaCalFF(calobsfile, phafile, eff_area_out, spectral_order=1, 
                 nsmo_calspec=1, nsmo_obsspec=41, calobs_err=2, apcorr_err=5,
         figno=15, frametime=None, option=2, fudgespec=1.75, sens_rate=0.01,
         ignore_pha_quality=False,wlshift=None,
         auto=False, autoexclude=[], db=None, chatter = 1, clobber=True):
   ''' prepare effective area output from PHA file 
   
   input parameters: 
   
       calobsfile: the calspec file after being uvotified
       phafile   : the extracted spectral file, e.g., sw00056762002ugu_1ord_5.pha
       eff_area_out: name for a file to write the effective area to
       spectral_order=1  - higher orders need code revision
       nsmo_calspec: n-point boxcar smoothing of the calobsfile data (input 'n')
       nsmo_obsspec: n-point boxcar smoothing of the phafile (input 'n')
       figno: figure number which will be opened for selecting bad regions
       calobs_err, apcorr_err - percent error in calibration flux, aperture correction
       sens_rate - rate of change of detector sensitivity per year
       clobber: True to overwrite older version of eff_area_out file
       auto: do not ask for input, read exclude ranges from autoexclude,
             i.e., autoexclude= [[20,40],[50,1200]]
       db the datebase record with the observation data and exclude areas to be edited       
       
   output written to eff_area_out file    
   2012-03-08 added coi-correction column
   2012-07-14 edits to accept database record ; automatic flagging bad regions in record
   2014-08-03 updated to use new coi correction based on box-coi area
   '''
   try:
      import pyfits
   except:
      from astropy.io import fits as pyfits   
   import numpy as np
   from uvotpy import uvotmisc
   import pylab as plt
   import os
   import uvotgrism
   from uvotpy import uvotgetspec 
   try:
      from stsci.convolve import boxcar
   except: 
      from convolve import boxcar
    
   if db != None: 
      not2include = db[1]['exclude']
      name = db[1]['name']
      if wlshift == None:
         wlshift = db[1]['wlshift']
      expo = db[1]['exposure']
   else: 
      not2include = []
      name = calobsfile  
      wlshift = 0 
      expo = -1
   print("applying initial wlshift to spectrum = ",wlshift)   
   if chatter > 2: 
      print("exclude: ",not2include)
      print("input argument db : ",db)
      
   try:
      hdu = pyfits.open(phafile)
      if chatter > 2: print(phafile+" opened successfully")
   except:
      raise IOError(" error opening "+phafile )

   wheelpos = hdu[1].header['wheelpos']  
   _good = np.isfinite(hdu[2].data['lambda'])    
   wavein = hdu[2].data['lambda'][_good] + wlshift
   flx = hdu[2].data['flux'][_good]
   if wheelpos < 500:
      qwav = np.where( (wavein > 1600) & (wavein < 7000) & np.isfinite(flx) & (flx >= 0.))   # limit to valid range
   else:   
      qwav = np.where( (wavein > 2400) & (wavein < 7000) & np.isfinite(flx) & (flx >= 0.))   # limit to valid range
   wavein = wavein[qwav]
   flx = flx[qwav]
   obs = hdu[2].data
   hdr2 = hdu[2].header
   pixno = obs['pixno'][_good][qwav]
   countrate = obs['netrate'][_good][qwav]
   quality = obs['quality'][_good][qwav]
   if ignore_pha_quality: quality = 0*quality
   exposure = obs['exposure'][_good][qwav]
   if expo == -1 : expo = exposure[0]
   if frametime == None: 
      frametime = hdu[2].header['framtime']
   print("frame time = ",frametime) 
   bkg1 = obs['bg_l'][_good][qwav] 
   bkg2 = obs['bg_r'][_good][qwav]
   if 'COSP1RAT' in obs.names: # get rates in coi-box
        obstot = obs['COSP1RAT']  # total rate in coi-box
        obsbck  = obs['COBG1RAT']  # background rate in coi-box
        coi_width = hdr2['coiwidth']
        obsnet = obstot-obsbck
        print("coincidence loss version 2, coi width=",coi_width)
   elif (option != 1):
       print("coincidence loss version 1 ")
   else:
       raise IOError("spectrum file does not include coincidence box data")

   if chatter > 2: print('@computing error')
   # percent error in the observation due to counting statistics only
   obs_err = np.sqrt(np.abs(countrate*exposure+bkg1+bkg2)) / exposure / countrate /100 

   # total error RMS of calibration spectrum and observed count rate
   tot_error = np.sqrt(obs_err**2+calobs_err**2+apcorr_err**2)
    
   hdr1 = hdu[1].header
   tstart = hdr1['tstart']
   sens_corr = sensitivityCorrection(tstart,sens_rate=sens_rate)
   hist = hdr1.get_history()
   
   C_1 = uvotmisc.get_dispersion_from_header(hdr1)

   SIG1 = uvotmisc.get_sigCoef(hdr1)  
   
   trackwidth = 2.5
   try:
      trackwidth = float(uvotmisc.get_keyword_from_history(hist,'TRACKWID'))
   except:
      pass   
   
   if chatter > 0:
      print("dispersion is ", C_1)
      print("trackwidth is ", trackwidth)
      print("SIG1       is ", SIG1)
   
   # get the coincidence loss 
   bkgrate = 0.5*(bkg1 + bkg2)
   
   if 'COSP1RAT' in obs.names:   
       print("coi version 2 computation")
       print("chatter = ",chatter)
       fcoi, vcoi = uvotgetspec.coi_func(pixno,wavein,
           obstot[_good][qwav],
           obsbck[_good][qwav],
           option=1,
           background=False,
           debug=False,chatter=chatter)       
       # background cr coi    
       fbgcoi = uvotgetspec.coi_func(pixno,wavein,
           obstot[_good][qwav],
           obsbck[_good][qwav],
           option=1,
           background=True,
           debug=False,chatter=chatter)   
   else:
       print("coi version 1 !!")
       # net spectrum cr coi
       fcoi = uvotgetspec.coi_func(pixno,wavein,countrate,bkgrate,
           sig1coef=SIG1,option=option,
           fudgespec=fudgespec,background=False,trackwidth=trackwidth,
           debug=False,chatter=chatter)
       
       # background cr coi    
       fbgcoi = uvotgetspec.coi_func(pixno,wavein,countrate,bkgrate,
           sig1coef=SIG1,option=option,
           fudgespec=fudgespec,background=True,trackwidth=trackwidth,
           debug=False,chatter=chatter)

   # make wavelength adjustment to match the observed wavelengths to the calibration file.
   
   tab1 = uvotmisc.rdTab(calobsfile)

   happy = False
   ans = "No"
   #fig1 = plt.figure(figno+1) ; fig1.clf()
   dw=0
   titl = 'uvotcal.EffAreaCalFF:  '+name+" t_exp="+str(expo)

   if wlshift == None:  
       print("applying wavelength shifts ")
       fig = plt.figure(figno)
       try:
          # correlate the spectra and shift 
          if chatter > 0:
              print("spectrumpixshift input:")
          dw = uvotgetspec.spectrumpixshift(wavein, flx, tab1[:,0].flatten(),
               tab1[:,1].flatten(), wmin=2400, wmax=4200, 
               spectrum=False, delwav=True, chatter=chatter )
          print("SUGGESTED WAVELENGTH SHIFT = ", dw) 
       except:
           print("WARNING: could not get wavelength shift:")
           dw = 0
           pass
   #
       while not happy:
           fig.clf()
           ax = fig.add_subplot(111)
           ax.plot(wavein, flx,label='observed')
           ax.plot(wavein, boxcar(flx,(20,)),'k',alpha=0.4,label='obs smoothed')
           ax.plot(tab1[:,0],tab1[:,1],label='calibration src')
           ax.set_xlim(1600,7000)
           ax.set_ylim(1e-16,1e-12)
           ax.set_title(titl)
           ax.legend(loc=0,markersize=6)
      
           if (not auto): ans = input(
              "See Plot. Happy with wavelength (no offset) ? (yes/no/abort): ")
           else: ans = 'Y'   
      
           try:
               if ans.upper()[0:4] == 'ABORT':
                   break
                   raise 
           except: 
               pass

           if (ans.upper()[0] == 'Y'):  
               print(" ") 
               happy = True
               break
           if (not auto): 
              ans=input("give number of angstroms to shift (cumulative) : ")
           else: ans = 0
      
           try:
               ans = float(ans)
               wavein += ans
               wlshift += ans
               clf()
           except:
               print("there seems to be a problem with the shift")
     
   if (not auto): 
           print("Now you have the opportunity to exclude ranges for the effective area.")
   
   wave, EffArea, wout, effa, expo, toterr, qual = EffAreaCal(
                        calobsfile,wavein,countrate, 
                        C_1, expo=exposure, toterr=tot_error, 
                        quality=quality, spectral_order=1,
                        nboxcar=nsmo_calspec, nsmooth=nsmo_obsspec)
   print("call EffAreaCal done")                 
   effa *= sens_corr            
   qw = np.ones(len(wout), dtype=bool)
   qw_ = np.ones(len(wout), dtype=bool)
   
   if len(not2include) > 0:
     if len(not2include[0]) > 0:
        for weg in not2include:
           wlo,wup = weg
           qw_ = qw_ & ((wout < wlo) | (wout > wup))

   fig1= plt.figure(figno+1)
   fig1.clf()
   ax1 = fig1.add_subplot(111)
   ax1.plot( wout,effa,'b' )
   if not qw_.all(): 
       ax1.plot(wout[qw_ == False],effa[qw_ == False],'.r',markersize=4)
   if wheelpos > 500:
       ax1.set_xlim(1600,7000)
       ax1.set_ylim(-10,35)
   else:
       ax1.set_xlim(2400,7000)   
       ax1.set_ylim(-10,50)
   ax1.set_title(titl)
   happy = False

   while not happy:
      print("corrected effective area OK?")
      if (not auto): 
         ans = input('See plot. Happy with effective area ? (yes/no)')
      else: 
      # --> code to exclude auto regions
         if len(autoexclude) == 0:   
            ans = 'Y'
         else:
            for xrange in autoexclude:
                try:
                   (wlo,wup) = xrange
                   qw = qw & ((wout < wlo) | (wout > wup)) 
                except:
                   print("autoexclude invalid range : ",xrange," --- skipping ")
                   pass 
            # completed automatic exclusions into qw                  
            happy = True
            break                  
      
      # the following lines interactively let you make modifications (auto = False)  
      if ans.upper()[0] == 'Y': 
         happy = True
         break
     
      ans = input('give wavelength range to exclude as lower,upper: ')
      
      try:
         wlo,wup = np.array(ans.split(','),dtype=float)
         ans = input('confirm exclude range = %6.1f,%6.1f  (Y): '%(wlo,wup)) 
         if ans.upper()[0] == 'Y': 
            qw = qw & ((wout < wlo) | (wout > wup))
            if not qw_.all(): plt.plot(wout[qw_ == False],effa[qw_ == False],'.r')
            plt.plot(wout[qw], effa[qw], 'b', alpha=0.8)
            if not qw.all(): plt.plot(wout[qw == False], effa[qw == False],'.k',)
            plt.title(name)
            plt.ylim(-10,42)
            print("next should be plot and then ask last time happy?")
            not2include.append([wlo,wup])
            #continue
      except:
         print("Try again: there seems to be a problem reading the answer")
         print("ans = ",ans)
         #print "wlo, wup = ",wlo,wup       

      print("qw_.all()=",qw_.all())
      ax1.cla()
      if not qw_.all(): ax1.plot(wout[qw_ == False],effa[qw_ == False],'.r',markersize=8)
      ax1.plot(wout[qw], effa[qw], 'b', alpha=0.8)
      if not qw.all(): ax1.plot(wout[qw == False], effa[qw == False],'.k',markersize=8)
      #plt.plot(eanorm[:,0],eanorm[:,1],'m',lw=0.5,alpha=0.5)
      ax1.set_xlim(1600,7000)
      ax1.set_ylim(-10,42)
      #plt.plot( wout[qw], effa[qw], )
      ax1.set_title(titl)
   print("finished filtering bad areas out of EA")    
      
   # happy! write output
   qw = qw & qw_  # merge old and new
   
   if  os.access(eff_area_out,os.F_OK) & (not clobber):
      print("output file is present and clobber is not set")
      eff_area_out = 'temp_effarea_outfile.txt'   
      print("writing instead to file (and clobbering if there): "+eff_area_out)
   
   try:
     coi = fcoi(wout)
   except:
     print(" problem with call to fcoi. wout = ",wout)
     print(" returning function fcoi(wave) ")
     print("argument: ",arange(1700,4500,100))
     print("fcoi: ",fcoi(arange(1700,4500,100)))
     raise RuntimeError("Help!")

   print("start write output file to ",eff_area_out)   
   f = open(eff_area_out,'w')
   qflags = uvotgetspec.quality_flags()
   N = len(wout)
   qual[(qw == False)] = qflags['bad']
   
   f.write("#effective area from observation %s\n" % (phafile))
   f.write("#calspec file used %s\n"%(calobsfile))
   f.write("#TSTART = %12.1f, sensitivity correction applied to effArea = %5.2f\n"%(tstart,sens_corr))
   f.write("#program uvotcal.EffAreaCalFF(), smoothing calspec:%s, obs:%s\n"%(nsmo_calspec,nsmo_obsspec))
   f.write("#  coi factor has not been applied \n")
   f.write("# background coi factor for selected\n# wave, coi-factor\n")
   
   for i in range(1,N,int(N/9)):
      f.write("#BGCOI %7.1f  %7.4f \n" % (wout[i], fbgcoi(wout[i])))
   f.write("# \n")
   f.write("#wave(A), effArea(cm2), exposure(s), total error (percent), quality flag, coi_factor\n")
   for i in range(1,N):    
      f.write("%6.1f  %7.3f  %10.3f  %8.3f  %2i  %8.4f\n" % (wout[i],effa[i],expo[i],toterr[i],qual[i],coi[i]) )        
   f.close()
   print("finished with writing the effective area to "+eff_area_out)
   if db != None: 
      db[1]['exclude']=not2include  # replace the exclusion 
      db[1]['wlshift']=wlshift
      db[1]['nproc'] += 1
      return db
   else: return None   
   

def sensitivityCorrection(swifttime,sens_rate=0.01):
   '''
   give the sensitivity correction factor to multiply to the rate/flux 
   
   swifttime : float
      time of observation since 2005-01-01 00:00:00  in seconds, usually TSTART 
      
   Notes
   -----
    A 1%/year decay rate since 2005-01-01 has been assumed and 
    the length of the mean Gregorian year was used
   '''
   sens_corr = 1.0/(1.0 - sens_rate*(swifttime-126230400.000)/31556952.0 )  
   return sens_corr
      

def mergeEffAreas(EAfiles, eff_area_out,figno=15,markers=None,wait=False,nsiglim=5,randerr=2):
   '''The effective areas are read in from a list EAfiles, and 
   merged to output file eff_area_out
   The input data are assumed to be on a whole angstroem grid
   
   nsiglim : standard deviation limit for data 
   randerr  : additional random unspecified UVOT detector stability? error 
   
   2012-03-08 added coi-correction support
   2012-07-30 renamed exp expo' divided plot into subplots
   '''
   from uvotpy.uvotmisc import rdTab
   import numpy as np
   from uvotgrism import quality_flags
   import pylab as plt
   
   qflags = quality_flags()
   wlmin = 7000
   wlmax = 1600

   # find the wavelength range by a quick browse through the files
   for fil in EAfiles:
      t = rdTab(fil)
      wlmin = np.min([t[0,0], wlmin])
      wlmax = np.max([t[-1,0], wlmax])
      
   wave = np.arange(wlmin, wlmax+1, 1)
   effArea = np.zeros(len(wave),dtype=float)
   quality = np.zeros(len(wave),dtype=int)
   exposure = np.zeros(len(wave),dtype=float)
   coi_corrected_ea = np.zeros(len(wave),dtype=float)
   sum1 = np.zeros(len(wave),dtype=float)  # weighted ea sum
   sum2 = np.zeros(len(wave),dtype=float)  # sum weights
   sum3 = np.zeros(len(wave),dtype=float)
   sum4 = np.zeros(len(wave),dtype=float)
   coi_err = np.zeros(len(wave),dtype=float)
   effArea_err = np.zeros(len(wave),dtype=float)+0.01
   ndata_err = np.zeros(len(wave),dtype=float)

   fig = plt.figure(figno)
   plt.clf()
   plt.rcParams['legend.fontsize'] = 'small'
   
   # add the good data to the arrays
   k = -1
   for fil in EAfiles:   
      t = rdTab(fil)
      k += 1
      w = t[:,0]         # wave
      ea = t[:,1]        # effective area (no coi)
      expo = t[:,2]      # exposure time
      terr = t[:,3]      # percent error
      qa = t[:,4]        # quality flag
      coi = t[:,5]       # coi-factor to multiply eff. area with
      
      good = (qa == qflags['good'])
         
      if markers == None:
         ax1 = plt.subplot(3,2,1)
         plt.plot(w[good],ea[good],'.',markersize=3)
         #plt.ylabel('Eff Area (cm$^2$) uncorrected')
         ax2 = plt.subplot(2,1,2)
         plt.plot(w[good],ea[good]*coi[good],'.',markersize=3)
         #plt.ylabel('Eff Area (cm$^2$) coi-corrected')
      else:
         ax1 = plt.subplot(3,2,1)
         plt.plot(w[good],ea[good],markers[k],markersize=4)     
         #plt.ylabel('Eff Area (cm$^2$) uncorrected')
         ax2 = plt.subplot(2,1,2) 
         plt.plot(w[good],ea[good]*coi[good],markers[k]+'--',markersize=4)      
         #plt.ylabel('Eff Area (cm$^2$) coi-corrected')
      if wait: 
        ans = input(fil+' > ')
                  
         
      w = w[good]
      ea = ea[good]
      expo = expo[good]
      terr = terr[good]
      coi = coi[good]
      
      coi_e = (coi-1.0) * 1.5 # placeholder -> 1.5% of coi-correction 
      # is error not included in error output
      coi_e[ coi_e < 0. ] = 0.
      wt = 1.0/(terr**2)
      
      # add the good data to the array (one wave at a time)
      for i in range(len(w)):
         p = np.where(w[i] == wave)
         effArea[p] += ea[i] * expo[i]
         exposure[p] += expo[i]
         sum1[p] += ea[i]*wt[i]
         sum2[p] += wt[i]             
         sum3[p] += terr[i]**2
         coi_corrected_ea[p] += ea[i]*wt[i]*coi[i]
         coi_err[p] += coi_e[i]
   
   LEG = []
   LEG1 = []
   for ll in EAfiles: 
      LEG.append(ll.split('/')[-1])   # list of filenames 
      LEG1.append(ll)   # list of fully qualified filenames
   plt.subplot(3,2,1)
   plt.legend(LEG,bbox_to_anchor=(1.02, 1.1), loc=2, borderaxespad=0.,numpoints=5,ncol=2)
   plt.ylabel('Eff Area (cm$^2$) uncorrected')
   plt.subplot(2,1,2)
   plt.ylabel('Eff Area (cm$^2$) coi-corrected')  
   plt.xlabel('$\lambda (\AA)$',fontsize=14)
     
   # now remove all points with no data
   p = (exposure > 1e-4)
   wave = wave[p]
   effArea = effArea[p] 
   exposure = exposure[p]
   sum1 = sum1[p]
   sum2 = sum2[p]
   sum3 = sum3[p]
   coi_corrected_ea = coi_corrected_ea[p]
   coi_err = coi_err[p]
   
   # normalise 
   effArea /= exposure
   sum1 /= sum2   # this is effArea weighted by the inverse variance of the relative errors
   coi_corrected_ea /= sum2 # this is effArea weighted by the inverse variance and coi-corrected
   coi_err /= len(EAfiles)  # mean coi-error (percentage)
   effArea_err1 = np.sqrt( (np.sqrt(sum3)/len(EAfiles))**2 + randerr**2)   # RMS error in 
   # effective area due to errors streamed down (percent)
   effArea_err2 = np.sqrt(1./sum2 + randerr**2)  # error based on weights
   
   # add to plot
   plt.subplot(3,2,1)
   plt.plot(wave, sum1,'k',lw=2.5)
   LEG.append('merged')
   plt.title('effective area - merging data')
   plt.legend(LEG,bbox_to_anchor=(1.02, 1.1), loc=2, borderaxespad=0.,numpoints=5,ncol=2)
   plt.subplot(2,1,2)
   plt.plot(wave, coi_corrected_ea,'k',lw=2.5)
   
   # write outout file
   
   f = open(eff_area_out,"w")
   f.write("#file1: effective area with errors determined from downflowing errors in input...\n")
   f.write("#effective area files used are:\n")
   for fil in EAfiles:  f.write("# +file: %s\n" % (fil))
   f.write("#program uvotcal.mergeEffAreas()\n")
   f.write("#        - error weighted - exposure weighted; need to add systematic uncertainty coi-correction.\n")
   f.write("#wave(A) - effArea(cm2)   -   effArea        - effArea_error(percent)  - coi-corrected effArea - alternative error \n")
   
   for i in range(len(wave)):
      f.write("%6.1f   %6.2f   %6.2f  %6.2f   %6.2f  %6.2f\n"%(wave[i],sum1[i], effArea[i], \
          effArea_err2[i],coi_corrected_ea[i],effArea_err1[i]))
   f.close()
   

   # second iteration: error based on spread individual eff. area curves:

   k = -1
   
   effArea_err[p] = effArea_err1/100.0  # start with error based on weights from individual spectra ; add RMS from mean
   
   for fil in EAfiles:   
      t = rdTab(fil)
      k += 1
      w = t[:,0]         # wave
      ea = t[:,1]        # effective area (no coi)
      terr = t[:,3]      # percent error
      qa = t[:,4]        # quality flag
      coif = t[:,5]      # coi factor
      
      good = (qa == qflags['good'])
      w = w[good]
      ea = ea[good]
      terr = terr[good]
      wt = 1.0/(terr**2)
      #
      # add the good data to the array (one wave at a time)
      dist = []
      
      for i in range(len(w)):
         p1 = np.where(w[i] == wave)
         d1 = (ea[i]*coif[i] - coi_corrected_ea[p1])         # difference from coi-corrected mean
         dist.append(np.abs(d1)/(terr[i]/100.0*sum1[p1]*nsiglim))
         effArea_err[p1] += d1**2 * wt[i]  # add to error
         ndata_err[p1] += 1
      d2 = np.array(dist)    # change list to numpy array
      if d2.mean() > 0.1:
         print("possibly large deviation from mean: "+fil+" has average d=",d2.mean())    
   # now remove all points with no data
   effArea_err = np.sqrt(effArea_err[p]/sum2)   # only points p were valid
   effArea_err *= 100.0/sum1  # percentage error 
   
   # write outout file
   
   f = open(eff_area_out+'.err',"w")
   
   f.write("#file2: effective area error from RMS differences observed EA curves...\n")
   f.write("#effective area files used are:\n")
   for fil in EAfiles:  f.write("# +file: %s\n" % (fil))
   f.write("#program uvotcal.mergeEffAreas()\n")
   f.write("#        - error weighted - exposure weighted; need to add systematic uncertainty coi-correction.\n")
   f.write("#wave(A) - effArea(cm2)   -   effArea        - effArea_error(percent)"+\
       "  - coi-corrected effArea - input error percent - ndata\n")
   
   for i in range(len(wave)):
      f.write("%6.1f   %6.2f   %6.2f  %6.2f   %6.2f  %6.2f  %4i\n"%(wave[i],sum1[i], effArea[i], \
          effArea_err2[i],coi_corrected_ea[i],effArea_err[i],ndata_err[i]))
   f.close()
   #
   print("completed merge effective area files ",EAfiles)


   
def coi_spectral_calibration(calobsfile,effective_area_file,phafile,coi_output_file,
    coi_length=26,frametime=None, fudgebg=1.0, fudgespec=1.55, option=1, chatter=1):
   '''
   use the calibration spectrum and effective area to determine incoming counts per frame 
   use the extracted spectrum and background to determine the total 
   measured counts in the spectrum per frame, and the measured background counts 
   in the spectrum per frame, averaged over the coi_length (in pixels).
      
   background coi set to classic theroy when fudgebg=1   
   output the coi_correction as a function of counts per frame
   
   2012-03-05 NPMK 
   
   '''   
   from uvotpy.uvotmisc import rdTab
   import numpy as np
   import uvotgrism
   from convolve import boxcar
   

   h_planck = 6.626e-27  # erg - s
   lightspeed = 2.9979e10  # cm/sec
   h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
   
   calobs = rdTab(calobsfile)
   wcal = calobs[:,0]
   fcal = calobs[:,1]
   
   effarea = rdTab(effective_area_file)
   
   wmin = effarea[:,0].min()
   wmax = effarea[:,0].max()

   hdu = pyfits.open(phafile)
 
   if frametime == None: frametime = hdu[1].header['FRAMTIME']
   alpha = hdu[1].header['deadc']
   pixno = hdu[2].data['pixno']
   wave = hdu[2].data['lambda']
   countrate = hdu[2].data['netrate']
   exposure = hdu[2].data['exposure']
   bkg1 = hdu[2].data['bg_l'] 
   bkg2 = hdu[2].data['bg_r']
   
   # limit to wavelengths in the effective area 
   v = (wave > wmin) & (wave < wmax)
   pixno = pixno[v]
   wave = wave[v]
   countrate = countrate[v]
   exposure = exposure[v]
   bkg1 = bkg1[v]
   bkg2 = bkg2[v]
 
   hdr1 = hdu[1].header
   hist = hdr1.get_history()
   
   C_1 = [uvotmisc.get_keyword_from_history(hist,'DISP1_0')]
   try:
      coef=uvotmisc.get_keyword_from_history(hist,'DISP1_1')
      if coef != None: C_1.append(coef)
      try:
         coef=uvotmisc.get_keyword_from_history(hist,'DISP1_2')
         if coef != None: C_1.append(coef)
         try:
            coef=uvotmisc.get_keyword_from_history(hist,'DISP1_3')
            if coef != None: C_1.append(coef)
            try:
               coef = uvotmisc.get_keyword_from_history(hist,'DISP1_4')
               if coef != None: C_1.append(coef)
            except:
               pass
         except:
            pass   
      except:
         pass    
   except:
      pass   
   C_1 = np.array(C_1,dtype=float)
   trackwid = float(uvotmisc.get_keyword_from_history(hist,'TRACKWID') )
   SIG1 = [uvotmisc.get_keyword_from_history(hist,'SIGCOEF1_0')]
   try:
      coef=uvotmisc.get_keyword_from_history(hist,'SIGCOEF1_1')
      if coef != None: SIG1.append(coef)
      try:
         coef=uvotmisc.get_keyword_from_history(hist,'SIGCOEF1_2')
         if coef != None: SIG1.append(coef)
         try:
            coef=uvotmisc.get_keyword_from_history(hist,'SIGCOEF1_3')
            if coef != None: SIG1.append(coef)
            try:
               coef = uvotmisc.get_keyword_from_history(hist,'SIGCOEF1_4')
               if coef != None: SIG1.append(coef)
            except:
               pass
         except:
            pass   
      except:
         pass    
   except:
      pass   
   SIG1 = np.array(SIG1,dtype=float)
   
   frames = exposure/frametime
   bkgrate = 0.5*(bkg1+bkg2)
   
   # average over the coi_length, and convert to counts per frame
   obs_countsperframe = boxcar((countrate + bkgrate) * frametime, (coi_length,))
   bkg_countsperframe = boxcar( bkgrate * frametime, (coi_length,) )

   # convert calibration spectrum flux to counts per frame
   hnu = h_c_ang / wave # ergs per photon   
   dAdpix = uvotgrism.dAngstrom_dpix_wave(wave,C_1,sp_order=1)  # width pixel in angstroms
   fcalbin = uvotgrism.rebin(wcal,fcal,wave, mode='interpolate',N=1) 
   EffArea = uvotgrism.rebin(effarea[:,0],effarea[:,1],wave, mode='interpolate',N=1)
   
   cal_countsperframe = EffArea * fcalbin * dAdpix / hnu * frametime 
   
   # using the point-source coi-correction formula for the background, 
   # find the predicted incident background 
   # correct to 314 subpixels in size. 
   sigma1 = np.polyval(SIG1, pixno)
   factor = 314./trackwid/sigma1/2
   # the fudge factor accounts for a larger effective lenght scale for the coi-loss 
   # for the spectrum, it depends on something else, maybe background, maybe brightness
   
   try:
      yy = (1.0 - fudgebg*factor*bkg_countsperframe)
      yy[ yy > 0.00001 ] = 0.00001
      bkg_cpf_incident = (-1.0/alpha) * np.log(yy)/(fudgebg*factor)
      if option == 1:
         if fudgespec == None: fudgespec=1.0+31.*bkg_cpf_incident
         yy = (1.0 - fudgespec*factor*obs_countsperframe)
         yy[ yy > 0.00001 ] = 0.00001
         obs_cpf_incident = (-1.0/alpha) * np.log(yy)/(fudgespec*factor)
         # tested this - works to within 10% 
      if option == 2:
         yy = (1.0 - factor*obs_countsperframe)
         yy[ yy > 0.00001 ] = 0.00001
         obs_cpf_incident = (-1.0/alpha) * np.log(yy)/factor
         # classic coi formula 
      if option == 3:
         if fudgespec == None: fudgespec = 0.5
         yy = (1.0 - factor*(obs_countsperframe+fudgespec*bkg_countsperframe))
         yy[ yy > 0.00001 ] = 0.00001    
         obs_cpf_incident = (-1.0/alpha) * np.log(yy)/(factor*(1+fudgespec*bkg_countsperframe/obs_countsperframe))
         # treating fudgespec as an additive component - ie. background contributes to larger area.
         # over estimate incident obs_cpf for increasing bkg_countsperframe; similar results as option 1
   except:
      print("ERROR: counts in background too high for point source coi formula.")
      return
   
   
   # now we compute the coi-correction \psi
      
   psi = (cal_countsperframe+bkg_cpf_incident) / (obs_countsperframe) 
               
   # perhaps recast coi(rate)? 
        
   f = open(coi_output_file,"w")
   f.write('# coi correction to the observed rate as function of input rate and observed rate\n')
   f.write('# input: calobsfile %s \n'% (calobsfile))
   f.write('# input: effective_area_file %s \n' % (effective_area_file))
   f.write('# input: observation phafile %s \n' % (phafile))
   f.write('# coi length (pix) %i \n' % (coi_length))
   f.write('# wave  observed_rate   incident_rate   COI_correction \n')
   for i in range(len(wave)):
      f.write("%7.1f  %10.5f  %10.5f  %10.5f\n" % (wave[i],\
          obs_countsperframe[i]/frametime,\
          (cal_countsperframe[i]+bkg_cpf_incident[i])/frametime,psi[i] ) )   
   f.close()
   #
   return wave, frametime, obs_countsperframe,cal_countsperframe,bkg_countsperframe,\
          bkg_cpf_incident, obs_cpf_incident, psi 


def _util_img_position_ea(angle,anchor_det,wave,C_1=None):
   ''' 
   2012-02-08 NPMK utility to calculate the approx. position of the image where 
   a certain wavelength falls of a spectrum with given anchor point, angle, and 
   dispersion C_1. Dispersion defaults to central dispersion if not passed. 
   
   example:
   
   position = _util_img_position_ea(35.1,array([1098.8,1025.7]),2700.0, )
   
   '''
   from uvotpy.uvotmisc import uvotrotvec
   from numpy import polyval,array, asarray
   import uvotgrism
   
   if C_1 == None: C_1 = array([4.1973e-10,-1.3010e-6,1.4366e-3,3.2537,2607.6])
   dx = (uvotgrism.pix_from_wave(C_1, wave))[0]
   angle = float(angle)
   theta = -(180.0-angle) # angle to rotate 
   ank_img = asarray(anchor_det) - array([104,78])
   pos = ank_img + array( uvotrotvec(dx,100.,theta) )
   return pos

def _util_fix_eff_area_stuff(file1='effective_area_160_table1.1.txt',file2='effective_area_160_table2.1.txt'):
   '''read the effective area tables and add the positions '''
   from uvotpy.uvotmisc import rdList
   import numpy as np
   
   f1 = rdList(file1)
   f2 = rdList(file2)
   Nrow = len(f1)
   cross_ref_no = [] ; zone = []
   anker_det_x = [] ;  anker_det_y = []
   exposure = [] ;     angle = []
   ea_17 = [] ; xy_17x = [] ; xy_17_y = []
   ea_19 = [] ; xy_19x = [] ; xy_19_y = []
   ea_21 = [] ; xy_21x = [] ; xy_21_y = []
   ea_23 = [] ; xy_23x = [] ; xy_23_y = []
   ea_25 = [] ; xy_25x = [] ; xy_25_y = []
   ea_27 = [] ; xy_27x = [] ; xy_27_y = []
   ea_28 = [] ; xy_28x = [] ; xy_28_y = []
   ea_32 = [] ; xy_32x = [] ; xy_32_y = []
   ea_36 = [] ; xy_36x = [] ; xy_36_y = []
   ea_40 = [] ; xy_40x = [] ; xy_40_y = []
   ea_45 = [] ; xy_45x = [] ; xy_45_y = []
   ea_50 = [] ; xy_50x = [] ; xy_50_y = []
   ea_60 = [] ; xy_60x = [] ; xy_60_y = []
   for i in range(Nrow): 
      cross_ref_no.append(f1[i][0] )
      zone.append(        f1[i][1] )
      anker_det_x.append( f1[i][2] )
      anker_det_y.append( f1[i][3] )
      exposure.append(    f1[i][4] )
      angle.append(       f1[i][5] )
      
      ea_17.append(     f1[i][6] )
      ea_19.append(     f1[i][7] )
      ea_21.append(     f1[i][8] )
      ea_23.append(     f1[i][9] )
      ea_25.append(     f1[i][10])
      ea_27.append(     f1[i][11] )
      
      ea_28.append(     f2[i][6] )
      ea_32.append(     f2[i][7] )
      ea_36.append(     f2[i][8] )
      ea_40.append(     f2[i][9] )
      ea_45.append(     f2[i][10])
      ea_50.append(     f2[i][11] )
      ea_60.append(     f2[i][12] )
      
      anchor_det = np.array([anker_det_x[-1],anker_det_y[-1]],dtype=float)
      
      pos = _util_img_position_ea(angle[-1],anchor_det,1700.) ; x,y = list(pos)
      xy_17x.append(x) ; xy_17_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,1900.) ; x,y = list(pos)
      xy_19x.append(x) ; xy_19_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,2100.) ; x,y = list(pos)
      xy_21x.append(x) ; xy_21_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,2300.) ; x,y = list(pos)
      xy_23x.append(x) ; xy_23_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,2500.) ; x,y = list(pos)
      xy_25x.append(x) ; xy_25_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,2700.) ; x,y = list(pos)
      xy_27x.append(x) ; xy_27_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,2800.) ; x,y = list(pos)
      xy_28x.append(x) ; xy_28_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,3200.) ; x,y = list(pos)
      xy_32x.append(x) ; xy_32_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,3600.) ; x,y = list(pos)
      xy_36x.append(x) ; xy_36_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,4000.) ; x,y = list(pos)
      xy_40x.append(x) ; xy_40_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,4500.) ; x,y = list(pos)
      xy_45x.append(x) ; xy_45_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,5000.) ; x,y = list(pos)
      xy_50x.append(x) ; xy_50_y.append(y)      
      pos = _util_img_position_ea(angle[-1],anchor_det,6000.) ; x,y = list(pos)
      xy_60x.append(x) ; xy_60_y.append(y)      

   wave = 100.*np.array([17,19,21,23,25,27,28,32,36,40,45,50,60],dtype=float)
   eff_area = np.zeros(len(wave)*Nrow).reshape(len(wave),Nrow)
   positions = np.zeros(len(wave)*Nrow*2).reshape(len(wave),2,Nrow)
   
   for i in range(Nrow):
      eff_area[:,i] = [ea_17[i],ea_19[i],ea_21[i],ea_23[i],ea_25[i],ea_27[i],\
              ea_28[i],ea_32[i],ea_36[i],ea_40[i],ea_45[i],ea_50[i],ea_60[i]]
      positions[:,0,i] = [xy_17x[i],xy_19x[i],xy_21x[i],xy_23x[i],xy_25x[i],\
                          xy_27x[i],xy_28x[i],xy_32x[i],xy_36x[i],xy_40x[i],\
                          xy_45x[i],xy_50x[i],xy_60x[i]]              
      positions[:,1,i] = [xy_17_y[i],xy_19_y[i],xy_21_y[i],xy_23_y[i],xy_25_y[i],\
                          xy_27_y[i],xy_28_y[i],xy_32_y[i],xy_36_y[i],xy_40_y[i],\
                          xy_45_y[i],xy_50_y[i],xy_60_y[i]]           
   return wave, eff_area, positions, cross_ref_no, zone, anker_det_x, anker_det_y, exposure, angle      


def getCal(wave,clockposition='200',all='False'):
      '''read the effective area ARF calibration file from the CALDB
      
      SUPERSEDED by routine uvotio.specResp()
      '''
      
      import os
      
      CALDB = os.getenv('CALDB')
      
      dir = CALDB+'/data/swift/uvota/cpf/arf/'

      if clockposition == '200':#
         file = dir+'/swugu0200_20041120v101.arf'
      if clockposition == '160': 
         file = dir+'/swugu0160_20041120v101.arf'
      if clockposition == '1000':
         file = dir+'/swugu1000_20041120v101.arf'
      if clockposition == '955': 
         file = dir+'/swugu0955_20041120v101.arf'
      f = pyfits.open(file)
      bintable = f[1]
      coldefs = bintable.get_coldefs()
      # enrergy in KeV, wave in A, SpecResp in cm2
      ENERGY_LO = bintable.data.field(0)
      ENERGY_HI = bintable.data.field(1)
      WAVE_MIN = bintable.data.field(2)
      WAVE_MAX = bintable.data.field(3)
      SPECRESP = bintable.data.field(4)
      f.close()
      wav = (WAVE_MIN + WAVE_MAX)*0.5
      allvals = None
      if all == 'True': allvals = ENERGY_LO,ENERGY_HI,WAVE_MIN,WAVE_MAX,SPECRESP
      sresp = interpol(wave, wav,SPECRESP)
      return sresp, allvals

# zeroth order calibration ###############################################################

def zo_ea(calfile,extimg,hdr,figno=11,makecalplot=False, mode=1, wheelpos=160):
   ''' rough estimate of zeroth order effective area, for a single extimg. 
   
   Parameters
   ==========
   calfile : path
      calibration ascii file (uvotified)
   extimg : numpy array
      the extracted image (third extension of the spectrum data file)
   hdr : fits header
      the header containing the exposure time (first extension fits file)
   figno : int
      number of first figure 
   mode : int   
      0 for background subtracted     
   
   Notes
   =====
   This analysis uses the fit to the dispersion 
   from zemax: p = 48-110715/(wave-1220) 

   where wave is in A, p wavelength bins.
    
   p = 0 for wave =  (4329 peak of the B-band).
   
   Typical offset needed = -8
   
   Method: the images is displayed. Zoom in. 
   start selection of points all along the zeroth order from left to right
   wait 
   while you wait the background is computed using the uvotgetspec background routine
      and the maximum netCR is computed
      wait till it asks for the offset from max(CRnet) and prints the maximum location,
      for example, 202
   Then enter a list of offsets (up to 4), like 200,202,205,210
   The results are plotted and written to an ascii files (one per offset).
   rename the files as you see fit.    
    '''
   from uvotpy import uvotmisc 
   import numpy as np
   import uvotgrism
   from stsci.convolve import boxcar
   from scipy import interpolate
   import pylab as plt
   import uvotphot
   from scipy.interpolate import splrep,splev
   
   exposure = hdr['exposure']

   h_planck = 6.626e-27  # erg - s
   lightspeed = 2.9979e10  # cm/sec
   h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
   
   # dispersion - use simple solution for now 
   # from wave to pix  
   #xwave = np.arange(1650,6800,1)
   #xp = 48.-110715/(xwave-1220)
   #del_lambda = 1./110715*(xwave-1220)**2   # delta wavelengths per 1 pixel = (dlam/dpix) term
   # from pix p to wave
   if wheelpos == 955:
      p = np.arange(-90,25)
      # the dispersion was derived using the inverse of the zemax w-data and the pix positions for this
      # position p=0 corresponds to 4265.15A
      C = array([ -1.38243428e-08,  -3.27891071e-06,   2.34458158e-04])
      w  = 1./ polyval(C,p)   
      w1 = 1./ polyval(C,p-0.5)
      w2 = 1./ polyval(C,p+0.5)
      crpix_cal = np.zeros(len(p))
      p0 = p[78]
      w0 = (w1[78]+w2[78])/2.0
      print("at reference pixel the wavelength boundaries are :(%6.1f,%6.1f)"%(w1[78],w2[78]))
      sf = 100
   elif wheelpos == 160:   
      p = np.arange(-255,18,1)
      w  = 1220. - 110715./(p-35.6)   
      w1 = 1220. - 110715./(p-36.1)
      w2 = 1220. - 110715./(p-35.1)
      crpix_cal = np.zeros(len(p))
      p0 = p[255]  # was at point 243!
      w0 = (w1[255]+w2[255])/2.0
      print("at reference pixel the wavelength boundaries are :(%6.1f,%6.1f)"%(w1[255],w2[255]))
      sf = 60
   
   # get calibration spectrum as count rate(per sec) per cm2 per angstrom for each pixel
   cal = uvotmisc.rdTab(calfile)   
   wcal = cal[:,0]
   crcal = boxcar(cal[:,1]/(h_c_ang/wcal), (3,))
   # linear interpolating spline which can be used for integrating
   tck = interpolate.splrep(wcal,crcal,k=3 , s=1.5e-6)
   for i in range(len(p)):
      crpix_cal[i] = interpolate.splint(w1[i],w2[i],tck,)
   
   # plot extimg
   plt.figure(figno)
   if makecalplot:
      pp = p-p0
      q = np.isfinite(w1) & np.isfinite(w2)
      pp = pp[q]
      w1 = w1[q]
      w2 = w2[q]
      cr = 0.5*(splev(w1,tck,) + splev(w2,tck,))*(w2-w1)
      x = 10./(pp+300)
      ii = list(range(len(x)-1,-1,-1))
      tck = splrep(x[ii],cr[ii],s=1e-5 )
      #plt.plot(pp, cr,'.',markersize=2 )
      plt.plot(pp[ii],splev(x,tck)[ii],'k',lw=2,label="reference spectrum")
      plt.xlabel("distance in spectral bins to reference at $\lambda$%i$\AA$"%(int(w0)))
      plt.ylabel("count rate per spectral bin")
      plt.text(-235,0.031,"reference spectrum WD1657+343",fontsize=12)
      plt.ylim(0,0.0465)
      plt.xlim(-245,40)
      return pp,w1,w2,cr, wcal,crcal 
   bgimg = findBackground(extimg)[4]
   if mode == 0:
      plt.imshow(extimg-bgimg,vmin=0,vmax=30)
   else:
      plt.imshow(extimg,vmin=0,vmax=30)
   plt.contour(extimg,levels=[30,40,50,70,90,130,200,350])
   
   print("pick the zeroth order locations from begin to end")
   y = input("ready ?")
   xy = plt.ginput(n=30,timeout=40)
   x = np.array(xy)[:,0]
   y = np.array(xy)[:,1]
   
   plt.plot(x,y,'+k')
   tck2 = interpolate.splrep(x,y,)
   
   # extract counts spectrum, background
   x1 = int(np.min(x)+0.5)
   x2 = int(np.max(x)+0.5)
   # redefine x,y 
   x = np.arange(x1,x2+1)
   y = interpolate.splev(x,tck2)
   iy = np.array( (y+0.5), dtype=int )
   crobs = np.zeros(len(x))
   bgobs = np.zeros(len(x))
   for i in range(len(x)):
      crobs[i] = extimg[iy[i]-10:iy[i]+10,x1+i].sum()
      bgobs[i] = 0.5*(bgimg[iy[i]-40:iy[i]-20,x1+i].sum()+bgimg[iy[i]+20:iy[i]+40,x1+i].sum())
      
   crobs = crobs/exposure
   bgobs = bgobs/exposure
   
   # compute coincidence correction factor "coif" of background based on 314.26 pixels 
   frametime=0.0110322
   deadc = 0.984
   bkg_rawcpf = bgimg[:,x1:x2].mean()/deadc * 314.26 * frametime/exposure
   bkg_coicpf = -np.log(1.-bkg_rawcpf)
   coif = bkg_coicpf/bkg_rawcpf
   print("coi factor adopted = ",coif)
   
   # call clip routine from uvotgrism 
   bgmask = uvotgrism.clipmask(bgobs,sigclip=2.5,fpos=False)
   bgmean = np.mean(bgobs[bgmask])
   bgobs[bgmask == False] = bgmean
   # boxcar smooth bgobs
   bg = boxcar( bgobs, (30,))
   crnet = (crobs-bg)*(2*coif) # approximate net coif is twice that in background for crnet << bg
   crnet[crnet <= 0.] = 0.01
   print("background counts per frame (20) and (315 subpixels) = ",bg.mean()*0.0110302*np.array([1.,315./20.]))

   wh_ea = uvotphot.rd_uvot_effearea('wh')

   p1 = np.arange(len(crobs))
   k = np.where(crnet == np.max(crnet))[0][0]
   ans = input("give a list of values, like 2,3,4 for the zeropoint (@lambda_4329). Maximum CRnet is at %3i >> "%(k))
   print(" new values k = ",ans,"  trying conversion... ")
   try:
      klist = np.array(ans.split(','),dtype=int)
      print(" new values k = ",klist)
   except: 
      print("expected  a list 200,210,205 etc.")
      klist = [k, k+5, k+10]
      print("continuing but with : ",klist)  
      pass 
   k_p = klist
   if (type(k_p) != list) | len(k_p) == 0:
      k_p = [int(k),] 
      print(" reset values k = ",k_p)
      
   fig2 = plt.figure()
   ax1 = fig2.add_axes([0.12,0.66,0.8,0.28])
   ax2 = fig2.add_axes([0.12,0.10,0.8,0.48])
   
   kleur = ['k','m','b','darkgreen','orange']
   ii = -1
   p1_ = p1
   crnet_ = crnet
   k_ = k
   for k in k_p:
       print("now computing ",k)
       ii = ii+1
       print(kleur[ii])
       p1 = p1_ - k
   
       k1 = p1[0]
       k2 = p1[-1]
       pp = (p >= k1) & (p < k2)
       print("k1,k2,len(pp)=",k1,k2,len(np.where(pp)[0]))
       if len(p1) > len(np.where(pp)[0]): 
          p1 = p1[0:len(np.where(pp)[0])]
          crnet = crnet_[0:len(np.where(pp)[0])]
       else:
          crnet = crnet_          
   
       # effective area = cr_obs/cr_cal
       eff_area = crnet/crpix_cal[pp]
       ea_wave = 1220. - 110715./(p1-35.6) 
       q = np.isfinite(eff_area)
       tck = splrep(ea_wave,eff_area,s=sf)       
   
       fout = open("zeroth_order_eff_area_k"+str(k_-k)+".txt","w")  
       for i in range(len(p1)): 
           fout.write( "%2i %3i %6.1f %7.4f %8.4f  %8.4f \n"%(i, 
              p1[i],ea_wave[i],eff_area[i],crnet[p1][i],crpix_cal[pp][i] ))
       fout.close()

       ax1.plot(p1,eff_area,'-',lw=1.5,color=kleur[ii],label="p=%i"%(k_-k))
       ax1.set_xlabel("distance p-p' (bins)")     
       #ax2.plot(ea_wave , eff_area,'.',markersize=2,color=kleur[ii])
       ax2.plot(ea_wave , splev(ea_wave,tck),ls='-',color=kleur[ii],lw=1.5,alpha=0.7,label="zeroth order EA p=%i"%(k_-k))
       ax2.set_xlabel('$\lambda (\AA)$',fontsize=14)
   ax2.plot(0.5*(np.array(wh_ea[0])+np.array(wh_ea[1])), wh_ea[2], 
            ls='--',lw=2,color='darkorange',label="white Effective Area")
   ax2.legend(loc=0)
   ax1.legend(loc=0)
   ax2.set_ylabel("Effective Area (cm$^2$)")
   ax1.set_ylabel("Effective Area (cm$^2$)")
   ax2.set_xlim(1600,7500)
   ax1.set_xlim(-200,40)
   
   # TBD adjust parameter k to match up cal spectrum trying to get p0 line up to k



   
# ======== flux calibration zemax model 

def makeNewFluxCalFile(wheelpos=None,spectralorder=1,outfile="newFluxCalFile.fits",
    Xank=1129, Yank=1022, chatter=0):
   """ 
   Combines the individual effective area, zemax flux model, and old cal file
   
   Find calibration files in local directory for incorporation.
   
   2013-05-05 NPM Kuin
   """
   import os
   #import pyfits
   
   CALDB = os.getenv('CALDB')
   if CALDB == '': 
      raise ioerror('please define environment variable CALDB before proceeding')
   UVOTPY = os.getenv('UVOTPY')
   if UVOTPY == '': 
      raise ioerror('please define environment variable UVOTPY before proceeding')
   
   if wheelpos == 160:
      oldcalfile = CALDB + '/data/swift/uvota/cpf/arf/swugu0160_20041120v101.arf'
      grism = 'UGRISM'
   if wheelpos == 200:
      oldcalfile = CALDB + '/data/swift/uvota/cpf/arf/swugu0200_20041120v101.arf'
      grism = 'UGRISM'
   if wheelpos == 955:
      oldcalfile = CALDB + '/data/swift/uvota/cpf/arf/swugv0955_20041120v101.arf'
      grism = 'VGRISM'
   if wheelpos == 1000:   
      oldcalfile = CALDB + '/data/swift/uvota/cpf/arf/swugv1000_20041120v101.arf'
      grism = 'VGRISM'
   calfiles = UVOTPY+"/calfiles/"
   # use local directory to seach for the calibration files
   status = os.system("ls -1  > tmp.1")
   if status != 0:
      raise runtimeerror( "FAIL: ls -1 "+calfiles+" > tmp.1")
   f = open("tmp.1")
   clist = f.readlines()
   f.close()      
   #status = os.system("rm tmp.1") 
   
   if len(clist) == 0: 
      raise runtimeerror( "WARNING XYSpecResp: calfiles directory seems empty" )
         
   status = os.system("cp -f "+oldcalfile+" "+outfile)   
   if status != 0:
      raise runtimeerror( "FAILED to copy oldcalfile")
         
   for cl in clist:
       cl1 = cl
       cl = cl.split("\n")[0] 
       print('processing: ',cl)

       try:       
          if int(cl[5:9]) == wheelpos:
             cl = cl.split("_")
             if (cl[3] == 'o'+str(spectralorder)):
                sporder = cl[3]
                anchor = [int(cl[1].split('ay')[0].split('ax')[1]),int(cl[1].split('ay')[1]) ]
                dxy = cl[2]
                ver = cl[4].split('.')[0].split('v') # [date,version]
                dist = (Xank-anchor[0])**2 + (Yank-anchor[1])**2
                print(sporder," - ",anchor," - ",dxy," - ",ver," - ",dist)
                extname = "SPECRESP"+grism+"%04d"%wheelpos+"_"+cl[1]+"_"+cl[3]
                extname = extname.upper()

                command = "fparkey '"+extname+"'  '"+cl1.split("\n")[0]+"[1]'  extname"
                print(command); status = os.system(command)
                if status != 0:
                   raise runtimeerror( "FAILED to modify extname in header of "+
                         calfiles+cl1.split("\n")[0]+"+1 to "+extname)

                command = "fchecksum "+cl1.split("\n")[0]+"  update+  datasum"
                print(command); status = os.system(command)
                if status != 0:
                   raise runtimeerror( "FAILED to update checksum "
                         +calfiles+cl1.split("\n")[0])

                command = "fappend "+cl1.split("\n")[0]+"+1 "+outfile
                print(command); status = os.system(command)
                if status != 0:
                   raise runtimeerror( "FAILED to append "
                         +calfiles+cl1+"+1 to "+outfile)

                if chatter > -1: 
                   print("adding \norder=%i\nanchor = %s\ndxy = %s\nversion = %s\ndistance = %7.2f\n" \
                          %(spectralorder,anchor, dxy, ver, dist,))

       except:
          print("problem at cl=",cl1)
          pass          
   
   # add model   
   modfile = '/calibration/grism/fluxcal2/simplified_zmxmodel_'+str(wheelpos)+'.fits'
   wr_simplified_zmxmodel_fits(wheelpos,clobber=True)
   status = os.system("fappend "+modfile+"+1 "+outfile)       
   if status != 0:
      raise RuntimeError( "FAILED to apend "+modfile+"+1 to "+outfile)
   # done   
   print("makeNewFluxCalFile done")


def anker2field(anker,modelfile='/Volumes/data/grism/fluxcal2/simplified_zmxmodel_160.fits'):
   from scipy import interpolate
   import numpy as np
   wave, xf, yf, xpix, ypix, flux, valid = rd_simplified_zmxmodel_fits(file=modelfile)
   q = np.isfinite(xpix[3,:]) & np.isfinite(ypix[3,:])
   tckx = interpolate.bisplrep(xpix[3,q],ypix[3,q],xf[q])
   tcky = interpolate.bisplrep(xpix[3,q],ypix[3,q],yf[q])
   return interpolate.bisplev(anker[0],anker[1],tckx),interpolate.bisplev(anker[0],anker[1],tcky)

def wr_simplified_zmxmodel_fits(wheelpos,spectralorder = 1, addfields=False, option=0, clobber=True):
   '''Write first order scaled simplified zemax model to a file 
   
   '''
   import pyfits
   import datetime
   import numpy as np
   from uvotpy import uvotio
   
   if addfields:
      outfile = '/calibration/grism/fluxcal2/simplified_zmxmodel_fld_'+str(wheelpos)+'.fits'
   else:
      outfile = '/calibration/grism/fluxcal2/simplified_zmxmodel_'+str(wheelpos)+'.fits'
      
   filtername = str(wheelpos)
   zmxwav,zmxorder, xfzmx, yfzmx, xpix, ypix, zmxflux, zmxvalid = zemax_simply_corrected(wheelpos,option=option)
   zmxwav = np.array(1000*zmxwav,dtype=int)
   nw = 16
   if wheelpos > 500: nw = 12
   wave = np.zeros(nw*784,dtype=int).reshape(nw,784)
   
   indices=np.zeros(nw*784,dtype=int).reshape(nw,28,28)
   for i in range(nw):
      indices[i,:,:] = i
   for i in range(28):
      indices[:,i,:] = i
      for j in range(28):
         indices[:,i,j] = j
            
   for i in range(784): wave[:,i] = zmxwav
   anow = datetime.date.today()
   datestring = anow.isoformat()[0:4]+"-"+anow.isoformat()[5:7]+"-"+anow.isoformat()[8:10]
   if wheelpos < 400: 
      filtername = 'UGRISM'
   else:
      filtername = 'VGRISM'   
   energy_lo = uvotio.angstrom2kev(7000.)
   energy_hi = uvotio.angstrom2kev(1900.)
   hdu = pyfits.PrimaryHDU()
   hdulist=pyfits.HDUList([hdu])
   hdulist[0].header.update('TELESCOP','SWIFT   ','Telescope (mission) name')                       
   hdulist[0].header.update('INSTRUME','UVOTA   ','Instrument Name')                                
   hdulist[0].header.update('COMMENT','SCALED ZEMAX NORMALIZED FLUX MODEL')
   hdulist[0].header.update('DATE',datestring,'creation date')                                
   col11 = pyfits.Column(name='WAVE',format='I',array=wave.flatten(),unit='A')
   col12 = pyfits.Column(name='XPIX',format='E',array=xpix[1,:,:].flatten(),unit='pix')
   col13 = pyfits.Column(name='YPIX',format='E',array=ypix[1,:,:].flatten(),unit='pix')
   col14 = pyfits.Column(name='FLUX',format='E',array=zmxflux[1,:,:].flatten(),)
   if addfields:
     col15 = pyfits.Column(name='XFIELD',format='E',array=xfzmx,unit='mm' )
     col16 = pyfits.Column(name='YFIELD',format='E',array=yfzmx,unit='mm' )
   col17 = pyfits.Column(name='VALID',format='B',array=zmxvalid[1,:,:].flatten(), )
   if addfields:
     cols1 = pyfits.ColDefs([col11,col12,col13,col14,col15,col16,col17,])
   else: cols1 = pyfits.ColDefs([col11,col12,col13,col14,col17,])
   tbhdu1 = pyfits.new_table(cols1)
   tbhdu1.header.update('EXTNAME','ZEMAXMODEL_'+str(wheelpos),'Name of this binary table extension')
   tbhdu1.header.update('TELESCOP','Swift','Telescope (mission) name')
   tbhdu1.header.update('INSTRUME','UVOTA','Instrument name')
   tbhdu1.header.update('COMMENT','POSITIONS HAVE BEEN APPROXIMATED - ')
   tbhdu1.header.update('COMMENT','DATA ONLY TO BE USED FOR FLUX CALIBRATION')
   tbhdu1.header.update('FILTER',filtername)
   tbhdu1.header.update('WHEELPOS',wheelpos)
   tbhdu1.header.update('ORIGIN','MSSL/UCL','source of FITS file')
   tbhdu1.header.update('CREATOR','uvotcal.py','MSSL uvotpy python library')
   tbhdu1.header.update('VERSION',1.0)
   tbhdu1.header.update('FILENAME',outfile)
   tbhdu1.header.update('SP_ORDER','ORDER('+str(1)+')','spectral order')  
   tbhdu1.header.update('DATE',datestring,'creation date')                                

   tbhdu1.header.update('HDUCLASS','OGIP','format conforms to OGIP standard')
   tbhdu1.header.update('HDUCLAS1','RESPONSE','MODEL RESPONSE DATA')
   tbhdu1.header.update('HDUCLAS2','SPECRESP','type of calibration data')   
   tbhdu1.header.update('CCLS0001','CPF','dataset is a calibration product file')
   tbhdu1.header.update('CCNM0001','SPECRESP','Type of calibration data')
   tbhdu1.header.update('CDES0001',filtername+' NORMALISED SPECTRAL RESPONSE MODEL','Description')
   tbhdu1.header.update('CDTP0001','DATA','Calibration file contains data')
   tbhdu1.header.update('CVSD0001','2004-11-20','UTC date when calibration should first be used')
   tbhdu1.header.update('CVST0001','00:00:00','UTC time when calibration should first be used')
   tbhdu1.header.update('CBD10001','FILTER('+filtername+')','Parameter boundary')
   tbhdu1.header.update('CBD20001','ENERG('+str(energy_lo)+'-'+str(energy_hi)+')keV','spectral range')
   tbhdu1.header.update('CBD30001','RADIUS(0-10)pixel','Parameter boundary')
   tbhdu1.header.update('CBD50001','WHEELPOS('+str(wheelpos)+')','Filter/Mode Selection')
   tbhdu1.header.update('CBD60001','ORDER('+str(spectralorder)+')','spectral order')  
   hdulist.append(tbhdu1)
   hdulist.writeto(outfile,clobber=clobber)
 
def rd_simplified_zmxmodel_fits(file='/calibration/grism/fluxcal2/simplified_zmxmodel_160.fits',wheelpos=160):
   import pyfits
   import numpy as np
   if wheelpos == 200:
      file = '/calibration/grism/fluxcal2/simplified_zmxmodel_'+str(wheelpos)+'.fits'
   xf, yf = None, None
   hdu = pyfits.open(file)
   t = hdu[1].data
   wave = t['WAVE']
   xpix = t['xpix'].reshape(16,784)
   ypix = t['ypix'].reshape(16,784)
   flux = t['flux'].reshape(16,784)
   try:
     xf = t['xfield']
     yf = t['yfield']
   except: pass  
   valid = t['valid'].reshape(16,784)
   return wave, xf, yf, xpix, ypix, flux, valid      

def zemax_simply_corrected(wheelpos,option=0, plot=False):
    '''Quick and dirty simple model refit.
     
       first prepare align_zemax and correctAnkPos,  then run
       
       option = 0 was used for initial flux calibration Feb-Mar 2013
       option = 1 ? 
       
       2012-01-17 NPMK
    '''
    from scipy import interpolate
    from . import zemax
    from uvotpy import uvotgetspec
    import numpy as np
    
    ## start with the simplified zemax model << not sure why I put these lines here (? during wavecal 2013-04-28)
    #
    #wr_simplified_zmxmodel_fits(wheelpos,spectralorder = 1, addfields=True, clobber=True)
    #
    #simzmx = '/Volumes/data/grism/fluxcal2/simplified_zmxmodel_fld_'+str(wheelpos)+'.fits'
    #
    #wave, xf, yf, xpix, ypix, flux, valid = rd_simplified_zmxmodel_fits(file=simzmx,wheelpos=wheelpos)
    #
    if wheelpos == 160:
       # set align_zemax to (-60,0)  
       # set correctAnkPos to scale relative to the boresight by 0.973
       # No further distortion correction made
       wave,orders,xf,yf,xc,yc,xp,yp,flx,vori = zemax.read_zemax_flux(wheelpos=160)
       xf2=xf[1,3,:]
       yf2=yf[1,3,:]
       xp2, yp2 = zemax._align_zemax(xp,yp,160, test="simple")
       xp2, yp2 = zemax.correctAnkPos(xp2.flatten(), yp2.flatten(), 160, test="simple")  # note: everything scaled by same factor
       xp2 = xp2.reshape(5,16,784)
       yp2 = yp2.reshape(5,16,784)
       xp2[np.isnan(xp2)] = -450
       yp2[np.isnan(yp2)] = -450
       flx[np.isnan(flx)] = 0.
       # xp2,yp2,flx is the new model on the detector
       # find the translation for the field coordinates
       q = np.where(xp2[1,3,:] > 310)
       bx,by = uvotgetspec.boresight('uc160')
       tck_x = interpolate.bisplrep(xp2[1,3,:].flatten()[q],yp2[1,3,:].flatten()[q],xf2.flatten()[q],)
       tck_y = interpolate.bisplrep(xp2[1,3,:].flatten()[q],yp2[1,3,:].flatten()[q],yf2.flatten()[q],)
       bfx = interpolate.bisplev(bx,by,tck_x)
       bfy = interpolate.bisplev(bx,by,tck_y)
       # found -0.02267, -0.0020953 deg [which comes to (*6550.4) -148.5, -13.7 pix]
       xf2 = xf2 - bfx
       yf2 = yf2 - bfy
       if plot:
         try:
           import pylab as plt
           plt.plot(xp2[1,3,:],yp2[1,3,:],'.b')
           plt.plot([597],[1292],'s',color='c',markersize=16)
           plt.plot([733],[1261],'s',color='orange',markersize=16)
           plt.plot([880],[1540],'s',color='orange',markersize=16)
           for i in range(784): plt.text(xp2[1,3,i],yp2[1,3,i],"%6.2f"%(flx[1,6,i]),fontsize=9,horizontalalignment='center')
           plt.title("simple model fit with flux level at 350nm")
         except: pass
       # note the flux is normalised to field position bfx,bfy
       return wave,orders, xf2, yf2, xp2, yp2, flx, vori 

    elif (wheelpos == 200) | (wheelpos == 1000):
       # for this model, use the distortion correction and alignment from the wavecal
       # **for all wavelength points!** ** This affects the dispersion, so cannot be used **
       # ** for the wavecal **
       wave,orders,xf,yf,xc,yc,xp,yp,flx,vori = zemax.read_zemax_flux(wheelpos=wheelpos)
       k = 3 
       if wheelpos == 1000: k = 5
       nw = len(wave)
       xf2=xf[1,k,:]
       yf2=yf[1,k,:]
       xp2, yp2 = zemax._align_zemax(xp,yp,wheelpos, ) #test="simple") would imply about (-140,0) ??? 
       xp2, yp2 = zemax.correctAnkPos(xp2.flatten(), yp2.flatten(),wheelpos=wheelpos, )  # note: old scale factor
       xp2 = xp2.reshape(5,nw,784)
       yp2 = yp2.reshape(5,nw,784)
       xp2[np.isnan(xp2)] = -450
       yp2[np.isnan(yp2)] = -450
       flx[np.isnan(flx)] = 0.
       # xp2,yp2,flx is the new model on the detector
       # find the translation for the field coordinates
       q = np.where(xp2[1,k,:] > -50)
       if wheelpos == 200:
          bx,by = uvotgetspec.boresight('ug200')
       elif wheelpos == 1000:
          bx,by = uvotgetspec.boresight('vg1000')                 
       tck_x = interpolate.bisplrep(xp2[1,k,:].flatten()[q],yp2[1,k,:].flatten()[q],xf2.flatten()[q],)
       tck_y = interpolate.bisplrep(xp2[1,k,:].flatten()[q],yp2[1,k,:].flatten()[q],yf2.flatten()[q],)
       bfx = interpolate.bisplev(bx,by,tck_x)
       bfy = interpolate.bisplev(bx,by,tck_y)
       # found -0.02267, -0.0020953 deg [which comes to (*6550.4) -148.5, -13.7 pix]
       xf2 = xf2 - bfx
       yf2 = yf2 - bfy
       if plot:
         try:
           import pylab as plt
           plt.plot(xp2[1,k,:],yp2[1,k,:],'.b')
           #plt.plot([597],[1292],'s',color='c',markersize=16)
           #plt.plot([733],[1261],'s',color='orange',markersize=16)
           #plt.plot([880],[1540],'s',color='orange',markersize=16)
           for i in range(784): plt.text(xp2[1,6,i],yp2[1,6,i],"%6.2f"%(flx[1,6,i]),fontsize=9,horizontalalignment='center')
           if wheelpos == 200: plt.title("simple model fit with flux level at 350nm")
           if wheelpos == 1000: plt.title("simple model fit with flux level at 460nm")
         except: pass
       # note the flux is normalised to field position bfx,bfy
       return wave,orders, xf2, yf2, xp2, yp2, flx, vori 
    
    elif wheelpos == 955:
       # for this model, use the distortion correction and alignment from the wavecal
       # **for all wavelength points!** ** This affects the dispersion, so cannot be used **
       # ** for the wavecal **
       wave,orders,xf,yf,xc,yc,xp,yp,flx,vori = zemax.read_zemax_flux(wheelpos=955)
       k = 5
       nw = len(wave)
       xf2=xf[1,k,:]
       yf2=yf[1,k,:]
       if option == 0:
          xp2, yp2 = zemax._align_zemax(xp,yp,wheelpos, ) #test="simple") would imply about (-140,0) ??? 
          xp2, yp2 = zemax.correctAnkPos(xp2.flatten(), yp2.flatten(), wheelpos=wheelpos, )  # note: old scale factor
       elif option == 1:
          xp2, yp2 = zemax._align_zemax(xp,yp,wheelpos, test="simple")  # would imply about (-140,0) ??? 
          xp2, yp2 = zemax.correctAnkPos(xp2.flatten(), yp2.flatten(), wheelpos=wheelpos, test="simple" )  # note: old scale factor     
       else:
          raise RuntimeError ("option parameter value in error")          
       xp2 = xp2.reshape(5,nw,784)
       yp2 = yp2.reshape(5,nw,784)
       xp2[np.isnan(xp2)] = -450
       yp2[np.isnan(yp2)] = -450
       flx[np.isnan(flx)] = 0.
       # xp2,yp2,flx is the new model on the detector
       # find the translation for the field coordinates
       q = np.where(xp2[1,k,:] > -50)
       bx,by = uvotgetspec.boresight('vc955')
       tck_x = interpolate.bisplrep(xp2[1,k,:].flatten()[q],yp2[1,k,:].flatten()[q],xf2.flatten()[q],)
       tck_y = interpolate.bisplrep(xp2[1,k,:].flatten()[q],yp2[1,k,:].flatten()[q],yf2.flatten()[q],)
       bfx = interpolate.bisplev(bx,by,tck_x)
       bfy = interpolate.bisplev(bx,by,tck_y)
       # found -0.02267, -0.0020953 deg [which comes to (*6550.4) -148.5, -13.7 pix]
       xf2 = xf2 - bfx
       yf2 = yf2 - bfy
       if plot:
         try:
           import pylab as plt
           plt.plot(xp2[1,k,:],yp2[1,k,:],'.b')
           #plt.plot([597],[1292],'s',color='c',markersize=16)
           #plt.plot([733],[1261],'s',color='orange',markersize=16)
           #plt.plot([880],[1540],'s',color='orange',markersize=16)
           for i in range(784): plt.text(xp2[1,6,i],yp2[1,6,i],"%6.2f"%(flx[1,6,i]),fontsize=9,horizontalalignment='center')
           plt.title("simple model fit with flux level at 460nm")
         except: pass
       # note the flux is normalised to field position bfx,bfy
       return wave,orders, xf2, yf2, xp2, yp2, flx, vori     
    else:
       raise IOError("zemax_simply_corrected: this wheelpos not supported")   
       
def _get_calibrated_zemax_flux(wheelpos, correct_field=False, chatter=0):
   ''' BROKEN 
   The anchor position adjustments to the zemax model are implemented 
       in zemax.correctAnkPos(xank,yank,wheelpos=160).  
   The wavelength calibration scaling is done in uvotcal.getZmx(xfield,yfield,160,test="flux")
   Get the curvature polynomial y(x) from uvotgrism.spec_curvature(160,anchor_in_img,)

   Method
   ------
   - get unscaled zemax flux model
   - determine correct anchor positions for the field points (align_zemax, field position ofset)
   - for each field point, get the dispersion, spectrum slope/angle, 
   - get the curvature for each of those
   - calculate the position in detector coordinates of the zemax model points 
         corrected for scale and curvature in Zxnew, Zynew, which are the 
         updated xpix, ypix
   
   '''
   import zemax,uvotgrism
   from uvotpy import uvotmisc
   import numpy as np
   
   if wheelpos < 500: refwav = 2600. 
   else: refwav = 4200.
   
   Zwave, Zorders, Zxfield, Zyfield, Zxcentroid, Zycentroid, Zxpix, Zypix, flux, \
   valid = zemax.read_zemax_flux(wheelpos=wheelpos,file=None, rotate=False,chatter=chatter, )
   # Note that the anchor cannot be gotten from Zxpix,Zypix without shifting the input field positions.
   # xpix, ypix has not been aligned either    

   xf = Zxfield[1,3,:] ; yf = Zyfield[1,3,:] 
   
   #specslope = np.zeros(28*28,dtype=float)
   #val  = np.ones(5*16*28*28,dtype=bool).reshape(5,16,28*28)
   
   for i in range(28*28):
      if chatter > 0: print("processing field point # ",i)
      if valid[1,3,i]:
         # get the zemax model for the selected field point
         C_zero, C_1, C_2, C_3, C_min1, xpix, ypix, zmxdis, zmxwav, wave, theta = getZmx( \
             xf[i], yf[i], wheelpos, mode='spline',chatter=chatter, \
             test="flux", kk=3,s=50,wpixscale=None,xlimits=[15,2100])
         # Here the zemax model has been corrected by aligning the zemax model boresight to
         # the observed boresight. The uncorrected anchor position can be found now. 
         # xpix and ypix are the unscaled(wpixscale) zemax positions adjusted for the 
         #   boresight, but with shift from align_zemax applied.
         # zmxdis,zmxwav are the dispersion and wavelengths after adjustment of wpixscale was applied.
         
         
         # ancher position adjustement (biquadratic fit between zemax and the observed positions):
         xpin = xpix[1,3]
         ypin = ypix[1,3]
         anchor = np.array(zemax.correctAnkPos(xpin,ypin,wheelpos=wheelpos,test='flux',chatter=0),dtype=float)

         # get the curvature coefficients of the first order (apply with dis)
         curvecoef = uvotgrism.spec_curvature(wheelpos,anchor,order=1,)
         
         # use pair (wav, dis) for valid points 
         dis = zmxdis[1]
         wav = zmxwav[1]
         if chatter > 1: print("check 1 file point #:",i,"  wav:",wav,"    dis:",dis)
         if (wav == None): 
            valid[1,:,i] = False 
            print('skipping field point ',i)
            continue
         if (len(wav) < 3): 
            valid[1,:,i] = False 
            print('skipping field point ',i)
            continue
         if not (zmxwav[1] == 2600.).any() :
         # some valid points but anchor not included
         # determine where the anchor should be and add that to dis and wav
            try:
               inv_disp = np.polyfit(wav,dis,2)
               xdis = np.polyval(inv_disp, refwav)
               k = xwav.searchsorted( refwav)
               dis = list(dis)
               dis.insert(k, xdis)
               dis = np.array(dis)
               wav = list(wav)
               wav.insert(k, refwav)
               wav = np.array(wav)
            except:
            # give up
               valid[1,:,i] = False 
               print('skipping field point ',i)
               continue
         # adjust yoff normal to the dis array, then rotate 
         if chatter > 1: print("check 2 file point #:",i,"  wav:",wav,"    dis:",dis)
         yoff = np.polyval(curvecoef, dis)
         xr,yr = uvotmisc.uvotrotvec(dis,yoff,-(180.-theta[1])) 
         # put the data into the correct positions in the result array & add anchor offset
         xnew = np.zeros(len(Zwave)) -999
         ynew = np.zeros(len(Zwave)) -999
         for j in range(len(Zwave)):
            w = Zwave[j]
            q = np.abs(w*1.0e4 - wav) < 0.03
            if q.any(): 
               xnew[j] = xr[q]+ anchor[0]
               ynew[j] = yr[q]+ anchor[1] 
               
         # fix 'valid' using zmxwav  
         valid[1,xnew == -999,i] = np.nan      
         xnew[xnew == -999] = np.nan
         ynew[ynew == -999] = np.nan
         
         # update
         Zxpix[1,:,i] = xnew
         Zypix[1,:,i] = ynew  
         dis = 0
         wav = 0
               
   if correct_field:
      Zxfield, Zyfield = align_zemax_field(Zxfield[1,3,:], Zyfield[1,3,:],wheelpos,)
   else:
      Zxfield, Zyfield = Zxfield[1,3,:], Zyfield[1,3,:] 
            
   return Zwave,Zorders, Zxfield,Zyfield,  Zxpix,Zypix, flux, valid
   
   
def normalize_spectrum(wave, spectrum, standard_wave, standard_spectrum, ):
   ''' 
   Take a spectrum and normalize it against another spectrum
   after correlating them. 
   This is useful for comparing to the zemax flux calculation 
   where the model throughput has been normalized to the 
   spectrum at the boresight position.
   
   2012-05-22 NPMK
   '''
   from uvotgrism import spectrumshift
   
   # correlate the spectra and shift 
   k, (w1,s2) = spectrumpixshift(wave, spectrum, standard_wave,
       standard_spectrum, wmin=None, wmax=None, spectrum=True )
   
   # ratio 
   f = spectrum/s2   
   return w1, f 
   

def _select_zmxpoints(wave, spectrum, zmx=[190,210,235,260,285,320,
     350,380,420,460,500,540,580,620,660,700], ):
   '''interpolate to the zemax wavelengths '''   
   from numpy import interp
   return zmx,interp(zmx*10.0,wave,spectrum)
   

def find_zmxflux_match(w_spect,norm_spect,w_model,modelflux):
    '''
    find least squares match of normalized spectrum to zemax_flux models
    returns index k of best fit in modelflux[1,:,k] 
    
    -- never used
    '''   
    import numpy as np
    # sanity check wavelengths match
    if w_spect != w_model: 
       print("find_zmxflux_match : wavelenths don't match ", w_spect, w_model)
       raise
       return None
    # check modelflux has 3 dimensions 
    if len(modelflux.shape) != 3: 
       print("modelflux should be the same shape as the flux array returned by read_zemax_flux()")
       return    
    
    N = (modelflux.shape)[-1]
    lsq_val = np.zeros(N)
    
    for i in range(N):
       lsq_val[i] = ((norm_spect - modelflux[1,:,i])**2).sum()
    
    return np.where(np.min(lsq_val) == lsq_val)   
 
 
def _helper_flx(xpix2,ypix2,xp,yp,leg=[]): 
    '''Find array element in fluxes to plot for given position '''
    leg.append('[%6.1i,%6.1i]'%(xp,yp))  # legend
    q = (xpix2[1,3,:] > (xp-d)) & (xpix2[1,3,:] < (xp+d)) & (ypix2[1,3,:] > (yp-d)) & (ypix2[1,3,:] < (yp+d))
    return where(q)[0],[0], leg
    
def _patch_arf_HDU(hdutopatch, hdutouse, fituptowave, scale):
   '''for arf files with missing part of the wave length coverage,
   patch the missing range with a scaled version of another EA curve.
   Parameters
   ---------- 
   hdutopatch : fits.hdu.table.BinTableHDU 
      hdu of arf to be patched
   hdutouse : fits.hdu.table.BinTableHDU
      hdu of arf to use for patching
   fituptowave : float
     below this wavelength fit the scaled hdutouse EA data to hudtobepatched 
   scale : float
     factor to multiply EA of hdutouse
   
   Returns
     patched hdu
     
   Notes:
     error is doubled in patched region   
   '''  
   import numpy as np
   from astropy.io import fits
   import pyfits, astropy
   
   isinstance(hdutopatch,(pyfits.hdu.table.BinTableHDU,astropy.io.fits.hdu.table.BinTableHDU))
   isinstance(hdutouse      ,(pyfits.hdu.table.BinTableHDU,astropy.io.fits.hdu.table.BinTableHDU))
   
   qin = hdutopatch.data['wave_min'] >  fituptowave
   qp  = hdutouse.data['wave_min'] <= fituptowave
   
   nin = len(np.where(qin)[0])
   np  = len(np.where(qp )[0])
   if isinstance(hdutopatch,astropy.io.fits.hdu.table.BinTableHDU):
      hdu = fits.new_table( hdutopatch.columns, header=hdutopatch.header, nrows=nin+np)
   elif isinstance(hdutopatch,pyfits.hdu.table.BinTableHDU) :   
      hdu = pyfits.new_table( hdutopatch.columns, header=hdutopatch.header, nrows=nin+np)
   else :
      print("problem with HDU definition")   

   for name in ['ENERG_LO', 'ENERG_HI', 'WAVE_MIN', 'WAVE_MAX']: 
      hdu.data[name][nin:] = hdutouse.data[name][qp] 
      hdu.data[name][:nin] = hdutopatch.data[name][qin]
   hdu.data['specresp'][nin:] = hdutouse.data['specresp'][qp] * scale
   hdu.data['specresp'][:nin] = hdutopatch.data['specresp'][qin]
   hdu.data['sprs_err'][nin:] = hdutouse.data['sprs_err'][qp] * 2.0 
   hdu.data['sprs_err'][:nin] = hdutopatch.data['sprs_err'][qin]

   hdu.header['COMMENT'] = "Patched effective area below "+str(fituptowave)+" A"      
             
   return hdu 
    
    
    
def sum_spectra( Ylist ,shfts=None,wheelpos=1000,chatter=0):
   ''' input a list of Y things where, for Y = [Y1 .. YN], for each Yi, 
   dis=Yi[0][0], spnet=Yi[0][1], bg=Yi[1][0], C_1 = Yi[2][0], hdr = Yi[3]
   and the indices for the effective limits for the spectrum are m1=Yi[4],m2=Yi[5]
   
   shifts is a list of the same length with pixel offsets to be applied 
   to perfectly line up the spectra.
   
   wheelpos is the position of the filterwheel
   
   This routine then generates a common wavelength scale, and resamples 
   all spectra on that scale, after dividing each wavelength bin by the total
   exposure time. 
   
   need to put in mask to multiply out -> circular aperture
   
   '''
   from numpy import arange,zeros, polyval,polyfit,min,max,array
   if wheelpos < 500: 
      w = arange(1700,7200,1) 
   else:
      w = arange(2600,7200,1)
   expmap=zeros(len(w))
   sprate=zeros(len(w))

   if shfts == None:
      ddis = zeros(len(Ylist)) 
   else:
      ddis = array(shfts)
   
   yn=0
   for Y in Ylist:
      dis = Y[0][0]+ddis[yn] ; spnet = Y[0][1] ; bg = Y[1][0]
      yn += 1
      C_1 = Y[2][0]
      hdr = Y[3] ; m1=Y[4] ; m2=Y[5]
      aa = arange(m1,m2)
      wav = polyval(C_1,dis[aa])        # wavelength middle bin
      wavlo = polyval(C_1,dis[aa]+0.5)  # lower wavelength boundary 
      wavup = polyval(C_1,dis[aa]-0.5)  # upper ...
      spnet = spnet[aa]
      w1 = min(wav) ; w2 = max(wav) 
      i1 = w.searchsorted(w1)-1 ; i2 = w.searchsorted(w2) 
      if i1 <= 0: i1 = 1
      if i2 >= len(w)-1: i2 = len(w)-2
      expmap[i1:i2] += hdr['exposure']
      for j in range(i1,i2):
         spbin = 0.
         wlo = 0.5*(w[j] + w[j-1]) 
         wup = 0.5*(w[j] + w[j+1])
         klo = wavlo.searchsorted(wlo)-1 ; kup=wavup.searchsorted(wup)
         if len(wavup) <= kup:
            break
         if chatter > 5: print('====  ',j,' - ',klo,kup,' - w: ',wlo,wup , ' - wavlo[klo]: ',\
            wavlo[klo],' - wavhi[kup]: ',wavup[kup],'  =======') 
         if abs(wavlo[klo] - wavup[kup-1]) < 1e-8: 
            spbin = spnet[klo]*(wup-wlo)
            if chatter > 5: print('A - ',spbin)
         elif abs(wavlo[klo] - wavup[kup-2]) < 1e-8:
            spbin  = spnet[klo]*(wavlo[klo+1]-wlo) 
            if chatter > 5: print('B1 - ',spbin, spnet[klo],wavlo[klo+1], wlo)
            spbin += spnet[kup]*(wup-wavup[kup-1])
            if chatter > 5: 
               print('B2 - ',spnet[kup]*(wup-wavup[kup-1]),spnet[kup], wup, wavup[kup-1])
         else:
            spbin = spnet[klo]*(wavlo[klo+1]-wlo) 
            if chatter > 5: print('C1 - ',spbin, spnet[klo+1],wavlo[klo], wlo)
            spbin += spnet[kup]*(wup-wavup[kup-1])
            if chatter > 5: 
               print('c2 - ',spnet[kup]*(wup-wavup[kup-1]),spnet[kup], wup, wavup[kup-1])
               print('kup, klo ', kup, klo)
            for k in arange(klo+1,kup-1): 
               spbin += spnet[k]*(wavup[k]-wavlo[k])
               if chatter > 5: print('C3 -  ',k,'   ', spnet[k]*(wavup[k]-wavlo[k]),\
                  spnet[k],  wavup[k], wavlo[k])
         sprate[j] += spbin
         if chatter > 4: print('sprate[',j,'] = ',sprate[j])
      
   sprate = sprate/expmap 
   return w, sprate   

###################################################################################
# photometry & coincidence loss ###################################################

# Region Boundaries
#------------------- 
# (from http://hea-www.harvard.edu/RD/funtools/regbounds.html)
#
# The golden rule for spatial region filtering was first enunciated by
# Leon VanSpeybroeck in 1986:
#
# Each photon will be counted once, and no photon will be counted more
# than once.
#
# Image boundaries : radially-symmetric shapes (circle, annuli, ellipse)
#------------------
#
# For image filtering, pixels whose center is inside the boundary are
# included. This also applies non-radially-symmetric shapes. When a
# pixel center is exactly on the boundary, the pixel assignment rule
# is:
#
#    * the outer boundary of a symmetric shape does not include such pixels
#    * the inner boundary of a symmetric shape (annulus) includes such pixels 
#
# In this way, an annulus with radius from 0 to 1, centered exactly on
# a pixel, includes the pixel on which it is centered, but none of its
# neighbors. These rules ensure that when defining concentric shapes,
# no pixels are omitted between concentric regions and no pixels are
# claimed by two regions. When applied to small symmetric shapes, the
# shape is less likely to be skewed, as would happen with
# non-radially-symmetric rules. These rules differ from the rules for
# box-like shapes, which are more likely to be positioned adjacent to
# one another.
#
# Image Boundaries : non-radially symmetric shapes (polygons, boxes)
#------------------
#
# For image filtering, pixels whose center is inside the boundary are
# included. This also applies radially-symmetric shapes. When a pixel
# center is exactly on the boundary of a non-radially symmetric
# region, the pixel is included in the right or upper region, but not
# the left or lower region. This ensures that geometrically adjoining
# regions touch but don't overlap.
#

def counts_in_circle(img, xy_c, radius):
    import numpy as np
    # NB xpix and ypix arrays run from 0 and have shape img.shape
    ypix, xpix = np.indices(img.shape)

    r_tmp = np.sqrt((xpix - xy_c[0])**2 + (ypix - xy_c[1])**2)

    mask = (r_tmp < radius)

    pix_in_reg = img[mask]

    npix = pix_in_reg.size

    counts = pix_in_reg.sum()

    return (counts, npix)

def counts_in_annulus(img, xy_c, r1, r2):
    import numpy as np
    # NB xpix and ypix arrays run from 0 and have shape img.shape
    ypix, xpix = np.indices(img.shape)

    r_tmp = np.sqrt((xpix - xy_c[0])**2 + (ypix - xy_c[1])**2)

    mask = (r_tmp >= r1) & (r_tmp < r2)

    pix_in_reg = img[mask]

    npix = pix_in_reg.size

    counts = pix_in_reg.sum() 

    return (counts, npix)

def process_galaxy_photometry(infile="M82_UVW2_sk.fits",outfile="raw_photometry.out",
   regionfile="M82_glx.reg",exposurefile="M82_UVW2_ex.fits", bkgannulus=[15.0,22.0], 
   ext=0,clobber=False, chatter=0):
   '''
   process the image and exposure file using the regions in the region file 
   the WCS info is read from the header
   
   Parameters
   ----------
   infile, exposurefile, outfile, regionfile : str
     file names of input uvot sky image, exposure, output file name, region file
     
   ext : int
     extension to use
       
   bkgannulus : list  
     inner and outer radius of background annulus in arcsec
     
   clobber : bool
     if set, the old file will be overwritten
       
   chatter : int
     verbosity
   
   '''
   import os
   import numpy as np
   from astropy.io import fits
   from astropy.wcs import wcs
  
   if not os.access(infile,os.F_OK):
      raise IOError ("image file not found")
   if not os.access(exposurefile,os.F_OK):
      raise IOError ("exposure file not found")
   if not os.access(regionfile,os.F_OK):
      raise IOError ("region file not found")
          
   pixscale = 0.5036   # pixel size in arcsec     
          
   (ds9version,orifile,epoch,xwcs),(signs,boxtype,position,box) = _parse_DS9regionfile(regionfile)
   nregion = len(position)

   img = fits.open(infile)
   w_img = wcs.WCS(img[ext].header) 
   # convert positions to pixel coordinates   
   pix_xy_img = w_img.wcs_world2pix(np.asarray(position,dtype=float),1)
      
   ex  = fits.open(exposurefile)
   #w_ex = wcs.WCS(ex[ext].header) # no WCS header
   # convert positions to pixel coordinates   
   # pix_xy_ex = w_ex.wcs_world2pix(np.asarray(position,dtype=float),1)
   
   r1,r2,r3 = 5.0/pixscale,bkgannulus[0]/pixscale,bkgannulus[1]/pixscale
   ri1 = int(r1+0.5) ;ri2 = int(r2+0.5) ; ri3 = int(r3+0.5)
   print("radii : ",r1,r2,r3)
   
   print("k, src_counts, src_npix, src_exposure, bkg_counts, bkg_npix, bkg_exposure  src_rate  bkg_rate, position")
   for k in range(nregion):
      # subimage
      simg = img[ext].data[pix_xy_img[k][0]-ri3-1:pix_xy_img[k][0]+ri3+2,pix_xy_img[k][1]-ri3-1:pix_xy_img[k][1]+ri3+2]
      sex = ex[ext].data[pix_xy_img[k][0]-ri3-1:pix_xy_img[k][0]+ri3+2,pix_xy_img[k][1]-ri3-1:pix_xy_img[k][1]+ri3+2]
      src_counts,src_npix = counts_in_circle(simg,[ri3+1,ri3+1],r1)
      src_exp, sexp_npix = counts_in_circle(sex,[ri3+1,ri3+1],r1)
      bkg_counts, bkg_npix = counts_in_annulus(simg,[ri3+1,ri3+1],r2,r3)
      bkg_exp, bexp_npix = counts_in_annulus(sex,[ri3+1,ri3+1],r2,r3)
      bkg_exposure = bkg_exp/bexp_npix
      src_exposure = src_exp/sexp_npix
      src_rate = src_counts / src_exposure
      bkg_rate = bkg_counts / bkg_exposure * src_npix / bkg_npix
      print("%03i %10.1f %7.1f %8.2f %10.1f %7.1f %8.2f - %8.3f  %8.3f - (%7.1f,%7.1f)"%\
       (k, src_counts, src_npix, src_exposure, bkg_counts, bkg_npix, bkg_exposure,\
       src_rate ,bkg_rate, pix_xy_img[k][0],pix_xy_img[k][1] ))
   img.close()
   ex.close()  


def process_galaxy_photometry_maghist(infile="galaxy_sum.img",outfile="maghist.out",
   regionfile="regions.reg", bkgannulus=[15.0,22.0], clobber=False, chatter=0):
   '''
   Use positions in a region file to process the photometry in the galaxy summed image
   using uvotmaghist repeatedly using an annulus background for each region. 
  
   Parameters
   ----------
   infile, outfile, regionfile : str
     file names of input uvot sky image, maghist output file, a region file
     
   bkgannulus : list  
     inner and outer radius of background annulus in arcsec
     
   clobber : bool
     if set, the old file will be overwritten
       
   chatter : int
     verbosity
  
   Notes
   -----
             
   '''
   import os
  
   if not os.access(infile,os.F_OK):
      raise IOError ("image file not found")
   if not os.access(regionfile,os.F_OK):
      raise IOError ("region file not found")
          
   (version,filename,epoch,wcs),(signs,boxtype,position,box) = _parse_DS9regionfile(regionfile)
   
   bkgr1 = str(bkgannulus[0])
   bkgr2 = str(bkgannulus[1])
   if clobber : clob = "yes"
   else: clob = "no"

   region_number = len(position)
   if chatter > 1 : print("number of positions = %i"%region_number)
   
   # check that the coordinate system/epoch is FK5 
   if epoch.upper() != "FK5": 
      print("found for system/epoch: ",epoch)
      raise RuntimeError("Fatal Error: region file coordinate system is not FK5")
   
   # main loop
   for k in range(0,region_number):
      if k > 0 : clob = "No"
      if chatter > 1: print(position[k])
      ra,dec = position[k][0],position[k][1]
      if chatter > 1: print("coordinate: ", ra,dec) 
      # write background annulus
      f = open("bkgreg.01234.reg","w")
      f.write('fk5\ncircle(%s,%s,%s")\n-circle(%s,%s,%s")\n' % (ra,dec,bkgr2,ra,dec,bkgr1))
      f.close()
      # write source circle
      f = open("srcreg.01234.reg","w")
      f.write('fk5\ncircle(%s,%s,5.0")\n' % (ra,dec))
      f.close()
      # run uvotmaghist
      command = "uvotmaghist infile="+infile+" outfile="+outfile+" plotfile=None"+\
      " srcreg=srcreg.01234.reg apercorr=None bkgreg= bkgreg.01234.reg clobber="+clob  
      print("processing region number : %i" % k)
      if chatter > 1 : print(command)          
      status = os.system(command)
      
   # Done  


def _parse_DS9regionfile(file,chatter=0):
   '''
   parse the region file
   
   Note
   ----
   return structure with data
   so far only for circle() 
   does not grab colour or annotation metadata
   '''
   F = open(file)
   f = F.readlines()
   F.close()
   
   signs = []
   position = []
   box = []
   boxtype = []
   
   try:
      if f[0].split(":")[1].split()[0] == "DS9":
         version=f[0].split(":")[1].split()[-1]
      else:
         version="0"
      filename = f[1].split(":")[1].split("\n")[0]
      epoch = f[3].split("\n")[0].split(":")[0] 
      wcs = 'wcs'
      r = f[3].split("\n")[0].split(":")
      if len(r) > 1:  # other coordinate system definition 
         wcs = r[1]                              
   except:
     print("Error reading region file : ",file)
   try: 
     r = f[4].split("\n")[0]
     if chatter > 3: 
        print("line# 4",r)
     if r[0:4].upper() == "WCS":
        wcs = r
     elif len(r) == 0:
        do_nothing = True       
     elif r.split("(")[0] == "circle" :
        signs.append("+")
        boxtype.append("circle")
        values = r.split("(")[0].split(")")[0].split(",")
        position.append(values[0:2])
        box.append(values[2:])
     elif r.split("(")[0] == "-circle" :
        signs.append("-")
        boxtype.append("circle")
        values = r.split("(")[0].split(")")[0].split(",")
        position.append(values[0:2])
        box.append(values[2:])  
     else:
        print("problem with unknown region type - update _parse_DS9regionfile() ")
                
   except:
     print("problem reading end header region file ")   
   
   for k in range(5,len(f)):
     try:
        r = f[k].split("\n")[0]
        if chatter > 3:
           print("line# ",k,' line=',r)
        elif r == "\n":
           continue
        elif r.split("(")[0] == "circle" :
           signs.append("+")
           boxtype.append("circle")
           values = r.split("(")[1].split(")")[0].split(",")
           position.append(values[0:2])
           box.append(values[2:])
        elif r.split("(")[0] == "-circle" :
           signs.append("-")
           boxtype.append("circle")
           values = r.split("(")[1].split(")")[0].split(",")
           position.append(values[0:2])
           box.append(values[2:])       
        else:
           print("problem with unknown region type - update _parse_DS9regionfile() ")
                
     except:
        print("problem reading region record number = ",k)
        
   return (version,filename,epoch,wcs),(signs,boxtype,position,box)


def uvot_alpha(frametime,readout=171e-6):
   """return the deadtime correction factor
   
   parameters
   ----------
   frametime : float, float array[:]
     frametime is seconds
   
   kwargs: dict
     optional keyword arguments, possible values are:  
     
     - **readout** : float
     readout time in seconds (rows*clockspeed, 285*600ns )  
   
   returns
   -------
   alpha : float or float array[:]
      the dead time correction factor 
   """
   return (frametime - readout)/frametime


def cr_per_arcsec_to_cpf_in_aperture(cr,frametime,aperture_radius):
    """
    convert from counts per second per arcsec to counts per frame (no deadtime correction)
    
    Parameters
    ----------
    cr : float, float array[:]
       count rate per arcsec**2

    kwargs : dict
    - **frametime** : float
       frame time in sec
       
    - **aperture_radius** : float 
       radius aperture in arcsec
       
    Returns
    -------
    cpf : float, float array[:]
       counts per frame   
    """
    from numpy import pi
    return cr * frametime * pi * aperture_radius**2 
    

def read_maghistout(fitsfile):
   """
   read the output from uvotmaghist 

   parameters
   ----------
   fitsfile : file path, file object or file-like object
      File to be opened.
      
   Returns
   -------
   rawdata : list
     a list of the raw data extracted using uvotmaghist    

   """
   from astropy.io import fits as pyfits
   hdu = pyfits.open(fitsfile)
   bintable = hdu["maghist"]
   tab = bintable.data
   rawdata = {}
   rawdata.update({"exposure":tab["exposure"]} )
   rawdata.update({"tmid_swift": 0.5*(tab["tstart"]+tab["tstop"])} ) # swift time mid-exposure
   rawdata.update({"src_area": tab["src_area"]})  # in arcsec**2
   rawdata.update({"std_area": tab["std_area"]})
   rawdata.update({"bkg_area": tab["bkg_area"]})
   rawdata.update({"raw_tot_cnts": tab["raw_tot_cnts"]})
   rawdata.update({"raw_std_cnts": tab["raw_std_cnts"]})
   rawdata.update({"raw_bkg_cnts": tab["raw_bkg_cnts"]})
   rawdata.update({"ft": tab["framtime"]})
   rawdata.update({"filt": tab["filter"]})
   rawdata.update({"apercorr": tab["ap_factor"]})    # aperture correction (multiply)
   rawdata.update({"lsscorr" : tab["lss_factor"]})  # large scale sensitivity (divide)
   rawdata.update({"senscorr": tab["senscorr_factor"]}) # sensitivity degradation (multiply)
   return rawdata 
   


def raw_to_cpf_to_rates(fitsfile,A=1.0):
   """
   reads raw data from the fitsfile, and computes counts per frame and count rates 
   
   Parameters
   ----------
   fitsfile : str, path
      file name, path of the output from uvotmaghist
   
   kwargs : dict
   - data : dict
   
   Returns
   -------
   data : dict
     raw data, count per frame data, count rates 
     **tmid_swift** : key
       swift time mid-time of exposure
     **exposure** : key
       exposure time (sec) corrected for dead time
     **src_area** : key
       area in arcsec**2 of source region
     **std_area** : key
       area in arcsec**2 of standard 5" region
     **bkg_area** : key
       area in arcsec**2 of background region
     **raw_tot_cnts** : key
       total counts measured in source region
     **raw_std_cnts** : key
       total counts measured in standard region
     **raw_bkg_cnts** : key
       total counts measured in background regions
     **ft** : key
       frame time                
     **filt** : key
       filter name
     **apercorr** : key 
       aperture correction to source region result the result needs to be 
       multiplied by
     **lsscorr** : key
       lss correction factor the result needs to be divided by    
     **senscorr** : key
       sensitivity degradation factor the result needs to be multiplied to
     **tot_cpf** : key
       total counts per frame in standard source region
     **bkg_cpf** : key
       background counts per frame in standard source region
     **bkg** : key
       measured background rate per arcsec^2     
     **bkg_coi_f** : key
       background coi-factor to multiply raw rate to
     **bkg_coi** : key
       background counts per frame corrected using standard coi-loss Poole et al. 2008
     **tot_coi** : key  
       total counts per frame corrected using standard coi-loss Poole et al. 2008
     **bkg_coi_n** : key  
       background   extended coi correction Kuin(...)  (5" background radius)
     **tot_coi_n** : key
       total cpf    extended coi correction Kuin(...)
     **bkg_coi_p** : key
       bkg_coi corrected with the polynomial Poole et al.
     **tot_coi_p** : key  
       tot_coi corrected with the polynomial Poole et al.
     **bkg_coi_np** : key
       bkg_coi_n corrected with the polynmial Poole et al.
     **tot_coi_np** : key  
       tot_coi_n corrected with the polynomial Poole et al.
     **net_coi_p** : key
       net source rate (standard coi) corrected with polynomial Poole et al.
     **net_coi_np** : key  
       net source rate (extended coi) corrected with polynomial Poole et al.
     **net_cpf_poserr : key
       net source count/frame, +positive error (Kuin & Rosen) 
     **net_cpf_negerr : key
       net source count/frame, -negative error (Kuin & Rosen) 
     **net_src** : key
       count rates (counts/sec in area) corrected for sensitivity and LSS, 
       but not apercorr
     **net_src_n** : key 
       count rates (counts/sec in area) corrected for sensitivity and LSS,
       but not apercorr 
     **bkg_out** : key
       count rates (counts/sec in area) corrected for sensitivity and 
       LSS, but not apercorr
     **tot_out_n** : key  
       count rates (counts/sec in area) corrected for sensitivity and LSS, 
       but not apercorr
     **net_src_poserr : key
       net source count rate, +positive error (Kuin & Rosen); LSS, sensitivity corrected 
     **net_src_negerr : key
       net source count rate, -negative error (Kuin & Rosen); LSS, sensitivity corrected 
     **bkg_poserr : key
       background count rate, +positive error (Kuin & Rosen); LSS, sensitivity corrected 
     **bkg_negerr : key
       background count rate, -negative error (Kuin & Rosen); LSS, sensitivity corrected 
     
   Notes
   -----
   The aperture used for the source in uvotmaghist has been ignored here. 
   Instead, the standard 5" aperture and data has been used.  
   The error in the LSS correction has been ignored.
   
   NPM Kuin, Mullard Space Science Laboratory, University College London, 1 March 2013   
   """
   from uvotpy import uvotmisc
   from numpy import pi, sqrt, abs
   from . import uvotcoi 
   
   background_radius = 6.5 # was 5.02 Poole et al. 
    # 5.02 arcsec would be correct, but adds to net counts in new 
   
   data = read_maghistout(fitsfile)
   ft = data["ft"]   # frame time
   lss = data["lsscorr"]  # inverse large scale sensitivity factor
   sen = data["senscorr"]   # sensitivity degradation factor
   alpha = uvotcoi.uvot_alpha(ft)  # dead time correction
   telap = data["exposure"] / alpha # elapsed time
   nf = telap/ft  # number of frames
   
   # compute counts measured in a frame (no dead time correction)
   # source+bkg aperture = std area
   tot_cpf = data["raw_std_cnts"]/nf
   # background aperture = source aperture 
   bkg_cpf = data["raw_bkg_cnts"]*data["std_area"]/data["bkg_area"]/nf     
   data.update({"tot_cpf": tot_cpf})
   data.update({"bkg_cpf": bkg_cpf})

   # background coi-factor   
   bkg = data["raw_bkg_cnts"]/telap/data["bkg_area"]  # background rate per arcsec^2
   bkg_coi_f = uvotcoi.uvot_bkg_coi_factor_1(bkg,frametime=ft,aperture_radius=background_radius)
   data.update({"bkg": bkg})
   data.update({"bkg_coi_f": bkg_coi_f})

   # apply the theoretical coi-correction (_n - new coi with neighbor)   
   tot_coi = uvotcoi.coi_point_source_1(tot_cpf,alpha=1)
   bkg_coi = bkg_coi_f * bkg_cpf
   data.update({"bkg_coi": bkg_coi})
   data.update({"tot_coi": tot_coi})

   #  extended coi correction Kuin(...) 
   tot_coi_n, bkg_coi_n = uvotcoi.coi_point_source_2(tot_cpf,bkg_cpf,A=A)
   data.update({"bkg_coi_n": bkg_coi_n})
   data.update({"tot_coi_n": tot_coi_n})
    
   # apply the polynomial correction (_p)
   tot_coi_p = tot_coi * uvotcoi.uvot_coi_correction_polynomial(tot_cpf)
   tot_coi_np = tot_coi_n * uvotcoi.uvot_coi_correction_polynomial(tot_cpf)
   bkg_coi_p = bkg_coi * uvotcoi.uvot_coi_correction_polynomial(bkg_cpf)
   bkg_coi_np = bkg_coi_n * uvotcoi.uvot_coi_correction_polynomial(bkg_cpf)
   data.update({"bkg_coi_np": bkg_coi_np})
   data.update({"tot_coi_np": tot_coi_np})
   data.update({"bkg_coi_p": bkg_coi_p})
   data.update({"tot_coi_p": tot_coi_p})
         
   # net counts per frame source (no dead time correction yet done)
   net_coi_p = tot_coi_p - bkg_coi_p  # Poole et al. standard point source but bkg coi for 6.5" rad
   net_coi_n = tot_coi_np -bkg_coi_np # New extended coi-correction (with polynomial correction)
   data.update({"net_coi_p": net_coi_p})
   data.update({"net_coi_n": net_coi_n})
   
   # count rates (counts/sec in area) corrected for sensitivity and LSS, but not apercorr
   net_src   = net_coi_p / alpha / ft /lss *sen
   net_src_n = net_coi_n / alpha / ft /lss *sen
   bkg_out   = bkg_coi_p / alpha / ft /lss *sen
   bkg_out_n = bkg_coi_n / alpha / ft /lss *sen 
   data.update({"net_src": net_src})
   data.update({"net_src_n": net_src_n})
   data.update({"bkg_out": bkg_out})
   data.update({"bkg_out_n": bkg_out_n})
   
   # errors 
   # - (binomial measurement error -> use coi_point_source_2( x+err) etc.
   tc = tot_cpf * nf  # measured counts total in 5" radius aperture
   bc = bkg_cpf * nf  # background counts total in 5" radius aperture
   tv = sqrt( tc*(nf-tc)/nf )  # binomial distribution sqrt(variance) source counts
   bv = sqrt( bc*(nf-bc)/nf )  # sqrt(variance) background counts
   tcpf1 = (tc + tv)/nf * tot_coi_n/tot_cpf  # total count per frame measured plus error
   tcpf2 = (tc - tv)/nf * tot_coi_n/tot_cpf  # total count per frame measured minus error
   bcpf1 = (bc + bv)/nf * bkg_coi_f # total count per frame coi-corrected plus error
   bcpf2 = (bc - bv)/nf * bkg_coi_f # total count per frame coi-corrected minus error
   ncpf1 = tcpf1 - bcpf2  # net cpf coi-corrected net +error
   ncpf2 = tcpf2 - bcpf1  # net cpf coi-corrected net -error 
   data.update({"net_cpf_poserr": abs( +ncpf1 - net_coi_n) })
   data.update({"net_cpf_negerr": abs( -ncpf2 + net_coi_n) })
   net_src_poserr =  ncpf1 /alpha /ft /lss *sen - net_src_n
   net_src_negerr = -ncpf2 /alpha /ft /lss *sen + net_src_n
   bkg_poserr =  bcpf1 /alpha /ft /lss *sen - bkg_out
   bkg_negerr = -bcpf2 /alpha /ft /lss *sen + bkg_out
   data.update({"net_src_poserr": abs(net_src_poserr)})
   data.update({"net_src_negerr": abs(net_src_negerr)})
   data.update({"bkg_poserr": abs(bkg_poserr)})
   data.update({"bkg_negerr": abs(bkg_negerr)})
   
   # - LSS error, ignored
   # - aperture error (smaller aperture not implemented)
   # sens. error negligible
   # post procesing: 
   # CR->flux error
   # zero point error CR->mags
   
   return data
##########################################################################################    
   
def ea_plot_weak():
   '''
   make plot of EA using weak stars
   
   
   import cal4
   from uvotpy.uvotmisc import rdTab
   cd ../fluxcal2/eff_area_200/ea_2.5_sums
   Z=cal4.plot_ea_200_default(s=15)
   wav_out,ea_out,error,weight,w,c,ea =Z
   cd eff_area_160/ea_2.5_sums/
   t = rdTab('eamerged_weak_160_groups_A.EA.err')
   cd ../../eff_area_1000/ea_2.5_sums/
   t1 = rdTab('eamerged_weak_1000_groups_X.EA.err')

   cd ../../eff_area_955/ea_2.5_sums/
   t2 = rdTab('eamerged_weak_955_groups_X.EA.err')

In [76]: fill_between(wav_out,ea_out*(1-0.01*error),ea_out*(1+0.01*error),color='blue',alpha=0.3)
Out[76]: <matplotlib.collections.PolyCollection at 0x11be76a10>

In [77]: plot(wav_out,ea_out,'b.-',label='UV nominal',lw=1,markersize=2.5)
Out[77]: [<matplotlib.lines.Line2D at 0x117d232d0>]

In [78]: fill_between(t[:,0],t[:,4]*(1-0.01*t[:,5]),t[:,4]*(1+0.01*t[:,5]),color='y',alpha=0.5)
Out[78]: <matplotlib.collections.PolyCollection at 0x117d229d0>

In [79]: plot(t[:,0],t[:,4],'.-',markersize=2.5,color='darkgreen',lw=1.,label='UV clocked')
Out[79]: [<matplotlib.lines.Line2D at 0x117d307d0>]

In [80]: legend()
Out[80]: <matplotlib.legend.Legend at 0x117d30450>

In [81]: title(u'Effective area at default position using weak sources')
Out[81]: <matplotlib.text.Text at 0x120aac6d0>

In [82]: xlabel(u'$\lambda (\AA)$ ',fontsize=16)
Out[82]: <matplotlib.text.Text at 0x120080290>

In [83]: ylabel(u'Effective Area ($cm^2$)')
Out[83]: <matplotlib.text.Text at 0x120aa0c50>

In [84]: ylabel(u'Effective Area ($cm^2$)',fontsize=14)
Out[84]: <matplotlib.text.Text at 0x120aa0c50>

In [85]: title(u'Effective Area at default position using weak sources')
Out[85]: <matplotlib.text.Text at 0x120aac6d0>

fill_between(t1[45:,0],t1[45:,4]*(1-0.01*t1[45:,5]),t1[45:,4]*(1+0.01*t1[45:,5]),color='pink',alpha=0.5)
Out[129]: <matplotlib.collections.PolyCollection at 0x120f4c650>

In [130]: plot(t1[45:,0],t1[45:,4],'.-',markersize=2.5,color='darkred',lw=1.,label='Visual nominal')
Out[130]: [<matplotlib.lines.Line2D at 0x120f4c510>]

plot(t2[:,0],t2[:,4],'.-',markersize=2.5,color='purple',lw=1.,label='Visual clocked')
Out[164]: [<matplotlib.lines.Line2D at 0x11aeaee10>]

In [165]: fill_between(t2[:,0],t2[:,4]*(1-0.01*t2[:,5]),t2[:,4]*(1+0.01*t2[:,5]),color='c',alpha=0.2)
Out[165]: <matplotlib.collections.PolyCollection at 0x11aeb1890>

In [166]: ylim(0,50)
Out[166]: (0, 50)

In [167]: legend()
Out[167]: <matplotlib.legend.Legend at 0x11aeae550>

In [168]: xlim(1600,6000)
Out[168]: (1600, 6000)


   '''
        
##########################################################################################    
def measlines2cal3(file,anchor=4200.,acc = 1e-3,listf=False, chatter=1):
   '''
   Input sw000xxx_meas.txt measured line positions.
   
   returns list of (p,w) for line positions, corrected to have 
   p=0 at anchor, and dispersion 
   2013-08-08 npmk while check uvotgrASPcorr error
   '''
   from astropy.table import Table
   
   if not listf:
     t = Table.read(file,comment='#',format='ascii')
     c1 = t['col1'].data  # dis
     c2 = t['col2'].data  # wave
   else:
     c1 = file[:,0]
     c2 = file[:,1]  

   
   C1 = polyfit(c1,c2,3)
   Cnew = polyfit(c1,c2,3)
   delta = (C1[-1]-anchor)/C1[-2]
   c1 = c1+delta
   n=10
   while (abs(delta) > acc) &(n>0):
      Cnew = polyfit(c1,c2,3)
      delta = (Cnew[-1]-anchor)/Cnew[-2]
      c1 = c1+delta
      n = n-1
      if chatter > 0 : print(n, Cnew[-1])
   
   C1 = Cnew      
   dis1 = []
   print("dis1=[")
   for i in range(len(c1)): 
      dis1.append([c1[i],c2[i]])
      print("[%7.1f, %7.1f],"%(c1[i],c2[i]))  
   print("],")
   print("dispers1st=[%12.5e,%12.5e,%10.4f,%8.2f],"%(C1[0],C1[1],C1[2],C1[3]))
   print("dis1corr=%8.2f,"%(C1[3]-anchor))
   return dis1, C1          

   
##########################################################################################    
    # VISIBLE GRISM 2ND ORDER Eff Area    
##########################################################################################    
#  Visible grism, second order   00031417004_1 (U Sco) overexposed first order, second order 
#  dist12 ~ 660 (model 546.3) Estimate second order flux from observation, first order from
#  photometry. H-alpha? visible in first order (from comparison to Fred Walters spectrum).
#  just before sensitivity drops off completely. Second order response below 3800A seems 
#  very low. Spectrum is for 29 Jan 2010. Found optical spectrum. Need to use photometry 
# to extend below 3400. Strong H+He lines + rising continuum. 
#
##########################################################################################    
# for reference : optical spectrum Anupama et al. J/A+A/559/A121                                                                   uvot photometry (vega mag) 
# file         date        JD (HJD)        T-T0  phase                  datetime                uvot spectrum                      uvw2      uvm2      uvw1       u         v 
#27may09.dat   27/May/2009 2454979.2974H         0.166 3600s quiescence 2009-05-27T19:08:15.36  
#29jan10.dat   28/Jan/2010 2455225.5212    0.831 0.256                  2010-01-29T00:30:31.68  31417004+1  2010-01-29T01:08:04  9.1\pm0.1  ---------  8.5\pm0.2 --------- 9.3\pm0.2
#   -f_lam :                                                                                                                     5.46e-15              3.66e-15           2.43e-15
#30jan10.dat   29/Jan/2010 2455226.5152    1.825 0.064                  2010-01-30T00:21:53.28  31417004+17 2010-01-30T00:52:21  9.3\pm0.2  9.7\pm0.2  9.0\pm0.1 ---------10.1\pm0.2
#   -f_lam                                                                                                                       5.58e-15   7.28e-15   3.87e-15           2.64e-15       
#31jan10.dat   30/Jan/2010 2455227.5110    2.821 0.874                  2010-01-31T00:15:50.40  
#
#


##########################################################################################    
#   check coincidence loss relation between predicted rate and observed rate 
#   units(cpf), including background 
##########################################################################################  
  
def convert_flux_to_rate(
    obsfile='/calibration/grism/fluxcal2/eff_area_200/all_2.5/sw00055500010ugu_1ord_1.pha',
    calspectrum='/calibration/grism/fluxcal2/uvotified/gd153_stisnic_003_uvotified.ascii',
    ):
    """
    parameters
    ----------
    obsfile : path
       observed spectrum
    calspectrum : path
       calibrated spectrum
    
    USE calibration.calcoi.normalize_spec()
       
    """
    from astropy.io import fits,ascii
    import numpy as np
    from uvotpy.uvotio import readFluxCalFile
    from scipy import interpolate
    from uvotpy import uvotmisc
    from uvotpy.uvotio import sensitivityCorrection
    
    # get total observed rate and background 
    f = fits.open(obsfile)
    obswave = f[2].data['lambda']
    sigcoef = uvotmisc.get_sigCoef(f[1].header)
    dis = f[2].data['pixno']
    sigma = np.polyval(sigcoef, dis)
    swifttime = f[2].header['tstart']
    if 'filevers' in f[0].header:
        if f[0].header['filevers'] == 2:
            obsrate = f[2].data['netrate']+f[2].data['bkgrate1']  #version 2 file 
            obsbkgrate = f[2].data['bkgrate1'] #version 2 file
        else:
            # version 1 file        
            obsrate = (f[1].data['counts']/f[1].data['exposure'] )  # version 1 file
            obsbkgrate = f[2].data['bg_L']          
    else:
        obsbkgrate =  f[2].data['bg_L'] # version 1 file
        obsrate = (f[1].data['counts']/f[1].data['exposure'] ) # version 1 file
        
    g = ascii.read(calspectrum)
    calwave = g['col1']
    calflux = g['col2']
    z =  readFluxCalFile(160,anchor=[1023,1085],spectralorder=1,chatter=5)
    hdu,fnorm = z
    q = np.arange(len(hdu.data['WAVE_MIN'])-1,-1,-1)
    w = (0.5*(hdu.data['WAVE_MIN']+hdu.data['WAVE_MAX']))[q]     
    r = (( hdu.data['SPECRESP'] )*fnorm(w))[q]
    specrespfunc = interpolate.interp1d(w, r,  bounds_error=False, fill_value=np.NaN)
    coef = uvotmisc.get_dispersion_from_header(f[1].header)
    binwidth = np.polyval(coef,dis+0.5) - np.polyval(coef,dis-0.5)      # width of each bin in A (= scale A/pix)
    fbinwidth=interpolate.interp1d(obswave,binwidth,bounds_error=False, fill_value=np.NaN)
    
    h_planck = 6.626e-27  # erg - s
    lightspeed = 2.9979e10  # cm/sec
    h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
    hnu = h_c_ang/calwave
    
    senscorr=1/sensitivityCorrection(swifttime)
    ft = 0.0110329  # frames/seconds
    factor = 314.26/(5*sigma)
    # flux = hnu*(sprate*senscorr)*fcoi(wave)/specresp1func(wave)/binwidth   # [erg/s/cm2/angstrom]  
    calrate = calflux*specrespfunc(calwave)*fbinwidth(calwave) /(hnu*senscorr)
    fcalrate = interpolate.interp1d(calwave,calrate,bounds_error=False, fill_value=np.NaN)
    calrate = (fcalrate(obswave)+obsbkgrate)*ft*factor 
    fcalrate = interpolate.interp1d(obswave,calrate,bounds_error=False, fill_value=np.NaN)
    fobsrate = interpolate.interp1d(obswave,obsrate*ft*factor,bounds_error=False, fill_value=np.NaN)
    return fcalrate, fobsrate, obswave 
    
    
########################################################    
# reprocess script  2014-06-05  NPMK 
# -  adding flag for auto processing to cal4.py tables    
########################################################    
    
# end uvotcal.py ============    
