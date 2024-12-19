'''zemax.py
These routines are for processing the zemax data
 
   1. read the zemax ascii files in and organize them
   2. scale, rotate and shift the zemax sky inputs to 
      get detector coordinates ignoring distortion
   3. Create a grid of detector coordinates and calculate
      the sky positions (using the inverse transform of 2.)   
   4. Difference the projected sky and calculated DET positions
      to get a distortion field.
   5. Interpolate the distortion   
   
Detailed instructions:

  some default parameters have been set in the calls, like scaling factor
  rotation and offsets. In the examples below these are not modified.
  
  First,  
  >>> import zemax

  The output from the zemax optical model are in ascii files which are
  read using a call 
  >>> (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, 
            ypix, valid) =  zemax.rdzemax()
   
  xpix,ypix are arrays with the detector coordinates in pixels with the centre at 
    (1023.5,1023.5) with dimensions for order, wavelength and point. To convert to
    DET coordinates, add [77,77]
  xfield,yfield are arrays with the input rays in degrees with the same dimensions 
     as xpix,ypix
  wave and order are arrays with the wavelengths (in micrometers) and order     

  To scale and map the sky coordinates to detector coordinates without 
  distortions (just a scaling, rotation, and offset), use:
  
  >>> (xscaled ,yscaled) = matchFieldToPix(xfield,yfield,
           order = 1, wave = 2, angle = -64.6, scale = 6550.4)

  To generate the input rays field coordinates for a regular grid of detector coordinates,
  run 
  
  >>> (xdet, ydet, xfield, yfield) = make_zemax_input()

  These are useful for finding the input field locations in a regular grid on the sky 
  that, in the absence of distortions from the grism, would yield a regular grid on 
  the detector.   
   
2009-04-06      
'''
# 2009-02-12 changed the xdetoff,ydetoff offsets in call rdzemax to 0.,0. by default
# 2009-02-16 debugged all the offsets globally, so matchDetToField quite different now.
# 2009-04-06 added some code to read_zemax4() to read s single first order; to rdzemax() 
#            to rotate/shift back the original calculation with the erronous 65 deg rotation.

import numpy
import numpy as N
from uvotpy import uvotmisc
#from uvotmisc import *
from scipy import interpolate
from os import getenv
from astropy.io import fits as pyfits
#import math

__version__ = '2009-02-11 NPMK (MSSL)'   

# 2009_02_11 GLOBAL change coordinate system to DET aby adding [77,77] to pix coordinates
#            [1023.5,1023.5] => [1100.5,1100.5]


def correctAnkPos(xpin,ypin,wheelpos=200,test=None,chatter=0):
  ''' The anker positions from the unscaled zemax model are 
   modified with the following correction.
           
   define  
          Z1 = -Xdet + 0.63462 Ydet (tan -32.4deg =  0.63462 )
          Z2 = -Xdet - 1.57575 Ydet (tan  57.6deg = -1.57575)
          LSQ fit to the observed wavelength and offset are
          
          dLam = -5.05095 -0.01850036*Z1
          doff = -48.454 - 0.0184636*Z2 (excluding 2 outlyers)
          
          The centre of the transform where dLam,doff = 0,0  is where
          Xdet-img,Ydet-img = 949.22, 1065.52
          
          The transform is nearly a linear scaling from that point 
          with a scale factor of about 0.98153 
                  
Revised method using dx=xobs-xank, dy=yobs-yank
   define  
          Z1 = -xobs + 0.63462 yobs (tan -32.4deg =  0.63462 )
          Z2 = -xobs - 1.57575 yobs (tan  57.6deg = -1.57575)
          LSQ fit to the observed wavelength and offset are
          
          dl = 0.53582*dx - 0.844323*dy
          do = 0.844323*dx + 0.53582*dy
          
          (qq is good points)
          array([ 0.00524, 1.7045]) = numpy.polyfit(z1[qq],dl[qq],1)
          array([ 0.00305, 3.7705]) = numpy.polyfit(z2[qq],do[qq],1)
          dLam = 1.70 + 0.00524*Z1
          doff = 3.77 + 0.00305*Z2 (excluding 2 outlyers)
          
          The centre of the transform where dLam,doff = 0,0  is where
          Z1 = -324.5 ,  Z2 = -1236.1
          Xdet,Ydet = 586, 412
          
          The transform is nearly a linear scaling from that point 
          with a scale factor of about 0.996 
         
   then the corrections to Xdet,Ydet are 
          alp = 32.4 * math.pi / 180.0
          Xdet = Xdet - cos(alp)*dLam -sin(alp)*doff
          Ydet = Ydet + sin(alp)*dLam -cos(alp)*doff
         
  '''
         
  if type(xpin) == float or type(xpin) == numpy.float64: 
     xpin = N.array([xpin])
  else:
     xpin = N.asarray([xpin]).flatten()   
  if type(ypin) == float or type(xpin) == numpy.float64: 
     ypin = N.array([ypin])
  else:
     ypin = N.asarray([ypin]).flatten()
        
  if chatter > 4: 
      print('correctAnkPos (',xpin,',',ypin,'wheelpos=',wheelpos,'test=',test,')')

  if wheelpos == 200: 
     if chatter > 4: print('UV Nominal')
     # convert to DET-image 
     xp = xpin - 77
     yp = ypin - 77
     Z1 = -xp + 0.63462*yp
     Z2 = -xp - 1.57575*yp
     dlam = -5.05095 -0.01850036*Z1
     doff = -48.454  -0.018436  *Z2
     alp = 32.4 * N.pi / 180.0
     xp = xpin + (-N.cos(alp)*dlam - N.sin(alp)*doff)
     yp = ypin + ( N.sin(alp)*dlam - N.cos(alp)*doff)
     return xp,yp

  elif wheelpos == 160:
      if test == "flux":
        # version 2 calibration using zemax + flux model
        # make sure X,Y are simple 1 dimensional arrays
        xpin = N.asarray(xpin)
        ypin = N.asarray(ypin)
        if (len(xpin.shape) > 1) | (len(ypin.shape) > 1):
           print("input array to correctAnkPos is wrong dimension:",type(xpin),xpin.shape,"  ",type(ypin),ypin.shape)
        # converto to image (physical) coordinate   
        X = xpin - 104  
        Y = ypin - 78
        # see flux_comparison_to_effective*.pages for notes, 
        # best fit to zemax_flux, week of 2013-01-01 
        # residuals in anchor (2.8,2.8) pixels RMS, 3 pixels in radius, so about 10A
        # 
        tckx = [N.array([-50.,-50.,-50.,1223.49126478, 2100., 2100. , 2100.]),\
        N.array([ -50., -50., -50. ,  1126.86855584, 1900. , 1900. , 1900.]),\
        N.array([ 161.95076371,  -78.89912033,   76.56623909, -176.42777675,\
        -42.54958083,   66.36035772,  -46.62329962,   82.7023882 ,\
         39.93029571,  -87.67033295,   29.18443101, -112.12711925,\
         45.07270813,   25.8779904 ,  -84.1278129 ,  -42.65666229]),  2,2]
        tcky = [N.array([ -50., -50., -50. ,1021.54170247, 2100., 2100.,2100.]),\
        N.array([ -50. , -50. , -50. , 1023.83018478, 1900., 1900. , 1900. ]),\
        N.array([-388.70677938,  153.37701596,  -70.26913226,  125.97729059,\
        127.48396157,  -28.20009961,   15.83861805,  -72.1504286 ,\
        -60.00561045,   35.30521939,  -22.16711781,    3.53266776,\
        287.86613624,  -49.73984126,  -15.55817971,  -53.68367463]), 2,2]
        xp = X + interpolate.bisplev(X, Y, tckx) 
        yp = Y + interpolate.bisplev(X, Y, tcky)                    
        return xp+104,yp+78
        
      elif test == "simple":    
        #    first guess;  simple scaling with 0.973 with respect to the boresight 
        X = N.asarray(xpin) - 1129.1
        Y = N.asarray(ypin) - 1022.3
        X = 0.973 * X 
        Y = 0.973 * Y
        return X + 1129.1, Y + 1022.3
        
      elif test == 'simple_corrected':
        # 2013-05-11 corrected anchor position using simple 
        if chatter > 4: 
           print('UV Clocked version 2 calibration 2013-05-10 NPMK')
        npnt = len(xpin)
        dx = N.zeros(npnt)
        dy = N.zeros(npnt)
        # assume the input has already been scaled by the simple model scaling:
        #   X = xpin - 1129.1
        #   Y = ypin - 1022.3
        #   X = 0.973 * X 
        #   Y = 0.973 * Y
        #   X = X + 1129.1 -104 -11
        #   Y = Y + 1022.3 -78  +5
        # 
        # xp,yp are the anchor in img coordinates
        xp = xpin - 104
        yp = ypin -  78
        # the bilinear forms from calibrating between the simple zemax 
        # model and the observations referenced to the observed anchors (img)
        #  Not forcing the boresight correction to be about zero
        tck_dx_11=[N.array([100.,100.,1156.554,2200.,2200.]),
                   N.array([-100.,-100.,868.3828,2000.,2000.]),
            N.array([-10.58800631,-1.65560843,13.39884122,-6.15361629,
            1.22695531,3.76402771,-3.06903265,4.10603678,-6.65929437]),1,1]
        tck_dx_21 = [N.array([100.,100.,100.,1145.19889314,2200.,2200.,2200.]),
            N.array([-100.,-100.,851.16216701,2000.,2000.]),
            N.array([-11.43878627,  1.00885436,16.88944849,-7.75064942,
                      -3.24734859,  8.3292758 ,-4.62759882, 6.11436238,
                      -1.57131022,-19.21339703, 0.51990376,-9.97359666]),2,1]
        tck_dy_12=[N.array([100.,100.,1172.275,2200.,2200.]),
            N.array([-100.,-100.,-100.,1075.197,2000.,2000.,2000.]),
            N.array([10.19347628,-17.26046227, 5.01814597,-8.08003115,
                      3.91293559,  8.86212983,-8.39262201, 1.11986956,
                    -18.76687171, -6.44986454,15.13866687, 2.27766693]),1,2]        
        #  put in weigths to pin down boresight
        #tck_dx_11 = [N.array([100.,100.,2200.,2200.]), N.array([-100.,-100.,2000.,2000.]),
        #    N.array([ -4.81093068,  -0.29231751,  -8.93141859,  15.59013236]),1,1]
        #tck_dy_11 = [N.array([100.,100.,2200.,2200.]),N.array([ -100.,-100.,2000.,2000.]),
        #    N.array([ -4.99692597,  -3.25416013,  19.47843113,  -9.14871286]),1,1]
        #tck_dx_12 = [N.array([100.,   100.,  2200.,  2200.]),
        #    N.array([ -100.,  -100.,  -100.,  2000.,  2000.,  2000.]),
        #    N.array([ 29.98237727, -55.81341527,  62.26860013, 
        #    3.65285845, 5.20894757,  11.06899946]),1,2]
        #tck_dy_12 = [N.array([  100.,   100.,  2200.,  2200.]),
        #    N.array([ -100.,  -100.,  -100.,  2000.,  2000.,  2000.]),
        #    N.array([-13.14184908,   8.60363512, -22.18512146,
        #    2.16964947, 20.4593104 , -20.44840889]),1,2]
        #tck_dy_22 = [N.array([  100.,   100.,   100.,  2200.,  2200.,  2200.]),
        #    N.array([ -100.,  -100.,  -100.,  2000.,  2000.,  2000.]),
        #    N.array([  11.13689074,  -74.16515394,   60.89105333,  -12.04378053,
        #    99.74817022, -123.31821113,   -3.38538815,  -56.71926279,
        #    75.65488355]),2,2]
        #tck_dx_22 = [N.array([  100.,   100.,   100.,  2200.,  2200.,  2200.]),
        #    N.array([ -100.,  -100.,  -100.,  2000.,  2000.,  2000.]),
        #    N.array([ -35.5519192 ,   74.7684609 ,   -1.26948618,   76.18803363,
        #    -152.14992411,   94.43166875,  -49.08716344,  128.37644766,
        #    -48.81027073]),2,2]
        # now load the dx,dy corrections to the model       
        for i in range(npnt):
            dx[i] = interpolate.bisplev(xp[i],yp[i],tck_dx_11)  
            #dx[i] = interpolate.bisplev(xp[i],yp[i],tck_dx_22)  
            dy[i] = interpolate.bisplev(xp[i],yp[i],tck_dy_12)
        #  shift in mean when not fixing the boresight anchor     
        dx = dx + 23.0
        dy = dy -  9.7
        # apply the correction to the simple zemax model     
        xp = xpin - dx
        yp = ypin - dy
        #  Fixing the correction to be zero at the original boresight anchor position leads to 
        #  no good correction (a correction worse than no correction.) Assuming that 
        #  the boresight anchor position was approximate and the adjustment from this solution
        #  gives a better anchor position, the new position is (1117.36,1027.22). The difference
        #  amounts to (-11.79,+4.96). 
        return xp, yp
                
      else:             
        if chatter > 4: print('UV Clocked version 1 calibration')
        X = xpin - 104
        Y = ypin - 78
        xp = X + 2.26132467e-05*X*X -1.17844816e-01*X +  90.1636849 + \
          1.16578350e-05*Y*Y  -6.78937589e-03*Y  -3.3658264
        yp = Y + 17.9031280 -1.58468408e-02*Y + 2.81777369e-03*X -2.8932900
        return xp+104,yp+78
        
  elif wheelpos == 1000:
     if chatter > 4: print('V nominal')
     X = xpin - 104
     Y = ypin - 78
     c1x = N.array([ -3.31706085e-06,  -2.12684998e-02,   2.30316863e+01])
     c2x = N.array([  4.12241297e-03,  -4.88680416e+00])
     c1y = N.array([ -3.06719634e-02,   3.32000692e+01]) # wrt yobs -3.15596678e-02,   3.41526717e+01])
     c2y = N.array([ -0.00877312,  7.7688439 ])
     xp = X + N.polyval(c1x,X)+N.polyval(c2x,Y)
     yp = Y + N.polyval(c2y,X)+N.polyval(c1y,Y)   
     return xp+104,yp+78
     
  elif wheelpos == 955:
     if chatter > 4: print('V clocked')
     X = xpin - 104
     Y = ypin - 78
     xp = xpin
     yp = ypin
     q = N.where((X < 2200.) & (Y < 2200.) & (X > -50.) & (Y > -50.))[0]
     tckx = [N.array([-50.,-50.,1100.825930,2200.,2200.]), \
             N.array([-50.,-50., 958.671647,2200.,2200.]), \
      N.array([-157.63428897,-130.94200779,-40.23248474,7.73189807,38.07534651, \
        61.38915254,184.63422561,176.9951507,200.7860891]), 1, 1] 
     tcky = [N.array([-50.,-50.,1131.6768029,2200.,2200.]), \
             N.array([-50.,-50.,1029.671610,2200.,2200.]), \
      N.array([-150.43784311,0.94833005,124.98862246,-130.98068849,9.05561256, \
       134.91166323,-144.02971497,24.23372216,112.36889942]), 1, 1]
     for i in q:  
       xp[i] = X[i] + interpolate.bisplev(X[i], Y[i], tckx) 
       yp[i] = Y[i] + interpolate.bisplev(X[i], Y[i], tcky)                 
     return xp+104,yp+78
  else:
     if chatter > 4: print('undefined wheelpos') 
     return xpin,ypin  

def _align_zemax(xpix,ypix,wheelpos,test=None,chatter=0):
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
            
      elif test == "simple":
            
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

def _align_zemax_field(xfield,yfield,wheelpos,test=None,chatter=0):
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
      return xfield, yfield


def prepDistortion(setno='nominal_3.8',order = 1, wav = 2, angle = -64.6):
   ''' prepare data at given order and wavelength for determining scale factor '''
   (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, ypix, 
         valid) = rdzemax(setno=setno,xdetoff=-49.0,ydetoff=+12.6)
   rx, ry = matchFieldToPix(xfield,yfield,order = order, wave = wave, angle = angle )
   wav = N.where(wave == 0.260)    
   dx = xpix[order,wav,:] - rx 
   dy = ypix[order,wav,:] - ry
   xp = (xpix[order,wav,:]).squeeze()
   yp = (ypix[order,wav,:]).squeeze()
   xf = (xfield[order,wav,:]).squeeze()
   yf = (yfield[order,wav,:]).squeeze()
   return xp, yps, dx, dy, xf, yf             

def interpolateDistortion(x,y):
   ''' 
   provides distortion map: 
   rectifies DET position to a uniform grid with a pixel scale of 0.5496"/pixel.
   
   input: unrectified DET position (x,y)
   output: rectified position (x,y)
   
   ''' 
   # get the calculated map 
   xpix, ypix, dx, dy, xf, yf = prepDistortion()
   # interpolate to find correction
   xrec, yrec = uvotmisc.interpgrid(x, y, xpix, ypix, dx, dy)
   # apply correction
   xrec = xrec + x
   yrec = yrec + y
   return xrec, yrec

#def DetToSky(None):
#   return None

def SkyToDet(ra_src,dec_src,ra_pnt,dec_pnt,pa_pnt):
   ''' convert SKY position to position on detector, 
   including distortion from UV grism. 
   
   input: source position ra_src, dec_src in decimal degrees
          pointing position ra_pnt, dec_pnt in decimal degrees
                   roll angle pa_pnt in degrees (decimal)
   output: detector position
   
   There is an additional rotation due to the angle XI 
   between north in the Satellite frame and north in the instrument, 
   which has not been determined yet. 
   '''
   # get the mapping from sky to det (zemax optical model; 2600, first order)
   xpix, ypix, dumx, dumy, xf, yf = prepDistortion()
   
   # distance on sky source to boresight
   ra_del  = ra_src  - ra_pnt
   dec_del = dec_src - dec_pnt
   ra_del  = (numpy.array(ra_del) ).squeeze()
   dec_del = (numpy.array(dec_dec)).squeeze()

   # rotate around the boresight to line up ra_del, dec_del with xsky, ysky (degrees)
   # which are the zemax inputs relative to the zemax boresight. 
   XI = -3.5 # rotation angle in degrees between Satellite and instrument.
   xi = (XI-pa_pnt)/180.*math.pi
   m11 = math.cos(xi)
   m12 = math.sin(xi)
   m21 = -math.sin(xi)
   m22 = math.cos(xi)
   matrix = numpy.array([[m11, m12],
                         [m21, m22]], dtype = numpy.float64)
   
   # define some arrays to fill:
   coord = numpy.zeros( 2,len(ra_del), dtype = numpy.float32)
   xsky = numpy.zeros( len(ra_del), dtype = numpy.float32)
   ysky = numpy.zeros( len(ra_del), dtype = numpy.float32)
   xdet = xsky
   ydet = ysky
   
   # find sky positions                  
   for k in range(len(ra_del)):                  
      coord[:,k] = numpy.dot(matrix, [ra_del[k],dec_del[k]])
      xsky[k] = coord[0,k]
      ysky[k] = coord[1,k]
      xdet[k], ydet[k] = uvotmisc.interpgrid(xsky[k],ysky[k],xf,yf,xpix,ypix)   
                         
   return xdet, ydet

def make_zemax_input():
   ''' 
   calculate a for rectangular grid of detector coordinates the 
   zemax input rays xsk, ysk for a given pixel scale factor.
   The detector coordinates are presumed to be the coordinates
   of the 260nm first order point in the grism spectra.
   
   This will give a grid of source field (sky) positions, that 
   would result in an undistorted DET grid when used on the 
   lenticular filters. So the sky positions will be a regular grid
   in degrees on the sky.
   
   Then the zemax model ray-tracing will predict where the regular
   spectra fall in DET coordinates, including therein the grism 
   distortion as predicted by Zemax 
   
   v1.0 20080822 NPMK(MSSL)
   '''   
   nx = 28 # 42
   ny = 28 # 42
   nn = nx*ny
   f = 75 # 50  STEP SIZE
   x1 = f*(numpy.arange(nx,dtype=numpy.float64))
   y1 = f*(numpy.arange(ny,dtype=numpy.float64))
   xdet = numpy.zeros( nn, dtype=numpy.float64).reshape(nx,ny) 
   ydet = numpy.zeros( nn, dtype=numpy.float64).reshape(nx,ny) 
   # the offset should place the centre of the calculated positions at the centre 
   # centre of the field positions (i.e., in the boresight ray) for a zemax scale 
   # factor of 0.963. The offset is opposite in sign to the shift that had to be
   # applied to the earlier zemax results. The expectation is that the boundaries
   # of the grid will now follow the edge of the detector
   for k in range(ny):
      xdet[k,:] = x1 + 49.0 + 77.
   for k in range(nx): 
      ydet[:,k] = y1 - 12.6 + 77.
   xsk, ysk = matchDetToField(xdet.flatten(),ydet.flatten())
   return xdet, ydet, xsk.reshape(nx,ny), ysk.reshape(nx,ny)  

   
def matchDetToField(xdet,ydet, wheelpos=200,order = 1, wave = 2, angle = 64.6, scale = 6550.4):
    '''
    Apply to the zemax det coordinates at 2600 A, first order, the coordinate 
    transformation to the zemax model input parameters xpix,ypix to overlay the 
    zemax field coordinates.    
    
    There is an additional distortion left which is minimal at the center. 
    The zemax field coordinates are in degrees on the sky.

    These coordinates are completely internal tot the zemax model. Since the zemax model is
    offset from the borepoint, a correction from actual DET-pix coordinates to the zemax model 
    coordinates needs to be made beforehand to the input xdet, ydet values. 
    Input values of xdet,ydet = 1092.5,1104.5 give the zero angle of the boresight.
    that means, that relative to the center the model boresight (xdet',ydet') = -8.0, +4.0, for 
    zemax offsets (xdetoff,ydetoff) set to zero. 
    
    '''
    # angle = 64.6
    # move offsets xdetoff, ydetoff from rdzemax 
    # the zemax X-position (uncorrected is 1023.5-977.5, while the observed position is 1023.5-928.5 )
    # so the correction needs to be made to line up the boresight. 
    # this depends on the wheel pos
    #if wheelpos == 200 : borex = 928.5 + 77 ; borey = 1002.7 + 77   
    #elif wheelpos == 160: 
    #   print ' matchDetToField. Error: no offset boresight set'  
    #   return
    #elif wheelpos == 0: 
    #   #borex = 977.53 + 77 ; borey = 990.09 + 77 
    borex = 1100.5 - 8.0 ; borey = 1100.5 + 4.0
    #   print 'origin: ',borex,borey       # for zemax model input    
    #else:
    #   print ' matchDetToField Error: no offset boresight set'    
    #   return   
    x = ( xdet.copy() - borex ) / scale
    y = ( ydet.copy() - borey ) / scale
    #print 'matchDetToField x,y = ',x,y
    #x = ( xdet.copy()  +8.0 - 1023.5 ) / scale   977.5- 928.5= 49.0 actual: -1015.5  x=0 for xdet=1015.5
    #y = ( ydet.copy()  -4.0 - 1023.5 ) / scale   990.0-1002.6=-12.6 actual: -1027.5  y=0 for ydet=1027.5
    # scale = 6550.4  which gives ~0.5496"/pixel
    # focallength = 3406 +/- 26 mm
    # focallength = scale * (pixel/mm) * (180/pi)
    # where (pix/mm) = 0.009075
    x = x.flatten()
    y = y.flatten()
    return uvotmisc.uvotrotvec(x, y, angle)
   
   
def matchFieldToPix(xfield,yfield,order = 1, wave = 2, angle = -64.6, scale = 6550.4):
    '''
    Apply to the zemax field coordinate a rotation, scaling, and shift 
    to give an approximate DET pixel coordinate at 2600, first order. 
    There is an additional distortion left which is minimal at the center. 
    The zemax field coordinates are in degrees on the sky.
    '''
    #angle = - 64.6
    # move offsets xdetoff, ydetoff from rdzemax to here?
    # xdetoff = - 49.0 ; ydetoff = 12.6
    xdetoff = 0 ; ydetoff= 0 
    rx , ry = uvotmisc.uvotrotvec(xfield[order,wave,:],yfield[order,wave,:], angle)
    # scale = 6550.4  which is 0.549585"/pixel converts between degrees and 
    #         pixels and is ~ 3600./0.549585
    # position [1023.5-8.0, 1023.5+4] is the centre of the zemax model distortion
    #         in a [2048,2048] grid.
    rx = scale*rx - 8.0 + 1100.5 - xdetoff
    ry = scale*ry + 4.0 + 1100.5 - ydetoff
    # in this case, the 
    # focallength = 3406 +/- 26 mm
    # focallength = scale * (pixel/mm) * (180/pi)
    # where (pix/mm) = 0.009075
    return rx, ry


def rdzemax(wheelpos, setno=None, scale=1.00, xdetoff=0.0,ydetoff=0.0, chatter=0):
   '''
   Read the zemax data 
   required input: wheelpos value 
     if not a legal wheelpos value, a default uv grism nominal 3.8 is selected.
      
   Call and output: 
   
   data = rdzemax( , scale=None)
   (wave, order, xfield, yfield, xcentroid, ycentroid,  xpix, ypix, valid) = data
   
   '''   
   __version__ = '2009-09-21 NPMK (MSSL)'
   # 20081121: changes to accomodate different wheelpos parameter
   # 20090921: added handling for Vgrism both modes
   
   import uvotgrism
   
   if chatter > 1: 
     print('zemax.  wheelpos =   ',wheelpos)
     print('zemax.  setno    =   ',setno)
     print('zemax.  scale    =   ',scale)
     print('zemax.  xdetoff  =   ',xdetoff)
     print('zemax.  ydetoff  =   ',ydetoff)
   zemaxmodel = getenv('ZEMAXMODEL')
   if zemaxmodel == '': zemaxmodel = getenv('PWD')
   #
   if setno == None: 
      if wheelpos == 200:
         setno = 'nominal_3.8_final'
         dir = zemaxmodel+'/nominal_3.8/'
         filename = 'new_fields.txt'
         
      elif wheelpos == 161: 
         setno = 'UGC_160' 
         dir = zemaxmodel+'/UGC_160/'
         filename =  'boresight_clocked160_4-1deg_corrected_all_fields.TXT' 
         
      elif wheelpos == 160: 
         setno = 'UGC_160' 
         dir = zemaxmodel+'/UGC_160/'
         filename = 'ug_clocked_4.1.txt'
         
      elif wheelpos == 1000:
         setno = 'VGN_1000' 
         dir = zemaxmodel+'/VGN_1000/'
         filename = 'rot6-9-allorders-allwaves-allfields.TXT'#'onaxis_grism_rot6-9_image65.TXT'  #'vg_nominal_4.1.txt'
         
      elif wheelpos == 955:
         setno = 'VGC_955'
         dir = zemaxmodel+'/VGC_955/'
         filename = 'rot6-9_clocked7-36_allorders-allwaves-allfields.TXT' #'onaxis_grism_rot6-9_image65_clocked7-36.TXT' #'vg_clocked_4.1.txt'
         
      else:
         print('Error in call zemax.rdzemax -  assuming nominal UV ')
         setno = 'nominal_3.8_final'
         
   elif setno == 'nominal_4.1':
         dir = zemaxmodel+'/nominal_4.1/'
         filename = 'new_fields.txt'
         
   elif setno == 'check160':
         dir = zemaxmodel+'/UGC_160/'
         filename = 'boresight_clocked160_4-1deg_corrected.TXT'
         
   elif setno == 'nominal_3.8':  
         dir = zemaxmodel+'/nominal_3.8/'
                 
   elif setno == 'VGN_1000_bs': 
         dir = zemaxmodel+'/VGN_1000/'
         filename = 'onaxis_grism_rot6-9_image65.TXT'  #'vg_nominal_4.1.txt'
         
   elif setno == 'VGC_955':
         dir = zemaxmodel+'/VGC_955/'
         filename = 'onaxis_grism_rot6-9_image65_clocked7-36.TXT' #'vg_clocked_4.1.txt'
                
   else: 
      print('setno = ',setno)
      print('valid values for setno are: nominal_3.8_final(default), test, UGC_160, VGN_1000,VGC_980')
      print('>>> setting to UV grism nominal 3.8deg ')
      setno = 'nominal_3.8_final'
      dir = zemaxmodel+'/nominal_3.8/'
   #   
   if chatter > 1: 
      print('zemax.  setno   = ',setno)
      print('zemax.    dir   = ',dir)
      print('zemax. filename = ',filename)
   # ================================================================================================
   if setno == 'nominal_3.8':
      dir = zemaxmodel+'/nominal_3.8/'
      suf = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V']
      NS = len(suf)
      nc=5
      nw=7
      nf=256
      n = nc*nw*nf
      
      xfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      yfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xpix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ypix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      valid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 

      for i in range(NS-2):
         file = dir+suf[i]+'.TXT'  
         (wave, order, txfield, tyfield, txcentroid, tycentroid,txpix, 
             typix, tvalid) = read_zemax4(file, nf=12, pixscale=scale)
         k1 = i*12
         k2 = (i+1)*12
         xfield   [:,:,k1:k2] = txfield[:,:,:]
         yfield   [:,:,k1:k2] = tyfield[:,:,:]
         xpix     [:,:,k1:k2] = txpix[:,:,:] + xdetoff
         ypix     [:,:,k1:k2] = typix[:,:,:] + ydetoff
         xcentroid[:,:,k1:k2] = txcentroid[:,:,:]
         ycentroid[:,:,k1:k2] = tycentroid[:,:,:]
         valid    [:,:,k1:k2] = tvalid[:,:,:]
   
      file = dir+suf[NS-1]+'.TXT'
      (wave, order, txfield, tyfield, txcentroid, tycentroid,txpix, 
            typix, tvalid) = read_zemax4(file, nf=4, pixscale=scale)
      i = NS-1   
      k1 = i*12
      k2 = i*12+4
      xfield   [:,:,k1:k2] = txfield
      yfield   [:,:,k1:k2] = tyfield
      xpix     [:,:,k1:k2] = txpix + xdetoff
      ypix     [:,:,k1:k2] = typix + ydetoff
      xcentroid[:,:,k1:k2] = txcentroid
      ycentroid[:,:,k1:k2] = tycentroid
      valid    [:,:,k1:k2] = tvalid
      data = (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, 
            ypix, valid)          
      return data     
      
   if setno == 'check160':
      nc=5
      nw=7
      nf=1
      n = nc*nw*nf
      
      xfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      yfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xpix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ypix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      valid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
   
      file = dir+filename
      data = (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, 
            ypix, valid) = read_zemax4(file, nf=nf, nw=nw, pixscale=scale) 
            #, wavset=numpy.array([0.19,0.21,0.235,0.26,0.285,0.32,0.35]) )
            # offset [-215, -41] to WR52 004 observation 
      return data     

   if ((setno == 'VGN_1000_bs') ^ (setno == 'VGC_955_bs')):
      nc=5
      nw=12
      nf=1
      n = nc*nw*nf
      
      xfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      yfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xpix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ypix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      valid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      file = dir+filename
      (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, 
            ypix, valid) = read_zemax4(file, nf=nf, nw=nw, pixscale=scale, chatter=chatter
            , wavset=numpy.array([0.27,0.29,0.32,0.35,0.39,0.42,0.46,0.50,0.55,0.60,0.65,0.70]) )
            # offset [-215, -41] to WR52 004 observation 
        
      if setno == 'VGN_1000':
         # rotate each calculation by -2.6 deg around the point at 4200A
         # loop over all points
         from uvotpy.uvotmisc import uvotrotvec
         from uvotpy import uvotspec
         for k in range(nf):
            xpiv = xpix[1,5,k]
            ypiv = ypix[1,5,k]
            x = xpix[:,:,k].reshape(nc*nw) - xpiv
            y = ypix[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-2.6)
            xp += xpiv  
            yp += ypiv  
            xpix[:,:,k] = xp.reshape(nc,nw)
            ypix[:,:,k] = yp.reshape(nc,nw)
         dxy =  uvotgrism.boresight('vg1000') - numpy.array([1098.78,1101.57])
         xpix = xpix + dxy[0] #  -98.781
         ypix = ypix + dxy[1] #   +5.330
      else:
         dxy = uvotgrism.boresight('vc955') - numpy.array([1204.79,1168.20])
         xpix = xpix + dxy[0] #  -61.292
         ypix = ypix + dxy[1] # -126.202
      return (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, ypix, valid)     

   if ((setno == 'VGN_1000') ^ (setno == 'VGC_955')):
      from uvotpy import uvotspec
      nc=5
      nw=12
      nf=784
      n = nc*nw*nf
      
      xfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      yfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xpix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ypix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      valid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      file = dir+filename
      (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, \
            ypix, valid) = read_zemax4(file, nf=nf, nw=nw, pixscale=scale, chatter=chatter, \
            wavset=numpy.array([0.27,0.29,0.32,0.35,0.39,0.42,0.46,0.50,0.55,0.60,0.65,0.70]) ) 
       
      if setno == 'VGN_1000':
         print('\n READING ZEMAX MODEL SETNO VGN_1000 2009/09/21 \n')
         # rotate each calculation by -2.6 deg around the point at 4200A
         # loop over all points
         from uvotpy.uvotmisc import uvotrotvec
         for k in range(nf):
            xpiv = xpix[1,5,k]
            ypiv = ypix[1,5,k]
            x = xpix[:,:,k].reshape(nc*nw) - xpiv
            y = ypix[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,0.0)
            xp += xpiv  
            yp += ypiv  
            xpix[:,:,k] = xp.reshape(nc,nw)
            ypix[:,:,k] = yp.reshape(nc,nw)
         #dxy = numpy.array([0,0])       
         dxy = uvotgrism.boresight('vg1000') - numpy.array([1098.78,1101.57]) 
         # offset derived from boresight model vgn1000_bs detector coordinate
         xpix = xpix + dxy[0] #  -98.781
         ypix = ypix + dxy[1] #   +5.330
      else:
         print('\n READING ZEMAX MODEL SETNO VGC_955 2009/09/21 \n')
         dxy = uvotgrism.boresight('vc955') - numpy.array([1204.79,1168.20]) 
         # offset derived from boresight model vgn1000_bs detector coordinate
         xpix = xpix + dxy[0] #  -61.292
         ypix = ypix + dxy[1] # -126.202
      return (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, ypix, valid)      
       
   else: 
      nc=5
      nw=16
      nf=784
      n = nc*nw*nf
      
      xfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      yfield = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xpix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ypix = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      valid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
      file = dir+filename  
      (wave, order, xfield, yfield, xcentroid, ycentroid,xpix, 
             ypix, valid) = read_zemax4(file, nf=nf, nw=nw, pixscale=scale,chatter=chatter)
      if wheelpos == 160:
         # rotate each calculation around the point at 2600A
         # loop over all points
         if chatter > 0: print("rotating the model points 65 deg.")
         from uvotpy.uvotmisc import uvotrotvec
         for k in range(nf):
            # pixels
            xpiv = xpix[1,3,k]
            ypiv = ypix[1,3,k]
            x = xpix[:,:,k].reshape(nc*nw) - xpiv
            y = ypix[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-65.0)
            xp += xpiv-18.63  # -175.9+157.27
            yp += ypiv+111.94 # -102.2+214.14
            xpix[:,:,k] = xp.reshape(nc,nw)
            ypix[:,:,k] = yp.reshape(nc,nw)
            
            # centroids
            xpiv = xcentroid[1,3,k]
            ypiv = ycentroid[1,3,k]
            x = xcentroid[:,:,k].reshape(nc*nw) - xpiv
            y = ycentroid[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-65.0)
            xp += xpiv  # no shift
            yp += ypiv  # no shift
            xcentroid[:,:,k] = xp.reshape(nc,nw)
            ycentroid[:,:,k] = yp.reshape(nc,nw)
            
      if wheelpos == 1000:
         # rotate each calculation by -2.6 deg around the point at 4200A
         # loop over all points
         from uvotpy.uvotmisc import uvotrotvec
         for k in range(nf):
            xpiv = xpix[1,5,k]
            ypiv = ypix[1,5,k]
            x = xpix[:,:,k].reshape(nc*nw) - xpiv
            y = ypix[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-2.6)
            xp += xpiv  
            yp += ypiv  
            xpix[:,:,k] = xp.reshape(nc,nw) - 98.78
            ypix[:,:,k] = yp.reshape(nc,nw) +  5.33
            
      notvalid = N.where( valid == -1 )
      xfield[notvalid] = None
      yfield[notvalid] = None
      xpix[notvalid] = None
      ypix[notvalid] = None
          
      return (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, ypix, valid)               
          
   
   

def read_zemax4(file, nf=9, nw=9, nc=5, pixscale=1.00, chatter=0, 
    setwav = False, nskip=5, wavset=None):
   ''' 
    This routine reads the zemax files with the fourth version of the 
    format. nf = number of field positions, nw = number of wavelengths 
    points, nc = number of orders (0 order = 0, first order = 1, 
    second order = 2, third order = 3, and minus-one order = 4 )
    
    The pixel scale is adjusted before output using the pixscale keyword.
    Default pixscale = [reset to 1.00] ]was 0.963; perhaps 0.960 is better
    
    output is 
     
   '''
   __version__ = '2009-02-24 NPMK (MSSL)'

   n = nf*nw*nc
      
   if wavset == None:
      if nw < 7.5:
         wave = numpy.array(([0.19,0.214,0.26,0.33, 0.4,0.45,0.55])[0:nw]) 
      elif nw == 16: 
         wave = numpy.array(([0.19,0.21,0.235,0.26,0.285,0.32,0.35,0.38,0.42,0.46,0.5,0.54,0.58,0.62,0.66,0.7]))
   else: wave = wavset 
   
   if chatter > 4: 
      print('wavelengths set = ',wave)
      print('dimensions fields nf=',nf,'  waves nw=',nw,' orders nc=',nc)  
       
   order = numpy.array( ([0,1,2,3,-1])[0:nc] )
   
   xfield    = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 100.
   yfield    = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 100.
   xcentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 100.
   ycentroid = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 100.
   xpix      = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 10000.
   ypix      = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) + 10000.
   valid     = numpy.zeros( n,dtype=float ).reshape(nc,nw,nf) 
   
   f = open(file)
   recs = f.readlines()
   f.close()
   
   if chatter > 2:
      for i in range(nskip):
         print(recs[i])
   
   i = nskip
   
   for iw in range(nw):  
      for ic in range(nc):
            i += 1
            card = recs[i]
            if chatter > 3: print('from input for iw=',iw,'  ic=',ic,'  i=',i,' card= ',card) 
            i += 1
            if chatter > 3: print('from input for iw=',iw,'  ic=',ic,'  i=',i,' next card= ',recs[i])
            dum1, w2, dum2, cn = (card.split())[0:4]
            
            w2 = float(w2)
            cn = float(cn)
            cn = int(cn)
            if chatter > 1: print('read in: w2, cn = ',w2, cn)   
            if cn == -9:
               return -1
            wn = -1
            if setwav: 
               wave[iw] = w2
               wn_ = wn = iw
            else:    
               wn_ = numpy.where( w2 == wave )
               if (w2 == wave[wn_]): 
                  wn = wn_[0]
               else: 
                  print('error w2 not valid')
                  break
               if wn == -1 : 
                  print('error finding wavelength in wave ( wn )')
                  break
            if chatter > 1:
              print('config number=',cn,' wave number=',wn,'  wave=', w2)
            
            if nc == 1:
               cn = 0
               if chatter > 2: print('only one order to be read in: reset cn =' , cn)
            else: 
               cn=int(cn-1)
                  
            i += 1      # line of ********** 
            for j in range(nf):
               i += 1   # blank
               rec = recs[i]  # field positions
               if chatter > 4: print('rec:',rec)
               i += 1
               x3 , x4 = rec.split() 
               if chatter > 1: print('indexes cn=',cn,'  wn=',wn,'  j= ',j,'  i=',i)           
               xfield[cn,wn,j] = float(x3)
               yfield[cn,wn,j] = float(x4)
               if chatter > 1: print('field x,y = ',x3,x4)
               rec = recs[i]  # det coords or "No rays made it .."
               i += 1     # next record (blank)
               if (rec.split())[0] != 'No':
                  #rec = recs(i)
                  x3, x4 = rec.split()
                  if chatter > 2: print('centroid x,y =', x3, x4)
                  xcentroid[cn,wn,j] = float(x3)
                  ycentroid[cn,wn,j] = float(x4)
                  xpix [cn,wn,j] = float(x3)/0.009075+1100.5
                  ypix [cn,wn,j] = float(x4)/0.009075+1100.5
               else:
                  valid[cn,wn,j] = -1
                  if chatter > 1: print(' No valid solution here')
            i += 2 
   i += 1
   rec = recs[i]
   if chatter > 2: print(rec)    
      
   xpix = xpix * pixscale
   ypix = ypix * pixscale
      
   return (wave, order, xfield, yfield, xcentroid, ycentroid, xpix, ypix, valid)
         
   
def read_zemax_flux(wheelpos=None,file=None, rotate=False,chatter=0, ):
   ''' 
    This routine reads the zemax flux file. 
    
    output: 
    wave is wavelength in micrometers.
    
    config relates to the order (0 order = 0, first order = 1, 
    second order = 2, third order = 3, and minus-one order = 4 )
        
    input rays are defined by the (xfield, yfield) coordinates [deg]
    the position on the detector is given by [xcentroid,ycentroid] [mm]
    with a certain offset and scale (see calibration). 
    
    160:/Volumes/users/Users/kuin/zemaxmodel/UGC_160/flux/fields_x_y_784_no_image_rotation_flux_results_nominal_v1.txt
    200:/Volumes/users/Users/kuin/zemaxmodel/nominal_3.8/flux/new_field_points2_flux_results_v5_nominal.txt) 
    955:/Volumes/users/Users/kuin/zemaxmodel/VGC_955/fields_x_y_784_Visible-grism_50-6_rotation_flux_results_v1.txt
   1000:/Volumes/users/Users/kuin/zemaxmodel/VGN_1000/fields_x_y_784_Visible-grism_57-9_rotation_flux_results_v2.txt

   '''
   __version__ = '2012-04-11 NPMK (MSSL)'

   from uvotpy.uvotmisc import uniq
   import numpy as np
   
   if (wheelpos == None) & (file == None): 
      print("Insufficient arguments: read_zemax_flux(wheelpos=None,file=None, rotate=False,chatter=0,")
   
   if (wheelpos != None):
      if   wheelpos ==  160: 
         #file='/Volumes/users/Users/kuin/zemaxmodel/UGC_160/flux/new_field_points2_flux_results_v3.txt' changed dec 20, 2012
         file='/Volumes/users/Users/kuin/zemaxmodel/UGC_160/flux/fields_x_y_784_no_image_rotation_flux_results_nominal_v1.txt'
      elif wheelpos ==  200: 
         file='/Volumes/users/Users/kuin/zemaxmodel/nominal_3.8/flux/new_field_points2_flux_results_v5_nominal.txt' 
      elif wheelpos ==  955: 
         file='/Volumes/users/Users/kuin/zemaxmodel/VGC_955/fields_x_y_784_Visible-grism_50-6_rotation_flux_results_v1.txt'
      elif wheelpos == 1000: 
         file='/Volumes/users/Users/kuin/zemaxmodel/VGN_1000/fields_x_y_784_Visible-grism_57-9_rotation_flux_results_v2.txt'

   # read header lines 
   print("opening "+file)
   f = open(file)
   
   orders = [0,1,2,3,-1]
   wave = []
   config = []
   nf = -1
   fieldx = []
   fieldy = []
   flux = []
   xcen = []
   ycen = []
   for rec in f:
      #if rec[0] == "#":
      if rec[1:5] == "wave":
            dum1, w, dum2, c = rec.split()
            wave.append(float(w))
            config.append(float(c))
   f.close()
   if chatter > 2:
      print("wave:",wave)
      print("config:",config)
      print(" ")
   
   uniq_wave = uniq(wave)
   nwav = len(uniq_wave)
   uniq_config = (uniq(config))
   uniq_config.sort()
   ncon = len(uniq_config)   
   if chatter > 1:
      print("uniq_wave:",uniq_wave)
      print("uniq_config:",uniq_config)
      print(" ")
   print("zemax calculation for %i wavelengths and %i orders \n" % (nwav,ncon))
   
   f = open(file)
   # loop over all config and all wavel and get [fieldx,fieldy, centroidx, centroidy, flux 
   # skip the header
   happy = True
   while happy:
      rec = f.readline()  
      happy = (rec[0] == '#') & (rec[1:5] != 'wave') 
      if chatter > 0: print('header: %s'%rec)

   happy = True 
   while happy: 
      if (rec[1:5] == 'wave'): 
          dum1, w, dum2, c = rec.split()
          if chatter > 0:  print('2,wave = ',w,' config = ',c)
      elif rec[0] != '#':
          happy = False
          oeps = True
          break
      rec = f.readline()
      if chatter>2 : print(rec)   
   
   # read first data block
   happy = True
   while happy:   
      if oeps: 
         oeps = False
         rec1 = rec
      else:    
         rec1=f.readline()
      if len(rec1) < 4:
         print("length record < 4")
         happy = False
         break  
      else:       
         rec2=f.readline()
         if chatter > 3: 
            print("1:%s2:%sssizes:%i,%i,field#=%i" % (rec1,rec2,len(rec1),len(rec2),nf+1))
         fx1,fy1 = rec1.split()
         xcen1,ycen1,flx1 = rec2.split()
         fieldx.append(float(fx1))
         fieldy.append(float(fy1))
         flux.append(float(flx1))
         nf += 1
         xcen.append(float(xcen1))
         ycen.append(float(ycen1))
      
   #    
   nf += 1
   print("zemax calculation for number of field points = %i [ wavelength %s and order %s]\n" % (nf,w,c))
   
   # initialize arrays
   
   n = nf*nwav*ncon
   
   xfield    = numpy.zeros( n,dtype=float ).reshape(ncon,nwav,nf) 
   yfield    = numpy.zeros( n,dtype=float ).reshape(ncon,nwav,nf) 
   xcentroid = numpy.zeros( n,dtype=float ).reshape(ncon,nwav,nf) + 1e4
   ycentroid = numpy.zeros( n,dtype=float ).reshape(ncon,nwav,nf) + 1e4
   fluxes    = numpy.zeros( n,dtype=float ).reshape(ncon,nwav,nf) 
   valid     = numpy.zeros( n,dtype=bool ).reshape(ncon,nwav,nf) 
         
   xfield[0,0,:] = np.array( fieldx )
   yfield[0,0,:] = np.array( fieldy )
   fluxes[0,0,:] = np.array( flux   )
   xcentroid[0,0,:] = np.array( xcen )
   ycentroid[0,0,:] = np.array( ycen )
   
   skip = f.readline()
   if chatter > 1: print("0:"+skip)

   # the rest 
   for j in range(1,nwav*ncon):
   #for j in range(1,5):
      if chatter > 1: print("block number :",j)
      # read header + 1
      rec = f.readline()
      if chatter > 1: print("1:"+rec)
      skip = f.readline()
      if chatter > 1: print("2:"+skip)
      skip = f.readline()
      if chatter > 1: print("3:"+skip)
      if chatter > 1:
         print(50*"-")
         print("rec:",rec.split())
         print("   :",skip)
      try: 
        dum1,w,dum3,dum4 = rec.split()
        c = float(dum4)
        if wave[j] != float(w): 
           print("header #%i not consistent with wave %s :\n\t%s" % (j,wave[j],rec))
      except:
         print("trouble1?\nrec:%sskip:%s"%(rec,skip))
         print(dum1, w, dum3, dum4, " >>> ", c)
         print("returning")
         return
        
      # reset lists
      fieldx = []
      fieldy = []
      flux = []
      xcen = []
      ycen = []
      happy = True
      
      while happy:   
        rec1=f.readline()
        rec2=f.readline()
        if len(rec1) > 4: 
           if chatter > 4: print("_1:"+rec1, end=' '); print("_2:"+rec2, end=' ');print("j = ",j,"  nf = ",nf)
           try:
              fx1,fy1 = rec1.split()
              xcen1,ycen1,flx1 = rec2.split()
              fieldx.append(float(fx1))
              fieldy.append(float(fy1))
              flux.append(float(flx1))
              nf += 1
              xcen.append(float(xcen1))
              ycen.append(float(ycen1))
           except:
              print("trouble2: "+20*"=")
              print("1:%s"%rec1)
              print("2:%s"%rec2)
              print(40*"=")
              return
        else: happy=False      
        
      # copy to the output arrays
      kw = np.searchsorted(uniq_wave, w) - 1
      kc = np.searchsorted(uniq_config,c) 
      if chatter > 1: 
         print("block %i, index(config,wave)=(%i,%i), sum(flux)=%f, number=%i " % (j,kc,kw,(np.array(flux)).sum(),nf))

      xfield[kc,kw,:] = np.array( fieldx )
      yfield[kc,kw,:] = np.array( fieldy )
      fluxes[kc,kw,:] = np.array( flux   )
      xcentroid[kc,kw,:] = np.array( xcen )
      ycentroid[kc,kw,:] = np.array( ycen )
      xpix = xcentroid/0.009075+1100.5
      ypix = ycentroid/0.009075+1100.5
      valid = (xcentroid < 3000.) & (ycentroid < 3000.)
         
   f.close()        
   #return uniq_wave, orders, xfield, yfield, xcentroid, ycentroid, fluxes, valid
   if rotate :
         # wheelpos == 160 ! only for model v3
         # rotate each position calculation around the point at 2600A
         # loop over all points
         nc,nw,nf = xcentroid.shape
         if chatter > 0: print("rotating the flux model points 65 deg.")
         from uvotpy.uvotmisc import uvotrotvec
         for k in range(nf):
            # pixels
            xpiv = xpix[1,3,k]
            ypiv = ypix[1,3,k]
            x = xpix[:,:,k].reshape(nc*nw) - xpiv
            y = ypix[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-65.0)
            xp += xpiv-18.63  # -175.9+157.27
            yp += ypiv+111.94 # -102.2+214.14
            xpix[:,:,k] = xp.reshape(nc,nw)
            ypix[:,:,k] = yp.reshape(nc,nw)
            
            # centroids
            xpiv = xcentroid[1,3,k]
            ypiv = ycentroid[1,3,k]
            x = xcentroid[:,:,k].reshape(nc*nw) - xpiv
            y = ycentroid[:,:,k].reshape(nc*nw) - ypiv
            xp, yp = uvotrotvec(x,y,-65.0)
            xp += xpiv  # no shift
            yp += ypiv  # no shift
            xcentroid[:,:,k] = xp.reshape(nc,nw)
            ycentroid[:,:,k] = yp.reshape(nc,nw)
         #xcentroid[not valid] = np.nan
         #ycentroid[not valid] = np.nan  
         return uniq_wave, orders, xfield, yfield, xcentroid, ycentroid, xpix, ypix, fluxes, valid    
   else:
         #xcentroid[not valid] = np.nan
         #ycentroid[not valid] = np.nan  
         #xpix[not valid] = np.nan
         #ypix[not valid] = np.nan  
         q = (xpix == 1100.5) & (ypix == 1100.5)
         xpix[q] = None
         ypix[q] = None
         fluxes [q] = None
         valid[q] = False
         uniq_wave = np.asarray(uniq_wave)
         orders    = np.asarray(orders)      
         return uniq_wave, orders, xfield, yfield, xcentroid, ycentroid, xpix, ypix, fluxes, valid
      
    
      
def simple_disp(wheelpos):
   ''' provides the dispersion coeficients from the original calibration by 
       Sally, Wayne 
       These are based on a distance in pixels based on the zeroth order position
       The distance of the zeroth order anchor to the anchor point in first order
       at 420nm is 669pix [wheelpos=1000] and ?? pix [wheelpos=955]
       '''
   if wheelpos == 1000: 
        return N.array([0.00130680,4.05833,898.253])
   elif wheelpos == 955:
        return N.array([0.00159532,3.85725,1016.56])      
   else: 
        return N.array([0,0,0])          
     
def plotcont(Z,img=None,dxy=[-800,+530],iord=0,figno=1,levels=[0.05,0.1,0.2,0.5,0.9]):
    ''' plot contours zeroth order '''
    import pylab as plt
    import numpy as np
    wave, nord, xf, yf, xcent, ycent, xpix, ypix, flux, valid = Z  
    plt.figure(figno)
    
    if img != None:
       x = img.mean()
       plt.imshow(img,vmin=0.001*x,vmax=x)
    kw = len(wave)
    for i in range(kw): 
       xp = np.zeros(28) - 999.9
       yp = np.zeros(28) - 999.9
       print("wave = ",wave[i])
       
       for i1 in range(28):
          xp1 = xpix[iord,i,:].reshape(28,28)[:,i1]
          xp1 = xp1.flatten()
          q2 = np.where(xp1 != 1100.5)
          if len(q2[0]) > 0: 
             xp[i1] = xp1[q2[0]].mean()
        
       q1 = np.where(xp != -999.9)
       coef = np.polyfit(q1[0],xp[q1[0]],1)
       q1 = np.where(xp == -999.9)
       xp[q1[0]] = np.polyval(coef, q1[0]) 
       
       for i1 in range(28):
          yp1 = ypix[iord,i,:].reshape(28,28)[i1,:]
          yp1 = yp1.flatten()
          q2 = np.where(yp1 != 1100.5)
          if len(q2[0]) > 0: 
             yp[i1] = yp1[q2[0]].mean()
       
       q1 = np.where(yp != -999.9)
       coef = np.polyfit(q1[0],yp[q1[0]],1)
       q1 = np.where(yp == -999.9)
       yp[q1[0]] = np.polyval(coef, q1[0]) +dxy[1]
       xp += dxy[0]
       yp += dxy[1]
       print("xp: ", xp)
       print("yp: ",yp)
       
       C = plt.contour(xp,yp,flux[iord,i,:].reshape(28,28), levels=levels)
   
       
def interp_model(wheelpos):
    '''get a fit to the zemax model '''      
    uniq_wave, orders, xfield, yfield, xcentroid, ycentroid, xpix, ypix, fluxes, valid\
    =  read_zemax_flux(wheelpos=wheelpos)
    xpix1 = 0.94* xpix[1,:,:] + 300.0
    ypix1 = 0.94* ypix[1,:,:] +   0.0
   
    
def fitflux(wheelpos, err = None,chatter=0):
   from . import zemax, uvotcal
   import numpy as np  
   
   Z = uvotcal.get_calibrated_zemax_flux(wheelpos,chatter=chatter)
   wave, orders, xf, yf, xc, yc, xpix, ypix, flux, v = Z 
   
   wave = 1.0e4*np.array(wave)
   A = []
   B = []
   w1 = []
   w2 = []
   for i in range(28*28):
      q = np.isfinite(flux[1,:,i])
      zmxmodel = flux[1,:,i]
      if err == None: err = 0.001
      y = _fitfluxsub0(wave[q], zmxmodel[q], err,chatter=chatter)
      A.append(y[0])
      B.append(y[1])
      w1.append(y[2])
      w2.append(y[3])
   return A,B,w1,w2, wave, xf[1,:,:],yf[1,:,:],xc[1,:,:],yc[1,:,:],xpix[1,:,:],ypix[1,:,:],flux[1,:,:],v[1,:,:]     


def _fitfluxsub0(wave, zmxmodel, err, chatter=0):
   '''  
   Fit the flux from the zmxmodel with a simple function defined as follows:
   
     flx = A  for w < w1
     flx = A *sf(w,w2-w1) for w1=< w < w2
     flx = 0 for w2 < w
     
   which is a smoothed step function.  The fit can also be done to a real model. 
   err = RMS (1-sigma) errors in zmxmodel.
      
   return the parameters (A,w1,w2)
   
   related functions
     fitted_flux(A,B,w2,w1) interpolating function (interpolate.interp1d) returns the approximate model flux(w).
     _fitfluxsub1() difference function called by Levenberg-Marquardt mpfit.mpfit() 
         
   '''     
   import mpfit
   from numpy import array, arange,transpose, where, abs, min, zeros, atleast_1d, atleast_2d, sqrt
   
   # initial condition
   A = 1.0
   B = 0.0
   w1 = 2000.
   w2 = 6500.
   # define the variables for the function 'myfunct'
   fa = {'x':wave,'y':zmxmodel,'err':err}
   # define the parameters 
   p0 = (A,B,w1,w2)
   parinfo=[
   {'limited': [1,1],   'limits' : [ 0.0,    2.0],'value':     A,   'parname':  'A'    },
   {'limited': [1,1],   'limits' : [ 0.0,    2.0],'value':     B,   'parname':  'B'    },
   {'limited': [1,1],   'limits' : [1650.,7000.0],'value':    w1,   'parname': 'w1'    },
   {'limited': [1,1],   'limits' : [1650.,7000.0],'value':    w2,   'parname': 'w2'    },
   ]
   # call non-linear least-squares fit program (mpfit)
   Z = mpfit.mpfit(_fitfluxsub1,p0,functkw=fa,parinfo=parinfo,quiet=True)
   if (Z.status <= 0): 
      print('uvotgrism.runfit3.mpfit error message = ', Z.errmsg)
      print("parinfo has been set to: ") 
      for par in parinfo: print(par)
      A, w1, w2 = None, None, None
   elif (chatter > 0):   
      print("\nparameters and errors : ")
      for i in range(len(Z.params)): print("%10.3e +/- %10.3e\n"%(Z.params[i],Z.perror[i]))   
   #return Z
   A,B,w1,w2 = Z.params[:4]           
   return A,B,w1,w2

def _fitfluxsub1(p,fjac=None, x=None, y=None, err=None):
   import numpy as np
   
   (A,B,w1,w2) = p 
   model = fitted_flux(x,A,B,w1,w2)                 
   status = 0
   return [status, (y-model)/err]   
   
def fitted_flux(w,A,B,w1,w2):
   ''' 
   return the model fit to the flux for all wave lengths w

   '''
   import numpy as np
   w_ = np.asarray(w)
   f = w_*0.   
   f[ w_ < w1 ] = A
   f[ w2 < w_ ] = B
   q = (w_ >= w1) & (w_ <= w2)
   wm = 0.5*(w1+w2)
   f[q] = 0.5*(A-B)*(1+np.sin( -(w_[q]-wm)/(w2-w1)*np.pi )) + B 
   return f
   
   
   
