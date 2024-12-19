# -*- coding: iso-8859-15 -*-
#
# This software was written by N.P.M. Kuin (Paul Kuin) 
# Copyright N.P.M. Kuin 
# All rights reserved
# This software is licenced under a 3-clause BSD style license
# 
#Redistribution and use in source and binary forms, with or without 
#modification, are permitted provided that the following conditions are met:
#
#Redistributions of source code must retain the above copyright notice, 
#this list of conditions and the following disclaimer.
#
#Redistributions in binary form must reproduce the above copyright notice, 
#this list of conditions and the following disclaimer in the documentation 
#and/or other materials provided with the distribution.
#
#Neither the name of the University College London nor the names 
#of the code contributors may be used to endorse or promote products 
#derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
#THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
#PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
#EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
#PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
#OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
#WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
#OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
#ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from __future__ import division
__version__ = "1.0.0"

###################
# ZEMAX PSF 160 order 1 results
#
# 1.0.0 2015-01-30 moved from uvotcal.py 
#            


class PSF:

 """
 combines the software related to the UVOT grism PSF
 
 The routine readZemaxPSF to read the Zemax PSF data created originally
 for the centre of the uv clocked grism position. 
 
  
 """

 def __init__(self):
    self.status = 0

 def _readZemaxPSF(mydir = '/Users/kuin/zemaxmodel/PSF_160/data/', 
      lsffile='zemaxlsf0160_003.fit', 
      slit=[3.0], sigma=2.5, 
      version = 0.2,
      smooth = False,
      kernel = (32,32),
      clobber=True):
   '''
   Read the zemax PSF files for the 160 wheel position / first order / centre detector 
   Write the LSF file by summing the PSF in the cross-dispersion direction 
   
   Parameters
   ==========
   mydir : path
      location of the zemax model output
   lsffile : path
      name of output lsf file 
   slit : array
      polynomial coefficients (sigcoef1_x) for width slit of the first order (pixels) 
   sigma : float
      how many slit widths to include in cross-dispersion sum 
   version : x.x.x
      version number of the LSF file   
   smooth : bool
      apply boxcar smoothing on initial PSF to simulate detector smoothing
      1 pixel kernel = (36,36)    
   clobber : int
      verbosity 0...5 (0: None; 5: maximal)   
      
      
   half pixel size Ebins, with a boundary at the centre of the PSF.   
   
   '''
   
   import numpy as np
   import scipy.ndimage as ndimage
   from scipy.ndimage import convolve as conv
   from astropy.io import fits
   from uvotio import angstrom2kev
   import datetime
   from stsci.convolve import boxcar
   from uvotgetspec import singlegaussian
   
   # expand gaussian kernel
   px = singlegaussian(np.arange(4*kernel[0])-2*kernel[0],1.0,0.,kernel[0])
   py = singlegaussian(np.arange(4*kernel[1])-2*kernel[1],1.0,0.,kernel[1])
   pg = np.outer(px,py) # normalised gaussian 
   pg = pg/pg.sum()
   
   files = ('0-210um-1storder.TXT','0-235um-1storder.TXT','0-260um-1storder.TXT',\
           '0-2850um-1storder.TXT','0-320um-1storder.TXT','0-350um-1storder.TXT',\
            '0-380um-1storder.TXT','0-420um-1storder.TXT','0-460um-1storder.TXT',\
            '0-500um-1storder.TXT','0-540um-1storder.TXT','0-580um-1storder.TXT',\
            '0-620um-1storder.TXT','0-660um-1storder.TXT','0-700um-1storder.TXT')
	    
   #filemetadata = {'wave':0.,'spacing':0.250,'area':[32,32],'strehl_ratio':0.87, 
   #    'pupil_grid_size':[256,256],'image_grid_size':[128,128],'centre':[65,65],
   #    'centre_coord':[0.8,0.6]}
   filemetadata = {'wave':0., 'image_grid_size':[128,128],
                   'centre':[65,65], 'centre_coord':[0.8,0.6]}
   
   # data spacing is the same in each file, but for the longer wavelengths the size of 
   # the arrays is larger.
       
   metadata = list()    
   mm = len(files)
   datalist = list()
   profile=list()
       
   # units: wave(um), spacing(um), area(um^2),     	    

   i = -1
   for file in files:
      f = open(mydir+file)
      i += 1	    
      # Read header
      for i in range(8):
         m = f.readline()
      m = f.readline() ; filemetadata['wave'] = float(m.split()[0])
      m = f.readline() #; filemetadata['spacing'] = float(m.split()[3])
      m = f.readline() #; filemetadata['area'] = [float(m.split()[3]),float(m.split()[5])]
      m = f.readline() #; filemetadata['strehl_ratio'] = float(m.split()[2])
      m = f.readline() #; filemetadata['pupil_grid_size'] = [int(m.split()[3]),int(m.split()[5])]
      m = f.readline() ; filemetadata['image_grid_size'] = [int(m.split()[3]),int(m.split()[5])]
      m = f.readline() ; filemetadata['centre'] = [int((m.split()[3]).split(',')[0]), int(m.split()[4])]
      m = f.readline() ; filemetadata['centre_coord'] = [float((m.split()[2]).split(',')[0]),float(m.split()[3])]
      metadata.append(filemetadata.copy())
      
      # skip 2 lines
      m = f.readline()
      m = f.readline()

      # read data
      sz1 = filemetadata['image_grid_size'][0]
      sz2 = filemetadata['image_grid_size'][1]
      data = np.zeros(sz1*sz2,dtype=float).reshape(sz1,sz2)
      for row in range(sz1):
         m = f.readline()
	 datrow = m.split()
         for k in range(sz2): data[row,k] = float(datrow[k])
      
      if smooth: 
         data = conv(data,pg)
      datalist.append(data.copy())

      f.close()
      
   # one sub-pixel = 0.009075 mm = 9.075 um
   pix = 9.075  # pixel size in um  
   dx = 0.25    # resolution model in um
   # => approximate a pixel by 9/0.25 = 36 model data points (good to ~ 1%) 
   # so one pixel = ~36 points  
   # typical slit width (sigma) ~ 3 pix 
   # 5 sigma width is used for count rate (default) which comes to 5x3x36 = 540 grid points 
   # which means sum ~270 grid points on each side of spectrum (or however many there are provided)
   # if slit coefficients are given, 
   #     sigma = polyfit(slit, pixdist)
   #     npoints_sideways = 5*sigma*36/2 
   # where
   #     pixdist = |centre_coord[i]-centre_coord[k]|*1000./9.075 (centre_coor in mm?)
   #     
      
   m = len(metadata)
   wav = np.zeros(m,dtype=float) 
   centre_coord = (np.zeros(2*m,dtype=float)).reshape(2,m)
   print centre_coord.shape
   for i in range(m): 
       wav[i] = metadata[i]['wave']
       centre_coord[:,i] = metadata[i]['centre_coord']
   try:    
      centre_coord -= centre_coord[:,wav == 0.260]    
   except:
      centre_coord -= centre_coord[:,2]    
      print "problem finding anchor - assuming 3rd position in array"
      
   # distance from anchor in pixels   
   dist = np.sqrt((centre_coord**2).sum(0))*1e3/pix 
   print "distances from anchor \n   w    dist  "
   for (w,d) in zip(wav,dist):
       print "%6.1f  %6.1f"%(w*1e4,d )
   
   # prepare linespread, defined as a normalized profile relative to the 
   # wavelength/energy for RMF file input. This will be saved to a file. 
   # units are half-pixels (within 1%) 
      
   # implicit channels are 1 pix wide 
   # rotate data on angle near 0.26um found from nearest neighbors (centre-coord) 
   coef_ang = np.array([ 24.19064932, -42.98613931,  27.24415296,  29.42746355])
   angle = np.polyval(coef_ang,0.26)
   
   nmax = 0
   for kk in range(m):    # loop over all wavelength data arrays 
      #  rotate image 
      im = ndimage.rotate(datalist[kk],angle-180.) 
      #,reshape = False,order = 1,mode = 'constant',cval = out_of_img_val)
      
      # sum data across dispersion direction (all across) 
      # (MAX FWHM 4 pixels = 36um /0.25um = 144 rows@7000A) 
      # => so, the broadening is from the telescope optics.
      
      # here we determine how much to include of the cross-dispersion (slit) part 
      # of the PSF
      profy = im.sum(0)     # sum over axes 0; variation with axes 1
      profymem = profy
      if slit != None:
          nbins = im.shape[0]
          wid = np.polyval(slit, dist[kk] + np.arange(-im.shape[1]/2,im.shape[1]/2,1))
          for k in range(0,nbins):              
              # M = one-sided width of slit in half-pixels, within the array
              M = int(np.min([wid[k]*18*sigma,nbins/2]))
              klow = im.shape[0]/2-M              
              kmax = im.shape[0]/2+M
              profy[k] = im[klow:kmax,k].sum(0)
              # sum over slit enclosure of axes 0; variation with axes 1
          print "%i, integrate  [%i:%i,%i]"%(kk,klow,kmax,k)
      print "len profy (all) = ",len(profymem)
      print "len profy(selection) = ",len(profy)
      print "im shape",im.shape
      
      normy = profy.sum()   # total sum 
      #peakx = np.where(profy.max() == profy)[0][0]  # find location(s) maximum
      centrex = im.shape[1]/2  # half the dimension of axes 1
      # the coordinate across dispersion centred on peak in units of pixels:
      #profx = (np.arange(im.shape[1]) - centrex)*dx/pix
      
      # now rebin the data into """1/2 pix""" bins by summing 18 (=dx/pix/2) points 
      # making sure a boundary falls at the centre.
      xtra = -(centrex-centrex/18*18) # centre shift in 1/18th pixel 
      nbins = int(centrex/18)*2      # even number nbins 
      if nbins > nmax: nmax = nbins  # make sure nbins is not too large
      
      # Note that the arrays are all symmetric, so the profile in the 
      # dispersion direction can use the same scaling as in the x-dispersion
      # direction.
          
      Ebins = np.zeros( nbins )    # the energy of the bin in the dispersion direction   
      Ybins = np.zeros( nbins )    # the amplitude of the PSF integrated across the 
                 # dispersion, either all data, or limited to the slit width
          
      # since the resolution is 1/36th bin, we need to integrate over a bin in
      # the dispersion dimension
          
      # first bin 
      Ebins[0] = -nbins/2 
      Ybins[0] = profy[0:18-xtra].sum() 
          
      for k in range(1,nbins-1):
          Ebins[k] = -nbins/2 + k  
          Ybins[k] = profy[18*k-xtra:18*(k+1)-xtra].sum()
      if len(Ebins) > nmax: nmax = len(Ebins)
          
      # last bin	 
      Ebins[nbins-1] = -nbins/2 + (nbins -1)
      Ybins[nbins-1] = profy[18*(nbins-1)-xtra:].sum()  
      
      # normalise
      Ybins = Ybins/Ybins.sum()
      profile.append( [Ebins,Ybins] )
             
   result = np.zeros(m*nmax).reshape(m,nmax) 
   Epix = np.arange(nmax)-nmax/2  
   print "result dimensions ",result.shape 
          
   for k in range(m):
       [Ebins,Ybins] = profile[k]
       nbins = len(Ebins)
       print "Ebins dimension, final nbins",Ebins.shape, nbins
       result[k,nmax/2-nbins/2:nmax/2+nbins/2] = Ybins  
   
   channels = angstrom2kev(wav*1e4)   
   atime = datetime.datetime.today()   
   now = atime.isoformat()[0:19]
   # create primary header 
   #
   hdu0 = fits.PrimaryHDU()
   hdu0.header.update('TELESCOP','SWIFT   ',comment='Telescope (mission) name')
   hdu0.header.update('INSTRUME','UVOTA   ',comment='Instrument Name')
   hdu0.header.update('ORIGIN','MSSL/UCL','source of FITS file')
   
   hdulist=fits.HDUList([hdu0])
   #
   # data
   #
   col1 = fits.Column(name='CHANNEL ',format='E',array=channels )
   col2 = fits.Column(name='EPIX    ',format=str(nmax)+'E',array=Epix)
   col3 = fits.Column(name='LSF     ',format=str(nmax)+'E',array=result)
   cols = fits.ColDefs([col1,col2,col3])
   #hdu1 = fits.new_table(cols)
   hdu1 = fits.BinTableHDU.from_columns(cols)

   hdu1.header['CREATED'] = 'written by readZemaxPSF (psf.py) '
   hdu1.header['DATE']    = (str(now),'creation date')
   hdu1.header['AUTHOR']  = ('NPM Kuin (UCL/MSSL)','made by')
   hdu1.header['VERSION'] = (version,'version number of this file')
   hdu1.header['WHEELPOS']= (160,'filter wheel position')
   hdu1.header['FILTER']  = ('UGRISM','UVOT filter')
   hdu1.header['ORDERS']  = ('1','applicable orders')
   hdu1.header['TELESCOP']= ('SWIFT   ','Telescope (mission) name')
   hdu1.header['INSTRUME']= ('UVOTA   ','Instrument Name')
   hdu1.header['ORIGIN']  = ('MSSL/UCL','source of FITS file')
   hdu1.header['CREATOR'] = ('uvotcal.py','MSSL/UCL Python program')
   hdu1.header['HDUCLASS']= ('OGIP','format attemts to follow OGIP standard')
   hdu1.header['HDUCLAS1']= ('RESPONSE','subset of RESPONSE DATA')
   hdu1.header['HDUCLAS2']= ('SPECRESP','type of calibration data')
   hdu1.header['CCLS0001']= ('BCF','dataset is a basic calibration file')
   hdu1.header['CCNM0001']= ('SPECRESP','Type of calibration data')
   hdu1.header['CDES0001']= ('UGRISM LINESPREAD FUNCTION','Description')
   hdu1.header['CDTP0001']= ('DATA','Calibration file contains data')
   hdu1.header['CVSD0001']= ('2004-11-20','UTC date when calibration should first be used')
   hdu1.header['CVST0001']= ('00:00:00','UTC time when calibration should first be used')
   hdu1.header['CBD10001']= ('FILTER(UGRISM)','Parameter boundary')
   hdu1.header['CBD20001']= ('ENERG('+str(np.min(channels))+'-'+str(np.max(channels))+')keV','Parameter boundary')
   hdu1.header.add_comment(' ')
   hdu1.header.add_comment('Zemax data line spread are given for 15 channels.  ')
   hdu1.header.add_comment('Data for each channel is in half-pixel increments,  ')
   hdu1.header.add_comment('with middle of array centered on the channel')
   hdu1.header.add_comment(' ')
   hdu1.header.add_comment('Interpolated LSF needs to be converted to E using dispersion.')
   hdulist.append(hdu1)
   hdulist.writeto(lsffile,clobber=clobber)

   return metadata, datalist, result
   
   
   
 def write_rmf_file (rmffilename, wave, wheelpos, disp, 
    flux = None,
    anchor=[1000,1000], # only that one is currently available 
    spectralorder = 1,  # not possible for second order yet
    effarea1=None, effarea2=None, 
    lsfVersion='001', msg="",
    chatter=1, clobber=False  ):
   '''
   Write the RMF file for the first order spectrum 
   
   Parameters
   ----------
      rmffile : path, str
         file name output file
	 
      wave : ndarray
         mid-wavelengths of the bins in the spectrum
         
      flux : ndarray [default None]
         used to omit channels of invalid (NaN) or negative 
         flux values from the response file   
	 
      wheelpos : int
         filter wheel position 	 
	 
      disp : ndarray
         dispersion coefficients 
         
      anchor : 2-element list
         The anchor position is used to select the correct 
         effective area (important for the clocked grism modes).
         
      spectralorder : 1 
         ** Do not change **
         Only for the first order the RMF can currently be 
         computed. 
         
      effarea1, effarea2 : hdu, interpolating function
         do not use unless you know what you're doing
           
      lsfVersion : ['001','003']
         version number of the LSF file to be used.
         
      chatter : int
         verbosity
	 
      clobber : bool
         if true overwrite output file if it already exists
	 
   Returns
   -------
   Writes the RMF file 
   
   Notes
   -----
   
   The count rate has been corrected for coincidence loss (version 2 SPECTRUM extension).
   The spectral response varies in the clocked grisms, is nearly constant in the 
   nominal grisms. Therefore, in the clocked grisms, the rmf file needs 
   to be created specifically for the specific spectrum. 
   
   The rmf files have also energy bins corresponding to the wavelength 
   bins in the spectrum. These also show some variation from spectrum to
   spectrum.       
         
   The line spread function from the uv grism at default position is
   currently used for all computations. Since the RMF file encodes also
   the effective area, this version presumes given anchor position. 
   
   2014-02-27 code cleaned up. Speed depends on number of points
   2015-02-02 versioning lsf introduced; changed instrument FWHM values
   2015-02-04 error found which affects the longer wavelengths (> 3500A)
   2015-02-04 verified the LSF with a calibration spectrum, which shows 
              instrumental broadening of 10+-1 Angstrom and at long 
              wavelengths (6560) the same broadening predicted by the 
              Zemax optical model. 
   2015-02-09 There was a major overhaul of this routine, which is now 
              much improved.           
   2015-02-13 remove trimming of channels            
   ''' 
   try:
      from astropy.io import fits
   except:   
      import pyfits as fits
   import numpy as np
   import os
   from scipy.ndimage import convolve
   import uvotio
   import uvotgetspec
   from scipy import interpolate
   import datetime
   
   version = '150208'
   
   if not ((lsfVersion == '001') | (lsfVersion == '002') | (lsfVersion == '003') ):
      raise IOError("please update the calfiles directory with new lsf file.")
   
   now = datetime.date.today()
   datestring = now.isoformat()[0:4]+now.isoformat()[5:7]+now.isoformat()[8:10]
   if chatter > 0: print "computing RMF file."
   
   # telescope and image intensifier broadening   
   if lsfVersion == '001':
       if wheelpos < 500:
           instrument_fwhm = 2.7 # Angstroms in units of pix
       else:     
           instrument_fwhm = 5.8 # Angstroms in units of pix

   # get the effective area for the grism mode, anchor position and order at each wavelength
   if effarea1 != None:
      if len(effarea1) == 2:
          hdu, fnorm = effarea1          
          w = 0.5*(hdu.data['WAVE_MIN']+hdu.data['WAVE_MAX']) 
	  fnorm = fnorm(w)    
      else: 
          hdu = effarea1
          w = 0.5*(hdu.data['WAVE_MIN']+hdu.data['WAVE_MAX']) 
	  fnorm = 1    
   else:
       hdu, fnorm, msg = uvotio.readFluxCalFile(wheelpos,anchor=anchor,
                    spectralorder=spectralorder,msg=msg, chatter=chatter)
       w = 0.5*(hdu.data['WAVE_MIN']+hdu.data['WAVE_MAX'])   
       fnorm = fnorm(w)  
       
   r = hdu.data['SPECRESP']
   ii = range(len(w)-1,-1,-1)
   w = w[ii]
   r = r[ii]
   r = r * fnorm
   specrespfunc = interpolate.interp1d( w, r, bounds_error=False, fill_value=0.0 )

   # exclude channels with no effective area
   #wave = wave[(wave >= np.min(w)) & (wave <= np.max(w))]
   
   # exclude channels that have bad data
   #if flux != None:
   #    wave = wave[(np.isfinite(flux) & (flux >= 0.))]            
                   
   NN = len(wave)  # number of channels in spectrum (backward since we will use energy)
   if NN < 20:
      print "write_rmf_file: not enough valid data points.\n"+\
      " No rmf file written for wheelpos=",wheelpos,", order=",spectralorder
      return
   
   iNL = np.arange(0,NN,1,dtype=int)  # index energy array spectrum
   NL = len(iNL)   # number of sample channels
   
   aa = uvotgetspec.pix_from_wave(disp, wave, spectralorder=spectralorder)  # slow !
   tck_Cinv = interpolate.splrep(wave,aa,)  # B-spline coefficients to look up pixel position (wave)
   
   channel = range(NN)
   aarev   = range(NN-1,-1,-1)  # reverse channel numbers (not 1-offset)
   channel = np.array(channel) + 1  # start channel numbers with 1 

   # spectral response as function of energy
   resp = specrespfunc(wave[aarev]) 

   # wavelengths bounding a pixel 
   wave_lo = np.polyval(disp,aa-0.5)  # increasing
   wave_hi = np.polyval(disp,aa+0.5)
   
   # corresponding energy channels
   energy_lo = uvotio.angstrom2kev(wave_hi[aarev]) # increasing energy channels (index reverse)
   energy_hi = uvotio.angstrom2kev(wave_lo[aarev])
   #energy_mid = uvotio.angstrom2kev(wave[aarev])
   e_mid = 0.5*(energy_lo+energy_hi)  # increasing energy 

   # channel (integer) pixel positions (spectrum)
   d_lo = np.array(interpolate.splev(wave_lo, tck_Cinv,) + 0.5,dtype=int)
   d_hi = np.array(interpolate.splev(wave_hi, tck_Cinv,) + 0.5,dtype=int)
   d_mid = np.array(interpolate.splev(wave[aarev], tck_Cinv,) + 0.5,dtype=int) # increasing energy
   
   # output arrays
   n_grp = np.ones(NN)
   f_chan = np.ones(NN)
   n_chan = np.ones(NN) * NN
   matrix = np.zeros( NN*NN, dtype=float).reshape(NN,NN)
   
   # (only for original version pre-2015-02-05) instrumental profile gaussian assuming first order
   # second order needs attention: instrument + LSF
   if lsfVersion == '001':
      ww = uvotgetspec.singlegaussian(np.arange(-12,12),1.0,0.,instrument_fwhm)
      ww = ww/ww.sum().flatten()  # normalised gaussian 
      
   # get the LSF data for selected wavelengths/energies   
   lsfwav,lsfepix,lsfdata,lsfener = _read_lsf_file(lsfVersion=lsfVersion,wheelpos=wheelpos,)
                  
   # provide a spectral response for NN energies 
         
   for k in iNL:   # increasing energy
              
         # find the LSF for the current energy channel     
         lsf = _interpolate_lsf(e_mid[k],lsfener,lsfdata,lsfepix,)

         # convolution lsf with instrument_fwhm 
         if lsfVersion == '001':
            lsf = convolve(lsf,ww.copy(),)
         
         # lsf should already be normalised normalised to one
         # assign wave to lsf array relative to w at index k in matrix (since on diagonal)   
         # rescale lsf from half-pixels (centre is on a boundary) to channels 

         # find the range in pixels around the e_mid pixel at index k
         pdel = len(lsfepix)/2/2  # one-way difference; half->whole pixels

         # where to put in matrix row :
         qrange = np.arange( np.max([k-pdel,0]),  np.min([k+pdel,NN]), 1,dtype=int)

         # skipping by lsfepix twos so that we can add neighboring half-pixel LSF 
         qpix = np.arange(pdel*4-2,-1,-2,dtype=int)
         
         
         matrixrow   = np.zeros(NN)	 
         matrixrow[qrange] = lsf[qpix]+lsf[qpix+1]
	 # centre is in middle of a channel
         matrix[k,:] = matrixrow*resp[k]

   # for output
   if wheelpos < 500: 
      filtername = "UGRISM"
   else:
      filtername = "VGRISM"

   if chatter > 0 : print "writing RMF file"

   hdu = fits.PrimaryHDU()
   hdulist=fits.HDUList([hdu])
   hdulist[0].header['TELESCOP']=('SWIFT   ','Telescope (mission) name')                       
   hdulist[0].header['INSTRUME']=('UVOTA   ','Instrument Name')  
   hdulist[0].header['COMMENT'] ="revision 2015-02-08, version 003" 
    
   col11 = fits.Column(name='ENERG_LO',format='E',array=energy_lo,unit='KeV')
   col12 = fits.Column(name='ENERG_HI',format='E',array=energy_hi,unit='KeV') 
   col13 = fits.Column(name='N_GRP',format='1I',array=n_grp,unit='None')
   col14 = fits.Column(name='F_CHAN',format='1I',array=f_chan,unit='None')
   col15 = fits.Column(name='N_CHAN',format='1I',array=n_chan,unit='None' )
   col16 = fits.Column(name='MATRIX',format='PE(NN)',array=matrix,unit='cm**2' )
   cols1 = fits.ColDefs([col11,col12,col13,col14,col15,col16])
   tbhdu1 = fits.BinTableHDU.from_columns(cols1)    
   tbhdu1.header['EXTNAME'] =('MATRIX','Name of this binary table extension')
   tbhdu1.header['TELESCOP']=('SWIFT','Telescope (mission) name')
   tbhdu1.header['INSTRUME']=('UVOTA','Instrument name')
   tbhdu1.header['FILTER']  =(filtername,'filter name')
   tbhdu1.header['CHANTYPE']=('PI', 'Type of channels (PHA, PI etc)')
   tbhdu1.header['HDUCLASS']=('OGIP','format conforms to OGIP standard')
   tbhdu1.header['HDUCLAS1']=('RESPONSE','RESPONSE DATA')
   tbhdu1.header['HDUCLAS2']=('RSP_MATRIX','contains response matrix')   
   tbhdu1.header['HDUCLAS3']=('FULL','type of stored matrix')   
   tbhdu1.header['HDUVERS'] =('1.3.0','version of the file format')      
   tbhdu1.header['ORIGIN']  =('UVOTPY revision 2015-02-08','source of FITS file')
   tbhdu1.header['TLMIN4']  =( 1, 'First legal channel number')                           
   tbhdu1.header['TLMAX4']  =(NN, 'Last legal channel number')                           
   tbhdu1.header['NUMGRP']  =(NN, 'Sum of the N_GRP column')                           
   tbhdu1.header['NUMELT']  =(NN, 'Sum of the N_CHAN column')                           
   tbhdu1.header['DETCHANS']=(NN, 'Number of raw detector channels')                           
   tbhdu1.header['LO_THRES']=(1.0E-10, 'Minimum value in MATRIX column to apply')                           
   tbhdu1.header['DATE']    =(now.isoformat(), 'File creation date') 
   hdulist.append(tbhdu1)
   
   col21 = fits.Column(name='CHANNEL',format='I',array=channel,unit='channel')
   col22 = fits.Column(name='E_MIN',format='E',array=energy_lo,unit='keV')
   col23 = fits.Column(name='E_MAX',format='E',array=energy_hi,unit='keV')
   cols2 = fits.ColDefs([col21,col22,col23])
   tbhdu2 = fits.BinTableHDU.from_columns(cols2)    
   tbhdu2.header['EXTNAME'] =('EBOUNDS','Name of this binary table extension')
   tbhdu2.header['TELESCOP']=('SWIFT','Telescope (mission) name')
   tbhdu2.header['INSTRUME']=('UVOTA','Instrument name')
   tbhdu2.header['FILTER']  =(filtername,'filter name')
   tbhdu2.header['CHANTYPE']=('PI', 'Type of channels (PHA, PI etc)')
   tbhdu2.header['HDUCLASS']=('OGIP','format conforms to OGIP standard')
   tbhdu2.header['HDUCLAS1']=('RESPONSE','RESPONSE DATA')
   tbhdu2.header['HDUCLAS2']=('EBOUNDS','type of stored matrix')   
   tbhdu2.header['HDUVERS'] =('1.2.0','version of the file format')      
   tbhdu2.header['DETCHANS']=(NN, 'Number of raw detector channels')                           
   tbhdu2.header['TLMIN1']  =( 1, 'First legal channel number')                           
   tbhdu2.header['TLMAX1']  =(NN, 'Last legal channel number')                              
   tbhdu2.header['DATE']    =(now.isoformat(), 'File creation date')                           
   hdulist.append(tbhdu2)     
   hdulist.writeto(rmffilename,clobber=clobber)



 def _interpolate_lsf(en,lsfener,lsfdata,lsfepix,):
    """
    interpolate the LSF data for a different energy 
    
    parameters
    ===========
    en : float (not a list/array)
       energy (keV) for wavelength for which LSF is desired
    lsfwav : numpy array
       list of wavelengths at which we have LSF
    lsfepix : numpy array      
       for given wavelength, give LSF value for a channel 
       the size of each channel is 0.5 pixel
    lsfdata : numpy array
       2-d array. first index relates to lsfwav, second to channels
       
    returns
    =======
    lsf[channels] for wavelength w
    
    method
    ========
    the LSF data near w are linearly interpolated for each channel
        
    """
    import numpy as np
    if  not ((type(en) == float) | (type(en) == np.float32)):  
        print "en = ",en
        print "type en = ", type(en)
        raise IOError("_interpolate_lsf only works on one *en* element at a time")
    #  find index of the nearest LSF
    indx = np.argsort(lsfener) # indices  
    jj = lsfener.searchsorted(en,sorter=indx)
    j = indx[jj-1] 
    k = lsfener.shape[0]
    if j == 0: 
            lsf = lsfdata[0,:].flatten()
    elif ((j > 0) & (j < k) ):
            e1 = lsfener[j-1]
 	    e2 = lsfener[j]
 	    frac = (en-e1)/(e2-e1)
            lsf1 = lsfdata[j-1,:].flatten()
	    lsf2 = lsfdata[j,:].flatten()
	    lsf = ((1-frac) * lsf1 + frac * lsf2) 	 
    else:
	    lsf = lsfdata[k-1,:].flatten()
            
    return lsf


def _read_lsf_file(lsfVersion='003',wheelpos=160,):
    """
    2015-08-19 error in lsfepix order of row and column? 
    """
    import os
    from astropy.io import fits
    import uvotio
    
    UVOTPY = os.getenv('UVOTPY')
    if UVOTPY == '': 
        raise IOError( 'The UVOTPY environment variable has not been set; aborting RMF generation ')

    if lsfVersion == '001':
        try:  
            lsffile = fits.open(  UVOTPY+'/calfiles/zemaxlsf0160_v001.fit' ) 
        except:
            print "WARNING: the oldest zemaxlsf calfile has been read in (= wheelpos 160; version 001)"
            lsffile = fits.open(  UVOTPY+'/calfiles/zemaxlsf.fit' ) 
    else:
          # in later versions the instrumental broadening is already included in the Line Spread Function
          lsffile = fits.open(  UVOTPY+'/calfiles/zemaxlsf0160_v'+lsfVersion+'.fit' ) 
               
    if wheelpos < 500: 
        lsfextension = 1
    else:
        print "using the LSF model of the UV grism for the Visible grism "+\
        "until such time as the LSF model for the visible grism can be incorporated" 
        
    lsfener = lsffile[1].data['channel'][0:15]   # energy value (keV)
    lsfepix = lsffile[1].data['epix'][0,:]         # 156 or 158 values - offset in half pixels (to be converted to wave(wave))
    lsfdata = lsffile[1].data['lsf'][:15,:]      # every half pixel a value - - to resolve at shortest wavelengths 
    lsfwav = uvotio.kev2angstrom(lsfener)        # LSF wavelength
    lsflen = lsfdata.shape[1]
    lsffile.close()
   
    return lsfwav,lsfepix,lsfdata,lsfener


