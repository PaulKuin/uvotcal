# -*- coding: iso-8859-15 -*-
# coincidence loss ######### general ##########################################
'''
*****************************************************************************************
Software for correcting the Coincidence-loss in the OM/UVOT data (photometry and spectra)
*****************************************************************************************

The coincidence-loss correction implies the conversion from raw (measured) counts (or rate)
to the expected incident (as if the detector did not have coincidence-loss) counts (or rate).

The detector dead-time correction must be applied as the final stage. 

2012,2013 Paul Kuin   
'''
__version__ = 1.0

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

def _coibg1(obs_cntrate=None,original_data=True):
  '''Fig 6 of Breeveld ea (2010) Eriks' background computations
  computed-incident vs computed-observed in counts per sec. per pix
  '''
  import numpy as np
  import rationalfit 
  
  incident=10**np.array([-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.50,-1.25,-1.0])
  ratio_observed=np.array([0.99,0.995,0.99,0.99,0.975,0.955,0.920,0.855,0.750]) 
  raw_observed = ratio_observed*incident 
  if original_data: return incident, ratio_observed
  #frametime = 0.0110322
  #alpha = (frametime - 171e-6)/frametime
  #incident *= frametime*315.
  #observed = raw_observed* frametime*315.
  #in_cpf = -1.0/alpha * log(1.0 - observed)
  # fit ratio as function of log10( ) observed rate per second for 5" aperture  
  #C = [1.1542, -19.219, -0.31648]
  #D = [-23.9826, -5.5413]
  #ratio = (rationalfit.ratfunct(obs_cntrate,(C,D))[0] )  # observed/incident
  # fit log10(ratio) as function of log10 observed for observed < 0.4 c/s/(5" aperture)
  C = [11.66281389, 8.53597265, 3.01518394]
  D = [-0.08496817,  -0.30125974]
  ratio = 10**(rationalfit.ratfunct(obs_cntrate,(C,D))[0] )  # observed/incident
  return obs_cntrate/ratio

def _coibg2():
   """
   Data read from fig.7 (using jpeg in GIMP) of 
   Fordham, Moorhead, Galbraith: measured coi-effect in flat illuminated MIC detector
   four curves for different frametimes 1.18, 3.11, 5.69, 12.2 ms.
   
   Parameters
   ----------
     none
     
   Returns
   -------
   graph data: dict
     - **"1_18ms"** : array[:,2] with X,Y coordinates  
     - **"3_11ms"** : array[:,2] with X,Y coordinates  
     - **"5_69ms"** : array[:,2] with X,Y coordinates  
     - **"12_2ms"** : array[:,2] with X,Y coordinates  
     
      X-axis CR_in (photons/pix/sec)
      Y-axis CR_obs (counts/pix/sec) 
      pix is scale of physical pixel (8x8 subpixels)
   
   Note
   ----
   blurred data at low x-values left out, while the read-off error is large for the low 
   x-values.
   """
   import numpy as np
   origen = np.array([174,765])
   frametimes = np.array([12.2,5.69,3.11,1.18]) * 0.001 
   dx = (1417.-174)/8.
   dy = (17. - 765.)/7.5
   ft12_2ms = np.array([[174,765],[275,701],[419,614],[563,548],[761,478],[849,444],[1200,363]],dtype=float)
   ft5_69ms = np.array([[174,765],[566,525],[761,431],[851,379],[1200,243]],dtype=float)
   ft3_11ms = np.array([[174,765],[419,598],[560,515],[756,413],[844,350],[1190,184]],dtype=float)
   ft1_18ms = np.array([[174,765],[276,698],[414,600],[550,516],[745,408],[832,345],[1172,160]],dtype=float)
   # in physical coordinates
   data = {"12.2ms" :  (ft12_2ms-origen)/np.array([dx,dy]) ,
           "5.69ms" :   (ft5_69ms-origen)/np.array([dx,dy]) ,
	   "3.11ms" :   (ft3_11ms-origen)/np.array([dx,dy]) ,
	   "1.18ms" :   (ft1_18ms-origen)/np.array([dx,dy]) }	  
   return data
   
def om_coi_correction_polynomial(x):
   """
   correction polynomial to the OM coincidence-loss  
   
   Parameters
   ----------
   x : float
     the counts per frame in 6" radius circular aperture

   returns 
   -------
   coincidence loss polynomial corrected input: float, float array
      the factor to multiply the coi-corrected counts per frame/count rate 
      
   Note
   ----
   the polynomial is a small correction above the main coincidence loss correction      
   """
   #P = 1+ (-0.076+( 0.144 +(-0.64 +0.57*x)*x)*x)*x   
   return  x*(x*(x*(x*(-0.076)+0.144)-0.64)+0.57) + 1
   
def uvot_coi_correction_polynomial(x):
   """
   correction polynomial to the UVOT coincidence loss
   
   Parameter
   ---------
   x : float, float array
      the counts per frame in a 5" radius circular aperture
   
   returns 
   coincidence loss polynomial corrected input : float, float array
      the factor to multiply the coi-corrected counts per frame/count rate 
      
   Note
   ----
   the polynomial is a small correction above the main coincidence loss correction      

   """   
   coef = [0.031,0.029,-0.091,0.066,1]
   return  x*(x*(x*(x*coef[0]+coef[1])+coef[2])+coef[3])+coef[4]  

def coi_point_source_1(x,alpha=1):
   """
   theoretical coincidence loss correction for point source with **low** background
   including dead time correction

   Parameters
   ----------   
   x : float, float array
      the counts per frame in an aperture (e.g., 6" for OM, 5" for UVOT)
   
   kwargs : dict
      optional keyword arguments, possible values are:

      - **alpha** : float
        dead time correction (frame time CCD minus read out time)/frame time  
      
   Returns
   -------
   cpf : float , float array
      the counts per frame corrected for coincidence-loss and CCD dead time 
      
   Note
   ----
   A further correction polynomial of a few percent is needed  
   e.g., uvot_coi_correction_polynomial(x)       
   """
   from numpy import log
   return -log(1.0 - x) / alpha
   
def coi_point_source_2(m,q,A=1.0):
   """
   Theoretical coincidence loss correction for point source with **neighbouring  
      background emission**, and **no** dead time correction
   
   Parameters
   ----------
   
   m : float, float array
      the source+background counts per frame in an aperture (e.g., 6" for OM, 5" for UVOT)
      
   q : float, float array
      the background counts per frame in an aperture (e.g., 6" for OM, 5" for UVOT)
      
   kwargs : dict   
   
      - **A** : float 
        calibration constant for the strength of the neighbor causing coincidence loss     
        not to be changed by normal users
	
   returns
   -------
   coi corrected input parameters : list 
      the list contains corrected counts per frame m,q 
      - **m** : float or array   
        the source counts per frame corrected for coincidence-loss and CCD dead time
      - **q** : float or array
        the corrected background counts per frame  
        
   Note
   ----
   The new coi-correction is higher for sources with a high background     
   
   compare to the fuller iterated routine 
    
   _solve_coi_eq2(r,q,a,err=1e-4,niter=20):
           
   """
   from numpy import log, exp
   m1 = -log(1.0 - m)
   q1 = -log(1.0 - q)
   q2 = -log((1.0-q) / (1.0+0.5*A*q1**2*(1.0+q1**2/4. ) ) )
   m2 = -log((1.0-m) / (1.0+0.5*A*m1*q2*(1.0+m1*q2/4. ) ) )
   return m2, q2 
   
def uvot_bkg_coi_factor_1(bkg, frametime=0.0110322, aperture_radius=6.5):
    """ 
    Compute the coi-factor for UVOT background *without dead time correction* using 
    the classical theoretical formula, see Poole et al. 2008, MNRAS. 
    
    Parameters
    ----------
    
    bkg : float, float array size N
      the measured background counts per second per arcsec**2 
      
    kwargs : dict  
       - **frametime** : float, float array size N
         the frame time in seconds  
          
       - **aperture_radius** : float
         radius in arcsec of the background aperture. 
         The original default is 5.0 arcsec but 
         6.5 fits the flat background self-coincidence-loss 
    
    returns
    -------
    
    coi_factor : float, float array  
     
      the factor that the count rate needs to be multiplied by 
      to correct the coincidence loss
    """
    from numpy import pi
    x = bkg * frametime * pi * aperture_radius**2  # measured counts per frame 
    x1 =  coi_point_source_1(x,alpha=1) 
    return x1/x
    
def uvot_bkg_coi_factor_2(bkg, frametime=0.0110322, aperture_radius=6.5, A=1.0):
    """ 
    compute the coi-factor for UVOT background *without dead time correction* using 
    the extended theoretical formula
    
    Parameters
    ----------
    bkg : float, float array size N
      the measured background counts per second per arcsec**2 
      
    kwargs : dict  
    - **frametime** : float, float array size N
      the frame time in seconds   
      
    - **aperture_radius** : float
       radius in arcsec of the background aperture. The original default is 5.0 arcsec but 
       6.5 fits the flat background self-coincidence-loss 
       
    - **A** : float 
      calibration constant for the strength of the neighbor causing coincidence loss     
    
    returns
    -------
    coi_factor : float, float array   
      the factor that the count rate needs to be multiplied by 
      to correct the coincidence loss
    """
    # measured counts per frame, x:
    x  = cr_per_arcsec_to_cpf_in_aperture(bkg,frametime,aperture_radius)  
    x1 = coi_point_source_2(x,x,A=A,alpha=1) 
    return x1/x

def cr_per_arcsec_to_cpf_in_aperture(cr,frametime,aperture_radius):
    """
    Convert from counts per second per arcsec to counts per frame (no deadtime correction)
    
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
       
def _coi_spectral_calibration(calobsfile,effective_area_file,phafile,coi_output_file,
    coi_length=26,frametime=None, fudgebg=1.685, fudgespec=1.55, option=1, chatter=1):
   ''' *Development* COI study routine
   
   Parameters
   ----------
   calobsfile : str
      path+filename of the reference spectrum (ascii)
   
   effective_area_file : str
      path+filename of the effective area file (ascii) 
      
   phafile : str
      path+filename of the extracted spectrum file (fits)
      
   coi_output_file : str
      path+filename output file (ascii)
      
   coi_length : int
      number of pixels along dispersion to average spectrum over before calculation 
      coi-loss at pixel 
      
   frametime : float
      frame time in seconds
      
   fudgebg : float
      adjust background area (Poole et al.: 1=>315 subpixels, 
      extended coi: 1.685 = 6.5" radius, 531 subpixels)                   
   
   
   Notes
   -----
   use the calibration spectrum and effective area to determine incoming counts per frame 
   use the extracted spectrum and background to determine the total 
   measured counts in the spectrum per frame, and the measured background counts 
   in the spectrum per frame, averaged over the coi_length (in pixels).
      
   background coi set to classic theory when fudgebg=1   
   output the coi_correction as a function of counts per frame
   
   2012-03-05 NPMK 
   
   '''   
   from uvotmisc import rdTab
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
      print "ERROR: counts in background too high for point source coi formula."
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

   
def coi_func(pixno,wave,countrate,bkgrate,sig1coef=[3.2],option=2,
   fudgespec=1.0,coi_length=26,frametime=0.0110302, background=False,
   sigma1_limits=[2.6,4.0], trackwidth = 2.5, ccc = [-1.5,+1.5,-1.5,+1.5,-1.5,+1.5,+0.995],
   ccb = [+0.72,-0.72,0.995], ca=[0,0,3.2],cb=[0,0,3.2],debug=False,chatter=1):
   #ccb = [+2.68,-2.68,-3.3,+3.3,0.995], debug=False,chatter=1): - proper background 
   '''Compute the coincidence loss correction factor to the (net) count rate 
   as a function of wavelength  EXPERIMENTAL
   
   Parameters
   ----------
   pixno : array-like
      pixel number with origen at anchor
   wave : array-like
      wavelength in A
   countrate : array-like
      input count net rate must be aperture corrected
   bkgrate : array-like
      background rate for trackwidth
   kwargs : dict
   
      - **sig1coef** : list
      
        polynomial coefficients
	
      - **frametime** : float
      
        CCD frame time in seconds
	
      - **trackwidth** : float
      
        width of the extraction in standard deviations of the profile matched across the spectrum
	
      - **option** : int 
      
        . option = 1 : classic coi-loss, but spectrum is box like 10x32 pix across spectrum
	
        . option = 2 : bkg classic coi-loss, total (spectrum*poly+bkg*poly) with 
                     polynomial corrections for extended coi-loss. 
		     classical limit for ccc= [0,0,1] ; ccb[0,0,1] 
		     
      - **background** : bool
      
        if the background is `True` an interpolated function for the coi 
	correction factor in the background count rate is returned
        
	if the background is `False` an interpolated function for the 
	coi correction factor in the net target count rate is returned
      		       
      
   
   Returns
   -------
       coi_func : scipy.interpolate.interpolate.interp1d
          if **background** is `True` an interpolated function for the coi correction 
          factor in the background count rate while if **background** is `False` an 
          interpolated function for the coi correction factor in the net target 
          count rate is returned 
   
   Notes
   -----   
   defaults to the background coincidence loss equivalent to an area of 
   315 sub-pixels (about pi*5"^2 on the detector) 
   
   Also see the discussion of coincidence loss in Breeveld et al. (2010).
   Her correction for high background + high source rate was used as inspiration.
   
   - 2012-03-21 NPMK initial version
   - 2012-07-04 NPMK added into option 1 the white-correction for high 
     background (photometry) from Alice Breeveld (2010) 
   - 2012-07-24 NPMK modified source area to be same as background area  
   - 2012-07-24 NPMK start modifications extended coi-model 
   - 2012-09-25 NPMK simplify. Add extended coi-loss as polynomial using classic coi as approx. 
                   coefficients TBD. 
   - 2012-10-10 NPMK temporary option 3 to address consistent approach with Breeveld et al. and 
       the coi-work on point sources. Basically, it is not a reduction in the background but 
       a lack of accounting for additional losses in the main peaks (due to surrounding 
       high background levels stealing counts from source). Option 2 has now been optimized 
       to work. Basically, there are multiple practical solutions to the problem, the third 
       option will be more in line with the theoretical model for coincidence loss in the 
       UVOT. 	   
   '''   
   import uvotmisc
   import numpy as np
   try:
     from uvotpy import uvotgetspec as uvotgrism
   except:  
     import uvotgrism
   try:
      from convolve import boxcar
   except:
      from stsci.convolve import boxcar   
   from scipy import interpolate
    
   if not do_coi_correction:   # global - use when old CALDB used for fluxes.
      # set factor to one:
      return interpolate.interp1d(wave,wave/wave,kind='nearest',bounds_error=False,fill_value=1.0 ) 
      
   if type(trackwidth) != float:
      raise TypeError ( "trackwidth is not of type float, trackwidth type: ",type(trackwidth) )   
  
   alpha = (frametime - 0.000171)/frametime
   
   # mask bad and problematic data 
   if background: v = np.isfinite(countrate) &  (bkgrate > 1e-8) 
   else: v = np.isfinite(countrate) & np.isfinite(bkgrate) & (countrate > 1e-8) & (bkgrate > 1e-8) 
   countrate = countrate[v]
   bkgrate   = bkgrate[v]
   pixno     = pixno[v]
   wave      = wave[v]
   
   # reset v
   v = np.ones(len(countrate),dtype=bool)
   
   # correct cpf to 550 subpixels in size, 5 sigma total width (~17.5). (about 6.5" circle on lenticular filter)
   # this initial setting was changed to 315 to match to the Poole method for photometry, but actually, may be 
   # the correct choice after all for the background-background coi-correction (high backgrounds), see Kuin (2013)
   # study on coincidence loss. 
   
   sigma1 = np.polyval(sig1coef, pixno)
   sigma1[ sigma1 > sigma1_limits[1] ] = sigma1_limits[1]
   sigma1[ sigma1 < sigma1_limits[0] ] = sigma1_limits[0]
   
   # scaling the counts per frame 
   #  - get the background counts per pixel by dividing through 2*sigma1*trackwidth
   #  - scale background to number of pixels used in photometry coi-background  correction
   bgareafactor = 315.0/(2 *sigma1*trackwidth)  # strip 'trackwidth' sigma halfwidth determined by input
   factor = 315.0/(2 *sigma1*trackwidth)        # strip 'trackwidth' sigma halfwidth determined by input
   specfactor = 315.0/(2.*sigma1*2.5)           # aperture correction was assumed done on the input rate to be consistent with the 2.5 sigma Eff. Area
   # coi-area spectrum in limit (net = zero) must be background one, so same factor
   # Very high backgrounds deviate (Breeveld et al. 2010, fig 6; ccb=[+2.68,-2.68,-3.3,+3.3,0.995] matches plot)
   
   # one pixel along the spectrum, 2 x sigma x trackwidth across, aperture corrected countrate (not bkgrate)
   # works for the lower count rates: total_cpf = boxcar( (countrate*fudgespec + bkgrate) * frametime  ,(coi_length,))
   if not background: 
      tot_cpf = obs_countsperframe = boxcar((countrate + bkgrate) * frametime, (coi_length,))
      net_cpf = boxcar( countrate * frametime, (coi_length,))
   bkg_cpf = bkg_countsperframe = boxcar( bkgrate * frametime, (coi_length,) )  
   # PROBLEM: boxcar smooth does not work for pixels on the array ends. downturn coi-correction. Need something better.
   
   if chatter > 3: 
	 print "alpha  = ",alpha
	 print "number of data points ",len(countrate),"  printing every 100th"
	 print " i    countrate    obs counts/frame "
	 for ix in range(0,len(countrate),10):
	   if background: print "%4i %12.5f %12.5f " % (ix, bkgrate[ix],bkg_cpf[ix])
	   else: print "%4i %12.5f  %12.5f" % (ix, countrate[ix],obs_countsperframe[ix])
	   
   try:
	 	 
      bkg_cpf_incident = (-1.0/alpha) * np.log(1.0 - bgareafactor*bkg_countsperframe)/(bgareafactor)
      
      if option == 1:   # default 
	 # classic coi formula
         yy = 1.0 - specfactor*obs_countsperframe
	 v[ yy < 1e-6 ] = False
	 yy[ yy < 1e-6 ] = 1e-6    # limit if yy gets very small or negative !!
         obs_cpf_incident = (-1.0/alpha) * np.log(yy)/specfactor
	 
      if option == 2: 	 
	 # new default reverts to classic coi-formula when all coef = 0 except the last one, which must be 1. 
         # extended coi-loss coefficients ccc, ccb 
	 if background: v[bkg_cpf*factor >= 0.9999] = False 
	 else: v[tot_cpf*factor >= 0.9999] = False	 
         ccc = np.asarray(ccc)
         ccb = np.asarray(ccb)
	 # extended coi-loss correction of counts per frame  - add polynomial corrections
         if not background: 
	    total_cpf  = obs_countsperframe = boxcar((countrate * np.polyval(ccc,tot_cpf*specfactor) + \
	                 bkgrate * np.polyval(ccb,bkg_cpf*factor)) * frametime , (coi_length,)) 
         bkg_countsperframe = boxcar( bkgrate * np.polyval(ccb,bkg_cpf*factor) * frametime , (coi_length,)) 	   
         bkg_cpf_incident = (-1.0/alpha) * np.log(1.0 - factor*bkg_countsperframe)/(bgareafactor)
         if not background:
	    yy = 1.0 - factor*total_cpf
	    v[ yy < 1e-4 ] = False
	    yy[ yy < 1e-4 ] = 1e-4    # limit if yy gets very small or negative !!
            obs_cpf_incident = (-1.0/alpha) * np.log(yy)/factor
	 	     	 
      if option == 3: 	 
	 # extension reverts to classic coi-formula  . 
         # extended coi-loss coefficients ccc, ccb acting on variable z = cpf * ( 1 - cpf )
	 # high background coi-loss correction fits FIG 6 in Breeveld et al. 
	 if background: v[bkg_cpf*factor >= 0.9999] = False 
	 else: v[tot_cpf*factor >= 0.9999] = False	 
         # convert to actual cpf:
	 ##CPFnet = net_cpf*specfactor
	 CPFtot = tot_cpf*specfactor
         CPFbkg = bkg_cpf*factor
	 z_tot = CPFtot * (1. - CPFtot)   # binomial type of variable
	 z_bkg = CPFbkg * (1. - CPFbkg)
	 
	 # extended coi-loss CPF correction of counts per frame  - correct with polynomial corrections in z
         if not background: 
	    CPFtot_corr = CPFnet*(1. + np.polyval(ca,z_tot)) + CPFbkg*(1. + np.polyval(cb,z_tot))		          
	 CPFbkg_corr = CPFbkg*(1 + np.polycal(cb,z_bkg))  
         CPFbkg_in   =  (-1.0/alpha) * np.log(1.0 - CPFbkg_corr)
	 bkg_cpf_incident = CPFbkg_in/factor
         if not background:
	    yy = 1.0 - CPFtot_corr
	    v[ yy < 1e-4 ] = False
	    yy[ yy < 1e-4 ] = 1e-4    # limit if yy gets very small or negative !!
	    CPFtot_in = (-1.0/alpha) * np.log(yy)
            obs_cpf_incident = CPFtot_in/specfactor
         

   except:
      print "ERROR: probably the effective counts per frame are > 1."
      print "WARNING: Continuing Setting COI factor = 1.0"
      obs_cpf_incident = obs_countsperframe
   
   # notify user that some points were flagged bad
   if v.all() != True: 
      ngood = len( np.where(v)[0] )
      print "WARNING uvotgetspec.coi_func(): Some data were ignored \n"+\
      "in the determination of the COI factor, since they exceeded the theoretical limit! "
      print "                              number of good points used = ",ngood
   
   # compute the coi-correction factor
   if not background: coi_factor    = (obs_cpf_incident - bkg_cpf_incident) / (obs_countsperframe - bkg_countsperframe)
   bg_coi_factor = (bkg_cpf_incident)/(bkg_countsperframe)
   
   # debug info
   if (chatter > 4) & (not background):
      print "bkg_countsperframe bkg_cpf_incident obs_countsperframe obs_cpf_incident bg_coi_factor coi_factor"
      for i in range(len(obs_cpf_incident)):
         print "%3i  %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f" % (i,bkg_countsperframe[i],bkg_cpf_incident[i],\
	 obs_countsperframe[i],obs_cpf_incident[i],bg_coi_factor[i],coi_factor[i])
   
   # calibrate
   if chatter > 0: 
      if not background: print " coi_factor stats (min, mean, max): ",np.min(coi_factor),np.mean(coi_factor),np.max(coi_factor)
      print " bgcoi_factor stats (min, mean, max): ",np.min(bg_coi_factor),np.mean(bg_coi_factor),np.max(bg_coi_factor)
   
   # assume wave is monotonically increasing: 
   if not background: coi_func = interpolate.interp1d(wave[v],coi_factor[v],kind='nearest',bounds_error=False,fill_value=1.0 ) 
   coi_bg_func = interpolate.interp1d(wave,bg_coi_factor,kind='nearest',bounds_error=False,fill_value=1.0 )
   if debug:   return coi_func, coi_bg_func, (coi_factor,coi_bg_factor,factor,obs_cpf_incident,bkg_cpf_incident)
   elif background: return coi_bg_func 
   elif (not background): return coi_func

  
def _coi_func2(pixno,wave,countrate,bkgrate,sig1coef=[3.2],option=2,
   fudgespec=1.0,coi_length=26,frametime=0.0110302, background=False,
   sigma1_limits=[2.6,4.0], trackwidth = 1.0, ccc = [0.,-0.,0.40],
   ccb = [0.,-0.67,1.0], debug=False,chatter=1):
   '''Second order version of coi_func
   
   Only difference from coi_func are preset defaults.
   '''
   return coi_func(pixno,wave,countrate,bkgrate,sig1coef=[4.5],option=2,
   fudgespec=1.0,coi_length=26,frametime=0.0110302, background=False,
   sigma1_limits=[3.0,7.5], trackwidth = 1.0, ccc = [0.,-0.,0.40],
   ccb = [0.,-0.67,1.0], debug=False,chatter=1)
   

#========= revised coi-calibration 1: figures for the paper 2013-12-11
#
#  As was pointed out by Mat Page, the coi-effect does not depend on the effective area
#  and mixing them together may confuse people about their separate nature.
#  Dividing the observed flux by the flux of the source spectrum gives a  measure of the 
#  coincidence loss independent of the effective area, since any discrepancies are 
#  due to the coi-loss. That presumes that we have the correct effective area 
#  (in the absence of coincidence loss) available to calibrate the spectrum. 
#  In order to simplify matters, we use the effective area once to derive the 
#  expected count rate/pix from the known source spectrum. Then we only have the 
#  count rate/pix to compare. 
#  We do this at the default position to get the coincidence loss correction right.
#  We do this for all modes. 
#  We should see the coincidence loss as a function of count rate in source and background.
#
#  function _coi_refspectrum2CRppix(calspecfile, C_1, ) -> CRppix_spectrum, random_error, sys_error
#  function _coi_normspectrum(obsspecfile,) -> CR, BG_CR, normalised_CR, error (==1 +noise if no coi)
#
#  figuur normalised_CR source only , different sources
#  figuur normalised_CR source, different BG
#

   




#========= revised coi-calibration 2: corrected coi-formulation 2013-12-11
#
#  The formula for the coincidence loss has been derived. There are free parameters.
#  The calibration determines the free parameters independently for each grism mode.
#  They should only be dependent on the MCP and detector, not on the grism. 
#
#  
#


#=========== determination of the coincidence loss area using the fact that 
# in the limit for very high illumination the count rate per frame 
# goes to one per coincidence loss area. 
# Due to the Mod-8 due to coincidence loss pattern, averaging is needed.
# In the grism (uv clocked mode) the pattern repeats close to 24/cos(angle) 
# where the angle is from the slope of the grism on the detector with repect 
# to the X-axis. 
# Given the measured average counts per pixel, the coi area is then found.

def _find_grism_coi_area_size(
    file1="/data/novae/TPyx/00031973019/uvot/image/sw00031973019ugu_1ord_1_f.pha",
    angle = 35.7,
    ):
    """
    determine the grism coi area size by studying options
    
    input parameters
    ----------------
    file : path
       input spectrum (output from uvotgetspec,getSpec()
    angle : float
       reduced angle of the spectrum on the detector (i.e., angle = 180-full_angle)   
    
    Note
    ----
     second image that can be used:
     ="/data/novae/TPyx/00031973019/uvot/image/sw00031973019ugu_1ord_2_f.pha",
    """
    import numpy as np   
    from astropy.io import fits
    import uvotmisc
    
    f1 = fits.open(file1)
    ankx, anky = f1[3].header['ankximg'], f1[3].header['ankyimg']
    img = f1[3].data
    C_1 = uvotmisc.get_dispersion_from_header(f1[1].header)
    sig = uvotmisc.get_sigCoef(f1[1].header)
    cur = uvotmisc.get_curvatureCoef(f1[1].header)
    pix = np.arange(-ankx,img.shape[1]-ankx,1)
    exposure = f1[1].header['exposure']
    telapse = f1[1].header['telapse']
    deadc = f1[1].header['deadc']
    loss = (telapse - exposure/deadc)/telapse
    if loss > 0.01:
       print "loss of data " 
    q = (pix > 0) & (pix < 500) 
    wav = np.polyval(C_1,pix[q])
    sp = np.zeros(101*q.sum(),dtype=float).reshape(101,q.sum())   
    for i in range(q.sum()):
       y  = np.polyval(cur,pix[q][i])
       y0 = int(anky +y   - 50+0.5) 
       y1 = int(anky +y   + 51+0.5)
       sp[:,i] = img[y0:y1,ankx+pix[q][i]]
    # now the spectrum is straight in sp centred on 50 and 101 pix wide
    m, p, smsp = _find_coi_smoothing_spectrum(pix[q],wav,sp, sig,)
    Nframes = exposure/deadc/0.0110329
    return p, smsp
    
    
def _find_coi_smoothing_spectrum(pix,wav,sp,sig):
    """
    smooth the spectrum for the coi using different ways
    run through a set of parameters and put them + result in a list
    
    input parameters
    ----------------
    pix : float numpy array
        x coordinate referenced to anchor position 
    wav : float numpy array
        wavelength in x in angstrom
    sp  : float 2D numpy array
       spectrum rectified and centred size [101,len(pix)]
    sig : list, array
       polynominal coefficient for the sigma parameter of the gaussian fit
    
    output 
    ------
    meth : list 
       method id
    par : list
       list of parameter list 
    result : list of float array spectrum counts/pixel
                     
    """
    import numpy as np
    from stsci.convolve import boxcar
    
    meth = []
    par = []
    result = []

    method = 0
    Xcoi = np.arange(27,31,1) # average over Lcoi
    Ycoi = np.array([7,10,16,22,28,32,]) # average over half-width from centre
    n = int((sp.shape[0]+0.5)/2)
    for kx in Xcoi:
        for ky in Ycoi:
	    smsp = (boxcar(sp[n-ky:n+ky+1,:],(kx,)) ).mean(axis=0)
            meth.append(method)
	    par.append([kx,ky])
	    result.append(smsp)
	    
    return meth,par,result
    
    
    
    
    	    
