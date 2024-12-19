''' wavelength fitting subroutines
    These are useful to assign improved pixel distance coordinates to 
    the spectral lines in the calibration source, plot the result, and 
    improve the fit. 
'''
from uvotpy.uvotmisc import rdTab
from uvotcal import findAnker,getZmx 
from pylab import polyfit,polyval,plot,subplot,text,ginput,xlim,ylim,\
   xlabel,ylabel,title,clf, searchsorted
from numpy import zeros,array,where,min,max,arange,copy
 
def approx_solution(wheelpos,disp_coeff, source=None):
   ''' 
   Read the wavelengths from spectral lines seen in the calibration sources.
   The source for Vgrism is WR121 by default
   For the UVgrism WR86 (nominal default) WR52 (clocked default) [not implemented ]
   The dispersion solution comes from unscaled zemax solution.
   
   Paul Kuin 2009 '''
   C = disp_coeff.copy()   
   if ((wheelpos == 1000) ^ (wheelpos == 955) ^ (wheelpos == 160) ^ (wheelpos == 200)):
      if ((source == None) ^ (source == 'WR121')): 
         source = 'WR121'
         file = '/Volumes/data/grism/WR121/lines.txt'
      elif source == 'WR52':
         file = '/Volumes/data/grism/WR52/lines.txt'
      elif source == 'WR1':
         file = '/Volumes/data/grism/WR1/lines.txt'     
      elif source == 'SMC':
         file = '/Volumes/data/grism/HD5980SMC/lines.txt'                
      elif source == None:      
         print("wheelpos is not valid for this function")
         return
   else:
      print("unknown wheelpos: using source as filename lines.txt")
      file = source              
   print('reading '+file)         
   T = rdTab(file)       
   w_line = T[:,0]  # wave lengths lines
   g_line = T[:,1]  # weight lines (1=low, 2=medium, 3=high, 4=very high) 
   xx = arange(-400,1200)
   ww = polyval(C,xx)
   if wheelpos > 500:
      q = where((ww > 2500) & (ww < 7000) )
   else:
      q = where((ww > 1700) & (ww < 6000) )
   D = polyfit(ww[q],xx[q],4)      
   d_line = polyval(D,w_line)  # find the pixel distance from the anker
   return  (d_line, w_line, g_line)   

def measure_position( dis, spectrum, disp_coeff, approx_solution, disoffset=0,rminmax=[-3,3]):
   ''' interactive point and click. '''
   (d_line, w_line, g_line) =  approx_solution
   C = disp_coeff
   print('=================================== ') 
   print('click the correct positions of the lines')
   print('right click cancels last input ')
   print('middle button click ends input')
   print('you have 5 minutes ')
   print('=================================== ')
   clf()
   plotspectrum(dis+disoffset, C, spectrum, approx_solution)
   xx = input('hit return . . . ')
   n=1
   pointed = ginput(n=n,timeout=30)
   # match pointed values with wavelengths 
   
   d_correct = d_line.copy() 
   p = zeros(n)
   for k in range(n):  p[k] = pointed[k][0]
   print('measured positions are ',p)
   kd = searchsorted(d_line,p)  # find left index in d_line 
   if kd[0] == len(d_line):
      print('correcting position for line at ',w_line[-1],'  with ',p[0])
      ans = int( input('to confirm type 0 (-99=abort): ') )
      if ans == 0: 
         d_correct[-1] = p[0]
         g_line[-1] += 5
   elif kd == len(d_line)-1: 
      print('correcting position for line at ',w_line[-1],'  with ',p[0])
      ans = int( input('to confirm type 0 (-99=abort): ') )
      if ans == 0: 
         d_correct[-1] = p[0] 
         g_line[-1] += 5  
   elif kd == 0:
      print('correcting position for line at ',w_line[0],'  with ',p[0])
      ans = int( input('to confirm type 0 (-99=abort): ') )
      if ans == 0: 
         d_correct[0] = p[0]
         g_line[0] += 5   
   else:
      for k in range(n): 
         if ( abs(p[k] - d_line[kd[k]]) < abs(p[k] - d_line[kd[k]+1]) ): 
            rmin=rminmax[0]
            rmax=rminmax[1]
            if rmax+kd[k] > len(d_line)-1: rmax = len(d_line)-1-kd[k] 
            if rmin+kd[k] < 0: rmin = 0
            print('kd: ',kd[k],'    rmin .. rmax:   ',list(range(rmin,rmax)))
            for kk in range(rmin,rmax): print('correcting position for line ',kd[k]+kk,' at ', w_line[kd[k]+kk])
            ans = int( input('give the number of the correct line (-99=abort): ') )
            if ans == -99: return (d_correct, w_line, g_line)
            #print 'ans: ',ans,' original d_line[ans]: ',d_line[ans],'  new  d_line[ans]: ', p[k]
            d_correct[ans] = p[k]
            g_line[ans] += 5
         else:
            print("cannot locate the line \n")                                       
   print('ans: ',ans,' original d_line[',ans,']: ',d_line[ans],'  new  d_line[',ans,']: ', p[k])
   return (d_correct, w_line, g_line)  
   
def plotspectrum(dis_,C, spnet,z,linewidth=15,disoffset=0,g_min=0,linecolor='b'):
   (d_line_, w_line_, g_line) = z 
   q = where( g_line >= g_min )
   d_line = d_line_[q]
   w_line = w_line_[q]
   dis = copy(dis_)+disoffset
   dis1 = dis[where( (dis > -400) & (dis < 1450) )]
   spectrum = spnet[where( (dis > -400) & (dis < 1450) )]
   wave = polyval(C,dis1)   
   wmin = min(min(w_line) - 3*linewidth,min(wave))
   wmax = max(max(w_line) + 3*linewidth,max(wave))
   xlim(wmin,wmax)
   valid = where( (wave > wmin) & (wave < wmax)) 
   fmax = max(spectrum[valid])
   fmin = min(spectrum[valid])
   ylim(fmin, fmax*1.5)
   stripe = array([0.02, 0.12])*(fmax - fmin)
   
   subplot(211)
   plot(wave[valid],spectrum[valid],linecolor,ls='steps',lw=1)
   n = len(w_line)
   ylim(fmin, fmax*1.5)
   xlim(wmin,wmax)
   xlabel('$\lambda(\AA)$' )
   title('Correct the line position')
   for k in range(n): 
      w1 = w_line[k] - linewidth
      w2 = w_line[k] + linewidth
      aa=where( (w1 < wave) & (w2 > wave) )
      if len(aa[0]) > 0: 
         f = max( spectrum[aa] ) 
         yt = f + 0.15*(fmax-fmin)
         xt = w_line[k]
         plot([xt,xt],f+stripe,'k',lw=1.3 )
         ylim(fmin, fmax*1.5)
         text(xt, yt, str(w_line[k]), horizontalalignment='center', rotation='vertical',fontsize=9 )

   subplot(212)
   xlabel('pixel distance')
   plot(dis1[valid],spectrum[valid],linecolor,ls='steps',lw=1)
   for k in range(n): 
      w1 = w_line[k] - linewidth
      w2 = w_line[k] + linewidth
      aa=where( (w1 < wave) & (w2 > wave) )
      if len(aa[0]) > 0:
         f = max( spectrum[aa] ) 
         yt = f + 0.15*(fmax-fmin)
         xt = d_line[k]
         plot([xt,xt],f+stripe,'k',lw=1.3 )
         ylim(fmin, fmax*1.5)
         text(xt, yt, str(w_line[k]), horizontalalignment='center', rotation='vertical',fontsize=9 )
      
      
def write_result(filename, xxx_todo_changeme ,note=None, clobber='No',disoffset=0, g_min=5 ):    
   ''' write the pixel distance of the lines to the anchor 
       for note suggestion is something like : filestub+"ugv ext="+str(ext) 
   '''
   (d_correct, w_line, g_line) = xxx_todo_changeme
   import os
   import time
   q = where( g_line >= g_min )
   dd = d_correct[q]
   ww = w_line[q]
   gg = g_line[q]
   path = os.getcwd() + '/'
   if (not os.access(path+filename,os.F_OK)):
      f = open(path+filename, 'w')  
      f.write( "#"+time.ctime()+"    dispersion measured \n" )
      f.write( "# dis pixel offset used = "+str(disoffset)+" \n" )
      f.write( "# directory path = "+path+" \n" )
      if note != None: f.write("# "+ note + "\n")
      for k in range(len(dd) ):
         f.write( str(dd[k])+' \t'+str(ww[k])+' \t'+str(gg[k]-g_min)+' \n') 
   else: 
      if clobber == 'Yes':
         os.remove(path+filename) 
         write_result(filename,  (d_correct, w_line, g_line),note=note )
      else:
         print("Error: file exists, set clobber=\'Yes\' to overwrite  ")
      return
         
def new_disp_coeff(xxx_todo_changeme1, order = 4, g_min=0 ):
   ''' compute dispersion coefficients. Select only points with 
   g_line values >= g_min '''
   (d_line, w_line, g_line) = xxx_todo_changeme1
   q = where(g_line >= g_min)
   disp_coeff = polyfit( d_line[q], w_line[q], order )
   return disp_coeff    
   
def read_result(filename):
   f = open(filename,'r')
   lines = f.readlines()
   disoffset = lines[1].split()[6]
   f.close()
   t = rdTab(filename,' ')
   (d_line, w_line, g_line) = (t[:,0],t[:,1],t[:,2])
   return   (d_line, w_line, g_line), disoffset  
         
def waveAccPlot(wave_obs,pix_obs, wave_zmx, pix_zmx, disp_coef, obsid=None, 
    acc=None, order=None, wheelpos=200, figureno=1,legloc=[1,2],chatter=0):
   ''' 
   plots of the accuracy of the wavelength solution from zemax compared to
   the observed wavelengths
     
     
   x-axis = pix - pixel number referenced to [260nm in first order]

   Top panel: 
   y-axis: lambda - lambda_linear
   
   lambda linear:
      fit a linear term to the wavelengths lambda_lin = coef[0]+coef[1]*pix
   wave_obs, pix_obs: observed wavelengths points (green circles)
   wave_zmx, pix_zmx: calculated zemax points (or the interpolated solution (red crosses) 
   wave, pix = zemax dispersion relation 
   
   Bottom panel@
   y-axis: residuals:  
      wave_obs, pix_obs - wave(pix_obs)  (green circles)
      wave_zmx, pix_zmx - wave(pix_zmx)  (red crosses)
      
   disp_coef = coefficients in reverse order: if p is of length N, this the polynomial
   is as follows for coeff named p:
          y(x) = p[0]*(x**N-1) + p[1]*(x**N-2) + ... + p[N-2]*x + p[N-1]
   acc = accuracy in wavelength   
   order = order of polynomial disp_coef (default len(coef) )
   obsid = if given, append to  title                      
   '''
   import numpy as N
   from pylab import ioff,ion,arange, plot, subplot, xlim, ylim, title, xlabel, \
     ylabel, polyval, figure, contour, plt, legend, polyval, polyfit, savefig, \
     text , grid, clf, gca
   if wheelpos < 500:
      ref_wave = 2600.
      titl = 'Wavelength accuracy UV grism - '
      textstart = 1600
   else:
      ref_wave = 4200.
      titl = 'Wavelength accuracy V grism - '
      textstart = 2700   
   # zero of pix_obs forced to ref_wave (2600.or 4200.) for initial plot
   if order == None:
      order = len(disp_coef)
   dcoef = polyfit(wave_obs,pix_obs,order)
   doff = polyval(dcoef,ref_wave)
   pix_obs = pix_obs - doff
   if chatter > 0:
      print("adopted ref_wave = ",ref_wave)
      print("fit through observations pixel position of anchor = ",doff)
   
   n1, n2 = len(pix_obs), len(pix_zmx)
   pix1 = N.zeros( (n1+n2) )
   pix1[0:n1,] = pix_obs
   pix1[n1:(n1+n2),] = pix_zmx
   pix2 = pix1.min()
   pix3 = pix1.max()
   pix = N.arange(pix2,pix3)
   wav = polyval(disp_coef, pix)
   #           wavlin = disp_coef[-1]+disp_coef[-2]*xxx_pix  
   #           linear term in dispersion:
   w_obs  = wave_obs - (disp_coef[-1]+disp_coef[-2]*pix_obs)
   w_zmx  = wave_zmx - (disp_coef[-1]+disp_coef[-2]*pix_zmx)
   wavlin = wav - (disp_coef[-1]+disp_coef[-2]*pix)
   zero_offset = (wave_obs-polyval(disp_coef, pix_obs+doff) ).mean()
   zo = zero_offset
   if acc == None:
      wave_off = (wave_obs-polyval(disp_coef, pix_obs+doff) )
      acc = wave_off.std()
      print(' first estimate')
      print(' initial accuracy (all points) = ',acc) 
      # remove outlyers
      q_in = N.where(abs(wave_off-zo) < 3.* acc)
      acc = (wave_off[q_in]).std()
      print(' after removing outliers:')
      print(' accuracy of the fit = ',acc, ' angstrom')
   stracc =    str(((10*acc+0.5).__int__())/10.) +'$\AA$'
   zero_offset = ((10*zero_offset+0.5).__int__())/10.
   txt = '<$\Delta\lambda$> = '+str(zero_offset)+'$\AA\ \ \ \sigma_{observed-model}$ = '+stracc

   figure( num=figureno )

   subplot(211)
   plot(pix, wavlin, '-')
   plot(pix_obs,w_obs,'ob')
   plot(pix_zmx,w_zmx,'+r')
   ylabel('$\lambda$ - $\lambda_{linear}$  ($\AA$)')
   xlabel('pixels')
   if order == 4: 
     sord = 'fourth '
   elif order == 3:
     sord = 'third '
   elif order == 2:
     sord = 'second '
   elif order == 1: 
     sord = 'first '
   else: 
     sord = 'unknown '        
   legend((sord+'order fit','observed data','model'),loc=legloc[1])
   if obsid == None: obsid=''
   title(titl+obsid)
   # a = getp( gca )
   # setp(a, xlim=(pix1,pix2), xticks=[])

   subplot(212)
   w1 = wave_obs-polyval(disp_coef, pix_obs+doff)
   w2 = wave_zmx-polyval(disp_coef, pix_zmx)
   plot(wave_obs,w1, 'ob',label='_nolegend_')
   plot(wave_zmx,w2, '+r',label='_nolegend_')
   p0 = pix*0.
   p1 = p0 - acc+zo
   p2 = p0 + acc+zo
   plot(wav,p0,'-r',label='_nolegend_')
   plot(wav, p1,'--b',label='1-$\sigma$ limits')
   plot(wav, p2,'--b',label='_nolegend_' )
   ylabel('$\Delta\lambda$ ($\AA$)')
   xlabel('$\lambda$ ($\AA$)')
   a = gca()
   ylim = a.get_ylim()
   #if (ylim[0] > -16.0): ylim=(-16.0,ylim[1])
   ylim=(zo-2.1*acc,zo+2.1*acc)
   a.set_ylim(ylim)
   legend(loc=legloc[0])
   text(textstart,ylim[0]*0.9,txt)
   #a = getp( gca )
   #lim1, lim2 = 1.1*max(w1), 1.1*min(w1)
   #setp(a, xlim=(lim1,lim2)) #, xticks=[])
   savefig('accuracy.png')
   return acc, zero_offset

def makePlots(ra,dec,filestub,wheelpos,ext=1,ext1=1,ext2=None,lfilt1='uvw1',lfilt2=None,
    objectname=None,tfile='_measlines.txt',source='WR121',specwidth=20):
  ''' Easy call  when made from the object directory 
     e.g., ~/grism/WR121/ 
     
  makes plots and writes hyperlink in files for web pages!!!  
  
  verify that the line ids in the spectrum (pixel) are correct, or an additional error is present 
  in the anchor position.
     '''
  import uvotgrism
  from pylab import polyval,plot,title,xlabel,ylabel,savefig,clf,polyfit,arange,array,figure
  import os
  HOME = os.getenv('HOME')
  if tfile == '_measlines.txt': tfile = filestub+tfile
  if wheelpos == 200: 
     modedir = 'UVnominal'
     wzmx = array([1800,2000,2300,2600,2900,3200,3500,3800,4000,4200,4500,5000,5500,6000,6500])
     ord1 = 4
  elif wheelpos == 160:
     modedir = 'UVclocked'
     wzmx = array([1800,2000,2300,2600,2900,3200,3500,3800,4000,4200,4500,5000,5500,6000,6500])
     ord1 = 4
  elif wheelpos == 955:      
     modedir='Vclocked'
     wzmx = array([2530,2800,3000,3300,3500,3800,4000,4200,4500,5000,5500,6000,6500]) 
     ord1 = 3
  elif wheelpos == 1000:
     modedir='Vnominal'
     wzmx = array([2530,2800,3000,3300,3500,3800,4000,4200,4500,5000,5500,6000,6500])
     ord1 = 3
  else:
     print('illegal wheelpos (160,200,955,1000)')
     return     
  print("modedir = ", modedir)     
  obsid = filestub[2:]
  if (os.getcwd().split('/')[-1] == 'image' ) :
     os.chdir('../../../'+obsid+'/uvot/image')
     print('processing in directory ', os.getcwd())
  else:
     os.chdir(obsid+'/uvot/image/')
  textfil1='acc_web.txt'
  textfil2='spectrum_web.txt'
  filt1= 'uvw1'
  filt2= 'uvw1'
  if wheelpos > 500:
    p = arange(-310,430,40)
  else:
    p = arange(-370,1000,40)  
  t = rdTab(tfile)
  print(ra,dec,filestub, ext, ext1, ext2) 
  Y = uvotgrism.getSpec(ra,dec,filestub,ext,lfilt1=lfilt1,lfilt1_ext=ext1,lfilt2=lfilt2,\
      lfilt2_ext=ext2,spextwidth=specwidth,chatter=5)
  ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
         (C_1,C_2,img,H_lines,WC_lines),  hdr,m1,m2,aa,wav1 ) = Y
  if hdr['wheelpos'] != wheelpos:
     print(" AARRGGG!!!!!!!  wheelpos in call must be WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")       
  figure(5)
  clf()
  #  
  #  first make sure that the anchor point of the observed scale is correct at C_1[-1]
  #
  tpix = t[:,0]
  twav = t[:,1]
  # find the wavelength at zero pixel
  tzwav = polyfit(tpix,twav,3)
  # find the correction needed
  t_wav_diff = C_1[-1] - tzwav[-1]
  # make the correction but as a pixel difference using the linear term  C_1[-2]
  tpixdiff = t_wav_diff/C_1[-2]
  print("correcting the measured positions to use the same anchor point as the dispersion relation. ")
  print("originally the wavelength of the zero pixel (measured) was: ",tzwav[-1]," while the dispersion gave ", C_1[-1])
  print("pixel coordinate shifted by subtracting ", tpixdiff)
  t[:,0] -= tpixdiff
  print("new wavelength of the zero pixel (measured) is: ",polyfit(tpix,twav,3)[-1])
  #
  #  
  #
  w = polyval(C_1,p) ; RC = polyfit(w,p,4) ; pzmx=polyval(RC,wzmx)

  figure(5)
  clf()
  good = False
  disoff = 0.
  while good == False :
     z = approx_solution(wheelpos,C_1,source=source)
     plotspectrum(dis,C_1, spnet/hdr['exposure'],z,linewidth=15,disoffset=disoff,g_min=3,linecolor='k')
     print("IMPORTANT: line up the lines with the peaks in the f(pix) [bottom] image")
     ans = float(input("give a pixel offset for Figure 5 whuch is needed to fix it or zero to continue [-1.6]: "))
     if int(10*ans) != 0: 
       disoff = ans
       clf()
     else:
       good = True
  
  figure(4)
  clf()
  waveAccPlot(t[:,1],t[:,0]+disoff,wzmx,pzmx,C_1,obsid=filestub+'['+str(ext)+'] '+objectname,figureno=4,\
     legloc=[1,2],order=ord1,wheelpos=wheelpos)
  dir1 = dir2 = HOME+'/Sites/Grism/'+modedir+'/image/'
  fig1=filestub+'acc_'+str(ext)+'.png'
  fig1a=filestub+'acc_'+str(ext)+'.eps'
  fig2=filestub+'spec_'+str(ext)+'.png'
  fig2a=filestub+'spec_'+str(ext)+'.eps'
  savefig(dir1+fig1)
  savefig(dir1+fig1a)
  
  figure(5)
  clf()
  z = approx_solution(wheelpos,C_1,source=source)
  plotspectrum(dis,C_1, spnet/hdr['exposure'],z,linewidth=15,disoffset=0,g_min=3,linecolor='k')
  subplot(211)
  #plot(polyval(C_1,dis[aa]),spnet[aa]/hdr['exposure'],'k',ls='steps')
  title(objectname+'  '+filestub+'_'+str(ext))
  #xlabel('$\lambda(\AA)$',fontsize=16)
  ylabel('count rate')
  subplot(212)
  ylabel('count rate')
  savefig(dir1+fig2)
  savefig(dir1+fig2a)
  coor = str(anker[0])+','+str(anker[1])
  commd1='  <area shape=\"circle\" coords=\"'+coor+',50\" href="../image/'+fig1+' \" alt=\"'+objectname+' '+filestub+'_'+str(ext)+'\" />'
  os.system('echo '+'"""'+commd1+'"""'+' >> '+dir1+textfil1)
  commd2='  <area shape=\"circle\" coords=\"'+coor+',50" href="../image/'+fig2+' \" alt=\"'+objectname+' '+filestub+'_'+str(ext)+'\" />'
  os.system('echo '+'"""'+commd2+'"""'+' >> '+dir2+textfil2)
  os.chdir( '../../..')
  print('next')
  
  
def niceplot():  
  '''  This makes a plot with contours and two inserts: 

  fig3 = figure(6)
  plot(obsdetx,obsdety,'bo')
  uvotplot.contourpk((obs42x[q]+104-1100.5)*0.009075,(obs42y[q]+78-1100.5)*0.009075,dlam1[q],levels=[0.-3,-6,-9,-12,-13,-14,-15,-17,-20])
  ylim(-8,14)
  xlabel('det-X (mm)')
  ylabel('det-Y (mm)')
  title('V grism nominal accuracy of wavelength scale ($\AA$)')
  fig3.add_subplot(3,3,1,axisbg='#00f0f3')
  hist(dlam1[q],bins=11)
  xlabel('$\lambda(\AA)$',fontsize=14,weight='bold')
  ylim(0,9.)
  text(-35,7.5,'error anchor position')
  fig3.add_subplot(3,3,3,axisbg='#00f0f3')
  hist(acc1[q],bins=11)
  ylim(0,11)
  xlabel('$\lambda(\AA)$',fontsize=14,weight='bold')
  text(6.3,8.5,'accuracy')
  savefig('/Volumes/misc/root/Sites/Grism/Vnominal/image/accuracy.png')
  savefig('/Volumes/misc/root/Sites/Grism/Vnominal/image/accuracy.eps')
  '''
  print("no data")
  
