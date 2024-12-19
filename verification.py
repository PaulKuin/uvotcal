from uvotcal import *
import os


# verification programs
# 2013-11-28 split from uvotcal 
#            need to add some kind of driver and store results to compare to
# 2013-12-04 created set of basic wavecal data.

dataroot='/calibration/grism'
caldir = os.getenv('UVOTPY')+'/calfiles/'

def make_test_data():
   """draft for testing the data. Need to add some thing to prevent asking for offset. """
   for n in range(25):
      try: 
          _test_get_wavcal_1(n,resultfile='wavecal_test_withlent_0160.txt',
	            wavecal='wavecal1/with_lent/',
		    fig_on=False, 
		    fancyplot=False,
		    chatter=1,
		    uselent=True)
      except:
         pass
      try: 
          _test_get_wavcal_1(n,resultfile='wavecal_test_nolent_0160.txt',
	            wavecal='wavecal1/with_lent/',
		    fig_on=False, 
		    fancyplot=False,
		    chatter=1,
		    uselent=False)
      except:
         pass

   for n in range(30):
      try: 
          _test_get_wavcal_2(n,resultfile='wavecal_test_withlent_0200.txt',
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=True)
      except:
         pass
      try: 
          _test_get_wavcal_2(n,resultfile='wavecal_test_nolent_0200.txt',
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=False)
      except:
         pass
	 
   for n in range(21):
      try: 
          _test_get_wavcal_3(n,resultfile='wavecal_test_withlent_0955.txt',
	            mo=0955,
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=True)
      except:
         pass
      try: 
          _test_get_wavcal_3(n,resultfile='wavecal_test_nolent_0955.txt',
	            mo=955,
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=False)
      except:
         pass
	 
   for n in range(34):
      try: 
          _test_get_wavcal_3(n,resultfile='wavecal_test_withlent_1000.txt',
	            mo=1000,
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=True)
      except:
         pass
      try: 
          _test_get_wavcal_3(n,resultfile='wavecal_test_nolent_1000.txt',
	            mo=1000,
	            wavecal='wavecal1/with_lent/', 
		    chatter=1,
		    lent=False)
      except:
         pass
	 
def compare_testdata(date='now'):
   import os   
   try:
      resultfile='wavecal_test_withlent_'+date+'_0160.txt'
      os.system("diff wavecal_withlent_4dec2013_0160.txt "+resultfile+" > diff_nolent_160.txt")      
   except: 
      pass	 
   try:
      resultfile='wavecal_test_withlent_'+date+'_0200.txt'
      os.system("diff wavecal_withlent_4dec2013_0200.txt "+resultfile+" > diff_nolent_200.txt")      
   except: 
      pass	 
   try:
      resultfile='wavecal_test_withlent_'+date+'_0955.txt'
      os.system("diff wavecal_withlent_4dec2013_0955.txt "+resultfile+" > diff_nolent_955.txt")      
   except: 
      pass	 
   try:
      resultfile='wavecal_test_withlent_'+date+'_1000.txt'
      os.system("diff wavecal_withlent_4dec2013_1000.txt "+resultfile+" > diff_nolent_1000.txt")      
   except: 
      pass	 
   # examine the diff files    
      
	 
def _test_get_wpixscale(n,chatter=0):
   """Testing the new dispersion for the uv clocked grism
      
      2013-05-15 NPMK
      
   """	 
   import cal3
   import numpy as np
   from pylab import figure,plot,legend,title,subplot,ylabel,xlabel,xlim,savefig
   
   dis1 = cal3.get_curvedata_aslist("dis1",curvedata=cal3.uv0160)
   dis1corr = (cal3.get_curvedata_asarray("dis1corr",curvedata=cal3.uv0160)).flatten()
   obsid = cal3.get_curvedata_aslist("obsid",curvedata=cal3.uv0160)
   anchor = cal3.get_curvedata_asarray("anchor",curvedata=cal3.uv0160)
   field = cal3.get_curvedata_asarray("field",curvedata=cal3.uv0160)
   dis = np.array(dis1[n])
   d = dis[:,0] + dis1corr[n]
   w = dis[:,1]
   ankx,anky = anchor[0,n],anchor[1,n]
   Xfield,Yfield=field[0,n],field[1,n]
   test2 = ['full','fix1','fix2','linear']
   disp1 = []
   for tt in test2:
       print "************        test2 %s    **************"%(tt)
       print "************        test2 %s    **************"%(tt)
       print "************        test2 %s    **************"%(tt)
       C_zero, C_1, C_2, C_3, flux, xpix_, ypix_ , dd, ww, zmxwav, theta, = getZmx(
           Xfield,Yfield,160,mode='spline',test='simple_corrected',test2=tt,
	   wpixscale=None, chatter=chatter)
       disp1.append(C_1)
   print disp1    	   
   figure()
   subplot(2,1,1)
   plot(d,w,'oy',markersize=3,label='observed')
   plot(d,polyval(disp1[0],d),'vg-',markersize=4,label='linear+2nd')
   plot(d,polyval(disp1[1],d),'^r-',markersize=4,label='fix 2nd')
   plot(d,polyval(disp1[2],d),'>b-',markersize=4,label='linear+fix 2nd')
   plot(d,polyval(disp1[3],d),'<c-',markersize=4,label='linear')
   legend(loc=4)
   title(u'dispersion scaled zemax model %s anchor = %7.1f,%7.1f'%(obsid[n],ankx,anky))
   ylabel(u'derived $\lambda(\AA)$')
   xlim(-400,1100)
   subplot(2,1,2)
   plot(d,polyval(disp1[0],d)-w,'vg-',markersize=4,label='linear+2nd')
   plot(d,polyval(disp1[1],d)-w,'^r-',markersize=4,label='fix 2nd')
   plot(d,polyval(disp1[2],d)-w,'>b-',markersize=4,label='linear+fix 2nd')
   plot(d,polyval(disp1[3],d)-w,'<c-',markersize=4,label='linear')
   ylabel(u'difference to observed $\lambda(\AA)$')
   legend(loc=0)
   xlim(-400,1100)
   savefig('/Volumes/users/Users/kuin/Desktop/calibration-recent/wavecal_160/disp_check_%4i_%4i.png'%(ankx,anky))
   
	 
def _test_get_wavcal_1(n,resultfile=None,reportlimit=None,wavecal='wavecal1',
   fig_on=False,fancyplot=False,uselent=None,radec=None,chatter=0):
   """
      Testing the new anchor and dispersion for the uv clocked grism (0160)
                  
      input parameters
      ==================
      n : int
         n is the number of the observation in the CAL3 list
      resultfile : path
        resultfile is the file where a summary will be appended to. 
	parameters dumped are:
        target,obsid,ext,wheelpos,anchor_obs,anchorX,anchorY,rms_blue,rms_centre,rms_red
      
      reportlimit : float
         reportlimit will flag observations with anchor error larger than the limit given (A)
      
      wavecal : str
         wavecal selects the tree of data (with/without lenticular files)
	 
      fig_on : bool
         make error plots
      
      fancyplot: bool
         make fancy publication style error plot
	 
      chatter: int, [0...5]
      
         verbosity
      
      output parameters
      ==================
      
      notes
      ======
      
      2013-05-15 NPMK
      2013-08-28 NPMK add output to file
      
   """	
   import os 
   import numpy as np
   from pylab import figure,plot,legend,title,subplot,ylabel,xlabel,xlim,ylim,\
        savefig,grid,text,errorbar,polyval,hlines
   from astropy.io import fits
   from astropy import coordinates as coord
   import cal3
   import uvotgetspec
   
   if uselent == None:
      uselent = rawinput("lenticular image parameter uselent needs to be set (True/False)")
      if chatter > 0: 
         print uselent
   msg = " "
   uvotgetspec.give_result=True
   #caldir = '/Volumes/users/Users/kuin/dev/uvotpy.latest/calfiles/'
   calfile = 'swugu0160wcal20041120v002.fits'
   dis1 = cal3.get_curvedata_aslist("dis1",curvedata=cal3.uv0160)
   dis1corr = (cal3.get_curvedata_asarray("dis1corr",curvedata=cal3.uv0160)).flatten()
   obsid = cal3.get_curvedata_aslist("obsid",curvedata=cal3.uv0160)
   name = cal3.get_curvedata_aslist("obsdir",curvedata=cal3.uv0160)
   ext = (cal3.get_curvedata_asarray("ext",curvedata=cal3.uv0160)).flatten()
   anchor = cal3.get_curvedata_asarray("anchor",curvedata=cal3.uv0160)
   obj = name[n]
   if radec == None:
      if obj.upper() == "WR52":
         ra,dec = cal3.WR52_position
      elif obj.upper() == "WR86":
         ra,dec = cal3.WR86_position
      else:	 	 
         if obj.upper() == "TPYX": 
            obj = "T Pyx"	 
         position = coord.ICRSCoordinates.from_name(obj)
         ra = position.ra.degrees
         dec = position.dec.degrees
         #ra,dec = uvotgetspec.get_radec(objectid=obj)
   else:
      ra,dec = radec   
   datadir = dataroot+"/"+wavecal+"/"+name[n]+"/"+obsid[n]+"/uvot/image/"
   dis1 = np.array(dis1[n])
   d = dis1[:,0] + dis1corr[n]
   w = dis1[:,1]
   ankx,anky = anchor[0,n],anchor[1,n]
   msg += "%12s %10s %02i 0160 ank_obs= %7.1f %7.1f "%(obj,obsid[n],ext[n], ankx,anky)
   print "input parameters to getSpec:"
   print "RA = %s, Dec = %s"%(str(ra),str(dec))
   print "observed anchor (%7.1f,%7.1f)"%(ankx,anky)
   print "obsid = "+obsid[n]
   print "ext = ",ext[n]
   print "datadir = "+datadir
   print "calfile = "+caldir+calfile
   #cmd = "fthedit '"+datadir+"sw"+obsid[n]+"ugv_dt.img["+str(ext[n])+"]'   keyword=ASPCORR operation=replace value=NONE"
   #os.system(cmd)
   fitsec = anky > 1250   
   if fitsec: print "fitting second order"
   interact = True
   if interact: print "interactive on"
   Z = uvotgetspec.getSpec(ra,dec,obsid[n],ext[n], fit_second=fitsec, 
        indir=datadir,
        calfile=caldir+calfile, interactive=interact, 
	clobber=True, plot_img=False,
	plot_spec=False, plot_raw=False, 
	use_lenticular_image=uselent,
	chatter=chatter)
   print "call made to uvotgetspec.getSpec()"	
   #(specfile, lfilt1_, lfilt1_ext_, lfilt2_, lfilt2_ext_, attfile), (method), \
   #    (Xphi, Yphi, date1), (dist12, ankerimg, ZOpos), expmap, bgimg, bg_limits_used, bgextra = Z[0]
   ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
       (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Z[1]
   msg += "ank= %7.1f %7.1f "%(anker[0],anker[1])    
   #fit,(coef0,coef1,coef2,coef3),(bg_zeroth,bg_first,bg_second,bg_third),\
   #    (borderup,borderdown),apercorr,expospec= Z[2]	
   #counts, variance, borderup, borderdown, (fractions,cnts,vars,newsigmas) = Z[3] 
   #wav2p, dis2p, flux2p, qual2p, dist12p = Z[4][0]        
   #(present0,present1,present2,present3),(q0,q1,q2,q3), \
   #(y0,dlim0L,dlim0U,sig0coef,sp_zeroth),(y1,dlim1L,dlim1U,sig1coef,sp_first),\
   #(y2,dlim2L,dlim2U,sig2coef,sp_second),(y3,dlim3L,dlim3U,sig3coef,sp_third),\
   #(x,xstart,xend,sp_all,quality)  = fit      
   
   dxy_ank = [anker[0]-104-ankx,anker[1]-78-anky]
   msg += " delta=%6.1f %6.1f "%(dxy_ank[0],dxy_ank[1])
   print "Anchor difference in pixels (new - observed) %s"%(dxy_ank)
   offsetout='160_anchor_off_uvotgraspcorr-wcs.2.txt'
   f = open(offsetout,'a')
   f.write("%6.1f %6.1f %6.1f %6.1f\n"%(ankx,anky,dxy_ank[0],dxy_ank[1]))
   f.close()
   if fig_on:
      figure()
      subplot(2,1,1)
      plot(d,w,'+b',markersize=8,label='observed')
      plot(d,polyval(C_1,d),'^r-',markersize=4,label='getSpec()')
      legend(loc=4)
      title(u'test wavecal on %s+%i anchor = %7.1f,%7.1f'%(obsid[n],ext[n],ankx,anky))
      text(70,2600,u'wavecal: '+calfile,fontsize=10)
      ylabel(u'derived $\lambda(\AA)$')
      grid()
      xlim(-400,1100)
      subplot(2,1,2)
      plot(d,polyval(C_1,d)-w,'+b-',markersize=8,label='observed-getSpec()')
      ylabel(u'difference to observed $\lambda(\AA)$')
      legend(loc=0)
      xlim(-400,1100)
      grid()
      savefig('uvotgraspcorr_check_%4i_%4i_20130612.png'%(ankx,anky))
   if fancyplot:
      figure()
      subplot(2,1,1)
      #
      #plot(d,w-C_1[-2]*d,'ob',markersize=7,label='observed')
      #plot(d,polyval(C_1,d)-C_1[-2]*d,'^r',markersize=8,label='model')
      #legend(loc=4)
      #ylabel(u'$\lambda(\AA)$ - %5.1f d '%(C_1[-2]))
      #xlim(-390,970)
      plot(w,d-(w-2600.)/C_1[-2],'ob',markersize=7,label='observed',alpha=0.7) 
      plot(polyval(C_1,d),d-(w-2600.)/C_1[-2] ,'^r',markersize=8,label='model',alpha=0.7)
      legend(loc=1)
      ylabel(u'distance - ($\lambda$-2600)/%5.1f (pixels)'%(C_1[-2]),fontsize=14)
      ylabel(u'$\Delta$d (pixels)',fontsize=14)
      xlim(1670,6120)
      ylim(-180,30)
      
      subplot(2,1,2)
      measurement_error = polyval(C_1,d+0.5)-polyval(C_1,d-0.5)   # assume 1 pixels measurement error
      err = np.ones(len(d),dtype=float) * measurement_error 
      err[-1] = 1.67*err[-1] # peak position of 5805 feature difficult to measure
      delta_w = polyval(C_1,d)-w
      errorbar(w,delta_w,yerr=err,fmt='ob',markersize=7,label='observed-model')
      hlines(delta_w[:-1].mean(),1675,6015,linestyle='dashed',lw=1.3,colors='purple',label='mean error')
      #errorbar(d,polyval(C_1,d)-w,yerr=err,fmt='ob',markersize=7,label='observed-model')
      ylabel(u'$\Delta\lambda(\AA)$',fontsize=14)
      legend(loc=0)
      xlim(1670,6120)
      xlabel(u'$\lambda (\AA)$',fontsize=14)
      #savefig('./fancy_dispersion_error_%4i_%4i_date.eps'%(ankx,anky),dpi=500,bbox_inches=3.506)  # 88 mm
      print "saved eps plot in current directory "
   dwv = polyval(C_1,d)-w # delta-wave for all observations
   q = [w<2000,(w>2000) & (w<4500), w>4500]
   rms_blue = dwv[q[0]].std()
   rms_mid = dwv[q[1]].std()
   rms_red = dwv[q[2]].std()
   msg += "%7.2f %7.2f %7.2f"%(rms_blue,rms_mid,rms_red)
   os.system('rm *.gat.fits attcorr.asp radec.txt detpix.out radec.usnofull radec.usno radec.sesame search.ub1')
   if resultfile != None:   
      os.system('echo "'+msg+'" >> '+resultfile)
   
def _test_get_wavcal_2(n,resultfile=None,reportlimit=None,lent=True,
     wavecal='wavecal1',plotit=False,radec=None,chatter=0):
   """Testing the new anchor and dispersion for the uv nominal grism
      
      n is the number of the observation in the CAL3 list
      
      resultfile is the file where a summary will be appended to
        target,obsid,ext,wheelpos,anchor_obs,anchorX,anchorY,rms_blue,rms_centre,rms_red
      
      reportlimit will flag observations with anchor error larger than the limit given (A)
      
      wavecal selects the tree of data (with/without lenticular files)
      
      radec = list of ra, dec in degrees
      
      2013-05-15 NPMK
      2013-08-07       
   """	
   import os 
   import numpy as np
   from pylab import figure,plot,legend,title,subplot,ylabel,xlabel,xlim,savefig,grid,text
   from astropy.io import fits
   from astropy import coordinates as coord
   import cal3
   import uvotgetspec
   from uvotio import fileinfo
   
   msg = ""
   logx = "= = = = = = = = = = = \nuvotcal._test_get_wavecal("
   logx += str(n)+",resultfile"+str(resultfile)+",\nreportlimit="+str(reportlimit)
   logx += ",\nwavecal="+wavecal
   logx += ",\nradec=%s\nlent=%s\n"%(radec,lent,)

   uvotgetspec.give_result=True
   #caldir = '/home/kuin/dev/uvotpy.garfield/calfiles/'
   calfile = 'swugu0200wcal20041120v001.fits'
   dis1 = cal3.get_curvedata_aslist("dis1",curvedata=cal3.uv0200)
   dis1corr = (cal3.get_curvedata_asarray("dis1corr",curvedata=cal3.uv0200)).flatten()
   obsid = cal3.get_curvedata_aslist("obsid",curvedata=cal3.uv0200)
   name = cal3.get_curvedata_aslist("obsdir",curvedata=cal3.uv0200)
   ext = (cal3.get_curvedata_asarray("ext",curvedata=cal3.uv0200)).flatten()
   anchor = cal3.get_curvedata_asarray("anchor",curvedata=cal3.uv0200)
   obj = name[n]
   print "OBJECT: "+obj
   print "dis1: ", dis1[n]
   if radec == None:
      if obj.upper() == "WR52":
         ra,dec = cal3.WR52_position
      elif obj.upper() == "WR86":
         ra,dec = cal3.WR86_position
      else:	 	 
         if obj.upper() == "TPYX": 
            obj = "T Pyx"	 
         position = coord.ICRSCoordinates.from_name(obj)
         ra = position.ra.degrees
         dec = position.dec.degrees
         #ra,dec = uvotgetspec.get_radec(objectid=obj)
   else:
      ra,dec = radec   
   print ra,dec
   datadir = dataroot+"/"+wavecal+"/"+name[n]+"/"+obsid[n]+"/uvot/image/"
   print "data dir: "+datadir
   specfile, lfilt1, lfilt1_ext, lfilt2, lfilt2_ext, attfile = fileinfo(
         "sw"+obsid[n],
	 ext[n],
	 directory=datadir,
	 wheelpos=200,
	 chatter=chatter)
   lfilt1_aspcorr = 0
   lfilt2_aspcorr = 0
   print "fileinfo: ",specfile, lfilt1, lfilt1_ext, lfilt2, lfilt2_ext, attfile 
   if lfilt1 != None:
      try: 
         print "try loading "+datadir+"/sw"+obsid[n]+"u"+lfilt1[-2:]+"_sk.img"
         hdu_1 = pyfits.getheader(datadir+"/sw"+obsid[n]+"u"+lfilt1[-2:]+"_sk.img",lfilt1_ext)
      except:
         print "try loading "+datadir+"/sw"+obsid[n]+"u"+lfilt1[-2:]+"_sk.img.gz"
         hdu_1 = pyfits.getheader(datadir+"/sw"+obsid[n]+"u"+lfilt1[-2:]+"_sk.img.gz",lfilt1_ext)
      if hdu_1["ASPCORR"] != "NONE": lfilt1_aspcorr = 1   
      print "lenticular filter 1 aspect correction is %s "%(hdu_1["ASPCORR"])	 
   if lfilt2 != None:
      try: 
         hdu_2 = pyfits.getheader(indir+"/sw"+obsid+"u"+lfilt2[-2:]+"_sk.img",lfilt2_ext)
      except:
         hdu_2 = pyfits.getheader(indir+"/sw"+obsid+"u"+lfilt2[-2:]+"_sk.img.gz",lfilt2_ext)	 
      print "lenticular filter2 aspect correction is %s "%(hdu_2["ASPCORR"])	 
      if hdu_2["ASPCORR"] != "NONE": lfilt2_aspcorr = 1   
   if dis1[n] == None:   # jury rig some kind of dispersion here NEED TO ADD TO CAL3
      print "DIS1 dummy values"
      dis1 = np.array([[-364, 1642], [-303.6, 1764], [-254.4, 1909], 
      [-101.1, 2295.0], [-66.3, 2400], [-25.1, 2525], [25.3, 2693], 
      [40.9, 2729], [89.5, 2900.0], [111.6, 2980], [174.1, 3200], 
      [230.4, 3406], [313.1, 3686], [565.5, 4649.0], [835.2, 5687], 
      [865.3, 5801], [877.4, 5860]])
      dis1cor = 0
   else:   
      dis1 = np.array(dis1[n])
      dis1cor = dis1corr[n]
   d = dis1[:,0] + dis1cor
   w = dis1[:,1]
   ankx,anky = anchor[0,n],anchor[1,n]
   msg += "%12s %10s %02i 0200 ank_obs= %7.1f %7.1f "%(obj,obsid[n],ext[n], ankx,anky)
   print "input parameters to getSpec:"
   print "RA = %s, Dec = %s"%(str(ra),str(dec))
   print "observed anchor (%7.1f,%7.1f)"%(ankx,anky)
   print "obsid = "+obsid[n]
   print "ext = ",ext[n]
   print "datadir = "+datadir
   print "calfile = "+caldir+calfile
   logx += "input parameters to getSpec:\n"
   logx += "RA = %s, Dec = %s\n"%(str(ra),str(dec))
   logx += "observed anchor (%7.1f,%7.1f)\n"%(ankx,anky)
   logx += "obsid = "+obsid[n]+"\n"
   logx += "ext = %i \n"%(ext[n])
   logx += "datadir = "+datadir+"\n"
   logx += "calfile = "+caldir+calfile+"\n"
   #cmd = "fthedit '"+datadir+"sw"+obsid[n]+"ugu_dt.img["+str(ext[n])+"]'  keyword=ASPCORR operation=replace value=NONE"
   #os.system(cmd)
   #logx += cmd+"\n"
   fitsec = anky > 1250   
   if fitsec: 
      print "fitting second order"
      logx +=  "fitting second order\n"
   interact = True
   if interact: 
      print "interactive on" 
      logx +="interactive on\n"
   Z = uvotgetspec.getSpec(ra,dec,obsid[n],ext[n], fit_second=fitsec, 
        indir=datadir,curved='no',
        calfile=caldir+calfile, interactive=interact, 
	clobber=True, plot_img=False,
	plot_spec=False, plot_raw=False, 
	use_lenticular_image=lent,
	chatter=chatter)
   print "call made to uvotgetspec.getSpec()"	
   ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
       (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Z[1]
   
   dxy_ank = [anker[0]-104-ankx,anker[1]-78-anky]
   print "Anchor difference in pixels (new - observed) %s"%(dxy_ank)
   logx += "Anchor difference in pixels (new - observed) %s\n"%(dxy_ank)
   msg += "ank= %7.1f %7.1f "%(anker[0],anker[1])    
   msg += " delta=%6.1f %6.1f "%(dxy_ank[0],dxy_ank[1])
   msg += " asp_lent= %i %i "%(lfilt1_aspcorr,lfilt2_aspcorr)
   if lent:
      offsetout='200_anchor_with_lent_uvotgraspcorr'
   else:
      offsetout='200_anchor_no_lent_uvotgraspcorr'
   
   f = open(offsetout,'a')
   f.write("%6.1f %6.1f %6.1f %6.1f\n"%(ankx,anky,dxy_ank[0],dxy_ank[1]))
   f.close()
   if plotit: 
      figure()
      subplot(2,1,1)
      plot(d,w,'+b',markersize=8,label='observed')
      plot(d,polyval(C_1,d),'^r-',markersize=4,label='getSpec()')
      legend(loc=4)
      title(u'test wavecal on %s+%i anchor = %7.1f,%7.1f'%(obsid[n],ext[n],ankx,anky))
      text(70,2600,u'wavecal: '+calfile,fontsize=10)
      ylabel(u'derived $\lambda(\AA)$')
      grid()
      xlim(-400,1100)
      subplot(2,1,2)
      plot(d,polyval(C_1,d)-w,'+b-',markersize=8,label='observed-getSpec()')
      ylabel(u'difference to observed $\lambda(\AA)$')
      legend(loc=0)
      xlim(-400,1100)
      grid()
      savefig('/Volumes/users/Users/kuin/Desktop/calibration-recent/wavecal_200/uvotgraspcorr_check_%4i_%4i_20130612.png'%(ankx,anky))

   dwv = polyval(C_1,d)-w # delta-wave for all observations
   q = [w<2000,(w>2000) & (w<4500), w>4500]
   rms_blue = dwv[q[0]].std()
   rms_mid = dwv[q[1]].std()
   rms_red = dwv[q[2]].std()
   msg += "%7.2f %7.2f %7.2f"%(rms_blue,rms_mid,rms_red)
   os.system('rm *.gat.fits attcorr.asp radec.txt detpix.out radec.usnofull radec.usno radec.sesame')
   if resultfile != None:   
      os.system('echo "'+msg+'" >> '+resultfile)
   f = open(obsid[n]+"_"+str(ext[n])+".log","a")
   f.write(logx)
   f.close()
   
   
def _test_get_wavcal_3(n,resultfile=None,reportlimit=None,wavecal='wavecal1',chatter=0,
    mo=955,plotit=False,withlent=True,radec=None):
   """Testing the new anchor and dispersion for the vis nominal & clocked grism
     
      n is the number of the observation in the CAL3 list
      
      resultfile is the file where a summary will be appended to
        target,obsid,ext,wheelpos,anchor_obs,anchorX,anchorY,rms_blue,rms_centre,rms_red
      
      reportlimit will flag observations with anchor error larger than the limit given (A)
      
      wavecal selects the tree of data (with/without lenticular files)
      
      2013-05-15 NPMK
      2013-08-07/28 NPMK add output       
   """	 
   import os
   import numpy as np
   from pylab import figure,plot,legend,title,subplot,ylabel,xlabel,xlim,savefig,grid,text
   from astropy.io import fits
   import cal3
   import uvotgetspec
   
   msg = " "
   uvotgetspec.give_result=True
   
   if mo == 1000:
      calfile = 'swugv1000wcal20041120v001.fits'
      dis1 = cal3.get_curvedata_aslist("dis1",curvedata=cal3.vis1000)
      dis1corr = (cal3.get_curvedata_asarray("dis1corr",curvedata=cal3.vis1000)).flatten()
      obsid = cal3.get_curvedata_aslist("obsid",curvedata=cal3.vis1000)
      name = cal3.get_curvedata_aslist("obsdir",curvedata=cal3.vis1000)
      ext = (cal3.get_curvedata_asarray("ext",curvedata=cal3.vis1000)).flatten()
      anchor = cal3.get_curvedata_asarray("anchor",curvedata=cal3.vis1000)
   elif mo == 955:   
      calfile = 'swugv0955wcal20041120v001.fits'
      dis1 = cal3.get_curvedata_aslist("dis1",curvedata=cal3.vis0955)
      dis1corr = (cal3.get_curvedata_asarray("dis1corr",curvedata=cal3.vis0955)).flatten()
      obsid = cal3.get_curvedata_aslist("obsid",curvedata=cal3.vis0955)
      name = cal3.get_curvedata_aslist("obsdir",curvedata=cal3.vis0955)
      ext = (cal3.get_curvedata_asarray("ext",curvedata=cal3.vis0955)).flatten()
      anchor = cal3.get_curvedata_asarray("anchor",curvedata=cal3.vis0955)
   obj = name[n]
   if radec == None:
      if obj.upper() == 'WR121':
         ra,dec = cal3.WR121_position
      else:
         ra,dec = uvotgetspec.get_radec(objectid=obj)
   else:
      ra,dec = radec   
   if withlent:
      datadir = dataroot+"/"+wavecal+"/"+name[n]+"/"+obsid[n]+"/uvot/image/"
   else:
      datadir = dataroot+"/"+wavecal+"/"+name[n]+"/"+obsid[n]+"/uvot/image/"
   if dis1 == None:   # jury rig some kind of dispersion here NEED TO ADD TO CAL3
      dis1=np.array([[ -102.3,  3609.0],[  -22.3,  4070.0],
      [   11.0,  4267.0],[   74.9,  4649.0],[  153.5,  5130.0],
      [  241.9,  5696.0],[  271.4,  5880.0],[  380.3,  6580.0],])
      dis1cor = 0
   else:   
      dis1 = np.array(dis1[n])
      dis1cor = dis1corr[n]
   d = dis1[:,0] + dis1cor
   w = dis1[:,1]
   ankx,anky = anchor[0,n],anchor[1,n]
   msg += "%12s %10s %02i %04i ank_obs= %7.1f %7.1f "%(obj,obsid[n],ext[n],mo, ankx,anky)
   print "input parameters to getSpec:"
   print "RA = %s, Dec = %s"%(str(ra),str(dec))
   print "observed anchor (%7.1f,%7.1f)"%(ankx,anky)
   print "obsid = "+obsid[n]
   print "ext = ",ext[n]
   print "datadir = "+datadir
   print "calfile = "+caldir+calfile
   f = fits.open(datadir+"sw"+obsid[n]+"ugv_dt.img",mode="update")
   f[int(ext[n])].header['aspcorr'] = 'None'
   f.flush()
   f.close()
   fitsec = False # anky > 1250   
   if fitsec: print "fitting second order"
   interact = True
   if interact: print "interactive on"
   Z = uvotgetspec.getSpec(ra,dec,obsid[n],ext[n], fit_second=fitsec, 
             indir=datadir,curved='no',
             calfile=caldir+calfile, 
	     interactive=interact, 
	     clobber=True, 
	     plot_img=plotit,
	     use_lenticular_image=withlent,
	     plot_spec=plotit, 
	     plot_raw=plotit, 
	     chatter=chatter, 
	     predict2nd=False)	
   print "call made to uvotgetspec.getSpec()"	
   ( (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
       (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Z[1]
   try:    
      hdr2 = fits.getheader(datadir +"sw"+obsid[n]+"ugv_dt.img",int(ext[n]))   
      print "(gr)aspcorr: ",hdr2['aspcorr']
      if hdr2['aspcorr'] == 'GRASPCORR' :
         asp = 1
      else:
         asp = 0
   except:
      asp = -1
      pass	 
   dxy_ank = [anker[0]-104-ankx,anker[1]-78-anky]
   msg += "ank= %7.1f %7.1f "%(anker[0],anker[1])    
   msg += " delta=%6.1f %6.1f "%(dxy_ank[0],dxy_ank[1])
   print "Anchor difference in pixels (new - observed) %s"%(dxy_ank)
   if mo == 1000:
      if withlent:
         offsetout ='1000_anchor_lent'
      else:
         offsetout='1000_anchor_off_uvotgraspcorr'
   elif mo == 955:   
      if withlent:
         offsetout = '955_anchor_lent'
      else:
         offsetout='955_anchor_off_uvotgraspcorr'
   f = open(offsetout,'a')
   f.write("%6.1f %6.1f %6.1f %6.1f  %i\n"%(ankx,anky,dxy_ank[0],dxy_ank[1], asp))
   f.close()
   if plotit:
      figure()
      subplot(2,1,1)
      plot(d,w,'+b',markersize=8,label='observed')
      plot(d,polyval(C_1,d),'^r-',markersize=4,label='getSpec()')
      legend(loc=4)
      title(u'test wavecal on %s+%i anchor = %7.1f,%7.1f'%(obsid[n],ext[n],ankx,anky))
      text(70,2600,u'wavecal: '+calfile,fontsize=10)
      ylabel(u'derived $\lambda(\AA)$')
      grid()
      xlim(-400,1100)
      subplot(2,1,2)
      plot(d,polyval(C_1,d)-w,'+b-',markersize=8,label='observed-getSpec()')
      ylabel(u'difference to observed $\lambda(\AA)$')
      legend(loc=0)
      xlim(-400,1100)
      grid()
   if mo == 1000:
      if not withlent and plotit:
         savefig('/Users/kuin/Desktop/calibration-recent/wavecal_1000/uvotgraspcorr_check_%4i_%4i_20130612.png'%(ankx,anky))
   elif mo == 955:  
      if not withlent and plotit: 
         savefig('/Users/kuin/Desktop/calibration-recent/wavecal_955/uvotgraspcorr_check_%4i_%4i_20130612.png'%(ankx,anky))
   dwv = polyval(C_1,d)-w # delta-wave for all observations
   q = [w<3100,(w>3100) & (w<5500), w>5500]
   rms_blue = dwv[q[0]].std()
   rms_mid = dwv[q[1]].std()
   rms_red = dwv[q[2]].std()
   msg += "%7.2f %7.2f %7.2f"%(rms_blue,rms_mid,rms_red)
   os.system('rm *.gat.fits attcorr.asp radec.txt detpix.out radec.usnofull radec.usno radec.sesame')
   if resultfile != None:   
      os.system('echo "'+msg+'" >> '+resultfile)

def _check_anchor_disp_zmx(mode=0,db=None):
   """ 
   compare the anchor and dispersion from the calibration observations to the model
   for mode = 0 prior to writing a calibration file
   for mode = 1 use the calibration file
   
   2013-06-12 updated to use calibration file 
   """
   import cal3
   import numpy as np
   import uvotgetspec
   
   z = [1100.5,1100.5]         # image center in det coordinates 
   o = [1100.5-104,1100.5-78]  # image center in physical image coordinates (same for all grism images)
   sc = []
   a = []
   # loop over wavecal observations
   if db == None: db = cal3.uv0160 
   for this in db:
      Xfield,Yfield = this['field']
      C_zero, C_1, C_2, C_3, flux, xpix_, ypix_ , dd, ww, zmxwav, theta = \
           getZmx(Xfield,Yfield,160,mode='spline',chatter=3,wpixscale=1., \
           wpixscale2=1., kk=3,s=None,test='simple',xlimits=[10,2100])
      k = np.where(zmxwav == 0.260)
      zmx_anchor = [ xpix_[1,k], ypix_[1,k] ]	
      if mode == 0:
         C_1 = this['dispers1st']
         anchor = this['anchor']
         # simple zemax model anchor, dispersion, wavelengths
         zmx_dis1   = dd[1]
         zmx_wav1st = ww[1]*1.0e4
         obs_dis1 = np.polyval(C_1,zmx_wav1st)
         # get list of dispersion scaling to zmx
         sc.append( np.polyfit(zmx_dis1,obs_dis1,2) )
         # get list of anchor relative to image center, first observed, then simple zmx model   
         a.append( [anchor[0]-o[0],anchor[1]-o[1], zmx_anchor[0]-z[0],zmx_anchor[1]]-z[1])
      elif mode == 1:
      	 name = this['obsdir']
         if name.upper() == "TPYX": name = "T Pyx"
	 ra,dec = uvotgetspec.get_radec(objectid=name)
	 obsid = this['obsid']
	 ext = this['ext']
	 uvotgetspec.give_result = True
	 indir = dataroot+'/wavecal/'+this['obsdir']+'/'+obsid+'/uvot/image/'
	 Z = uvotgetspec.getSpec(ra,dec,obsid,ext,indir=indir,)
	 # list for each calibration observation the dispersion in sc and anchor in a
	 # for each observation the list order is calfile value, getZmx value, observed value
	 sc.append( [Z[1][2][0],C_1,this['dispers1st'] ] )
	 a.append( [Z[1][0][3],zmx_anchor,this['anchor'] ] )
   return a, sc      

# end verification programs 
