# -*- coding: iso-8859-1 -*-
'''
cal4.py   Flux Calibration : observations  

This module contains the files and parameters used for the calibration. 
There are dependencies on the file system.

calibration file data are entered and updates here. Putting the updates 
in another file did not work very well (had to restart, recompiling did 
not work).

2021-01-06 added code to read uvot obsids from local archive

'''
__version__ = 1.4


''' Here I dump some higher level description of the first reduction  
   
   * initial reduction 
   (1) list file info: ftlist sw00055200028ugu_dt.img hk include=exposure,aspcorr,pa_pnt,date-obs
   (2) reduce the data using default (TRACKWID= 1.0) 
   (3) enter the observation in cal4.py and grism_info_***.*** 
   (4) copy the *.pha file to eff_area_***/all_1.0
   (5) set uvotgrism.trackwidth = 2.5
   (6) reduce the data again
   (7) copy the *.pha file to eff_area_***/all_2.5 
   
   (1) stage data using _stage_from_archive
   (2) process to spectrum files, e.g., 
       cd /calibration/grism/fluxcal3/pha15
       cal4.main_ea(wheelpos,make_pha_files=True,coi_width=15)
   (3) run plot_all_spectra() to find which ones failed automatic processing
   (4) adjust wavelength scale using uvotspec.adjust_wavelength_manually() 
   (5) make individual ea files ; set 
       
    
   
'''


def main_ea(wheelpos,key='groups',keyvalue='A',debug=0,phafiledir=None, 
    eff_area_outdir=None, calobs_err=2, apercorr_err=5, xyrange=[100,100],
    ign_q=False,auto=False, use_obs=None, wlshift=0.,
    phatype='_f',
    make_pha_files=False, coi_width=None, use_lenticular_image=True,
    chatter=1, clobber=True):
   '''main calls for doing coi/effective area analysis 
   
      
      what is does:
     
      when auto = False:  run the effective area computation interactively
      else, do it automatically using the exclude regions. 
       
      debug =1: list the calibration file parameters, and parts 
      of the spectrum to exclude for the first order flux calibration.
      calfildir,calfiles,calfilID,obs_160,obs_200,obs_955,obs_1000,group, caldict, obsdict
      
      Parameter: 
         select sets to process from list based on key, and keyvalue
         
      If keyvalue='name', select multiple names like
         keyvalue = 'name1' or 'name2' or 'name3'

      ign_q: ignore the data quality from the PHA file
      
      wlshift : 
        when 0. set all wavelength offsets to zero
        when None use the correction in the database. 
   
   '''   
   from . import uvotcal
   import os
   import time
   from . import cal6
   from sgwgs import retrieve_obsid_path_from_local_archive
   
   obs_955 = cal6.obs_955
   obs_1000 = cal6.obs_1000
   
   if chatter > 0: 
      print("cal4.main_ea(",wheelpos,", key=",key,", keyvalue=",keyvalue,", debug=",debug)
      print("   , phafiledir=",phafiledir,", calobs_err=",calobs_err,", apercorr_err=",apercorr_err)
      print("   , xyrange=",xyrange,", auto=",auto,", chatter=",chatter," clobber=",clobber)
      
   tt = time.struct_time(time.gmtime())
   date = "%4i%2s%2sT%2s%2s"%(tt.tm_year,str(tt.tm_mon).zfill(2),str(tt.tm_mday).zfill(2),str(tt.tm_hour).zfill(2),str(tt.tm_min).zfill(2))
    
   datadir = '/calibration/grism/fluxcal2/'
   #datadir2 = retrieve_obsid_path_from_local_archive(obsid)   # local archive 
   obsdirs = ['P041C','P330E','P177D','GD153','WD0320-539','WD1057+719',
              'WD1657+343','WD1121+145','G63-26',"BPM16274","LTT9491",
              'GD50','GD108','AGK+81_266','BD+25_4655','BD+33_2642'] # directory names in datadir 
   
   calfildir = '/calibration/grism/fluxcal2/uvotified/'
   calfiles = [ 'p041c_stisnic_003_uvotified.ascii',
                'p330e_stisnic_003_uvotified.ascii',
                'p177d_stisnic_003_uvotified.ascii',
                'gd153_stisnic_003_uvotified.ascii', 
                'wd0320_stisnic_003_uvotified.ascii' ,
                'wd1057_719_stisnic_003_uvotified.ascii',
                'wd1657_343_stisnic_003_uvotified.ascii',
                'wd1121_stis_preview_uvotified.ascii',
                'g63-28_stis_ngsl_v2_uvotified.ascii',
                '../calspectra/bpm16274_f_eso.dat',  
                '../calspectra/ltt9491_002.ascii',  
                '../calspectra/gd50_004.ascii',  
                '../calspectra/gd108_005.ascii',  
                '../calspectra/agk_81d266_stisnic_003.ascii',
                '../calspectra/bd_25d4655_002.ascii',
                '../calspectra/bd_33d2642_fos_003.ascii',
                ]
   calfilID = ['p041c','p330e','p177d','gd153','wd0320-539','wd1057+719',
               'wd1657+343','wd1121+145','g63-26',"bpm16274","ltt9491",
               'gd50','gd108','agk+81_266','bd+25_4655','bd+33_2642']  # names
   caldict = dict([calfilID[i],calfiles[i]] for i in range(len(calfiles)))
   obsdict = dict([calfilID[i],obsdirs[i]] for i in range(len(calfiles)))
   
   if chatter > 2: 
      print('caldict = ', caldict)
      print('obsdict = ', obsdict)
   
   if phafiledir == None   : phafiledir = datadir+'eff_area_'+str(wheelpos)+'/all_2.5/'
   if eff_area_outdir == None : eff_area_out = datadir+'eff_area_'+str(wheelpos)+'/ea_2.5/'
   
   tt = time.struct_time(time.gmtime())
   date = "%4i%2s%2sT%2s%2s"%(tt.tm_year,str(tt.tm_mon).zfill(2),str(tt.tm_mday).zfill(2),str(tt.tm_hour).zfill(2),str(tt.tm_min).zfill(2))
   logfile = '/calibration/grism/fluxcal2/obs_py/log_obs_'+str(wheelpos)+'.'+date+'.py'
   log = open(logfile,'w') # clobber 
   log.write('__version__=%s'%(date))
   log.write('def obs_%s:'%(str(wheelpos)))
   log.write('   obs_%s=dict(\n'%(str(wheelpos)))
   obs_160 = ''
   obs_200 = ''
      
   # data base of the uv clocked grism flux calibration observations
# the wlshift values correspond to the original pha files [before 2013]   
   obs_160 = dict(
#removed 1680174 since the spectrum is contaminated with the first order of a bright background star   
set01=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203011',wheelpos=160, ext=1, anchor=[694.,1369.],dateobs='2005-10-07',exposure=650.,roll=145.1,groups=['b'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2700,7000]],notes="smooth, but very low effective area .. exposure problem or too near the edge ?"),     
set02=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203011',wheelpos=160, ext=2, anchor=[743.,1368.],dateobs='2005-10-07',exposure=715.,roll=145.1,groups=['b'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2700,7000]],notes='dip at -175pix'),
set03=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904002',wheelpos=160, ext=1, anchor=[868.,1560.],dateobs='2011-11-19',exposure=1215,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set04=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904002',wheelpos=160, ext=2, anchor=[873.,1585.],dateobs='2011-11-19',exposure=1215,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set05=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904004',wheelpos=160, ext=1, anchor=[867.,1542.],dateobs='2011-11-20',exposure=1207,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set06=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904004',wheelpos=160, ext=2, anchor=[870.,1546.],dateobs='2011-11-20',exposure=1207,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set07=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904004',wheelpos=160, ext=3, anchor=[871.,1569.],dateobs='2011-11-21',exposure=1207,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set08=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904004',wheelpos=160, ext=4, anchor=[876.,1571.],dateobs='2011-11-21',exposure=1207,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set09=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055904004',wheelpos=160, ext=5, anchor=[888.,1566.],dateobs='2011-11-21',exposure=1207,roll=207.0,groups=['b'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1716],[2950,7000]],notes=''),
set10=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=2, anchor=[820.,1497.],dateobs='2011-11-22',exposure=1215,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=-3,exclude=[[1600,3200],[5200,7000]],notes=''),
set11=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=1, anchor=[917.,1573.],dateobs='2011-11-22',exposure=1215,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=10,exclude=[[1600,3200],[5200,7000]],notes=''),
set12=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=3, anchor=[903.,1574.],dateobs='2011-11-22',exposure=1215,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=0,exclude=[[1600,3200],[5200,7000]],notes=''),
set13=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=4, anchor=[903.,1563.],dateobs='2011-11-22',exposure=1215,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=15,exclude=[[1600,3200],[5200,7000]],notes=''),
set14=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=5, anchor=[908.,1552.],dateobs='2011-11-22',exposure=1215,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=5,exclude=[[1600,3200],[5200,7000]],notes=''),
set15=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763004',wheelpos=160, ext=1, anchor=[912.,1533.],dateobs='2011-11-29',exposure=1207,roll= None,groups=['b'],name='P177D',nproc=0,wlshift=-13,exclude=[[1600,3200],[5200,7000]],notes=''),
set16=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203007',wheelpos=160, ext=1, anchor=[832.,1311.],dateobs='2005-09-30',exposure=965.,roll=152.3,groups=['c'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]] ,notes='second order starts ar +27pix=2700A; deep dip'),
set17=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203007',wheelpos=160, ext=3, anchor=[833.,1306.],dateobs='2005-09-30',exposure=1716,roll= None,groups=['c'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]] ,notes='one deep dip'),
set18=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203008',wheelpos=160, ext=1, anchor=[830.,1297.],dateobs='2005-10-04',exposure=1944,roll=146.3,groups=['c'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2500,7000]] ,notes='dips, contamination > 70pix.(2500A)'),
set19=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203007',wheelpos=160, ext=2, anchor=[850.,1327.],dateobs='2005-09-30',exposure=1774,roll= None,groups=['c'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]] ,notes='3 deep dips'),
set20=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055203007',wheelpos=160, ext=4, anchor=[879.,1330.],dateobs='2005-09-30',exposure=1775,roll= None,groups=['c'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]] ,notes='3 deep dips'),
set21=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056763002',wheelpos=160, ext=6, anchor=[989.,1502.],dateobs='2011-11-22',exposure=949.,roll= None,groups=['c'],name='P177D',nproc=0,wlshift=-21,exclude=[[1600,3100],[5200,7000]] ,notes=''),
set22=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055903002',wheelpos=160, ext=1, anchor=[1637.,1596.],dateobs='2011-11-17',exposure=250.,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1920],[2960,3100],[3350,7000]],notes=''),
set23=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055903002',wheelpos=160, ext=2, anchor=[1622.,1573.],dateobs='2011-11-17',exposure= 67.,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1920],[2960,3100],[3350,7000]],notes=''),
set24=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903002',wheelpos=160, ext=3, anchor=[1621.,1599.],dateobs='2011-11-17',exposure=1018,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1820],[2960,3100],[3350,7000]],notes=''),
set25=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903004',wheelpos=160, ext=1, anchor=[1600.,1582.],dateobs='2011-11-18',exposure=1215,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1920],[2960,3100],[3350,7000]],notes=''),
set26=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903004',wheelpos=160, ext=2, anchor=[1618.,1596.],dateobs='2011-11-18',exposure=1215,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1920],[2960,3100],[3350,7000]],notes=''),
set27=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903004',wheelpos=160, ext=3, anchor=[1617.,1581.],dateobs='2011-11-19',exposure=1215,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1920],[2960,3100],[3350,7000]],notes=''),
set28=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903004',wheelpos=160, ext=4, anchor=[1629.,1590.],dateobs='2011-11-19',exposure=1215,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1820],[2960,3100],[3350,7000]],notes=''),
set29=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903006',wheelpos=160, ext=2, anchor=[1664.,1608.],dateobs='2011-11-23',exposure=1453,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1820],[2960,3100],[3350,7000]],notes=''),
set30=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903006',wheelpos=160, ext=3, anchor=[1666.,1555.],dateobs='2011-11-23',exposure=518.,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1820],[2960,3100],[3350,7000]],notes=''),
set31=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055901002',wheelpos=160, ext=1, anchor=[1829.,1717.],dateobs='2011-08-13',exposure=1264,roll= None,groups=['d'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,1820],[2730,7000]],notes=''),
set32=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=5, anchor=[1620.,1574.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=0,exclude=[[1600,3200]], notes='(investigate second order  > 4900)' ), 
set33=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=7, anchor=[1606.,1573.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-33,exclude=[[1600,3200],[5800,7000]], notes='(investigate second order  > 4900)' ),
set34=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00056762002',wheelpos=160, ext=1, anchor=[1656.,1589.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-35,exclude=[[1600,3200],[5800,7000]], notes='(investigate second order  > 4900)' ),
set35=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=2, anchor=[1635.,1554.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-33,exclude=[[1600,3200],[5800,7000]], notes='(investigate second order  > 4900)' ),
set36=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=3, anchor=[1645.,1562.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-30,exclude=[[1600,3200],[5800,7000]], notes='(investigate second order  > 4900)' ),
set37=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=6, anchor=[1651.,1559.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-40,exclude=[[1600,3200],[5800,7000]], notes='(investigate second order  > 4900)' ),
set38=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056761002',wheelpos=160, ext=1, anchor=[1856.,1703.],dateobs='2011-08-12',exposure=1264.5,roll= None,groups=['d'],name='P177D',nproc=0,wlshift=-50,exclude=[[1600,3208],[6000,7000]], notes='(investigate second order  > 4900)'),
#set39=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00056771002',wheelpos=160, ext=1, anchor=[1829.,1703.],dateobs='2011-08-12',exposure=1264,roll= None,groups=['d'],name='P330E',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]], notes='(*not included) - first order overlap- make another check)'),
set40=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057955002',wheelpos=160, ext=1, anchor=[1814.,1736.],dateobs='2011-08-11',exposure=1264,roll= None,groups=['d'],name='P041C',nproc=0,wlshift=-60,exclude=[[1600,3200],[5800,7000]], notes='verified - no apparent coi effect, but very noisy'),
set41=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055903006',wheelpos=160, ext=1, anchor=[1794.,1472.],dateobs='2011-11-23',exposure=1055,roll= None,groups=['h'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[3100,7000]],notes=''),
set42=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055900050',wheelpos=160, ext=1, anchor=[1334.,943.],dateobs='2008-04-22',exposure=772.,roll=054.5,groups=['g'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,2010],[2220,2365],[2670,7000]],notes=''),
set43=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055900050',wheelpos=160, ext=3, anchor=[1289.,891.],dateobs='2008-04-22',exposure=772.,roll=054.5,groups=['g'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,2010],[2220,2365],[2670,7000]],notes=''),
set44=dict(auto=True,fit2nd=False,offsetlimit=[100,3],phatype="_g",obsid='00056760002',wheelpos=160, ext=1, anchor=[1225.,899.],dateobs='2006-03-20',exposure=578.4,roll=059.1,groups=['g'],name='P177D',nproc=0,wlshift=0,exclude=[[1600,3600],[5200,7000]] ,notes='deadc prob.'),
set45=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055202018',wheelpos=160, ext=1, anchor=[1479.,0746.],dateobs='2005-10-04',exposure=1944,roll= None,groups=['g'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2745,7000]],notes=''),
#set46=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057952001',wheelpos=160, ext=1, anchor=[1216.,0831.],dateobs='2009-03-23',exposure=723.,roll= None,groups=['g'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,7000]],notes='overlap strong source'),
set47=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055201011',wheelpos=160, ext=1, anchor=[1444.,1329.],dateobs='2005-09-29',exposure=922.,roll=152.4,groups=['f'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]],notes='second order merge starts at 400pix,4000A[1],2nd order 70% weaker than 1600,1600'),
set48=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055201011',wheelpos=160, ext=2, anchor=[1447.,1308.],dateobs='2005-09-29',exposure=1631,roll=152.4,groups=['f'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]],notes=''),
set49=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055201011',wheelpos=160, ext=3, anchor=[1494.,1257.],dateobs='2005-09-29',exposure=923.,roll=152.4,groups=['f'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3100,7000]],notes=''),
set50=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200014',wheelpos=160, ext=1, anchor=[1130.,1003.],dateobs='2005-09-29',exposure=1395,roll= 151.7,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set51=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200014',wheelpos=160, ext=2, anchor=[1120.,1006.],dateobs='2005-09-29',exposure=1454,roll= None,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set52=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200014',wheelpos=160, ext=3, anchor=[1134.,1007.],dateobs='2005-09-29',exposure=805.,roll= None,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set53=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200014',wheelpos=160, ext=4, anchor=[1127.,1005.],dateobs='2005-09-29',exposure=1217,roll= None,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set54=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200012',wheelpos=160, ext=1, anchor=[1176.,1073.],dateobs='2005-05-25',exposure=1277,roll=288.3,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set55=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00056800006',wheelpos=160, ext=1, anchor=[1109.,1025.],dateobs='2006-06-17',exposure=1572,roll= None,groups=['aA'],name='P041C',nproc=0,wlshift=3,exclude=[[1600,3000],[5200,7000]],notes=''),
set56=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00056800006',wheelpos=160, ext=2, anchor=[1122.,1021.],dateobs='2006-06-17',exposure=1750,roll= None,groups=['aA'],name='P041C',nproc=0,wlshift=20,exclude=[[1600,3000],[5200,7000]],notes=''),
#set57=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00056800010',wheelpos=160, ext=1, anchor=[1111.,1010.],dateobs='2006-06-27',exposure=1749,roll= None,groups=['a'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3350],[5200,7000]],notes='overlap strong source'),
#set58=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00056800010',wheelpos=160, ext=2, anchor=[1115.,1007.],dateobs='2006-06-27',exposure=1750,roll= None,groups=['a'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3350],[5200,7000]],notes='overlap strong source'),
set59=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055900050',wheelpos=160, ext=2, anchor=[1300.,1078.],dateobs='2008-04-22',exposure=772.,roll=054.4,groups=['i'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[1600,2010],[2040,2150],[2220,2365],[2670,7000]],notes=''),
#set60=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057954001',wheelpos=160, ext=1, anchor=[1327.,1107.],dateobs='2009-03-23',exposure=723.,roll= None,groups=['i'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set61=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055201007',wheelpos=160, ext=1, anchor=[1157.,622.],dateobs='2005-05-25',exposure=1276,roll= None,groups=['j'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2700,7000]],notes=''),
set62=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055204009',wheelpos=160, ext=1, anchor=[788.,759.],dateobs='2005-10-03',exposure=812.,roll=147.0,groups=['e'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2600,7000]],notes=''),
set63=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055204011',wheelpos=160, ext=1, anchor=[755.,757.],dateobs='2005-10-07',exposure=692.,roll=143.0,groups=['e'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2600,7000]],notes='smooth; 2nd/3rd below 1st, overlap all'),
set64=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055204009',wheelpos=160, ext=2, anchor=[863.,676.],dateobs='2005-10-03',exposure=1826,roll= None,groups=['e'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2600,7000]],notes=''),
set65=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055204009',wheelpos=160, ext=3, anchor=[860.,674.],dateobs='2005-10-03',exposure=1826,roll= None,groups=['e'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2600,7000]],notes='slight dips in sum'),
set66=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200015',wheelpos=160, ext=1, anchor=[1098.,1025.],dateobs='2006-03-11',exposure=1734,roll=000.9,groups=['k'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='deadc problem'),
set67=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055200015',wheelpos=160, ext=2, anchor=[1098.,1022.],dateobs='2006-03-11',exposure=1734,roll=000.9,groups=['k'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='deadc problem'),
set68=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055200022',wheelpos=160, ext=1, anchor=[1078.,1025.],dateobs='2010-07-25',exposure=1298,roll=220.0,groups=['k'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='spectrum fainter than earlier spectrum; using this for eff. area.'),
set69=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056762002',wheelpos=160, ext=4, anchor=[1582.,1762.],dateobs='2011-11-20',exposure=1207,roll= None,groups=['k'],name='P177D',nproc=0,wlshift=-45,exclude=[[1600,3100],[5200,7000]],notes=''),
#set70=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00056770015',wheelpos=160, ext=1, anchor=[1211.,855.],dateobs='2006-03-23',exposure=435.,roll=067.6,groups=[''],name='P330E',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes='DOMINATED by neighbor'),
#set71=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057951001',wheelpos=160, ext=1, anchor=[1015.,1192.],dateobs='2009-03-23',exposure=723.,roll= None,groups=['k'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes=''),
#set72=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057953001',wheelpos=160, ext=1, anchor=[0903.,0917.],dateobs='2009-03-23',exposure=723.,roll= None,groups=['l'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes=''),
set73=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055200028',wheelpos=160, ext=1, anchor=[1086.,1073.],dateobs='2012-06-01',exposure=1158.,roll= 292.0,groups=['a'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='exclude from 2700 for trackwidth=2.5'),
set74=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055200028',wheelpos=160, ext=2, anchor=[1100.,1088.],dateobs='2012-06-01',exposure=814.,roll= 292.0,groups=['aA'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='exclude from 2700 for trackwidth=2.5'),
set75=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055200028',wheelpos=160, ext=3, anchor=[1116.,1060.],dateobs='2012-06-01',exposure=863.,roll= 292.0,groups=['aA'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2720,7000]],notes='exclude from 2700 for trackwidth=2.5'),
set76=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250021',wheelpos=160, ext=1, anchor=[1127.,1037.],dateobs='2012-06-05',exposure=1109.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set77=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250021',wheelpos=160, ext=2, anchor=[1129.,1026.],dateobs='2012-06-05',exposure=866.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set78=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250021',wheelpos=160, ext=3, anchor=[1130.,1039.],dateobs='2012-06-05',exposure=971.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set79=dict(auto=True,fit2nd=False,offsetlimit=[100,0],phatype="_f",obsid='00055900073',wheelpos=160, ext=1, anchor=[1144.,1003.],dateobs='2012-08-06',exposure=1207.,roll= 294.0,groups=['aA'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[2136,2333],[2660,7000]],notes=''),
set80=dict(auto=True,fit2nd=False,offsetlimit=[99,0],phatype="_f",obsid='00055900073',wheelpos=160, ext=2, anchor=[1166., 982.],dateobs='2012-08-06',exposure=1207.,roll= 294.0,groups=['aA'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[2130,2330],[2660,7000]],notes=''),
set81=dict(auto=True,fit2nd=False,offsetlimit=[99,0],phatype="_f",obsid='00055900073',wheelpos=160, ext=3, anchor=[1145.,1003.],dateobs='2012-08-06',exposure=1207.,roll= 294.0,groups=['aA'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[2135,2345],[2660,7000]],notes=''),
set82=dict(auto=True,fit2nd=False,offsetlimit=[98,0],phatype="_f",obsid='00055900073',wheelpos=160, ext=4, anchor=[1153., 999.],dateobs='2012-08-06',exposure=1158.,roll= 294.0,groups=['aA'],name='WD1657+343',nproc=0,wlshift=0,exclude=[[2135,2315],[2660,7000]],notes=''),
set83=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056760024',wheelpos=160, ext=1, anchor=[1133.,1025.],dateobs='2012-08-10',exposure=1010.6,roll= 280.0,groups=['aA'],name='P177D',nproc=0,wlshift=0,exclude=[[1600,3500]], notes='' ),
set84=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056760024',wheelpos=160, ext=2, anchor=[1143.,1028.],dateobs='2012-08-10',exposure=1010.6,roll= 280.0,groups=['aA'],name='P177D',nproc=0,wlshift=-9.0,exclude=[[1600,3500]], notes='' ),
set85=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00056760024',wheelpos=160, ext=3, anchor=[1138.,1028.],dateobs='2012-08-10',exposure=1007.2,roll= 280.0,groups=['aA'],name='P177D',nproc=0,wlshift=-10,exclude=[[1600,3500]], notes='' ),
set86=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250023',wheelpos=160, ext=1, anchor=[1131.,1039.],dateobs='2012-06-06',exposure=1109.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set87=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250023',wheelpos=160, ext=2, anchor=[1134.,1033.],dateobs='2012-06-06',exposure= 863.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2070,2177],[2720,7000]],notes='dip - data dropout'),
set88=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250023',wheelpos=160, ext=3, anchor=[1135.,1038.],dateobs='2012-06-06',exposure= 961.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set89=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00054250023',wheelpos=160, ext=4, anchor=[1130.,1027.],dateobs='2012-06-06',exposure=1060.,roll= 16.0,groups=['a'],name='WD0320-539',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set120=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600020',wheelpos=160, ext=1, anchor=[1206.,1018.],dateobs='2005-07-15',exposure=148.4,roll=291.3,groups=['Aa'],name='G63-26',nproc=0,wlshift=13,exclude=[[1600,2200],[2340,2640],[2850,3065],[3570,7000]],notes=''),
set121=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600020',wheelpos=160, ext=2, anchor=[1218.,1021.],dateobs='2005-07-15',exposure=123.2,roll=291.3,groups=['Aa'],name='G63-26',nproc=0,wlshift=16,exclude=[[1600,2200],[2850,3065],[3211,3360],[3570,7000]],notes=''),
set122=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055601010',wheelpos=160, ext=1, anchor=[1179., 579.],dateobs='2005-07-15',exposure=156.3,roll=291.4,groups=['C'],name='G63-26',nproc=0,wlshift=4,exclude=[[2911,3065],[1600,2720],[3570,7000]],notes=''),
set123=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055601010',wheelpos=160, ext=2, anchor=[1164., 561.],dateobs='2005-07-15',exposure=130.6,roll=291.4,groups=['C'],name='G63-26',nproc=0,wlshift=4,exclude=[[2920,3065],[1600,2720],[3570,7000]],notes=''),
set124=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055601010',wheelpos=160, ext=3, anchor=[1158., 559.],dateobs='2005-07-15',exposure=166.3,roll=291.4,groups=['C'],name='G63-26',nproc=0,wlshift=4,exclude=[[2860,3065],[1600,2720],[3570,7000]],notes=''),
set125=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055500014',wheelpos=160, ext=1, anchor=[1079.4,1005.7],dateobs='2005-08-15',exposure=663.,roll=268.8,groups=[''],name='GD153',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set126=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055500014',wheelpos=160, ext=2, anchor=[1084.4,1006.1],dateobs='2005-08-15',exposure=628.,roll=268.8,groups=[''],name='GD153',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set127=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055505002',wheelpos=160, ext=1, anchor=[1783.4,1671.3],dateobs='2011-08-11',exposure=1264.,roll=272.0,groups=['d'],name='GD153',nproc=0,wlshift=0,exclude=[[2720,7000]],notes=''),
set128=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057958001',wheelpos=160, ext=1, anchor=[1505.5,1470.6],dateobs='2012-10-09',exposure=1384.6,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-51,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set129=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057958001',wheelpos=160, ext=2, anchor=[1454.5,1495.6],dateobs='2012-10-09',exposure=1384.6,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-50,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set130=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057959001',wheelpos=160, ext=1, anchor=[1264.5,1303.0],dateobs='2012-10-09',exposure=1384.6,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-33,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set131=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057959001',wheelpos=160, ext=2, anchor=[1261.7,1304.0],dateobs='2012-10-09',exposure=1384.6,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-30,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set132=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057960002',wheelpos=160, ext=1, anchor=[1358.1,1059.8],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-33,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set133=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057960002',wheelpos=160, ext=2, anchor=[1362.0,1050.0],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-29,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set134=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057960002',wheelpos=160, ext=3, anchor=[1367.6,1062.8],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-40,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set135=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057963002',wheelpos=160, ext=1, anchor=[996.3,1425.1],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-20,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set136=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057963002',wheelpos=160, ext=4, anchor=[986.,1432.30],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-20,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set137=dict(auto=True,fit2nd=False,offsetlimit=[95,1.5],phatype="_f",obsid='00057964002',wheelpos=160, ext=1, anchor=[728.9,1265.2],dateobs='2012-10-06',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=4,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set138=dict(auto=True,fit2nd=False,offsetlimit=[101,0],phatype="_g",obsid='00057964002',wheelpos=160, ext=2, anchor=[641.2,1278.4],dateobs='2012-10-06',exposure=218.0,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes='fainter overlap'),
set139=dict(auto=True,fit2nd=False,offsetlimit=[101,0],phatype="_f",obsid='00057964002',wheelpos=160, ext=3, anchor=[729.8,1256.3],dateobs='2012-10-06',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3100],[5200,7000]],notes='fainter overlap'),
set140=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00057965002',wheelpos=160, ext=1, anchor=[780.9,1463.8],dateobs='2012-10-09',exposure=898.7,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=7,exclude=[[1600,3100],[5200,7000]],notes='faint, overlap'),
set141=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00057965002',wheelpos=160, ext=2, anchor=[780.8,1465.0],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=3,exclude=[[1600,3100],[5200,7000]],notes='faint, overlap'),
set142=dict(auto=True,fit2nd=False,offsetlimit=[99,0],phatype="_f",obsid='00057966002',wheelpos=160, ext=1, anchor=[596.1,1292.7],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=24,exclude=[[1600,3100],[5200,7000]],notes='faint, overlap'),
set143=dict(auto=True,fit2nd=False,offsetlimit=[99,0],phatype="_f",obsid='00057966002',wheelpos=160, ext=2, anchor=[596.3,1285.6],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=25,exclude=[[1600,3100],[5200,7000]],notes='faint,overlap'),
set144=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00057967002',wheelpos=160, ext=1, anchor=[709.4,652.6],dateobs='2012-10-09',exposure=1277.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-50,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set145=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00057967002',wheelpos=160, ext=2, anchor=[710.8,653.4],dateobs='2012-10-09',exposure=1277.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-51,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set146=dict(auto=True,fit2nd=False,offsetlimit=[98,2],phatype="_f",obsid='00057967002',wheelpos=160, ext=3, anchor=[708.5,640.6],dateobs='2012-10-09',exposure=901.3,roll=218.0,groups=[''],name='P041C',nproc=0,wlshift=-50,exclude=[[1600,3100],[5200,7000]],notes='overlap'),
set147=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057968002',wheelpos=160, ext=1, anchor=[1547.,1549.],dateobs='2012-10-23',exposure=740.9,roll=198.0,groups=['d'],name='P041C',nproc=0,wlshift=-57,exclude=[[1600,3100],[5200,7000]],notes=''),
set148=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057968002',wheelpos=160, ext=2, anchor=[1571.9,1565.6],dateobs='2012-10-23',exposure=1062.7,roll=198.0,groups=['d'],name='P041C',nproc=0,wlshift=-49,exclude=[[1600,3100],[5200,7000]],notes=''),
set149=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057968002',wheelpos=160, ext=3, anchor=[1576.0,1554.3],dateobs='2012-10-23',exposure=901.3,roll=198.0,groups=['d'],name='P041C',nproc=0,wlshift=-47,exclude=[[1600,3100],[5200,7000]],notes=''),
set150=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057969001',wheelpos=160, ext=1, anchor=[1541.2,1301.1],dateobs='2012-10-24',exposure=1008.6,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-45,exclude=[[1600,3100],[5200,7000]],notes=''),
set151=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057969001',wheelpos=160, ext=2, anchor=[1536.1,1301.7],dateobs='2012-10-24',exposure=1008.6,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-46,exclude=[[1600,3100],[5200,7000]],notes=''),
set152=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00057969001',wheelpos=160, ext=3, anchor=[1575.4,1155.5],dateobs='2012-10-24',exposure=413.5,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-80,exclude=[[1600,3100],[5200,7000]],notes=''),
set153=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057970002',wheelpos=160, ext=1, anchor=[1572.8,737.4],dateobs='2012-10-25',exposure=1107.,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-28,exclude=[[1600,3100],[5200,7000]],notes=''),
set154=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057970002',wheelpos=160, ext=2, anchor=[1572.4,736.9],dateobs='2012-10-25',exposure=1107.,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-39,exclude=[[1600,3100],[5200,7000]],notes=''),
set155=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057970002',wheelpos=160, ext=3, anchor=[1589.3,743.2],dateobs='2012-10-25',exposure=892.5,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-35,exclude=[[1600,3100],[5200,7000]],notes=''),
set156=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057971002',wheelpos=160, ext=1, anchor=[1154.9,618.8],dateobs='2012-10-26',exposure=1107.0,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-33,exclude=[[1600,3100],[5200,7000]],notes=''),
set157=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057971002',wheelpos=160, ext=2, anchor=[1113.2,604.1],dateobs='2012-10-26',exposure=839.3,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=-25,exclude=[[1600,3100],[5200,7000]],notes=''),
set158=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057972002',wheelpos=160, ext=1, anchor=[940.6,1560.5],dateobs='2012-10-27',exposure=1161.2,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3200],[5200,7000]],notes='bit fainter'),
set159=dict(auto=True,fit2nd=False,offsetlimit=[98,0],phatype="_f",obsid='00057972002',wheelpos=160, ext=2, anchor=[898.2,1529.7],dateobs='2012-10-27',exposure=1161.2,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=0,exclude=[[1600,3200],[5200,7000]],notes='bit fainter'),
set160=dict(auto=True,fit2nd=False,offsetlimit=[97,0],phatype="_f",obsid='00057972002',wheelpos=160, ext=3, anchor=[879.7,1524.5],dateobs='2012-10-27',exposure=677.9,roll=198.0,groups=[''],name='P041C',nproc=0,wlshift=9,exclude=[[1600,3200],[5200,7000]],notes='faint'),
set161=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055205001',wheelpos=160, ext=1, anchor=[1568.,1603.],dateobs='2012-11-04',exposure=1116.9,roll=123.0,groups=['d'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[3080,3250],[3500,7000]],notes=''),
set162=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055205001',wheelpos=160, ext=2, anchor=[1563.,1585.],dateobs='2012-11-04',exposure=663.6,roll=123.0,groups=['d'],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2976,3252],[3500,7000]],notes=''),
set163=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=1, anchor=[1216.,1310.],dateobs='2012-11-07',exposure=472.2,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1775,1905],[3000,7000] ],notes='DIP'),
set164=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=2, anchor=[1214.,1315.],dateobs='2012-11-07',exposure=686.8,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1771,1898],[2550,2750],[3200,7000]],notes='DIP'),
set165=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=3, anchor=[1228.,1295.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1800,2010],[2685,2870],[3100,7000]],notes='DIP'),
set166=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=4, anchor=[1228.,1294.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1850,2010],[2685,7000]],notes='DIPs'),
set167=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055207002',wheelpos=160, ext=1, anchor=[1546.,1299.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2954,3245],[3400,7000]],notes=''),
#set168=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055606002',wheelpos=160, ext=1, anchor=[1575.1,1564.5],dateobs='2012-11-19',exposure=1384.6,roll=138.0,groups=['d'],name='G63-26',nproc=0,wlshift=-37,exclude=[[1600,2600],[3700,7000]],notes='1st order cont.'),
#set169=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055607002',wheelpos=160, ext=1, anchor=[1221.5,1287.6],dateobs='2012-11-20',exposure=1277.3,roll=138.0,groups=[''],name='G63-26',nproc=0,wlshift=-22,exclude=[[1600,2600],[3700,7000]],notes='1st order cont.'),
#set170=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055608002',wheelpos=160, ext=1, anchor=[1559.2,1300.2],dateobs='2012-11-21',exposure=1224.2,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-41,exclude=[[1600,2740],[3700,7000]],notes='1st order cont.'),
#set171=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055609002',wheelpos=160, ext=1, anchor=[1546.8,727.9],dateobs='2012-11-22',exposure=1384.6,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-30,exclude=[[1600,2730],[3700,7000]],notes='1st order cont.'),
#set172=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055610002',wheelpos=160, ext=1, anchor=[962.7,1377.9],dateobs='2012-11-23',exposure=969.,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-70,exclude=[[1600,2475],[3700,7000]],notes='no lenticular'),
#set173=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055611002',wheelpos=160, ext=1, anchor=[742.1,1258.0],dateobs='2012-11-24',exposure=1322.6,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-22,exclude=[[1600,2350],[3700,7000]],notes='1st order cont. '),
#set174=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055612004',wheelpos=160, ext=1, anchor=[864.5,1510.4],dateobs='2012-11-28',exposure=337.3,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-62,exclude=[[1600,2900],[3700,7000]],notes='1st order cont.no lenticular; too short; contaminated?'),
#set175=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055215001',wheelpos=160, ext=1, anchor=[1626.4,678.4],dateobs='2012-12-20',exposure=50.8,roll=87.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2700,7000]],notes='nearby bright spectrum contaminates 2.5sig slit'),
#set176=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055215001',wheelpos=160, ext=2, anchor=[1547.4,707.2],dateobs='2012-12-21',exposure=579.5,roll=87.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2800,7000]],notes='nearby bright spectrum contaminates 2.5sig slit'),
#set177=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055215001',wheelpos=160, ext=3, anchor=[1545.3,709.7],dateobs='2012-12-21',exposure=579.5,roll=87.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2800,7000]],notes='nearby bright spectrum contaminates 2.5sig slit'),
#set178=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055216001',wheelpos=160, ext=1, anchor=[1181.4,644.0],dateobs='2012-12-22',exposure=1384.6,roll=87.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2330,2500],[2730,7000]],notes='next to bright cont spectrum'),
#set179=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055216001',wheelpos=160, ext=2, anchor=[1140.3,567.8],dateobs='2012-12-22',exposure=1008.6,roll=87.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[2350,2490],[2760,7000]],notes='next to bright cont. spectrum COI problem!'),
set180=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057530004',wheelpos=160, ext=1, anchor=[1111.0,1049.6],dateobs='2014-07-04',exposure=892.5,roll=228.0,groups=[''],name='AGK+81_266',nproc=0,wlshift=0,exclude=[[2800,7000]],notes='coi'),
set181=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057538004',wheelpos=160, ext=1, anchor=[1115.0,1067.0],dateobs='2014-07-05',exposure=892.5,roll=45.0,groups=[''],name='BD+25_4655',nproc=0,wlshift=0,exclude=[[1790,7000]],notes='coi'),
set182=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057530006',wheelpos=160, ext=1, anchor=[1153.8,1042.0],dateobs='2014-07-05',exposure=892.5,roll=226.0,groups=[''],name='AGK+81_266',nproc=0,wlshift=0,exclude=[[2800,7000]],notes='coi'),
#set18x=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=3, anchor=[1228.,1295.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1800,2010],[2685,2870],[3100,7000]],notes='DIP'),
#set18x=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=3, anchor=[1228.,1295.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1800,2010],[2685,2870],[3100,7000]],notes='DIP'),
#set18x=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=3, anchor=[1228.,1295.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1800,2010],[2685,2870],[3100,7000]],notes='DIP'),
#set18x=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055206001',wheelpos=160, ext=3, anchor=[1228.,1295.],dateobs='2012-11-07',exposure=633.6,roll=123.0,groups=[''],name='WD1057+719',nproc=0,wlshift=0,exclude=[[1800,2010],[2685,2870],[3100,7000]],notes='DIP'),
#setxxx=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00055614002',wheelpos=160, ext=1, anchor=[1158., 559.],dateobs='2012-11-27',exposure=1490.9,roll=138.,groups=[''],name='G63-26',nproc=0,wlshift=-24,exclude=[[1600,2200],[3700,7000]],notes=''),
)   

   obs_200 = dict(
#__version__= '20120820T1619'
# removed set86 which seems extra low and noisy
# 2012-08-10 removed& put back set92 and set 141. +set104, set93, set66 had DEADC wrong - reprocessed from Level 0 data
# 2012-08-20 removed again set66 : flux way too high CNTEXP suggest exposure ~ 182.9 instead of 893.6 but should be ~ 160s. 
#  
set01={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 131.75, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 1, 'wlshift': 0, 'groups': [''], 'nproc': 1, 'obsid': '00055500016', 'exclude': [[2730, 7000]], 'anchor': [1024.0, 1085.0], 'exposure': 1033.0,},
set02={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 131.75, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 2, 'wlshift': 0, 'groups': [''], 'nproc': 1, 'obsid': '00055500016', 'exclude': [[2730, 7000]], 'anchor': [1028.0, 1080.0], 'exposure': 1057.0,},
set03={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': 'weak blue star < 1736 ', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': [''], 'nproc': 1, 'obsid': '00055500010', 'exclude': [[2730, 7000]], 'anchor': [1066.0, 964.0], 'exposure': 195.0,},
set04={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250001', 'exclude': [[2850, 7000]], 'anchor': [992.0, 1086.0], 'exposure': 154.0,},
set05={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250001', 'exclude': [[2850, 7000]], 'anchor': [993.0, 1089.0], 'exposure': 104.0,},
set06={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 3, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250001', 'exclude': [[2850, 7000], [3000.0, 3200.0]], 'anchor': [996.0, 1083.0], 'exposure': 865.0,},
set07={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 4, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250001', 'exclude': [[2850, 7000]], 'anchor': [1000.0, 1039.0], 'exposure': 491.0,},
set08={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 5, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250001', 'exclude': [[2850, 7000], [3000.0, 3200.0]], 'anchor': [1021.0, 1066.0], 'exposure': 807.0,},
set09={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-23', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250019', 'exclude': [[2850, 7000]], 'anchor': [1037.0, 1093.0], 'exposure': 1217.0,},
set10={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-23', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250019', 'exclude': [[2850, 7000]], 'anchor': [1027.0, 1113.0], 'exposure': 1218.0,},
set11={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-23', 'ext': 3, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00054250019', 'exclude': [[2850, 7000]], 'anchor': [1021.0, 1095.0], 'exposure': 981.0,},
set12={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200001', 'exclude': [[2850, 7000]], 'anchor': [1001.0, 1117.0], 'exposure': 1456.0,},
set13={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200001', 'exclude': [[2850, 7000]], 'anchor': [999.0, 1082.0], 'exposure': 1042.0,},
set14={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200005', 'exclude': [[2850, 7000]], 'anchor': [1008.0, 1122.0], 'exposure': 1605.0,},
set15={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-30', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200010', 'exclude': [[2820, 7000]], 'anchor': [983.0, 1039.0], 'exposure': 978.0,},
set16={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-30', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200010', 'exclude': [[2704, 7000]], 'anchor': [986.0, 1045.0], 'exposure': 689.0,},
set17={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 359.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-03-12', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200016', 'exclude': [[2850, 7000]], 'anchor': [1001.0, 1078.0], 'exposure': 1734.0,},
set18={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 359.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-03-12', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055200016', 'exclude': [[1950, 2100], [2850, 7000]], 'anchor': [1004.0, 1080.0], 'exposure': 1700.0,},
set19={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 17.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-05-22', 'ext': 1, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055900056', 'exclude': [[2700, 7000]], 'anchor': [1004.0, 1069.0], 'exposure': 289.0,},
set20={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 17.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-05-22', 'ext': 2, 'wlshift': 0, 'groups': ['Aa'], 'nproc': 1, 'obsid': '00055900056', 'exclude': [[3100, 7000]], 'anchor': [1003.0, 1066.0], 'exposure': 289.0,},
set21={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'P041C', 'roll': 305.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-06-23', 'ext': 1, 'wlshift': -6.0, 'groups': ['AaB'], 'nproc': 1, 'obsid': '00056800008', 'exclude': [[4900, 7000], [1600.0, 2934.0], [4000.0, 4300.0]], 'anchor': [999.0, 1067.0], 'exposure': 1218.0,},
set22={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'P041C', 'roll': 305.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['AaB'], 'nproc': 1, 'obsid': '00056800011', 'exclude': [[1600.0, 2930.0],[3400,3500],[3970, 4300],[5100, 7000]], 'anchor': [994.0, 1061.0], 'exposure': 1749.0,},
set23={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'P041C', 'roll': 305.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-06-27', 'ext': 2, 'wlshift': 6, 'groups': ['AaB'], 'nproc': 1, 'obsid': '00056800011', 'exclude': [[4070, 4300], [5030, 7000], [1600.0, 2940.0]], 'anchor': [998.0, 1076.0], 'exposure': 1750.0,},
set24={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 320.0, 'notes': 'offsetlimit', 'wheelpos': 200, 'dateobs': '2008-07-23', 'ext': 2, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055900057', 'exclude': [[1600, 1760], [2660, 2800], [2700, 7000]], 'anchor': [1055.0, 1085.0], 'exposure': 624.0,},
set25={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 320.0, 'notes': 'offsetlimit', 'wheelpos': 200, 'dateobs': '2008-07-23', 'ext': 5, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055900057', 'exclude': [[1600, 1760], [3200, 7000], [2670.0, 2890.0]], 'anchor': [960.0, 1138.0], 'exposure': 624.0,},
set26={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 320.0, 'notes': 'offsetlimit', 'wheelpos': 200, 'dateobs': '2008-07-23', 'ext': 8, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055900057', 'exclude': [[1600, 1775], [2130, 2240], [2650, 2865], [3170, 7000]], 'anchor': [967.0, 1146.0], 'exposure': 624.0,},
set27={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD1657+343', 'roll': 320.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-07-23', 'ext': 11, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055900057', 'exclude': [[1600, 1775], [2650, 2915], [3120, 7000]], 'anchor': [935.0, 1171.0], 'exposure': 624.0,},
set28={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P041C', 'roll': 36.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2009-03-24', 'ext': 1, 'wlshift': 18.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00057950003', 'exclude': [[1600, 2740], [4700, 7000], [3030.0, 3200.0], [3486.0, 3990.0]], 'anchor': [1007.0, 1117.0], 'exposure': 955.0,},
set29={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-01', 'ext': 1, 'wlshift': -9, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056800014', 'exclude': [[3475, 3678], [4454, 7000], [1600.0, 2966.0]], 'anchor': [977.0, 1051.0], 'exposure': 863.0,},
set30={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-01', 'ext': 2, 'wlshift': 0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056800014', 'exclude': [[1600, 2650], [3475, 3700], [4454, 7000]], 'anchor': [975.0, 1125.0], 'exposure': 1060.0,},
set31={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-01', 'ext': 3, 'wlshift': 2, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056800014', 'exclude': [[3475, 3701], [4454, 7000], [1600.0, 2720.0]], 'anchor': [995.0, 1085.0], 'exposure': 1060.0,},
set32={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-01', 'ext': 4, 'wlshift': 5, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056800014', 'exclude': [[1600, 2650], [3475, 3690], [4454, 7000]], 'anchor': [999.0, 1074.0], 'exposure': 1060.0,},
set33={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 1, 'wlshift': -20.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600.0, 2945.0], [3550.0, 3800.0],[4050, 4450], [4900, 7000], ], 'anchor': [949.0, 1085.0], 'exposure': 863.0,},
set34={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 2, 'wlshift': -32.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600, 2970], [3550, 3800], [4050, 4450], [4900, 7000]], 'anchor': [986.0, 1089.0], 'exposure': 568.0,},
set35={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 3, 'wlshift': -10.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600, 2920], [2630, 2845], [3550, 3800], [4050, 4450], [4900, 7000], [4120.0, 4300.0]], 'anchor': [977.0, 1096.0], 'exposure': 1207.0,},
set36={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 4, 'wlshift': -14.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600, 2920], [3550, 3800], [4050, 4450], [4900, 7000]], 'anchor': [976.0, 1093.0], 'exposure': 1207.0,},
set37={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 5, 'wlshift': 0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600.0, 2930.0],[3550, 3800], [4050, 4500], [4900, 7000]], 'anchor': [999.0, 1077.0], 'exposure': 666.0,},
set38={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-10', 'ext': 6, 'wlshift': -19.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600.0, 2920.0], [3550.0, 3800.0],[4050, 4450], [4900, 7000]], 'anchor': [981.0, 1079.0], 'exposure': 863.0,},
set39={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD0320-539', 'roll': 359.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 1, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00054252003', 'exclude': [[3350, 7000], [1600.0, 1840.0]], 'anchor': [858.0, 1474.0], 'exposure': 593.0,},
set40={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD0320-539', 'roll': 359.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 2, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00054252003', 'exclude': [[1600.0, 1837.0],[3400, 7000]], 'anchor': [847.0, 1492.0], 'exposure': 629.0,},
set41={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-11', 'ext': 7, 'wlshift': -7.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [ [1600.0, 2935.0],[3550, 3800], [4050, 4450], [4900, 7000]],  'anchor': [1041.0, 1068.0], 'exposure': 912.0,},
set42={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-11', 'ext': 8, 'wlshift': -22, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[1600.0, 2935.0],[3550, 3790], [4050, 4450], [4900, 7000]],  'anchor': [1008.0, 1069.0], 'exposure': 617.0,},
set43={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P177D', 'roll': 335.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-11', 'ext': 9, 'wlshift': -2, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760020', 'exclude': [[3550, 3780], [4050, 4450], [4900, 7000], [1600.0, 2935.0]], 'anchor': [1020.0, 1057.0], 'exposure': 1010.0,},
set44={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 1, 'wlshift': -25, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [4900, 7000], [1600.0, 3000.0]], 'anchor': [1002.0, 1067.0], 'exposure': 469.0,},
set45={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 2, 'wlshift': -38.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [4900, 7000], [1600.0, 2930.0]], 'anchor': [1044.0, 1061.0], 'exposure': 715.0,},
set46={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 3, 'wlshift': -19.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [4900, 7000], [1600.0, 2920.0]], 'anchor': [1016.0, 1065.0], 'exposure': 814.0,},
set47={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 4, 'wlshift': -18.0, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [4900, 7000], [1600.0, 2945.0]], 'anchor': [1011.0, 1066.0], 'exposure': 715.0,},
set48={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 5, 'wlshift': 0, 'groups': ['FB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [1600,3650],[4900, 7000], [1600.0, 2905.0]], 'anchor': [849.0, 1323.0], 'exposure': 1100.0,},
set49={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P177D', 'roll': 330.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-12', 'ext': 6, 'wlshift': -33, 'groups': ['aB'], 'nproc': 1, 'obsid': '00056760022', 'exclude': [ [4900, 7000], [1600.0, 2950.0]], 'anchor': [1024.0, 1064.0], 'exposure': 813.0,},
set50={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'P177D', 'roll': 57.9, 'notes': '', 'wheelpos': 200, 'dateobs': '2006-03-20', 'ext': 1, 'wlshift': 47.5, 'groups': ['WaB'], 'nproc': 1, 'obsid': '00056760004', 'exclude': [ [1600,3600],[4900, 7000], [1600.0, 2935.0]], 'anchor': [1094.0, 984.0], 'exposure': 594.0,},
set51={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 151.2, 'notes': 'second order starts early', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 2, 'wlshift': 0, 'groups': ['G'], 'nproc': 1, 'obsid': '00054254008', 'exclude': [[2600, 7000], [1600.0, 1716.0]], 'anchor': [562.0, 665.0], 'exposure': 1099.0,},
set52={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 355.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['H'], 'nproc': 1, 'obsid': '00054251019', 'exclude': [ [2650.0, 7000.0]], 'anchor': [572.0, 916.0], 'exposure': 925.0,},
set53={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 355.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['H'], 'nproc': 1, 'obsid': '00054251021', 'exclude': [ [2723.0, 7000.0]], 'anchor': [572.0, 904.0], 'exposure': 665.0,},
set54={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 355.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['H'], 'nproc': 1, 'obsid': '00054251021', 'exclude': [ [2510.0, 2620.0], [2730.0, 7000.0]], 'anchor': [572.0, 904.0], 'exposure': 665.0,},
set55={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 355.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['HZ'], 'nproc': 1, 'obsid': '00054251023', 'exclude': [[2447.0, 2634.0], [2730.0, 7000.0]], 'anchor': [592.0, 878.0], 'exposure': 730.0,},
set56={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 264.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-26', 'ext': 1, 'wlshift': 0, 'groups': ['H'], 'nproc': 1, 'obsid': '00055202012', 'exclude': [[2185, 2320], [2521, 7000], [2233.0, 2373.0]], 'anchor': [576.0, 921.0], 'exposure': 1115.0,},
set57={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055501007', 'exclude': [[2733, 7000]], 'anchor': [683.0, 734.0], 'exposure': 168.0,},
set58={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 2, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055501007', 'exclude': [[2730, 7000]], 'anchor': [690.0, 734.0], 'exposure': 82.0,},
set59={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055501008', 'exclude': [[2730, 7000]], 'anchor': [695.0, 601.0], 'exposure': 171.0,},
set60={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055501010', 'exclude': [[2650.0, 7000.0]], 'anchor': [678.0, 572.0], 'exposure': 220.0,},
set61={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055501011', 'exclude': [[2632.0, 7000.0]], 'anchor': [634.0, 865.0], 'exposure': 1337.0,},
set62={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-30', 'ext': 1, 'wlshift': 0, 'groups': ['C'], 'nproc': 1, 'obsid': '00055201005', 'exclude': [[2720, 7000]], 'anchor': [645.0, 762.0], 'exposure': 1358.0,},
set63={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 348.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00054251008', 'exclude': [[2730, 7000]], 'anchor': [615.0, 870.0], 'exposure': 603.0,},
set64={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_g",'name': 'WD0320-539', 'roll': 348.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00054251010', 'exclude': [[2730, 7000]], 'anchor': [613.0, 874.0], 'exposure': 746.0,},
set65={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 348.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-08', 'ext': 1, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00054251012', 'exclude': [[2730, 7000]], 'anchor': [614.0, 872.0], 'exposure': 595.0,},
#set66={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1057+719', 'roll': 350.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 1, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055201001', 'exclude': [[2750, 7000]], 'anchor': [599.0, 880.0], 'exposure': 893.6,},
set67={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 2, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055201001', 'exclude': [[2750, 7000]], 'anchor': [604.0, 881.0], 'exposure': 1454.0,},
set68={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 3, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055201001', 'exclude': [[2750, 7000]], 'anchor': [602.0, 879.0], 'exposure': 1456.0,},
set69={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 4, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055201001', 'exclude': [[2750, 7000]], 'anchor': [603.0, 883.0], 'exposure': 1456.0,},
set70={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.7, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-14', 'ext': 5, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055201001', 'exclude': [[2750, 7000], [2675.0, 3000.0]], 'anchor': [605.0, 880.0], 'exposure': 1043.0,},
set71={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 264.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['I'], 'nproc': 1, 'obsid': '00055202014', 'exclude': [[2750, 7000]], 'anchor': [634.0, 956.0], 'exposure': 1170.0,},
set72={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 119.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-09-26', 'ext': 1, 'wlshift': 0, 'groups': ['J'], 'nproc': 1, 'obsid': '00054250017', 'exclude': [[2730, 7000]], 'anchor': [719.0, 958.0], 'exposure': 893.0,},
set73={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 119.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-09-26', 'ext': 2, 'wlshift': 0, 'groups': ['J'], 'nproc': 1, 'obsid': '00054250017', 'exclude': [[2730, 7000]], 'anchor': [727.0, 943.0], 'exposure': 1285.0,},
set74={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 37.0, 'notes': '2 dips, low', 'wheelpos': 200, 'dateobs': '2012-06-18', 'ext': 1, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00054256002', 'exclude': [[1600,1750],[2175, 2333], [2450, 2750], [2900, 7000]], 'anchor': [682.0, 1434.0], 'exposure': 1362.9,},
set75={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 37.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-18', 'ext': 2, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00054256002', 'exclude': [[1600,1750],[2071, 2235], [2387, 2548], [2900, 7000]], 'anchor': [639.0, 1434.0], 'exposure': 1166.1,},
set76={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-30', 'ext': 1, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00055202006', 'exclude': [[1855, 1970], [2349, 2428], [2800, 7000]], 'anchor': [712.0, 1396.0], 'exposure': 1338.0,},
set77={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 130.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 1, 'wlshift': 16.0, 'groups': ['N'], 'nproc': 1, 'obsid': '00055501019', 'exclude': [[3400, 7000]], 'anchor': [1223.0, 1499.0], 'exposure': 1185.0,},
set78={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 129.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 2, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00055504016', 'exclude': [[2730, 7000]], 'anchor': [843.0, 663.0], 'exposure': 699.0,},
set79={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054253020', 'exclude': [[1600,1690],[2730, 7000]], 'anchor': [797.0, 677.0], 'exposure': 1027.0,},
set80={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054253021', 'exclude': [[2730, 7000]], 'anchor': [796.0, 675.0], 'exposure': 1164.0,},
set81={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054253022', 'exclude': [[2300, 2355], [2730, 7000]], 'anchor': [796.0, 685.0], 'exposure': 571.9,},
set82={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 2, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054253022', 'exclude': [[1600,1825],[2540, 7000], [2232.0, 2376.0]], 'anchor': [812.0, 675.0], 'exposure': 285.0,},
set83={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 151.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 1, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054254008', 'exclude': [[1600,1700],[1930, 1985], [2730, 7000]], 'anchor': [756.0, 735.0], 'exposure': 1100.4,},
set84={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 151.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 3, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054254008', 'exclude': [[1600,1700],[2730, 7000]], 'anchor': [705.0, 767.0], 'exposure': 1041.0,},
set85={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 151.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 4, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054254008', 'exclude': [[1600,1700],[2730, 7000]], 'anchor': [734.0, 742.0], 'exposure': 921.0,},
#set86={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 3, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00054253022', 'exclude': [[2750, 7000]], 'anchor': [811.0, 677.0], 'exposure': 201.0,},
set87={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 195.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-12-09', 'ext': 1, 'wlshift': 0, 'groups': ['M'], 'nproc': 1, 'obsid': '00054254010', 'exclude': [[1911, 1960], [2200, 2600], [2750, 7000]], 'anchor': [593.0, 1056.0], 'exposure': 745.0,},
set88={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 129.9, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 1, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00055503026', 'exclude': [[2730, 7000]], 'anchor': [628.0, 1262.0], 'exposure': 1052.8,},
set89={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 129.9, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 2, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00055503026', 'exclude': [[1930, 1955], [2730, 7000]], 'anchor': [631.0, 1259.0], 'exposure': 1059.1,},
set90={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 1, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00054251045', 'exclude': [[1900, 1940], [2730, 7000]], 'anchor': [612.0, 1295.0], 'exposure': 1103.5,},
set91={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 40.8, 'notes': 'ragged spectrum', 'wheelpos': 200, 'dateobs': '2005-06-27', 'ext': 2, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00054251045', 'exclude': [[2730, 7000], [2200.0, 2310.0]], 'anchor': [628.0, 1299.0], 'exposure': 1284.8,},
set92={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 315.7, 'notes': 'bad background', 'wheelpos': 200, 'dateobs': '2005-03-31', 'ext': 1, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00054252002', 'exclude': [[2433, 2470], [2730, 7000]], 'anchor': [652.0, 1270.0], 'exposure': 585.9,},
set93={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': 'low bkg rate 0.0138/pix', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00055502009', 'exclude': [[1600, 1792], [1925, 2210], [2730, 7000]], 'anchor': [817.0, 1479.0], 'exposure': 1249.3,},
set94={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': 'low bkg rate 0.0139/pix', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00055502010', 'exclude': [[1600, 1785], [1960, 2230], [2730, 7000]], 'anchor': [813.0, 1478.0], 'exposure': 1033.4,},
set95={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 36.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2009-03-26', 'ext': 1, 'wlshift': 0, 'groups': ['OB'], 'nproc': 1, 'obsid': '00057951006', 'exclude': [[4735, 7000], [1600.0, 2700.0], [3770.0, 3990.0]], 'anchor': [920.0, 1273.0], 'exposure': 633.6,},
set96={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 36.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2009-03-26', 'ext': 1, 'wlshift': 15.0, 'groups': ['QB'], 'nproc': 1, 'obsid': '00057954004', 'exclude': [[4650, 7000],  [3783.0, 3981.0], [1600.0, 2895.0]], 'anchor': [1213.0, 1175.0], 'exposure': 632.6,},
set97={'auto':True,'fit2nd':False,'offsetlimit':[100,3],'phatype':"_f",'name': 'P041C', 'roll': 36.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2009-03-24', 'ext': 1, 'wlshift': 0, 'groups': ['KB'], 'nproc': 1, 'obsid': '00057953003', 'exclude': [ [4750, 7000], [1600.0, 2730.0], [3010.0, 3200.0], [3745.0, 3990.0]], 'anchor': [797.0, 994.0], 'exposure': 794.0,},
set98={'auto':True,'fit2nd':False,'offsetlimit':[102,3],'phatype':"_f",'name': 'P041C', 'roll': 36.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2009-03-26', 'ext': 1, 'wlshift': 0, 'groups': ['KB'], 'nproc': 1, 'obsid': '00057953006', 'exclude': [ [4000, 7000], [1600.0, 2710.0], [3030.0, 3190.0], [3765.0, 4100.0]], 'anchor': [800.0, 978.0], 'exposure': 633.6,},
set99={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 1, 'wlshift': 0, 'groups': ['L'], 'nproc': 1, 'obsid': '00054251001', 'exclude': [[2700, 7000]], 'anchor': [927.0, 590.0], 'exposure': 164.7,},
set100={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 2, 'wlshift': 20.0, 'groups': ['L'], 'nproc': 1, 'obsid': '00054251001', 'exclude': [[2700, 7000]], 'anchor': [925.0, 605.0], 'exposure': 340.1,},
set101={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 4, 'wlshift': 0, 'groups': ['L'], 'nproc': 1, 'obsid': '00054251001', 'exclude': [[1895, 1955], [2720, 7000]], 'anchor': [940.0, 622.0], 'exposure': 984.9,},
set102={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 5, 'wlshift': 0, 'groups': ['L'], 'nproc': 1, 'obsid': '00054251001', 'exclude': [[1880, 1960], [2700, 7000], [1900.0, 1999.0]], 'anchor': [917.0, 637.0], 'exposure': 1102.7,},
set103={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 296.5, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 6, 'wlshift': 0, 'groups': ['L'], 'nproc': 1, 'obsid': '00054251001', 'exclude': [[1880, 1960], [2700, 7000], [1900.0, 1999.0]], 'anchor': [904.0, 623.0], 'exposure': 983.9,},
set104={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 315.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-31', 'ext': 1, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00054254002', 'exclude': [[3239.0, 7000]], 'anchor': [1257.0, 1448.0], 'exposure': 971.4,},
set105={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': 'next to other spectrum', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 1, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00055204001', 'exclude': [[1600, 2450], [3300, 7000]], 'anchor': [1398.0, 1290.0], 'exposure': 2047.2,},
set106={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': 'next to other spectrum', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 2, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00055204001', 'exclude': [[1600, 2450], [2800, 7000]], 'anchor': [1399.0, 1287.0], 'exposure': 2047.6,},
set107={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': 'next to other spectrum', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 3, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00055204001', 'exclude': [[1600, 2450], [2800, 7000]], 'anchor': [1402.0, 1287.0], 'exposure': 1043.0,},
set108={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-30', 'ext': 1, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00055204007', 'exclude': [[2910, 7000]], 'anchor': [1334.0, 1324.0], 'exposure': 1338.9,},
set109={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 199.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-08-18', 'ext': 1, 'wlshift': 0, 'groups': ['S'], 'nproc': 1, 'obsid': '00055201009', 'exclude': [[2800, 7000]], 'anchor': [1488.0, 1092.0], 'exposure': 1188.4,},
set110={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 199.3, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-08-18', 'ext': 2, 'wlshift': 0, 'groups': ['S'], 'nproc': 1, 'obsid': '00055201009', 'exclude': [[2800, 7000]], 'anchor': [1496.0, 1051.0], 'exposure': 2073.6,},
set111={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 37.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-18', 'ext': 3, 'wlshift': 0, 'groups': ['F'], 'nproc': 1, 'obsid': '00054256002', 'exclude': [[1600,1750],[2054, 2160], [2900, 7000]], 'anchor': [620.0, 1425.0], 'exposure': 1363,},
set112={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 37.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-19', 'ext': 1, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054257002', 'exclude': [[2730, 7000]], 'anchor': [1349.0, 725.0], 'exposure': 1116.9,},
set113={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD0320-539', 'roll': 37.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-19', 'ext': 2, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054257002', 'exclude': [[2730, 7000]], 'anchor': [1337.0, 784.0], 'exposure': 1412.1,},
set114={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 310.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-20', 'ext': 1, 'wlshift': 0, 'groups': ['YB'], 'nproc': 1, 'obsid': '00057956002', 'exclude': [[1600, 2860], [3044, 3145], [4730, 7000]], 'anchor': [1332.0, 698.0], 'exposure': 1264.5,},
set115={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 308.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-21', 'ext': 1, 'wlshift': 0, 'groups': ['YB'], 'nproc': 1, 'obsid': '00057956004', 'exclude': [[1600, 2812], [3044, 3145], [5200, 7000], [2634.0, 2893.0]], 'anchor': [1363.0, 755.0], 'exposure': 1453.5,},
set116={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'P041C', 'roll': 308.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2012-06-20', 'ext': 1, 'wlshift': 7.0, 'groups': ['FB'], 'nproc': 1, 'obsid': '00057957002', 'exclude': [ [2650, 2960], [4745, 7000], [1600.0, 2850.0]], 'anchor': [627.0, 1386.0], 'exposure': 1264.5,},
set117={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-23', 'ext': 1, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054252008', 'exclude': [[2730, 7000]], 'anchor': [1371.0, 784.0], 'exposure': 981.9,},
set118={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-23', 'ext': 2, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054252008', 'exclude': [[2730, 7000]], 'anchor': [1358.0, 764.0], 'exposure': 745.1,},
set119={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 3, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054252008', 'exclude': [[2730, 7000]], 'anchor': [1356.0, 803.0], 'exposure': 1100.8,},
set120={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 4, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054252008', 'exclude': [[2730, 7000]], 'anchor': [1342.0, 784.0], 'exposure': 1041.5,},
set121={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 147.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-24', 'ext': 5, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00054252008', 'exclude': [[2730, 7000]], 'anchor': [1359.0, 816.0], 'exposure': 1100.3,},
set122={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 350.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-03', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00054253004', 'exclude': [[1906, 1945], [2670.0, 7000.0]], 'anchor': [1179.0, 678.0], 'exposure': 511.2,},
set123={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 350.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-03', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00054253006', 'exclude': [[1985, 2085], [2730, 7000]], 'anchor': [1223.0, 679.0], 'exposure': 423.6,},
set124={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 338.1, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-30', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00055203005', 'exclude': [[2800.0, 7000.0]], 'anchor': [1266.0, 693.0], 'exposure': 1338.8,},
set125={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 2, 'wlshift': 0, 'groups': ['X'], 'nproc': 1, 'obsid': '00055203003', 'exclude': [[1970, 1995], [2620, 2650], [2900, 7000]], 'anchor': [1360.0, 452.0], 'exposure': 106.8,},
set126={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-12', 'ext': 3, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00055203003', 'exclude': [[2750, 7000]], 'anchor': [1304.0, 697.0], 'exposure': 85.5,},
set127={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD1057+719', 'roll': 350.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-03-13', 'ext': 5, 'wlshift': 0, 'groups': ['Y'], 'nproc': 1, 'obsid': '00055203003', 'exclude': [[1823, 1890], [2750, 7000]], 'anchor': [1313.0, 702.0], 'exposure': 85.0,},
set128={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 358.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-05', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00055503010', 'exclude': [[2730, 7000]], 'anchor': [1148.0, 688.0], 'exposure': 1154.2,},
set129={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 358.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-05', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00055503012', 'exclude': [[2318, 2480], [2730, 7000]], 'anchor': [1143.0, 685.0], 'exposure': 1154.1,},
set130={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00055503013', 'exclude': [[2730, 7000]], 'anchor': [1249.0, 693.0], 'exposure': 1493.1,},
set131={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 346.6, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-04-12', 'ext': 1, 'wlshift': 0, 'groups': ['B'], 'nproc': 1, 'obsid': '00055503014', 'exclude': [[2730, 7000]], 'anchor': [1247.0, 701.0], 'exposure': 682.8,},
set132={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 196.8, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-12-09', 'ext': 1, 'wlshift': 0, 'groups': ['V'], 'nproc': 1, 'obsid': '00054252010', 'exclude': [[2730, 7000]], 'anchor': [1036.0, 618.0], 'exposure': 745.4,},
set133={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 359.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 1, 'wlshift': 0, 'groups': ['V'], 'nproc': 1, 'obsid': '00054253015', 'exclude': [[1600, 1758], [2730, 7000]], 'anchor': [1115.0, 613.0], 'exposure': 592.5,},
set134={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 359.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 2, 'wlshift': 0, 'groups': ['V'], 'nproc': 1, 'obsid': '00054253015', 'exclude': [[1690, 1765], [2730, 7000]], 'anchor': [1120.0, 609.0], 'exposure': 487.9,},
set140={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'GD153', 'roll': 129.4, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-11-09', 'ext': 1, 'wlshift': 0, 'groups': ['D'], 'nproc': 1, 'obsid': '00055504016', 'exclude': [[2250, 2370], [2730, 7000]], 'anchor': [863.0, 649.0], 'exposure': 790.5,},
set141={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 315.7, 'notes': 'background way high', 'wheelpos': 200, 'dateobs': '2005-03-31', 'ext': 2, 'wlshift': 0, 'groups': ['E'], 'nproc': 1, 'obsid': '00054252002', 'exclude': [[2850, 7000]], 'anchor': [652.0, 1269.0], 'exposure': 614.4,},
set142={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 146.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-10-22', 'ext': 1, 'wlshift': 0, 'groups': ['R'], 'nproc': 1, 'obsid': '00054251047', 'exclude': [[3500, 7000]], 'anchor': [1310.0, 1449.0], 'exposure': 1336.5,},
set143={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 359.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 1, 'wlshift': 0, 'groups': ['V'], 'nproc': 1, 'obsid': '00054253016', 'exclude': [[1600, 1775], [2730, 7000]], 'anchor': [1125.0, 638.0], 'exposure': 677.4,},
set144={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_g",'name': 'WD0320-539', 'roll': 359.2, 'notes': '', 'wheelpos': 200, 'dateobs': '2005-05-12', 'ext': 2, 'wlshift': 0, 'groups': ['V'], 'nproc': 1, 'obsid': '00054253016', 'exclude': [[1695, 1765], [2730, 7000]], 'anchor': [1136.0, 641.0], 'exposure': 429.7,},
set145={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1121+145', 'roll': 290.0, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2012-07-15', 'ext': 1, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055250081', 'exclude': [[1600,2000],[2900, 7000]], 'anchor': [1004.0, 1094.0], 'exposure': 768.0},
set146={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1121+145', 'roll': 290.0, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2012-07-15', 'ext': 2, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055250081', 'exclude': [[1600,2000],[2900, 7000]], 'anchor': [1012.0, 1093.0], 'exposure': 748.0},
set147={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1121+145', 'roll': 290.0, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2012-07-15', 'ext': 3, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055250081', 'exclude': [[1600,2000],[2900, 7000]], 'anchor': [1035.0, 1108.0], 'exposure': 749.0},
set148={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1121+145', 'roll': 290.0, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2012-07-15', 'ext': 4, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055250081', 'exclude': [[1600,2000],[2900, 7000]], 'anchor': [1014.0, 1086.0], 'exposure': 736.0},
set149={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1121+145', 'roll': 290.0, 'notes': 'stalloss', 'wheelpos': 200, 'dateobs': '2012-07-15', 'ext': 5, 'wlshift': 0, 'groups': ['a'], 'nproc': 1, 'obsid': '00055250081', 'exclude': [[1600,2000],[2900, 7000]], 'anchor': [1021.0, 1088.0], 'exposure': 761.0},
#set??=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='000',wheelpos=200, ext=, anchor=[0.,0.],dateobs='2000-00-00',exposure=0.,roll= None,groups=[''],name='',nproc=0,wlshift=0,exclude=[[6900,7000]],notes='')                ,
set150={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1657+343', 'roll': 311.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-07-31', 'ext': 1, 'wlshift': 0, 'groups': ['a'], 'nproc': 0, 'obsid': '00055900071', 'exclude': [[2691,3105],[3243,7000]], 'anchor': [1001.7, 1070.9], 'exposure': 1115.9,},
set151={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1657+343', 'roll': 311.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-07-31', 'ext': 2, 'wlshift': 0, 'groups': ['a'], 'nproc': 0, 'obsid': '00055900071', 'exclude': [[2750,7000]], 'anchor': [1066.9, 1070.1], 'exposure': 1460.4,},
set152={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1657+343', 'roll': 311.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-07-31', 'ext': 3, 'wlshift': 0, 'groups': ['a'], 'nproc': 0, 'obsid': '00055900071', 'exclude': [[2689,7000]], 'anchor': [1040.4, 1042.1], 'exposure': 1461.3,},
set153={'auto':True,'fit2nd':False,'offsetlimit':None,'phatype':"_f",'name': 'WD1657+343', 'roll': 311.0, 'notes': '', 'wheelpos': 200, 'dateobs': '2008-07-31', 'ext': 4, 'wlshift': 0, 'groups': ['a'], 'nproc': 0, 'obsid': '00055900071', 'exclude': [[2692,3060],[3250,7000]], 'anchor': [1044.3, 1074.7], 'exposure': 869.8,},
set154=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600003',wheelpos=200, ext=1, anchor=[835.4,1078.4],dateobs='2005-02-14',exposure=418.1,roll=89.8,groups=[''],name='G63-26',nproc=0,wlshift=26,exclude=[[1600,2400],[3030,3220],[3400,7000]],notes=''),
set155=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600018',wheelpos=200, ext=1, anchor=[1039.4,1073.3],dateobs='2005-04-04',exposure=302.4,roll=19.5,groups=['a'],name='G63-26',nproc=0,wlshift=12,exclude=[[1600,2160],[2180,2430],[3400,7000]],notes=''),
set156=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600018',wheelpos=200, ext=2, anchor=[766.4,1143.5],dateobs='2005-04-04',exposure=604.5,roll=19.5,groups=[''],name='G63-26',nproc=0,wlshift=8,exclude=[[1600,2100],[3400,7000]],notes=''),
set157=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00055600021',wheelpos=200, ext=1, anchor=[1016.6,1128.8],dateobs='2010-07-26',exposure=1460.4,roll=294.0,groups=['a'],name='G63-26',nproc=9,wlshift=0,exclude=[[1600,2250],[3400,7000]],notes=''),
set158=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_g",obsid='00055600014',wheelpos=200, ext=1, anchor=[1002.4,1115.3],dateobs='2005-03-28',exposure=1127.4,roll=35.8,groups=['a'],name='G63-26',nproc=0,wlshift=0,exclude=[[1600,2000],[2180,2430],[3300,7000]],notes=''),
set159=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00052973002',wheelpos=200, ext=1, anchor=[1313.0,1374.1],dateobs='2013-03-20',exposure=1331.4,roll=43.0,groups=['RB'],name='P041C',nproc=0,wlshift=10,exclude=[[1600,2865],[4085,4460],[4950,7000]],notes=''),
set160=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00052973002',wheelpos=200, ext=2, anchor=[1303.3,1381.3],dateobs='2013-03-20',exposure=1331.4,roll=43.0,groups=['RB'],name='P041C',nproc=0,wlshift=10,exclude=[[1600,2765],[4100,4480],[5000,7000]],notes=''),
set161=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057974002',wheelpos=200, ext=1, anchor=[617.5,915.7],dateobs='2013-03-21',exposure=1331.4,roll=43.0,groups=['HB'],name='P041C',nproc=0,wlshift=3,exclude=[[1600,2650],[4100,4480],[4536,7000]],notes=''),
set162=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057974002',wheelpos=200, ext=2, anchor=[617.4,915.7],dateobs='2013-03-21',exposure=1062.7,roll=43.0,groups=['HB'],name='P041C',nproc=0,wlshift=-5,exclude=[[1600,2650],[4120,4500],[4610,7000]],notes=''),
set163=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00057975002',wheelpos=200, ext=1, anchor=[940.7,487.1],dateobs='2013-03-22',exposure=1107.0,roll=43.0,groups=['LB'],name='P041C',nproc=0,wlshift=0,exclude=[[1600,2860],[4100,4460],[4990,7000]],notes=''),
set164=dict(auto=True,fit2nd=False,offsetlimit=None,phatype="_f",obsid='00442039000',wheelpos=200, ext=1, anchor=[1152.5,787.3],dateobs='2011-01-12',exposure=49.0,roll=None,groups=[''],name='BD+25_4655',nproc=0,wlshift=0,exclude=[[1715,1730],[1740,4000],[2750,7000]],notes='very bright v=9.61 blue calibration source for coi-loss'),
set165=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057530012',wheelpos=200, ext=1, anchor=[ 982.3,1142.2],dateobs='2014-07-05',exposure=763.1,roll=228.0,groups=[''],name='AGK+81_266',nproc=0,wlshift=0,exclude=[[2000,2150],[2980,7000]],notes='coi'),
set166=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_f",obsid='00057530014',wheelpos=200, ext=1, anchor=[ 987.0,1118.9],dateobs='2014-07-05',exposure=823.7,roll=226.0,groups=[''],name='AGK+81_266',nproc=0,wlshift=0,exclude=[[2960,7000]],notes='coi'),
)   
   obs_G202m65=dict(
set90=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850012',wheelpos=160, ext=1, anchor=[1101.,1032.],dateobs='2005-05-18',exposure=1458.3,roll= 13.5,groups=['aA'],name='G202-65',nproc=0,wlshift=53,exclude=[[1600,2500],[3600,7000]],notes=''),
set91=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850018',wheelpos=160, ext=1, anchor=[1080., 990.],dateobs='2005-09-08',exposure=318.4,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=20,exclude=[[1600,2500],[3600,7000]],notes=''),          
set92=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850020',wheelpos=160, ext=1, anchor=[1082.,1009.],dateobs='2005-09-08',exposure=381.7,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3600,7000]],notes=''),          
set93=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850022',wheelpos=160, ext=1, anchor=[1160.,1038.],dateobs='2005-10-02',exposure=1201.5,roll=237.2,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3600,7000]],notes=''),         
set94=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850023',wheelpos=160, ext=1, anchor=[1130.,1027.],dateobs='2005-10-05',exposure=1572.8,roll=233.9,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3600,7000]],notes=''),         
set95=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850026',wheelpos=160, ext=1, anchor=[1135.,1043.],dateobs='2005-10-07',exposure=1159.2,roll=247.0,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3600,7000]],notes=''),         
set96=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850030',wheelpos=160, ext=1, anchor=[1154., 990.],dateobs='2010-07-26',exposure=1461.4,roll=315.0,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3600,7000]],notes=''),         
set97=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851010',wheelpos=160, ext=1, anchor=[673., 998.],dateobs='2005-05-18',exposure=986.8,roll=13.6,groups=[''],name='G202-65',nproc=0,wlshift=15,exclude=[[1600,2500],[3700,7000]],notes=''),            
set98=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851012',wheelpos=160, ext=1, anchor=[1199., 569.],dateobs='2005-09-08',exposure=745.6,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),           
set99=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851012',wheelpos=160, ext=2, anchor=[1204., 572.],dateobs='2005-09-08',exposure=802.5,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),           
set100=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851013',wheelpos=160, ext=1, anchor=[1451., 754.],dateobs='2005-10-03',exposure=1132.2,roll=236.8,groups=[''],name='G202-65',nproc=0,wlshift=-9,exclude=[[1600,2500],[3700,7000]],notes=''),                
set101=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851015',wheelpos=160, ext=1, anchor=[1362., 830.],dateobs='2005-10-07',exposure=1021.5,roll=234.5,groups=[''],name='G202-65',nproc=0,wlshift=19,exclude=[[1600,2500],[3700,7000]],notes=''),                
set102=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055851015',wheelpos=160, ext=2, anchor=[1382., 767.],dateobs='2005-10-07',exposure=1055.9,roll=234.5,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),         
set103=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055852008',wheelpos=160, ext=1, anchor=[1086.,1469.],dateobs='2005-05-18',exposure=1456.6,roll=13.5,groups=[''],name='G202-65',nproc=0,wlshift=25,exclude=[[1600,2500],[3700,7000]],notes=''),         
set104=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055852010',wheelpos=160, ext=1, anchor=[ 675., 883.],dateobs='2005-09-08',exposure= 744.7,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=4,exclude=[[1600,2500],[3700,7000]],notes=''),         
set105=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055852010',wheelpos=160, ext=2, anchor=[ 697., 896.],dateobs='2005-09-08',exposure= 746.3,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=3,exclude=[[1600,2500],[3700,7000]],notes=''),         
set106=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055852012',wheelpos=160, ext=1, anchor=[ 843., 684.],dateobs='2005-10-03',exposure=1138.9,roll=237.5,groups=[''],name='G202-65',nproc=0,wlshift=-10,exclude=[[1600,2500],[3700,7000]],notes=''),               
set107=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055852014',wheelpos=160, ext=1, anchor=[ 772., 683.],dateobs='2005-10-12',exposure=1574.2,roll= 87.0,groups=[''],name='G202-65',nproc=0,wlshift=-8,exclude=[[1600,2500],[3700,7000]],notes=''),                
set108=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055853006',wheelpos=160, ext=1, anchor=[1120., 577.],dateobs='2005-05-18',exposure=1456.5,roll= 13.6,groups=[''],name='G202-65',nproc=0,wlshift=14,exclude=[[1600,2500],[3700,7000]],notes=''),                
set109=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055853008',wheelpos=160, ext=1, anchor=[1515.,1122.],dateobs='2005-09-08',exposure= 745.8,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=10,exclude=[[1600,2500],[3700,7000]],notes=''),                
set110=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055853008',wheelpos=160, ext=2, anchor=[1512.,1122.],dateobs='2005-09-08',exposure= 804.2,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=6,exclude=[[1600,2500],[3700,7000]],notes=''),         
set111=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055853009',wheelpos=160, ext=1, anchor=[1403.,1391.],dateobs='2005-10-03',exposure=1130.7,roll=235.9,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),         
set112=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055853011',wheelpos=160, ext=1, anchor=[1541.,1055.],dateobs='2005-10-07',exposure=1063.7,roll=230.8,groups=[''],name='G202-65',nproc=0,wlshift=-7,exclude=[[1600,2500],[3700,7000]],notes=''),                
set113=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854008',wheelpos=160, ext=1, anchor=[1541.,1055.],dateobs='2005-05-18',exposure=1456.7,roll= 13.5,groups=[''],name='G202-65',nproc=0,wlshift=21,exclude=[[1600,2500],[3700,7000]],notes=''),                
set114=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854010',wheelpos=160, ext=1, anchor=[ 971.,1426.],dateobs='2005-09-08',exposure= 650.7,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=8,exclude=[[1600,2500],[3700,7000]],notes=''),         
set115=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854010',wheelpos=160, ext=2, anchor=[ 995.,1432.],dateobs='2005-09-08',exposure= 805.7,roll=271.7,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),         
set116=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854012',wheelpos=160, ext=1, anchor=[ 843.,1391.],dateobs='2005-10-01',exposure= 417.8,roll=238.9,groups=[''],name='G202-65',nproc=0,wlshift=-3,exclude=[[1600,2500],[3700,7000]],notes=''),                
set117=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854012',wheelpos=160, ext=2, anchor=[ 843.,1391.],dateobs='2005-10-02',exposure=1271.2,roll=238.9,groups=[''],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2500],[3700,7000]],notes=''),         
set118=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854012',wheelpos=160, ext=3, anchor=[ 843.,1391.],dateobs='2005-10-02',exposure=1271.8,roll=238.9,groups=[''],name='G202-65',nproc=0,wlshift=-2,exclude=[[1600,2500],[3700,7000]],notes=''),                
set119=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055854014',wheelpos=160, ext=1, anchor=[ 715.,1315.],dateobs='2005-10-07',exposure=1038.9,roll=231.4,groups=[''],name='G202-65',nproc=0,wlshift=4,exclude=[[1600,1800],[3700,7000]],notes='was 1800-3800'),
set135=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850008',wheelpos=200, ext=1, anchor=[950.,1131.],dateobs='2005-03-25',exposure=926.0,roll=72.3,groups=['a'],name='G202-65',nproc=1,wlshift=18.4,exclude=[[1600,1800],[4000,7000]],notes=''),
set136=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850008',wheelpos=200, ext=2, anchor=[980.,1103.],dateobs='2005-03-25',exposure=127.7,roll=72.3,groups=['a'],name='G202-65',nproc=1,wlshift=18.5,exclude=[[1600,1800],[4000,7000]],notes=''),
set137=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850008',wheelpos=200, ext=3, anchor=[979.,1103.],dateobs='2005-03-25',exposure=379.3,roll=72.3,groups=['a'],name='G202-65',nproc=1,wlshift=17.2,exclude=[[1600,1800],[4000,7000]],notes=''),
set138=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850009',wheelpos=200, ext=1, anchor=[955.,1101.],dateobs='2005-03-17',exposure=98.9,roll=83.1,groups=['a'],name='G202-65',nproc=1,wlshift=18.2,exclude=[[1600,1800],[4000,7000]],notes=''),
set139=dict(auto=True,fit2nd=True,offsetlimit=None,phatype="_g",obsid='00055850028',wheelpos=200, ext=1, anchor=[1015.,1101.],dateobs='2005-11-07',exposure=1369.1,roll=199.4,groups=['a'],name='G202-65',nproc=0,wlshift=0,exclude=[[1600,2360],[4000,7000]],notes=''),
)
   if key == 'all': 
      key = 'wheelpos'
      keyvalue = wheelpos
      
   if debug == 3: 
      return calfildir, caldict     
   if debug == 2:
      return obs_G202m65 

   if ((key == 'groups') | (key == 'anchor') | (key == 'name') | 
       (key == 'obsid') | (key == 'obsid') | (key == 'wheelpos') ):
   
      if use_obs == None:
         if wheelpos == 160:
            datb = obs_160.copy() 
            group = _select_group(obs_160, key, keyvalue, xyrange=xyrange, ) 
         if wheelpos == 200: 
            datb = obs_200.copy()
            group = _select_group(obs_200, key, keyvalue, xyrange=xyrange, )
         if wheelpos == 955: 
            datb = obs_955.copy()
            group = _select_group(obs_955, key, keyvalue, xyrange=xyrange, )
         if wheelpos == 1000: 
            datb = obs_1000.copy()
            group = _select_group(obs_1000, key, keyvalue, xyrange=xyrange, )
      else: 
         try:
            os.system('cp /calibration/grism/fluxcal2/obs_py/'+str(use_obs)+
                      ' '+os.getenv('HOME')+'/pymodules/use_obs.py')
            from . import use_obs
            print("using obs database: ", use_obs.__version__)
                  
            if wheelpos == 160:
               datb = use_obs.obs_160
               group = _select_group(use_obs.obs_160, key, keyvalue, xyrange=xyrange, ) 
            if wheelpos == 200: 
               datb = use_obs.obs_200.copy()
               group = _select_group(use_obs.obs_200, key, keyvalue, xyrange=xyrange, )
         except:
             print("FATAL ERROR in main_ea: cannot read the database")
             raise 
 
   else: 
      group = []
      print("No groups selected for processing")
   
   if not _verify_pha(datb,phafiledir,chatter=1): 
      print(" ") 
         
   if debug == 1:    
      return calfildir,calfiles,calfilID,obs_160,obs_200,obs_955,obs_1000,group, caldict, obsdict
   
   if make_pha_files:
      import uvotgetspec
      if coi_width != None:
          uvotgetspec.coi_half_width = coi_width*0.5
      print("WARNING: WRITING PHA FILES TO CURRENT DIRECTORY")
      if wheelpos == 160:
         group = dict(_select_group(obs_160, key, keyvalue, xyrange=xyrange, ) )
         for setkey in list(group.keys()):
             obskey = group[setkey]
             if obskey['auto']: 
                 ra,dec = radeclist()[obskey['name']]
                 obsid = obskey['obsid']                 
                 obsid = obskey['obsid']                 
                 local_mounted, status, datadir2 = retrieve_obsid_path_from_local_archive(obsid)
                 # local archive ?
                 if local_mounted and status: 
                     indir = datadir2+obsid+'/uvot/image/'
                 else:
                     indir='/data/various/caltemp/'+obsid+'/uvot/image/'
                 uvotgetspec.getSpec(ra,dec,obsid,obskey['ext'],
                     indir=indir,
                     use_lenticular_image=use_lenticular_image,
                     fit_second=obskey['fit2nd'],
                     plot_img=False, plot_raw=False, plot_spec=False,
                     offsetlimit=obskey['offsetlimit'],
                     clobber=True,chatter=2,)
             else: 
                 print("skipping ",obskey)
                 #ra,dec = radeclist()[key['name']
                 #obsid = key['obsid']
                 #uvotgetspec(ra,dec,key['obsid'],obsid,indir='/data/various/caltemp/'+obsid+'/uvot/image/',
                 #    use_lenticular_image=True,fit_second=key['fit2nd'],clobber=True,chatter=2,offsetlimit=1)
                 
                 print("manual extraction") 
      elif wheelpos == 200:
         group = dict(_select_group(obs_200, key, keyvalue, xyrange=xyrange, ) )
         for setkey in list(group.keys()):
             obskey = group[setkey]
             if obskey['auto']: 
                 ra,dec = radeclist()[obskey['name']]
                 obsid = obskey['obsid']                 
                 local_mounted, status, datadir2 = retrieve_obsid_path_from_local_archive(obsid)
                 # local archive ?
                 if local_mounted and status: 
                     indir = datadir2+obsid+'/uvot/image/'
                 else:
                     indir='/data/various/caltemp/'+obsid+'/uvot/image/'
                 uvotgetspec.getSpec(ra,dec,obsid,obskey['ext'],
                     indir=indir,
                     use_lenticular_image=use_lenticular_image,
                     fit_second=obskey['fit2nd'],
                     plot_img=False, plot_raw=False, plot_spec=False,
                     offsetlimit=obskey['offsetlimit'], 
                     clobber=True,chatter=2,)
      elif wheelpos == 955:
         group = dict(_select_group(obs_955, key, keyvalue, xyrange=xyrange, ) )
         outlog = open('auto_log_955.txt','a')
         for setkey in list(group.keys()):
             obskey = group[setkey]
             if obskey['auto']: 
                 ra,dec = radeclist()[obskey['name']]
                 obsid = obskey['obsid']
                 if (obskey['offsetlimit'] != None) & (type(obskey['offsetlimit']) != list): 
                     offsetlimit = [obskey['offsetlimit']+100,2]
                 else: offsetlimit=16    # all observations for offsetlimit 22
                 print("next observation key="+setkey+"  obsid+ext="+obsid+"+"+str(obskey['ext']))
                 outlog.write( "next observation key="+setkey+"  obsid+ext="+obsid+"+"+str(obskey['ext'])+"\n")
                 local_mounted, status, datadir2 = retrieve_obsid_path_from_local_archive(obsid)
                 # local archive ?
                 if local_mounted and status: 
                     indir = datadir2+obsid+'/uvot/image/'
                 else:
                     indir='/data/various/caltemp/'+obsid+'/uvot/image/'
                 uvotgetspec.getSpec(ra,dec,obsid,obskey['ext'],
                     indir=indir,
                     use_lenticular_image=use_lenticular_image,
                     fit_second=obskey['fit2nd'],
                     plot_img=False, plot_raw=False, plot_spec=False,
                     offsetlimit=offsetlimit,
                     clobber=True,chatter=2,)
      elif wheelpos == 1000:
         from pylab import figure,draw
         group = dict(_select_group(obs_1000, key, keyvalue, xyrange=xyrange, ) )
         outlog = open('auto_log_1000.txt','a')
         for setkey in list(group.keys()):
             obskey = group[setkey]
             if obskey['auto']: 
                 ra,dec = radeclist()[obskey['name']]
                 obsid = obskey['obsid']
                 if (obskey['offsetlimit'] != None) & (type(obskey['offsetlimit']) != list): 
                     offsetlimit = [obskey['offsetlimit']+100,2]
                 else: offsetlimit=13    # all observations for offsetlimit 22
                 print("next observation key="+setkey+"  obsid+ext="+obsid+"+"+str(obskey['ext']))
                 outlog.write( "next observation key="+setkey+"  obsid+ext="+obsid+"+"+str(obskey['ext'])+"\n")
                 local_mounted, status, datadir2 = retrieve_obsid_path_from_local_archive(obsid)
                 # local archive ?
                 if local_mounted and status: 
                     indir = datadir2+obsid+'/uvot/image/'
                 else:
                     indir='/data/various/caltemp/'+obsid+'/uvot/image/'
                 uvotgetspec.getSpec(ra,dec,obsid,obskey['ext'],
                     indir=indir,
                     use_lenticular_image=use_lenticular_image,
                     fit_second=obskey['fit2nd'],
                     plot_img=debug==2, plot_raw=debug==2, plot_spec=debug==2,
                     offsetlimit=offsetlimit,
                     clobber=True,chatter=2,)
                 if debug == 2: 
                     print("last observation key="+setkey+"  obsid+ext="+obsid+"+"+str(obskey['ext']))
                     figure(2)
                     draw()
                     ans = input("continue")    
             else: 
                 print("skipping ",obskey)
                 outlog.write("skipping obskey = "+setkey+"\n")
                 #ra,dec = radeclist()[key['name']
                 #obsid = key['obsid']
                 #uvotgetspec(ra,dec,key['obsid'],obsid,indir='/data/various/caltemp/'+obsid+'/uvot/image/',
                 #    use_lenticular_image=True,fit_second=key['fit2nd'],clobber=True,chatter=2,offsetlimit=1)
                 
                 print("manual extraction") 
         outlog.close()          
      else:
         print("Nothing done. I need to verify that automated processing is possible for this mode.")  
      print("completed auto items")
      print("coi_width used = ",coi_width) 
      print("confirm half width",uvotgetspec.coi_half_width)      
      return   
   
   for g in group: 
      if g[1]['nproc'] > 10: 
         print("processed before ", g[1]['nproc'], '  times.')
         
      try: 
         print("\n\nprocessing ", g[0],'\n')
         calobsfile = calfildir+caldict[g[1]['name'].lower()]
         obsid = g[1]['obsid']
         gri = 'ugv_1ord_'
         ext = g[1]['ext']
         if ((wheelpos == 160) | (wheelpos == 200)): gri = 'ugu_1ord_'
         phafile = phafiledir+'sw'+obsid+gri+str(ext)+g[1]['phatype']+'.pha'      
         eff_area_out = eff_area_outdir+'/sw'+obsid+gri+str(ext)+g[1]['phatype']+'.EA'
         niet = g[1]['exclude']
      except:
         print("problem with g=",g)
         raise
      
      
      if chatter > 0: 
         print("parameters in call calobsfile:")
         print("  calobsfile   = ",calobsfile)
         print("  obsid        = ",obsid)
         print("  phafile      = ",phafile)
         print("  eff_area_out = ",eff_area_out)
         print("  exclude      = ",niet)
      
      try:
        g_update = g
        # for getting results with the old coi-formula, use option=2
        g_update = uvotcal.EffAreaCalFF(calobsfile, phafile, eff_area_out, spectral_order=1, 
                 nsmo_calspec=21, nsmo_obsspec=31, calobs_err=calobs_err, apcorr_err=apercorr_err,
                 figno=15, frametime=0.0110329, option=1, db=g, ignore_pha_quality=ign_q,
                 auto=auto, autoexclude=niet, wlshift=wlshift, chatter=chatter, clobber=clobber) 
        if type (g_update[1]) != dict:
           print("problem case")
           raise                 
      except: 
         print("########## WARNING: EffAreaCalFF failed for ",g) 
         pass
         
      if not auto: log.write(str(g[0])+'='+str(g[0])+',\n'+str(g_update)+',\n\n') 

   log.close()
   
         
   
   
   
def _select_group(dictio, keyname, keyvalue, xyrange=[100,100], chatter=0 ):
   ''' helper routine to main_ea() : 
   
   extract a sub-list of dictionaries in dictio that contains keyvalue in keyname
   
   special: for 'anchor' a range can be given as a list in xyrange (2 values)
   
   '''
   import numpy as np  
   group = []
   
   if type( dictio ) == dict:
      items = list(dictio.keys())
      if chatter > 1: print('input dictionary')
   elif type( dictio ) == list:
      items = list(range(len(dictio)))
      if chatter > 1: print('input list')
   else:
      print("wrong type input for dictio")
      raise      
    
   for k in items:
      if type(dictio) == dict:
         set1 = dictio[k]
         setname = k
      else: 
         set1 = dictio[k][1]
         setname = dictio[k][0]  
      key = set1[keyname]
      if keyname == 'groups':
         for m in key[0]:
            if chatter > 2: print(m, keyvalue, m == keyvalue)
            if m == keyvalue:
               group.append((setname,set1))   
      elif keyname == 'anchor':
         dx = xyrange[0]
         dy = xyrange[1]
         x = np.abs(key[0] - keyvalue[0])
         y = np.abs(key[1] - keyvalue[1])
         if (x <= dx) & (y <= dy) :
            group.append((setname,set1))  
      elif keyname == 'name':
         kval = np.array(keyvalue)
         if (key == kval).any(): group.append( (setname,set1))                         
      else: 
         if key == keyvalue: group.append( (setname,set1))             
   return group   


def _verify_pha(db,phafiledir,get_keyword=None,chatter=0):
   '''
   Check that the phafile is present when data in db
   
   
   '''
   import os
   status = True
   if get_keyword != None: 
      os.system('touch keyword.lis')
   for k in list(db.keys()):
      g = db[k]
      obsid = g['obsid']
      wheelpos = g['wheelpos']
      gri = 'ugv_1ord_'
      ext = g['ext']
      phatype= g['phatype']
      if ((wheelpos == 160) | (wheelpos == 200)): gri = 'ugu_1ord_'
      phafile = phafiledir+'/sw'+obsid+gri+str(ext)+phatype+'.pha'
    
      if (get_keyword == None) & (not os.access(phafile,os.F_OK)): 
         status = False
         print("WARNING : "+k+'  '+g['name']+'  '+phafile+" is missing!") 
      elif get_keyword != None:
         imgfile = '/calibration/grism/fluxcal/'+g['name']+'/'+obsid+'/uvot/image/sw'+obsid
         if wheelpos < 500: imgfile +='ugu_dt.img+'+str(ext)
         else: imgfile +='ugv_dt.img+'+str(ext)
         command = "echo "+imgfile+': >> keyword.lis'
         if chatter > 0: print(command)
         os.system(command)
         command = 'fkeyprint '+imgfile+' '+get_keyword+' >> keyword.lis'
         os.system(command)      
   if chatter > 0: print("verification completed")           
   return status


def EAmerge(wheelpos, database=None, name=None, key='groups',
    keyvalue='F', debug=0, eff_area_dir=None, phatype='_f',
    EAmerged_dir=None, xyrange=[100,100], markers=False, fig=10, chatter=1, clobber=True):   
    '''
    select and merge effective area files 
    
    example call: 
    
    cal4.EAmerge(200,database='obs_200_nproc_1.py',name='WD1657+343',
         key='groups',keyvalue='a',debug=0,
         eff_area_dir='/Volumes/data/grism/fluxcal2/eff_area_200/ea_2.5/',
         EAmerged_dir='/Volumes/data/grism/fluxcal2/eff_area_200/ea_2.5_sums/',)
    
    '''
    from . import uvotcal
    import os
    
    if (database == None):
       print("using database in cal4.main_ea()")
       datb = get_obs_xxx(wheelpos)    
       group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter ) 
    else:   
      try:
         os.system('cp /calibration/grism/fluxcal2/obs_py/'+str(database)+' '+os.getenv('HOME')+'/pymodules/use_obs.py')
         from . import use_obs
         print("using obs database: ", use_obs.__version__)
                  
         if wheelpos == 160:
               datb = use_obs.obs_160 
               group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter ) 
         if wheelpos == 200: 
               datb = use_obs.obs_200.copy()
               group = _select_group(datb, key, keyvalue, xyrange=xyrange, )
      except:
             print("FATAL ERROR in EAmerge: cannot read the database")
             print("input parameters used : ",database, key, keyvalue, xyrange)
             print("datb : ",datb)
             raise  
    
    EAfiles = []
    
    if name.lower() == 'all':
       name = 'all'
    elif name != None:
       group = _select_group(group, 'name', name, )
    else: name = 'all'   
       
    if key.lower() == 'all':
       key = 'wheelpos'
       keyvalue = wheelpos   
    
    if chatter > 4: print("printing group:  ", group)
    
    for g in group: 

      try: 
         name = g[1]['name']
         obsid = g[1]['obsid']
         gri = 'ugv_1ord_'
         ext = g[1]['ext']
         if ((wheelpos == 160) | (wheelpos == 200)): gri = 'ugu_1ord_'
         eafile = eff_area_dir+'/sw'+obsid+gri+str(ext)+g[1]['phatype']+'.EA'       
         EAfiles.append(eafile)   
         
      except:
         print("problem with ",g)
         raise
             
    if key != 'anchor':      
         eff_area_out = EAmerged_dir+'/eamerged_'+name+'_'+str(wheelpos)+'_'+key+'_'+str(keyvalue)+'.EA'
    else:
         eff_area_out = EAmerged_dir+'/eamerged_'+name+'_'+str(wheelpos)+'_'+key+str(keyvalue[0])+'_'+str(keyvalue[1])+'_rxy'+str(xyrange[0])+'_'+str(xyrange[1])+'.EA'
    
    m = None
    if markers: m = ['<','>','+','x','^','v','s','o','*','1','2','3','4','p','d','h',',','H','D','d','-','|']
    if chatter > 4: 
       print("EAmerge: EAfiles selected are:",EAfiles)
    uvotcal.mergeEffAreas(EAfiles, eff_area_out,figno=fig,markers=m,wait=False)
    # done

def plot_ea_200_default(s=20):    
   '''Working on producing a single curve out of this multitude ... '''
   from pylab import where,polyfit,clf,plot,legend,where,searchsorted,polyval,array,isfinite,xlim,ylim
   from scipy import interpolate
   from uvotmisc import rdTab 
   import numpy as np
   
   t_1657 = rdTab('eamerged_WD1657+343_200_anchor1000_1080_rxy70_70.EA.err')
   t_202  = rdTab('eamerged_G202-65_200_anchor1000_1080_rxy70_70.EA.err')
   t_63   = rdTab('eamerged_G63-26_200_groups_a.EA.err')
   t_1057 = rdTab('eamerged_WD1057+719_200_anchor1000_1080_rxy70_70.EA.err')
   t_0320 = rdTab('eamerged_WD0320-539_200_anchor1000_1080_rxy70_70.EA.err')
   t_177  = rdTab('eamerged_P177D_200_anchor1000_1080_rxy70_70.EA.err')   
   t_041  = rdTab('eamerged_P041C_200_anchor1000_1080_rxy70_70.EA.err') 
   
   q = where(t_1657[:,] > 1900.)

   q2 = where(t_1657[:,0] <= 1900.)

   c2_3 = polyfit(t_1657[q[0],0],t_1657[q[0],4],3)

   q1 = where(t_041[:,0] > 2800)

   clf()

   plot(t_1657[q[0],0],t_1657[q[0],4],'g',lw=1,label='WD1657+343')
   plot(t_1657[q2[0],0],t_1657[q2[0],4],'g',lw=2,label='WD1657+343')
#   plot(t_1657[q[0],0],polyval(c2_3,t_1657[q[0],0]),'g',lw=1,label='POLY 3 fit')
   plot(t_177[:,0],t_177[:,4],'r',lw=2,label='P177D')
   plot(t_0320[:,0],t_0320[:,4],'m',lw=1.5,label='WD0320-539')
   plot(t_1057[:,0],t_1057[:,4],'g',lw=1.5,label='WD1057+719')
   plot(t_041[q1[0],0],t_041[q1[0],4],'b',lw=1.5,label='P041C')
   plot(t_63[:,0],t_63[:,4]*1.,'r',lw=1.5,alpha=0.6,label='G63-26')
   #plot(t_202[:,0],t_202[:,4]*1.,'y',lw=1,label='G202-65')
   #plot(t_202[:,0],t_202[:,4]*1.05,'y',lw=3,label='1.05 x G202-65')
   #plot(t_202[:,0],t_202[:,4]*1.05,'k',lw=1,label='NOLABEL')

   legend(loc=0) 

   # derive weighted mean 
   
   # up to 2133A  WD1657+343, WD0320-539, WD1057+719 
   q1 = where(t_1657[:,0] < 2133.0)
   q2 = where(t_0320[:,0] < 2133.0)
   q3 = where(t_1057[:,0] < 2133.0)
   wmin = min([min(t_1657[q1[0],0]),min(t_0320[q2[0],0]),min(t_1057[q3[0],0])])
   ea1 = t_1657[q1[0],1]
   ea2 = t_0320[q2[0],1]
   ea3 = t_1057[q3[0],1]
   ea_out = []
   wav_out = []
   weight = []
   for i in range(int(wmin),2133,1):
      w = 1.0*i    
      i1 = searchsorted(t_1657[q1[0],0],w)
      i2 = searchsorted(t_0320[q2[0],0],w)
      i3 = searchsorted(t_1057[q3[0],0],w)
      ea = []
      we = []
      if t_1657[q1[0],0][i1] == w: 
         ea.append(   t_1657[q1[0],4][i1]/max([t_1657[q1[0],5][i1],3.0]) )
         we.append(1./max([3.,t_1657[q1[0],5][i1]]) )
      if t_0320[q2[0],0][i2] == w: 
         ea.append(   t_0320[q2[0],4][i2]/max([3.,t_0320[q2[0],5][i2]]) )
         we.append(1./max([3.,t_0320[q2[0],5][i2]]) )
      if t_1057[q3[0],0][i3] == w: 
         ea.append(   t_1057[q3[0],4][i3]/max([3.,t_1057[q3[0],5][i3]]) )
         we.append(1./max([3.,t_1057[q3[0],5][i3]]) )
      if len(ea) > 0: 
         wav_out.append(w)
         ea_out.append(sum(ea))
         weight.append(sum(we))

   # from 2134-2750 A    
   # use polynomial fit for WD1657 
   #q1 = where((t_1657[:,0] > 2133.0) & (t_1657[:,0] < 2750.0))
   q2 = where((t_0320[:,0] > 2133.0) & (t_0320[:,0] < 2750.0))
   q3 = where((t_1057[:,0] > 2133.0) & (t_1057[:,0] < 2750.0))
   q4 = where((t_63  [:,0] > 2133.0) * (t_63  [:,0] < 2750.0))
   wmin = 2134.0
   #ea1 = t_1657[q1[0],1]
   ea2 = t_0320[q2[0],1]
   ea3 = t_1057[q3[0],1]
   ea4 = t_63  [q4[0],1]
   for i in range(2134,2750,1):
      w = 1.0*i
      i1 = searchsorted(t_1657[:,0],w)
      i2 = searchsorted(t_0320[q2[0],0],w)
      i3 = searchsorted(t_1057[q3[0],0],w)
      i4 = searchsorted(t_63  [q4[0],0],w)
      
      ea = [polyval(c2_3,w)/t_1657[i1,5]]
      we = [1./t_1657[i1,5]]
      if t_0320[q2[0],0][i2] == w: 
         ea.append(   t_0320[q2[0],4][i2]/max([3.,t_0320[q2[0],5][i2]]) )
         we.append(1./max([3.,t_0320[q2[0],5][i2]]) )
      if t_1057[q3[0],0][i3] == w: 
         ea.append(   t_1057[q3[0],4][i3]/max([3.,t_1057[q3[0],5][i3]]) )
         we.append(1./max([3.,t_1057[q3[0],5][i3]]) )
      if t_63[q4[0],0][i4] == w: 
         ea.append(   t_63[q4[0],4][i4]/max([3.,t_63[q4[0],5][i4]]) )
         we.append(1./max([3.,t_63[q4[0],5][i4]]) )
      if len(ea) > 1: 
         wav_out.append(w)
         ea_out.append(sum(ea))
         weight.append(sum(we))
         
   # from 2750 - 2940 use P63-26
   q1 = where((t_63  [:,0] > 2750.0) * (t_63  [:,0] < 2940.0))
   #wmin = 2750.0
   ea1 = t_63[q1[0],1]
   for i in range(2750,2940,1):
      w = 1.0*i 
      i1 = searchsorted(t_63[q1[0],0],w)
      ea = []
      we = []
      if t_63 [q1[0],0][i1] == w: 
         ea.append(   t_63[q1[0],4][i1]/max([3.,t_63[q1[0],5][i1]]) )
         we.append(1./max([3.,t_63[q1[0],5][i1]]) )
      if len(ea) > 0: 
         wav_out.append(w)
         ea_out.append(sum(ea))
         weight.append(sum(we))

   # from 2940 onward combine P177D and P041C
   q1 = where((t_177[:,0] > 2940.0) & (t_177[:,0] < 6000.0))
   q2 = where((t_041[:,0] > 2940.0) & (t_041[:,0] < 6000.0))
   q3 = where((t_63 [:,0] > 2940.0) & (t_63 [:,0] < 6000.0))
   wmax = np.min( [ np.max(t_177[q1[0],0]), np.max(t_041[q2[0],0]), np.max(t_63[q3[0],0]) ] )
   ea1 = t_177[q1[0],1]
   ea2 = t_041[q2[0],1]
   ea3 = t_63 [q3[0],1]
   for i in range(2940,int(wmax),1):
      w = 1.0*i
      i1 = searchsorted(t_177[q1[0],0],w)
      i2 = searchsorted(t_041[q2[0],0],w)
      i3 = searchsorted(t_63 [q3[0],0],w)
      ea = []
      we = []
      if t_177[q1[0],0][i1] == w: 
         ea.append(   t_177[q1[0],4][i1]/max([3.,t_177[q1[0],5][i1]]) )
         we.append(1./max([3.,t_177[q1[0],5][i1]]) )
      if t_041[q2[0],0][i2] == w: 
         ea.append(   t_041[q2[0],4][i2]/max([3.,t_041[q2[0],5][i2]]) )
         we.append(1./max([3.,t_041[q2[0],5][i2]]) )
      if t_63 [q3[0],0][i3] == w: 
         ea.append(   t_63[q3[0],4][i3]/max([3.,t_63[q3[0],5][i3]]) )
         we.append(1./max([3.,t_63[q3[0],5][i3]]) )
      if len(ea) > 0: 
         wav_out.append(w)
         ea_out.append(sum(ea))
         weight.append(sum(we))
   for i in range(int(wmax),int(max(t_041[q2[0],0]))):  
      w = 1.0*i
      i2 = searchsorted(t_041[q2[0],0],w)
      ea = []
      we = []
      if t_041[q2[0],0][i2] == w: 
         ea.append(   t_041[q2[0],4][i2]/max([3.,t_041[q2[0],5][i2]]) )
         we.append(1./max([3.,t_041[q2[0],5][i2]]) )
      if len(ea) > 0: 
         wav_out.append(w)
         ea_out.append(sum(ea))
         weight.append(sum(we))

   # weight problem 2270-2310 - filter out
   wav_out= array(wav_out)
   ea_out = array(ea_out)
   weight = array(weight)
   v = isfinite(wav_out) & isfinite(ea_out) & isfinite(weight) & ((wav_out < 2270) | (wav_out > 2310))
   ea_out = ea_out[v]/weight[v]
   error  = 1/weight[v]
   wav_out = wav_out[v]
   weight = weight[v]
   plot(wav_out,ea_out,'k.',label='Best Fit',lw=1,markersize=3.0)
   
   c = interpolate.splrep(wav_out,ea_out,s=s,w=weight,xb=1600.,xe=wav_out[-1])
   w = array(list(range(1600,5000,1)))
   ea = interpolate.splev(w,c)
   v = where(ea > 0.)
   w = w[v[0]]
   ea = ea[v[0]]
   plot(w,ea,'k',lw=2,alpha=0.7,label='smoothing spline fit')
   legend(loc=0)
   xlim(1650,5000);ylim(0,22)   
   return wav_out,ea_out,error,weight,w,c,ea

def _plot_position(db,axes=None,title=None):
   '''plot the locations  '''
   import numpy as np
   import pylab as plt
   name = []
   ank = []
   M= len(db)
   for i in range(M): name.append(db[list(db.keys())[i]]['name'])
   for i in range(M):  ank.append(db[list(db.keys())[i]]['anchor'])
   name = np.array(name)
   ank  = np.array(ank)
   if axes == None:
      axes = plt.subplot(111)
   q=np.where(name == 'P041C') 
   if len(q[0]) > 0: axes.plot(ank[q[0],0],ank[q[0],1],'*',color='gold',markersize=9,label='P041C')   
   q=np.where(name == 'P177D') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'r^',markersize=7,label='P177D')
   q=np.where(name == 'WD1057+719') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'yo',markersize=5,label='WD1057+719')
   q=np.where(name == 'WD1657+343') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'bo',markersize=5,label='WD1657+343')
   q=np.where(name == 'WD1121+145') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'s',color='gold',markersize=5,label='WD1121+145')
   q=np.where(name == 'WD0320-539') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'or',markersize=5,label='WD0320-539')
   q=np.where(name == 'G202-65') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'b*',markersize=9,alpha=0.5,label='G202-65')
   q=np.where(name == 'G63-26') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'*m',markersize=9,label='G63-26')
   q=np.where(name == 'GD153') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'sb',markersize=5,label='GD153')
   q=np.where(name == 'GD50') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'<g',markersize=7,label='GD50')
   q=np.where(name == 'GD108') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'>r',markersize=7,label='GD108')
   q=np.where(name == 'LTT9491') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'h',color='pink',markersize=7,label='LTT9491')
   q=np.where(name == 'BPM16274') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'h',color='darkorange',markersize=7,label='BPM16274')
   q=np.where(name == 'AGK+81_266') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'*',color='purple',markersize=9,label='AGK+81_266')
   q=np.where(name == 'BD+25_4655') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'s',color='darkblue',markersize=7,label='BD+25_4655')
   q=np.where(name == 'BD+33_2642') 
   if len(q[0]) > 0:  axes.plot(ank[q[0],0],ank[q[0],1],'^',color='g',markersize=7,label='BD+33_2642')
   axes.legend(loc=0,title=title)
   
   
def _plot_box(x,y,dx,dy,axes=None):
   import pylab as plt
   if axes == None: axes=plt.subplot(111)
   axes.plot([x-dx,x+dx],[y-dy,y-dy],'k')
   axes.plot([x+dx,x+dx],[y-dy,y+dy],'k')
   axes.plot([x+dx,x-dx],[y+dy,y+dy],'k')
   axes.plot([x-dx,x-dx],[y+dy,y-dy],'k')   

def _plot_spectra_with_cal(db,allfiledir,calfile,arf1=None,tw=None,framtime=0.0110322,
     titel=' ',xlim=[],wvals=[],fudge=1.0,option=1,phatype='_f',chatter=0):
   '''
   To compare the spectra in DB with the standard spectrum
   
   Example:
    X = cal4.main_ea(200,debug=1,key='name',keyvalue='GD153',phafiledir='/calibration/grism/fluxcal2/eff_area_200/all_2.5/',
        eff_area_outdir='/calibration/grism/fluxcal2/eff_area_200/ea_2.5/',calobs_err=1.,apercorr_err=0.5,ign_q=True,auto=True,phatype="_f",)
    db200 = X[4]
    agroup = cal4._select_group(db200,'groups','a',chatter=1)
    g = cal4._select_group(agroup, 'name', 'GD153')
    allfiledir = '/Volumes/data/grism/fluxcal2/eff_area_200/all_2.5/'
    calfile = '/Volumes/data/grism/fluxcal2/uvotified/gd153_stisnic_003_uvotified.ascii'
    figure(5)
    cal4._plot_spectra_with_cal(g, allfiledir, calfile, chatter=2)
    
    xlim : [] plot range x-axis
    wvals : [] wavelength range for making statistics
       
   '''
   import numpy as np
   import os
   import pylab as plt
   import pyfits
   import uvotio
   import uvotmisc
   from stsci.convolve import boxcar
   from scipy import interpolate
   
   status = True
   anker = np.zeros(2)
   if type(db) == list: db = dict(db)

   _wave = []
   _flux = []
   _label = []
   #_err  = []

   plt.subplot(2,1,1)
   
   for k in list(db.keys()):
      if chatter > 0: print('processing ',db[k])
      g = db[k]
      obsid = g['obsid']
      wheelpos = g['wheelpos']
      gri = 'ugv_1ord_'
      ext = g['ext']
      anker = g['anchor']
      if ((wheelpos == 160) | (wheelpos == 200)): gri = 'ugu_1ord_'
      phafile1 = obsid+gri+str(ext)
      phafile = allfiledir+'/sw'+obsid+gri+str(ext)+phatype+'.pha'
    
      if (not os.access(calfile,os.F_OK)) : 
         status = False
         print("WARNING : "+'  '+calfile+" is missing!") 
      elif (not os.access(phafile,os.F_OK)): 
         status = False
         print("WARNING : "+k+'  '+g['name']+'  '+phafile+" is missing!") 
      else:
         f = pyfits.open(phafile)
         f.info()
         wave = f[2].data['lambda']
         pixno = f[2].data['pixno']
         rate = f[2].data['netrate']
         bkgrate = f[2].data['bg_r']
         cosprate = f[2].data['cosp1rat']
         cobgrate = f[2].data['cobg1rat']
         H =f[1].header
         tstart = H['tstart']
         hist = H.get_history()
         ank = uvotmisc.get_keyword_from_history(hist,'anchor1')
         if tw == None: 
            tw = uvotmisc.get_keyword_from_history(hist,'trackwidth')
            if type(tw) == str: tw = float(tw)
            if tw == None: tw = 2.5
         #print "trackwidth["+phafile+"] = ",tw   
         if framtime == None: frametime=H['framtime']
         if np.abs(H['framtime'] - 0.0110322)/0.0110322 > 0.01:
            print("WARNING: HEADER FRAME TIME is ",H['framtime'],' s')
         #anker = float(ank.split(',')[0].split(' ')[1]) , float(ank.split(',')[1].split(')')[0])
         anker = np.array(anker)
         fnew = uvotio.rate2flux(wave,rate,H['wheelpos'],
                bkgrate=bkgrate,
                co_sprate = cosprate,
                co_bgrate = cobgrate,
                pixno=pixno,
                anker=anker,
                #trackwidth=tw,
                frametime=framtime,
                fudgespec=fudge,
                swifttime=tstart,
                #option=option, 
                arf1=arf1,
                chatter=chatter,)
         _wave.append(wave)
         _flux.append(fnew)
         #plt.plot(wave,boxcar(fnew,(5,)),label=phafile1,lw=1.2)
         plt.plot(wave,fnew,label=phafile1,lw=1.) 
         _label.append(phafile1)   
         f.close()
   fcal = uvotmisc.rdTab(calfile)        
   plt.plot(fcal[:,0],fcal[:,1],'k--',lw=1.7,label='CALSPEC')   
   plt.legend(loc=0)
   plt.xlabel('$\lambda(\AA)$',fontsize=15)
   plt.ylabel('flux (erg /cm$^2$ /s /$\AA$)',fontsize=14)
   plt.ylim(1e-15,2e-13)
   if len(xlim) > 0: plt.xlim(xlim)
   plt.title(titel)
   plt.subplot(2,1,2)
   _fcal = interpolate.interp1d(fcal[:,0],fcal[:,1],)
   for i in range(len(_wave)):
      q = (_wave[i] > fcal[:,0].min() ) & (_wave[i] < fcal[:,0].max()) & np.isfinite(_flux[i])
      if wvals == []: wvals=[_wave[i][q].min(), _wave[i][q].max()]
      val = (_wave[i] > wvals[0] ) & (_wave[i] < wvals[1]) & np.isfinite(_flux[i])
      #print  "len q true=",np.array(q,dtype=int).sum()
      #print _wave[i][q]
      #print "  "
      #print _flux[i][q]
      ydif = (_flux[i][q]/_fcal(_wave[i][q])-1.0)
      vdif = (_flux[i][val]/_fcal(_wave[i][val])-1.0)
      print("%s : percent offset (mean)=%6.2f percent noise (std_dev)=%6.2f"%(phafile1,100*vdif.mean(),100*vdif.std())) 
      plt.plot(_wave[i][q],100.*ydif, label=_label[i] )
   plt.legend(loc=0)
   plt.xlabel('$\lambda(\AA)$',fontsize=15)
   plt.ylabel('percent flux measured off ',fontsize=14)
   plt.ylim(-50,50)
   plt.axhline(color='k')
   if len(xlim) > 0: 
      plt.xlim(xlim)
      #plt.hlines(0,xlim[0],xlim[1])
   if chatter > 0: print("done")     
   return status
    

def radeclist():
   '''dictionary of ra,dec in decimal degrees, ... for all calibration targets'''
   return {
   "P177D"     :[239.8065542,  47.6116139,],
   "P041C"     :[222.9915708,  71.721500, ],
   "P330E"     :[247.8909333,  30.1462528,],
   "GSPC P177-D"     :[239.8065542,  47.6116139,],
   "GSPC P041-C"     :[222.9915708,  71.721500, ],
   "GSPC P330-E"     :[247.8909333,  30.1462528,],
   "WD1657+343":[254.7129625,  34.3148667,],
   "WD1057+719":[165.1426000,  71.6341639,],
   "WD0320-539":[ 50.5617250, -53.7546333,],
   "WD1026+453":[157.4384458,  45.1179167,],
   "WD1121+145":[171.0632042,  14.2295889,],
   "WD0308-565":[ 47.4497000, -56.3970278,],
   "BPM16274"  :[ 12.5151542, -52.1376333,],
   "LTT9491"   :[349.8965125, -17.0918444,],
   "GD50"      :[ 57.2085042,  -0.9744806,],
   "GD108"     :[150.1968583,  -7.5585583,],
   "GD153"     :[194.2597125,  22.0313667,],
   "G63-26"    :[201.1275417,  20.4561333,],
   "WR86"      :[259.5960875, -34.4085080,],
   "WR52"      :[199.6166375, -58.1371111,],
   "WR121"     :[281.0548042,  -3.7993778,],
   "WR1"       :[ 10.8683208,  64.7598250,],
   "WR4"       :[ 40.2986292,  56.7304861,],
   "BD+25_4655":[329.92490  ,  26.43261  ,],
   "AGK+81_266":[140.33000  ,  81.724333 ,],
   "BD+33_2642":[237.99952  ,  32.94842  ,],
   }
  
def _stage_from_archive(wheelpos,fromdir,todir,):  
    """stage the flux calibration data """
    import os
    calfildir,calfiles,calfilID,obs_160,obs_200,obs_955,obs_1000,group, caldict, obsdict = main_ea(wheelpos,debug=1)
    if wheelpos == 160: 
       obs_set = obs_160
    elif wheelpos == 200:
       obs_set = obs_200
    elif wheelpos == 955:
       obs_set = obs_955
    elif wheelpos == 1000:
       obs_set = obs_1000
    else:
       raise IOError("wrong wheelpos entered" )           
    for key in list(obs_set.keys()):
        x = obs_set[key]
        print(os.system('cp -r '+fromdir+'/'+x['obsid']+'  '+todir))


def _find_wave_offsets(wheelpos,name=None,phadir='./'):
    """
    After shifting the spectra to match the reference using uvotspec.adjust_wavelength_manually()
    read the pha headers and collect the wavelength shift. 
    
    I need to be in the directory with the pha files
    """
    from astropy.io import fits
    if name == None:
        key = 'wheelpos'
        value = wheelpos
    else:
        key = 'name'
        value = name    
    obs_xxx = get_obs_xxx(wheelpos,key=key,keyvalue=value,phadir=phadir)
    if wheelpos > 500: ug = 'ugv_1ord_'
    else: ug = 'ugu_1ord_'
    woffsets = []
    obj = []
    ank = []
    expo = []
    fn = []
    for key in obs_xxx:
       x = obs_xxx[key] 
       filename = 'sw'+x['obsid']+ug+str(x['ext'])+x['phatype']+'.pha' 
       print('working on '+key+' : '+filename)
       hdr = fits.getheader(filename,2)
       if 'WAVSHFT' in hdr: 
           woffsets.append(hdr['WAVSHFT'])
           obj.append(x['name'])
           ank.append(x['anchor'])
           expo.append(x['exposure'])
           fn.append(filename)     
    return {'woffsets':woffsets,'name':obj,'anker':ank,'exposure':expo,'filename':fn}

def get_obs_xxx(wheelpos,
                phadir='./',
                key='wheelpos',keyvalue=None,
                xyrange=[2100,2100],
                chatter=0):
    X = main_ea(wheelpos,debug=1,key='all',
                chatter=0, clobber=False, 
                use_obs=None, phafiledir=phadir) 
    if keyvalue == None: 
        key = 'wheelpos'
        keyvalue = wheelpos             
    if wheelpos == 160: 
        datb = X[3]
        group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter )
    elif wheelpos == 200: 
        datb = X[4]
        group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter  )
    elif wheelpos == 955: 
        datb = X[5]
        group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter  )
    elif wheelpos == 1000: 
        datb = X[6]
        group = _select_group(datb, key, keyvalue, xyrange=xyrange, chatter=chatter )
    else:
        raise IOError("problem reading database cal4.main_ea(); wheelpos invalid")
    if key == 'wheelpos':       
       return datb        
    else:
       return dict(group) 


def plot_all_spectra(wheelpos,key='wheelpos',keyvalue=None,phatype=None):
    from pylab import figure,subplot,plot,title,xlim,ylim
    from astropy.io import ascii,fits
    obs_xxx = get_obs_xxx(wheelpos,key=key,keyvalue=keyvalue)
    calfiledir, caldict = main_ea(wheelpos, debug=3,)
    count = 0
    for key in obs_xxx:
       print(key)
       x = obs_xxx[key] 
       print(x)
       phat = phatype 
       if phat == None: phat=x['phatype']
       if wheelpos > 500: ug = 'ugv_1ord_'
       else: ug = 'ugu_1ord_'
       filename = 'sw'+x['obsid']+ug+str(x['ext'])+phat+'.pha' 
       print('working on '+key+' : '+filename)
       if count == 0:
          fig = figure()            
       subplot(3,3,count)     
       #subplot(count/3,count-count/3*3,count)
       # get, plot reference spectrum
       sp = ascii.read(calfiledir+caldict[x['name'].lower()])
       plot(sp['col1'],sp['col2'],'k',label='reference')
       f = fits.open(filename)
       d = f[2].data
       plot(d['lambda'],d['flux'],'b',lw=1.5,label='uvot',alpha=0.4)
       xlim( d['lambda'].min(), d['lambda'].max())
       q = (sp['col1'] > d['lambda'].min() ) & (sp['col1'] < d['lambda'].max())
       ylim(1e-14, sp['col2'][q].max())
       title(key+':'+x['name']+': '+filename ,fontsize=10)
       f.close()
       count += 1
       if count == 10: count = 0    
       
       
def _adjust_wavelengths_pha_spectra(wheelpos, name,ylim=[0.2e-14,2.0e-14],
        phatype=None,printall=True):
    """
    Short script to run through adjusting all pha spectra of object 
    
    Parameters
    ==========
    wheelpos : int
    name : str
       calibration object name as in caldict
    ylim : list
       range for plot
    phatype : [None, '_g','_f']      
       select pha file type:
       when 'None', the type is read from the cal4/6 database
       when 'g', the pha files of type 'g' (used graspcorr aspect) are done
       when 'f', the pha files of type 'f' (used a lenticular file for aspect) are done  
       
    returns
    =======
    list of files affected
    then interactive adjustment process
    
    use   _find_wave_offsets() to get statistics 
    
    must run from correct pha directory!  
    """
    import os
    from astropy.io import ascii,fits
    import uvotspec
    calfildir, caldict = main_ea(wheelpos,debug=3)
    obs_xxx = get_obs_xxx(wheelpos,phadir='./',)
    this = _select_group(obs_xxx,'name',name.upper(),xyrange=[3000,3000])       
    sp = ascii.read(calfildir+caldict[name.lower()])
    ug = 'ugu_1ord_'
    if wheelpos > 500: ug = 'ugv_1ord_'
    if printall:
       for s,x in this: 
          print(s, x['name']+':  '+x['obsid']+'  sw'+x['obsid']+ug+str(x['ext'])+"_f.pha  roll=",x['roll'],'  wlshift=',x['wlshift'])
    for s,x in this:
        if phatype == None:
           file='sw'+x['obsid']+ug+str(x['ext'])+x['phatype']+'.pha'
        else:
           file='sw'+x['obsid']+ug+str(x['ext'])+phatype+'.pha'
        if os.access(file,os.F_OK):   
              
           try:
               uvotspec.adjust_wavelength_manually(file=file,ylim=ylim,reference_spectrum=sp,recalculate=True)
           except: pass

def calibrate_delta_sensitivity():
    """
    * make list of files to process
    * obtain extracted spectra
    * extract in 20A bins from 1700A to 1900A the mean over the bin of the 
      net count rate+ error, net flux+error, coincidence loss factor, time of exposure, 
      source id, obsid, anchor position
    * put it all in a table
    *** OFF LINE:  
    * for each source & bin plot net count rate + err (time)
      and flux+err(time)
      and flux/reference spectrum(time)  
      
    """