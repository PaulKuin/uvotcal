# cal3.py
#
# calibration of the order layout uv clocked grism
# 
# 2013-04-10 add measured anchor position (anchor) where missing
#            add dispersion solution (dispers1st) where missing 
#            use scaled "simple" zemax model
#            add GSPC P041-C, GSPC P177D, and G63-26 spectra to fill in gaps at large wavelengths but
#            we use a different anchor at the Ca II lines anchor_397  
# 
# Define the data structure for each calibration spectrum "cur##" with ## from excel spreadsheet
# followed by utilities 
#
# in 2012 some of the upper left corner observations were repeated with slightly different 
# anchors to weigh the curvature solution. 
#
# copyright N.P.M. Kuin (2008-2016) No distribution is allowed without consent from the 
# copyright owner.
#
from numpy import array, zeros, asarray, atleast_1d, polyfit, polyval,\
 arange, where, linspace, meshgrid

#====================================================================
#
'''
the measured parameters for the anchor and dispersion of the calibration spectra plus 
utilities for basic operations

Each observation has a data structure labeled using the same number as in the excel file.
UV clocked in range 0-99, uv nominal 100-199.

Parameters in the data structures are:
  
  obsdir : str
    directory name of the object (e.g.,/Volumes/data/grism/`obsid`)
  obsid : str
    swift observation id
  ext : int
    extension number in fits file 
  *extname : str
    name extension 
  anchor : ndarray, list
    coordiante anchor in  ? image/ ? detector coordinates
  angle : float
    angle used to rotate the spectrum to determine the fitted coefficients  
  zmxangle : float
    angle from calibration file (zemax interpolated)          
  flatangle :float
    the angle that lines up similar wavelengths at the same Y value in orders(0,1)
  present0,present1,present2,present3 : bool
    indicate if order is present on detector.  
  dist12 : float
    measured pixel distance first and second order anchor both at 260nm
  *e_dist12 : float
    error in dist12    
  ratio12 : float
    the ratio between measured dist12 and the zemax model distance
  dispers1st : list, ndarray
    measured dispersion coefficients for the first order
  dispers2nd : list, ndarray
    measured dispersion coefficients for the second order 
  *obslines : list    
    identified line wavelength + pixel positions
  disp2nd_2000 : list, ndarray
    dispersion using 200nm as anchor 
  ratio12_2000 : float 
    ratio of observed anchor distance to zemax model distance 
  dist12_2000 : float
    observed pixel distance first order anchor at 260nm to second order anchor at 200nm    
  coef0,coef1,coef2,coef3 : list, ndarray
    measured coefficients of the Y-offset of the spectral track as a function 
    of wavelength for the zeroth, first, second and third orders
  dlim0L,dlim0U,dlim1L,dlim1U,dlim2L,dlim2U,dlim3L,dlim3U: int
    measured limits for each order in pixel coordinate (measured from anchor)
    L for lower, U for upper limit, order 0,1,2,3.
  sig0coef,sig1coef,sig2coef,sig3coef : list, ndarray
    coefficients for the width (gaussian fit sigma parameter) of the 
    spectral track, expressed as a polynomial as function of pixel distance
    to the anchor.  
  *v2_ratio12_200 : float
    ratio of the dist12_2000 pixel distance of first anchor at 260nm and 
    second order anchor at 200nm to the simplified zemax model value. 
 
Notes
-----
Transforming from the first to the second order dispersion is 
roughly given by coordinates::
  
   dis2 = dis1*1.75+620 pix
   lambda2 = polyval(coef1,dis1) angstrom
   
adopting dist12=620pix gives coef2=[-5.067e-11,-2.919e-08,3.76e-04,1.805,2600]
and dist12_2000=359.7 gives coef2_2000=[-5.067e-11,4.372e-08,3.682e-04,1.532,2000]               
'''
# read the table of coefficients/get the coefficients of the Y(dis) offsets and limits[]
# stored with array of angles used. 
# angle is the angle used to rotate the spectrum to determine the fitted coefficients
# the coefficients are valid between upper and lower limits
# the x-coordinate is distance to the first order anchor in pix 
# flatangle is the angle that lines up similar wavelengths at the same Y value in orders(0,1)
# zmxangle is angle from calibration file (zemax interpolated)
# present logicals indicate if order is present on detector.
# dist12 and ratio12 are the distance between the first order 
#   anchor and that of 2nd and 3rd orders.
# NPMK 2010-08-01  
# updated with second order parameters and dispersion 
# 2011-09-19 added dist12_2000, etc. ; commented out dist12 without second order disp.
#    those need to be redone by hand for 2000. 
# ===========================================
# source positions
# ===========================================
WR52_position = [199.61664,-58.13711]
WR86_position = [259.59608905, -34.408508567]
WR121_position = [281.05481,-3.79939]
TPYX_position = []
# ==================================================================== 
#  UV Clocked Grism
# ====================================================================
cur25=dict(
obsdir='WR52',
obsid="00057925001",
ext=1,
anchor=[1068.3,1733.6], #[1070.9,1733.6]),
field=[0.00596,   0.124779],
anchor1_2530=[1085.9,1721.4], 
anchor2_2000=[855,1870],
angle=34.7398,
zmxangle=34.7398,
flatangle=34.4,
dist12=564.6,
ratio12=0.917,
dis1= [[-342.9,1721.3],[-290.8,1814.1],[-239.5,1908.],[-199.6,2011.7],[-154.4,2106.5],\
    [-88.7,2298.],[-54.15,2405],[-14.3,2530],[6.4,2595],[62.3,2787],[98.5,2906],[141.5,3066.3],\
    [160.0,3137.4],[228.2,3409],[304.7,3718.],[334.7,3818],[364.7,3943],[458.9,4335],\
    [482.2,4436],[532,4649]],
dis2= [[128.7,1814.],[194,1902],[204,1915],[262.1,2000.],[433.3,2298],[491.2,2405]],
dis1corr=-7.9,
# ... but subtract 7.9 pix from dis1,dis2 coordinate to line up anchor     
dispers1st=[ 4.186865e-10, -1.423e-06, 1.7638e-03, 3.310, 2600.0],
dispers2nd=[  9.555e-04, 2.1619, 2600.00],
disp2nd_2000 = [ 9.555e-04, 1.543, 2000.0], #disp2nd_2000=array([  9.59503329e-04,   1.48874806e+00,   1.98815011e+03]), 
ratio12_2000=None, 
dist12_2000=254.22, # was dist12=245.25,
#
coef0 =  array([ -2.76254168e-04,  -5.25010466e-01,  -2.22160081e+02]) ,
dlim0L = -779,
dlim0U = -554,
present0 = True,
#
coef1=array([ -4.89951770e-08,   6.06448849e-05,  -2.55912368e-02, 1.71815034e-01]),
dlim1L = -367,
dlim1U = +529,
present1 = True,
#
coef2 = array([  3.95831960e-05,  -5.95688197e-02,   2.38158301e+01]),
dlim2L =  +31,
dlim2U = +540,
present2 = True,
coef3 = array([ -0.04777007,  45.70397429]),
dlim3L = +402,
dlim3U = +561,
present3 = True,
#
# the coefficients for the width of gaussian profile across the spectrum
# f = amplitude * exp( - ((y-y0)/sig)**2 ) 
# sig? = polyval( sig?coef, dis ) 
sig0coef = None,
sig1coef = array([3.2876e-6,-1.8132e-3,3.1770]),
sig2coef = array([6.9066e-6, -8.588e-3,6.751 ]),
sig3coef = array([5.6]),  # array([-0.0458,26.4] 
extname = "gu256770200I",
e_dist12 = 3,
) # revised 2013-04-11
#===================================================================
cur19=dict(
obsdir='WR52',
obsid="00057919001",
ext=1,
anchor=[1364.1,1506.6], #array([1365.3,1510.5]),
field=[0.055049, 0.088426],
anchor1_2530=[1381.5,1494.4],
anchor2_2000=[1137,1652],
angle=35.0667,
zmxangle=35.0667,
flatangle=34.7,
dist12=614.1, #613.1,
ratio12=0.963,
dis1=[[-348.1,1721.3],[-295.3,1814.1],[-245.5,1908],[-200.3,2000],[-168,2106],[-93.2,2297],[-56.5,2405],\
  [-16.2,2530],[6.3,2595],[63.,2787],[98.6,2906],[141.15,3066.3],[161.1,3137],[233.2,3405],[313.5,3718],\
  [370.6,3943],[469.8,4335],[496.96,4436],[549.9,4649],[827.5,5802],],
dis2=[[70,1721],[142,1814],[211,1908],[277,2000],[455.6,2297],[516.3,2405],[584.8,2530],[783.6,2906]],
# subtract 5.7 pixels for dis1,
dis1corr=-5.7, 
dispers1st=[  8.19421469e-10,  -1.73331974e-06,   1.60767791e-03, 3.28775111e+00,   2.60029921e+03],
dispers2nd=[  5.14933554e-04,   1.87652018e+00,   2.60011349e+03], # array([5.869e-4, 1.922, 2598.6]),
disp2nd_2000=[  5.14933554e-04,   1.51174125e+00,   2.00005239e+03], #array([  5.86900000e-04,   1.50695051e+00,   1.99237059e+03]), 
ratio12_2000=None,
dist12_2000=271.3, #259.50,
#
coef0=array([  4.14200681e-04,   4.14507076e-01,   9.21058394e+01]),
dlim0L= -745.5,
dlim0U= -554.0,
present0=True,
#
coef1=array([ -2.63072412e-08,   5.01517982e-05,  -2.33636992e-02, 2.53606187]),
dlim1L= -347.6,
dlim1U=  845.0,
present1= True,
#
coef2=array([  3.22303679e-05,  -5.21558612e-02,   2.60063811e+01]),
dlim2L=  31,
dlim2U= 958,
present2= True,
   #
coef3=array([ -3.01440747e-02,   4.04424451e+01]),
dlim3L = 456,
dlim3U = 967,
present3= True,
#
sig0coef=array([  4.83862904e-03,   7.859]),
sig1coef=array([  1.38586828e-06,  -7.29718907e-04,   3.199]),
sig2coef=array([  3.50974004e-06,  -3.60756604e-03,   4.957]),
sig3coef=array([5.88]),
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "gu256759694I",
e_dist12 = 3,
) # revised 2013-04-12
#===================================================================
cur09=dict(
obsdir='WR52',
obsid="00057909001",
ext=1,
anchor=[656,1568.3],#array([658.1,1568.5]),
field=[-0.0618443071809, 0.102016038569 ],
anchor1_2530=[673.1,1556.6],
anchor2_2000=[458,1702],
angle=34.278,
zmxangle=34.278,
flatangle=34.9,
dist12=551.9, # 545.2,
ratio12=0.9,
dis1=[[-347.1,1721.3],[-290.9,1814.1],[-241.8,1908],[-200.,2000],[-163.6,2106],[-91.2,2297],\
  [-56.4,2405],[-17.4,2530],[4.2,2595],[57.5,2787],[90.7,2906],[135.5,3066.3],[152.4,3137],\
  [221.5,3405],[303.8,3718]],
dis2=[[114.6,1814],[181.6,1908],[305.1,2105],],
# need offset dis1 of -4.1 
dis1corr=-4.1,
dispers1st=[ -1.12639582e-09,  -2.02256263e-06,   1.90626873e-03,3.37538180e+00,   2.60015417e+03],
dispers2nd=[  1.01653104e-03,   2.23077856e+00,   2.60071089e+03],  #array([8.8538e-4,2.0722,2600.1]),
disp2nd_2000=[  1.01653104e-03,   1.59239707e+00,   2.00047231e+03], #array([8.85380000e-04,   1.45266348e+00,   1.98347729e+03]), 
ratio12_2000=None,
dist12_2000=242.0,
#
# flatanglerr =0.3,
#
# minus first order present with lines
# 
#presentm1=True,
#
coef0=array([ -9.16716103e-04,  -1.22254483e+00,  -3.96309856e+02]),
dlim0L=-649.,
dlim0U=-597.6,
present0=True,
#
coef1=array([ -2.01059660e-08,   3.69589489e-05,  -3.12823047e-02, -4.85219123e-01]),
dlim1L=-376,
dlim1U=380,
present1=True,
#
coef2=array([  2.55116403e-05,  -4.95103458e-02,   11.684548]),
dlim2L=25,
dlim2U=380,
present2=True,
#
coef3=None,
dlim3L=None,
dlim3U=None,
present3=False,
#
sig0coef = None,
sig1coef=array([  6.63819996e-07,  -1.44455715e-03,   3.31091665e+00]),
sig2coef=array([ -1.93569712e-05,   1.37682064e-03,   4.81921157e+00]),
sig3coef=None,
#  *2013 NEW*
extname = "GU256752319I",
e_dist12 = 10,
) # revised 2013-04-12
#===================================================================
cur24=dict(
obsdir='WR52',
obsid="00057924001",
ext=1,
anchor=[1243.1,464.5],#array([1247.6, 464.3]),
field=[0.0382805345877, -0.075911242007],
anchor1_2530=[1261,452],
anchor2_2000=None,
angle= 34.91,
zmxangle= 34.91,
flatangle=35.3,
dist12=642.0,
ratio12=None,
dis1=[[-345.9,1721.1],[-294.7,1814.3],[-246.6,1908],[-203.2,2000],[-95.35,2297],\
      [-58.15,2405],[-17.9,2530],[6.5,2595],[95.5,2906],[164.0,3137],[235.8,3405],\
      [289.,3607],[322,3718],[348.5,3818],[377.4,3943],\
      [562.9,4649.],[867.1,5802] ],
# subtract 5.45 pix from dis1 
dis1corr=-5.45,       
dis2=[[477,2297],[607.8,2530]],  # possibly also [1100,3405???],[1250,????]],
dispers1st=[ 4.34826e-10, -1.4150e-06, 1.4792e-03, 3.214, 2600.146],
dispers2nd= [0, 1.781, 2600.7],
ratio12_2000=None,
dist12_2000=None,
disp2nd_2000=None, 
#  
coef0=array([-0.01087432,  6.74158999]),
dlim0L=-769,
dlim0U=-628,
present0=True,
#
#coef1=array([  1.67458724e-09,   5.13544267e-07,  -1.32243008e-02, 6.651060]),
coef1=array([ -5.08707297e-10,   2.40295893e-06,  -1.10023457e-02, 5.738]),
dlim1L=-359,
dlim1U=1243,
present1=True,
#  complete overlap
#coef2=array([  2.70722196e-06,  -1.33500943e-02,   6.4353966]),
coef2=array([  2.04088236e-06,  -1.10276823e-02,   6.758135]),
dlim2L=30,
dlim2U=1243,
present2=True,
#
#coef3=array([-0.01087432,  6.74158999]),
coef3=array([-0.01009546,  6.92509556]),
dlim3L=456,
dlim3U=1243,
present3=True,
#
sig0coef = None,
sig1coef=array([ -3.590306e-04,   3.33]),
sig2coef=array([5.8]),
sig3coef=None,
#  *2013 NEW*
extname = "GU256769652I",
e_dist12 = 10,
) # revised 2013-04-12
#===================================================================
cur02=dict(
obsdir='WR52',
obsid="00057902001",
ext=1,
anchor=[629.6,168.5],#array([ 631.3, 169.2]),
field=[-0.0631367326869, -0.120972531223],
anchor1_2530 = [647.1,155.7],
anchor2_2000 = [411,335],
angle=34.198412,
zmxangle=34.198412,
flatangle=36.35,
dist12=616.5,# 606.1,
ratio12=0.863,
dis1=[[-249.2,1908],[-205.5,2000],[-169.7,2105],[-95.3,2297],[-58.9,2405],\
      [-19.35,2530.],[0.9,2595],[57.3,2787],[94.2,2906],[136.1,3066.3],\
      [159.6,3137.4],[229.5,3409],[317.8,3718],[338.7,3818],[367.5,3943],\
      [463.9,4335],[492.5,4436],[545.3,4649]],
# subtract 2.35 pix from dis1 to linu up anchor
dis1corr=-2.35,
dis2=[[68.3,1721.1],[139.2,1814.3],[206.2,1908],[267.3,2000],[336.2,2105],\
      [452.5,2297],[511.6,2405],[582.5,2530],[618.0,2595],[669.2,2701],[729.3,2817]],
# overlap 1-2 above 2000A      
dispers1st=[  2.98314028e-10, -1.1390e-06, 1.55167e-03, 3.2275, 2600.058],
dispers2nd=[  5.67573556e-04,  1.90953, 2600.08], #array([7.355e-4,1.997,2597.8]),
ratio12_2000=None,
dist12_2000=268.08, #255.64,
disp2nd_2000=[ 5.67573556e-04, 1.511, 2000.0], #array([  7.35500000e-04,   1.48148653e+00,   1.98828040e+03]), 
#
#
coef0=None,
dlim0L=None,
dlim0U=None,
present0=False,
#
#coef1=array([ 9.73937625e-09,  -3.06535001e-05,  -1.51164505e-02, 4.073095]),
coef1=array([  1.40258121e-08,  -3.22748984e-05,  -1.57062672e-02, 1.3187861]),
dlim1L=-282,
dlim1U=+744,
present1=True,
#
coef2=array([-4.33414935e-05,   5.86902498e-03,  -1.41328757e+01]),
#coef2=array([ -4.16072066e-05,   4.17831266e-03,  -1.36628858e+01]),
dlim2L=35,
dlim2U=743,
present2=True,
#
coef3=array([ -9.45432166e-05,   9.52706536e-02,  -5.32881106e+01]) ,
#coef3=array([ -1.00806838e-04,   1.02751806e-01,  -6.25044780e+01]),
dlim3L=443,
dlim3U=744,
present3=True,
#
sig0coef = None,
sig1coef=array([  4.36285417e-07,  -2.08967811e-05,   3.23610095e+00]),
sig2coef=array([ 4.12331303]),
sig3coef=array([ 6.0]) ,
#  *2013 NEW*
extname = "gu256735030I",
e_dist12 = 2,
) # revised 2013-04-13
#===================================================================
cur04=dict(
obsdir='WR52',
obsid="00057904001",
ext=1,
anchor=[1288.2, 143.8],
field=[0.0439273850762, -0.125378743337],
anchor1_2530 = [1306.3,131.1],
angle= 34.87377,
zmxangle= 34.87377,
flatangle=35.1,
dist12=612,
ratio12=0.872,
dis1=[[-206,2000],[-171.1,2105],[-99.9,2297],[-64.5,2405],[-23.2,2530],[1.2,2595],\
      [57.5,2787],[93.1,2906],[136.2,3066.3],\
      [156.9,3137.4],[232.0,3409],[319.5,3718],[349.5,3818],[376.7,3943],\
      [478.7,4335],[508.6,4436],[560.3,4649],[857.2,5802]],
#
dis1corr=0.,
dis2=None,
#[[93.3,1721],[173,1814],[232.3,1908],[318.7,2000],[380.5,2105],[480.1,2297],\
# [532,2405],[610.5,2530],[810,?2906],[1270,?3724]], all averlap first order lines except at 1270pix
dispers1st=[ -4.75252197e-10, -1.56420580e-07, 1.15156689e-03, 3.1587, 2600.244],
dispers2nd=None,#array([  5.07415467e-05, 1.615, 2596.0]),
#
ratio12_2000=None,
dist12_2000=None,#241.0,
disp2nd_2000=None,#array([  5.07415467e-05,   1.57735106e+00,   2.00383916e+03]), 
#
coef0=None ,
dlim0L=None,
dlim0U=None,
present0=False,
#
coef1=array([  6.43274797e-09,  -1.18226914e-05,  -5.05960471e-03, 2.118450]),
dlim1L=-239,
dlim1U=+1280,
present1=True,
#
coef2=array([-0.01066998,  2.61601817]),
dlim2L=33,
dlim2U=1280,
present2=True,
#
coef3=array([-0.01066998,  2.61601817]),
dlim3L=443,
dlim3U=1280,
present3=True,
#
sig0coef = None,
sig1coef=array([  7.02988318e-04,   3.26598603e+00]),
sig2coef=array([4.2]),
sig3coef=array([6.]),
#  *2013 NEW*
extname = "gu256740791I",
e_dist12 = 0,
) # revised 2013-04-13
#===================================================================
cur23=dict(
obsdir='WR52',
obsid="00057923001",
ext=1,
anchor=[1040.3,204.2],#array([1042, 204.4]),
field=[0.00476658880196, -0.116650481629],
anchor1_2530=[1058.2,191.5],
angle=34.6939315,
zmxangle=34.6939315,
flatangle=35.3,
dist12=609.5,
ratio12=0.878,
dis1=[[-292.5,1814.3],[-245.3,1908],[-202.3,2000],[-166.0,2105],[-93.5,2297],\
     [-58.0,2405],[-16.7,2530],[6.9,2595],[64.3,2787],[99.7,2906],[142.9,3066],\
     [162.6,3137],[237.8,3409],[322.2,3718],[349.4,3818],[381,3943],\
     [478,4335],[508.4,4436.],[560.7,4649],[857.,5802]],
# adjust dis1 with -6.4 pix 
dis1corr=-6.4,    
dis2=None,
dispers1st=[ 3.13774e-10, -1.21289e-06, 1.4724e-03, 3.1966, 2600.15],
dispers2nd=[ 7.81485803e-04,   1.978,   2590.5],
#
ratio12_2000=None,
dist12_2000=264.8,
disp2nd_2000=[  7.81485803e-04,   1.42928020e+00,   1.99239505e+03], 
#
#
coef0=None,
dlim0L=None,
dlim0U=None,
present0=False,
#
coef1=array([  3.76962839e-09,  -1.22922325e-05,  -9.67998755e-03, 5.87880839e-01]),
dlim1L=-279,
dlim1U=1035,
present1=True,
#
coef2=array([ -6.33127213e-06,  -1.20625857e-02,   6.98890743e-01]),
dlim2L=30,
dlim2U=1035,
present2=True,
#
coef3=array([ -1.61665584e-05,   3.40421971e-03,  -5.15138796e+00]),
dlim3L=443,
dlim3U=1035,
present3=True,
#
sig0coef = None,
sig1coef=array([ -9.54079769e-07,   5.60333133e-04,   3.40536955e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256765515I",
e_dist12 = 0,
) # revised 2013-04-14
#===================================================================
cur21=dict(
obsdir='WR52',
obsid="00057921001",
ext=1,
anchor=[1820.1,1180.5],#array([1822.5,1179.6]),
field=[0.128655347199, 0.0360392892075],
anchor1_2530=[1838.,1167.7],
angle=35.595617,
zmxangle=35.595617,
flatangle=35.5,
dist12=650.2, #646.3,
ratio12=0.983,
dis1=[[-168.3,2105],[-96.3,2297],[-61.4,2405],[-19.7,2530],[2.8,2595],\
     [59.1,2787],[95.1,2906],[142.4,3066.],[160.4,3137],[237.,3409],\
     [320.4,3718],[350.,3818],[378.9,3943],[484.,4335],[508.8,4436],\
     [565.5,4649],[864.5,5802]],
# subtract 2.8 pix from dis1 
dis1corr=-2.8,    
dis2=[[83.6,1721.],[155.8,1814],[227.9,1908],[297.4,2000],[353.5,2105],\
      [476.0,2297],[533.7,2405],[608.3,2530],[652.8,2596],[762.4,2787],\
      [821.2,2906],[1091.1,3409],[1262.2,3718]],
dispers1st=[ 4.9335e-10, -1.3615e-06, 1.42188e-03, 3.1857, 2600.25],
dispers2nd=[  2.42182825e-04, 1.711, 2600.17],
# (dist12=653.4) 
# 3rd order: (dist12=652.4) array([-2.4681e-07,2.50951e-04,1.7732,2600.4])
# older array([  8.11418627e-04,   2.00874708e+00,   2.60796861e+03]),
# scaling dis = dis1*1.9+650, lamb = polyval(coef1,dis1)
ratio12_2000=None,
dist12_2000=283.5,#278.66,
disp2nd_2000=[ 2.42183e-04, 1.5316, 2000.52],
#array([  8.11418627e-04,   1.41212944e+00,   1.97914546e+03]), 
#
coef0=None,
dlim0L=None,
dlim0U=None,
present0=False,
#
coef1=array([ -2.53261346e-08,   4.47845184e-05,  -1.42643437e-02, -4.74910711e-01]),
dlim1L=-202,
dlim1U=1060,
present1=True,
#
coef2=array([  2.88509275e-05,  -4.00613883e-02,   2.36377682e+01]),
dlim2L=22,
dlim2U=1430,
present2=True,
#
coef3=array([  3.98991469e-05,  -7.92809918e-02,   5.84115121e+01]),
dlim3L=474,
dlim3U=1430,
present3=True,
#
sig0coef = None,
sig1coef=array([  5.02282900e-04,   3.21851425e+00]),
sig2coef=array([  2.42089354e-04,   4.02555684e+00]),
sig3coef=array([ 6.36145995]),
#  *2013 NEW*
extname = "gu256763893I",
e_dist12 = 0.7,
) # revised 2013-04-15
#===================================================================
cur16=dict(
obsdir='WR52',
obsid="00057916001",
ext=1,
anchor=[ 917.1, 454.3],
field=[-0.014884132513, -0.0776419671396],
anchor1_2530=[934.7,441.5],
angle=34.626588,
zmxangle=34.626588,
flatangle=35.5,
dist12=624.8,#623.7,
ratio12=0.899,
dis1=[[-346.5,1721.1],[-295.2,1814.3],[-246.4,1908],[-202.4,2000],[-163.4,2105],\
     [-92.2,2297],[-54.6,2405],[-15.7,2530],[8.4,2595],[66.1,2787],[98.2,2906],\
     [143.4,3066],[163.7,3137],[236.7,3409],[316.2,3718],[348.2,3818],[374.7,3943],\
     [555.6,4649],[842.5,5802]],
# subtract from dis1 7.7 pix to align anchor 
dis1corr=-7.7,    
dis2=[[180,1814],[239,1908],[466.5,2297],[522.6,2405],[593.4,2530],[796.4,2906]],
dispers1st=[ 5.010813e-10, -1.431e-06, 1.56989e-03, 3.230, 2600.2],
dispers2nd=None,#[2.076574e-04,1.83598,2600.31], -- dist12_2000 is not correct for the region (2013)
#array([ -1.284e-04,   1.9472,   2599.6]),
#
ratio12_2000=None,
dist12_2000=None,#292.6, #321.79, removed dist12_2000 and dispersion since does not fit with rest.
disp2nd_2000=None,#[2.076575e-04,1.6948,2000.25],#array([ -1.28400000e-04,   2.02472936e+00,   2.00002610e+03]), 
#
#
coef0=array([ -3.28134037e-02,  -4.61337148e+01,  -1.60993187e+04]),
dlim0L=-767,
dlim0U=-628,
present0=True,
#
coef1=array([ -1.34704540e-09,  -4.72530959e-06,  -1.47201609e-02, 5.46726954e-01]),
dlim1L=-280,
dlim1U=1065,
present1=True,
#            overlap -6 offset
coef2=array([ -6.98472808e-06,  -1.36095324e-02,  -5.61609248e+00]),
dlim2L=30,
dlim2U=1065,
present2=True,
#
coef3=array([ -1.12206880e-05,  -6.96155276e-03,  -8.10469178e+00]),
dlim3L=450,
dlim3U=1065,
present3=True,
#  based on lamda < 463.0 nm) : (increase in sigma afterwards )
sig0coef = None,
sig1coef=array([  1.83620154e-04,   3.31467964e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256753933I",
e_dist12 = 5,
) # revised 2013-04-15
#===================================================================
cur06=dict(
obsdir='WR52',
obsid="00057906001",
ext=1,
anchor=[1853.1,870.0],#array([1854.2, 871.4]),
field=[ 0.134159926708, -0.0140817821638 ],
anchor1_2530=[1870.9,857.3],
angle=35.51560,
zmxangle=35.51560,
flatangle=35.50,
dist12=659.2,#635.1,
ratio12=0.945,
dis1=[[-89.6,2297],[-58,2405],[-14.7,2530],[5.9,2595],[63.4,2787],\
     [100.8,2906],[146.,3066],[165.4,3137],[242.9,3409],[327.7,3718],\
     [359.5,3818],[388.8,3943],[570.8,4649],[878.3,5802]],
# offset dis1 subtract 5.95     
dis1corr=-5.95,
dis2=[[88.9,1721],[160.6,1814.3],[232.3,1908],[299,2000],\
     [360.8,2105],[489.3,2297],[620.3,2530],[837.8,2906],\
     [1102,3409],[1282,3718]],
dispers1st=[-8.83176e-10,6.1341e-07,6.6315e-04, 3.212, 2600.02],
dispers2nd=[ 2.44796e-04,1.6911, 2600.2],#array([7.0901e-4,1.9234,2590.1]),
#
ratio12_2000=None,
dist12_2000=290,#279.36,
disp2nd_2000=[2.4480e-04, 1.5074, 2000.1],#array([  7.09010000e-04,   1.41895647e+00,   1.99559848e+03]), 
#
#
coef0=None,
dlim0L=-767,
dlim0U=-628,
present0=False,
#
#coef1=array([  1.52030789e-05,  -7.30778821e-03,   1.30283732e+00]),
coef1=array([ -9.63575223e-09,   2.86353e-05,  -1.2420e-02, 1.340]),
dlim1L=-160,
dlim1U=1150,
present1=True,
#
#coef2=array([  3.67695938e-05,  -4.97219065e-02,   2.15179408e+01]),
coef2=array([  9.25233877e-06,  -1.44224685e-02,   1.58146100e+01]),
dlim2L=55,
dlim2U=1805,
present2=True,
#
coef3=array([  7.49075985e-05,  -1.28919235e-01,   7.12598551e+01]),
dlim3L=487,
dlim3U=1805,
present3=True,
#  this sig1 models the increase in width for lambda > 460nm:
sig0coef = None,
sig1coef=array([  4.44832428e-06,   3.36912419e-04,   3.14779625e+00]),
sig2coef=array([  4.70192703e-04,   3.81943153e+00]),
sig3coef=array([6.0]),
#  *2013 NEW*
extname = "gu256741875I",
e_dist12 = 2,
) # revised 2013-04-15
#===================================================================
cur22=dict(
obsdir='WR52',
obsid="00057922001",
ext=1,
anchor=[573.4,480.6],#array([ 571.9, 482.3]),
field=[-0.0723310346299, -0.0729594469731],
anchor1_2530=[590.9,468.1],
angle=34.1167,
zmxangle=34.1167,
flatangle=36.0,
dist12=604.,#604.5,
ratio12=0.873,
dis1=[[-343.8,1722],[-290.8,1814.3],[-243.2,1908],[-202.5,2000],[-163.6,2105],[-91.1,2297],\
      [-54,2405],[-13.4,2530],[7.8,2596],[64.3,2787],[96.1,2906],[142.1,3066],[161.5,3137],
      [232.2,3409],[312,3718],[341.7,3818],[373.5,3943],[543.2,4649],],
# shift dis1 by subtracting  7.2   
dis1corr=-7.2,  
dis2=[[135,1814],[204,1908],[262,2000],[446,2297],[509.6,2405],[573.2,2530],],
dispers1st=[  7.3045e-10, -1.41898e-06,1.5911e-03, 3.261265, 2600.184],
dispers2nd=[6.1117e-04,1.94256,2600.035], #array([  7.78031925e-04,   2.01359,   2597.69]),
#
ratio12_2000=None,
dist12_2000=264.5,#253.61,
disp2nd_2000=[6.11173e-04,1.5190,2000.31],#array([  7.78031925e-04,   1.46759141e+00,   1.98694381e+03]), 
#
#
coef0=array([ -3.73898072e-04,  -5.08747869e-01,  -1.48976203e+02]),
dlim0L=-786,
dlim0U=-582,
present0=True,
#
coef1=array([ -1.58199526e-09,  -1.57494264e-05,  -2.21808407e-02, 9.29756586e-01]),
dlim1L=-371,
dlim1U=692,
present1=True,
#
coef2=array([ -3.97530342e-05,   5.78477987e-06,  -9.08714894e+00]),
dlim2L=34,
dlim2U=680,
present2=True,
#
coef3=array([ -0.03094054, -12.07297117]),
dlim3L=434,
dlim3U=684,
present3=True,
#
sig0coef = None,
sig1coef=array([  2.65477889e-06,  -2.92350162e-04,   3.13985017e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256764973I",
e_dist12 = 2,
) # revised 2013-04-15
#===================================================================
cur08=dict(
obsdir='WR52',
obsid="00057908001",
ext=1,
anchor=[445.3,1305.0],#array([445,1305.6]),
field=[-0.0963809509685, 0.0603018133789],
anchor1_2530=[462.2,1292.9],
angle=33.979113,
zmxangle=33.979113,
flatangle=35.6,
dist12=None,
ratio12=None,
dis1=[[-343,1721],[-292.9,1814.3],[-245.2,1908],[-201.5,2000],[-164.5,2105],\
     [-92.2,2297],[-57.9,2405],[-19.5,2530],[-2.7,2596],[55.8,2787],[90.0,2906],\
     [131.9,3066],[151.1,3137],[218.0,3409],[299.2,3718]],
# subtract 1.8 pix from dis1  
dis1corr=-1.8,   
dispers1st=[ -1.40138e-09, -1.77316e-06, 1.9632489e-03, 3.3759,2600.09],
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
# minus first order present
#coefm1   = array([  8.30361876e-05,   1.98283659e-01,   1.56547428e+02])
#dlimm1L  = -1535
#dlimm1U  = -1391
#presentm1=True
#
coef0=array([  1.90094705e-04,   1.64251327e-01,   4.79414509e+01]),
dlim0L=-647,
dlim0U=-575,
present0=True,
#
coef1=array([ -2.98123636e-08,   2.05743852e-05,  -3.35452140e-02, 3.65714432e-01]),
dlim1L=-364,
dlim1U=337,
present1=True,
#  complete overlap:
coef2=array([  3.32365933e-05,  -4.19775359e-02,   8.75375074e-01]),
dlim2L=30,
dlim2U=335,
present2=True,
#
coef3=None,
dlim3L=None,
dlim3U=None,
present3=False,
#
sig0coef = None,
sig1coef=array([  4.48734468e-07,  -1.34603647e-03,   3.191425]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256746559I",
e_dist12 = 0,
) # revised 2013-04-15
#===================================================================
cur01=dict(
obsdir='WR52',
obsid="00057901001",
ext=1,
anchor=[307.5,467.9],#array([ 308.7, 466.1]),
field=[-0.11690201536, -0.0744620938528],
anchor1_2530=[324.4,455.3],
angle=33.538281,
zmxangle=33.538281,
flatangle=36.6,
dist12=582.2,#600,
ratio12=0.86,
dis1=[[-340.9,1721.1],[-288.3,1814.3],[-240.8,1908],[-195.2,2000],\
     [-160.0,2105],[-87.3,2297],[-52.7,2405],[-13.3,2530],[8.9,2596],\
     [64.2,2787],[96.8,2906],[141.0,3066.3],[160.3,3137],[229.4,3409],\
     [309,3718],[334,3818]],
# subtract 9.3 pix from dis1 
dis1corr=-9.3,    
dis2=[[66,1721],[200,1908],[260,2000]], # plotted netrates 1st & 2nd order together to find lines
dispers1st=[ -4.7336e-10, -1.5283e-06, 1.8517e-03, 3.3255495, 2600.799],
dispers2nd=[ 7.103657e-04, 2.04693, 2600.5], #array([7.74e-4,1.98,2600.3]),
#
ratio12_2000=None,
dist12_2000=260.0,
disp2nd_2000=[7.10366e-04, 1.5760, 2000.],#array([  7.74000000e-04,   1.42045740e+00,   1.98573257e+03]),
#
#
coef0= array([ -4.20056445e-05,  -3.80105540e-02,   2.58125943e+01]),
dlim0L=-767,
dlim0U=-556,
present0=True,
#
#coef1=array([ -3.58637922e-05,  -3.22124298e-02,  -8.22214695e-01]),
coef1=array([ -9.13946526e-09,  -3.89764344e-05,  -3.15512e-02, -1.00]),
dlim1L=-358,
dlim1U=+354,
present1=True,
#  
coef2=array([ -5.21060371e-05,  -8.47746943e-03,  -1.48003556e+01]),
dlim2L=+36,
dlim2U=+354,
present2=True,
#
coef3=None,
dlim3L=None,
dlim3U=None,
present3=False,
#
sig0coef = None,
sig1coef=array([3.25]),
sig2coef=array([4.16]),
sig3coef=None,
#  *2013 NEW*
extname = "gu256729814I",
e_dist12 = 20,
) # revised 2013-04-15
#===================================================================
cur12=dict(
obsdir='WR52',
obsid="00057912001",
ext=1,
anchor=[899.9,751.4],#array([901,752.9]),
field=[ -0.0180714720318, -0.0310439999923],
anchor1_2530=[917.4,738.8],
angle=34.64427914,
zmxangle=34.64427914,
flatangle=35.6 ,
dist12=616.2,
ratio12=None,
dispers1st=[ 5.91161e-10, -1.4740e-06, 1.5644277e-03, 3.239, 2600.43],
dispers2nd=[ 3.11830827e-04, 1.80138, 2600.399],
#
dis1=[[-350.8,1721.1],[-297.,1814.3],[-249.8,1908],[-207.5,2000],\
      [-166.6,2105],[-93.8,2297],[-59.3,2405],[-19.6,2530],[3.4,2596],\
      [60.9,2787],[96.0,2906],[140.7,3066],[157.8,3137],[233.3,3409],\
      [313,3718],[341.4,3818],[370,3943],[551.2,4649],[833.1,5802]],
# subtract 3.8 pix from dis1
dis1corr=-3.8,
dis2=[[445,2297],[583.5,2530],[783,2906],[1039,3409]],
# uncertaind due to overlap only somewhat good above 2300A.
ratio12_2000=None,
dist12_2000 = 265, 
disp2nd_2000=[ 3.1183e-04, 1.57998, 2000.21], 
#
coef0=array([-0.01787328,  2.43629045]),
dlim0L=-671,
dlim0U=580,
present0=True,
#
coef1=array([  2.3297e-09,  -3.1525e-06,  -1.79960e-02, 1.78545]),
#        array([-0.01787328,  2.43629045]),
dlim1L=-361,
dlim1U=1150,
present1=True,
#  complete overlap with first order
coef2=array([ -9.97995219e-07,  -1.79190236e-02,   1.6425]),
#        array([-0.01787328,  2.43629045]),
dlim2L=30,
dlim2U=1275,
present2=True,
#  complete overlap with first order
coef3=array([-0.01837, 1.568]),
dlim3L=470,
dlim3U=1275,
present3=True,
#
sig0coef = None,
sig1coef=array([  5.13706980e-06,  -6.16055737e-04,   3.18421373e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256747634I",
e_dist12 = 30,
) # revised 2013-04-16
#===================================================================
cur18=dict(
obsdir='WR52',
obsid="00057918001",
ext=1,
anchor=[928.9,1347.5],#array([929.7,1349]),
field=[ -0.0151417130615, 0.0646047743725],
anchor1_2530=[946.4,1335.2],
angle=34.6715,
zmxangle=34.6715,
flatangle=35.2,
dist12=602.2,# 590.8,
ratio12=0.904,
dis1=[[-343,1721.1],[-295.1,1814.3],[-246.5,1908],[-202.7,2000],[-163.3,2105],\
     [-91.8,2297],[-57.3,2405],[-17.8,2530],[5.5,2596],[61.1,2787],[97.4,2906],\
     [139,3066],[158.4,3137],[228.7,3409],[311.3,3718],[367.4,3943],\
     [536.5,4649],[810.9,5802]],
# subtract 5.8 pix from dis1 for anchor for proper distance
dis1corr=-5.8,
dis2=[[126,1814],[198,1908],[249.7,2000],[451,2297],[568.2,2530],[755.4,2906]],
dispers1st=[ 5.41124e-10, -1.49456e-06, 1.7228e-03,3.2772,2600.094],
dispers2nd=[ 7.11049e-04, 1.9693, 2600.73],#array([8.0856e-4,2.0182,2588.9]),
#
ratio12_2000=None,
dist12_2000=253.4,
disp2nd_2000=[ 7.11049e-04,1.4730,2000.051],#array([  8.08560000e-04,   1.46413725e+00,   1.99233532e+03]), 
#
coef0=array([ -2.64721798e-03,  -3.46671875e+00,  -1.12759301e+03]),
dlim0L=-672,
dlim0U=-570,
present0=True,
#
coef1=array([ -1.42034748e-08,   2.96951189e-05,  -2.51513207e-02, 1.441824]),
dlim1L=-374,
dlim1U=842,
present1=True,
#
coef2=array([  3.20586115e-05,  -4.98734058e-02,   1.47199039e+01]),
dlim2L=22,
dlim2U=842,
present2=True,
#
coef3=array([  8.80095234e-05,  -1.38944333e-01,   5.34797638e+01]),
dlim3L=428,
dlim3U=842,
present3=True,
#
sig0coef = None,
sig1coef=array([  3.97435858e-06,  -1.35937305e-03,   3.47172475e+00]),
sig2coef=array([ -5.49141614e-04,   4.508189]),
sig3coef=array([6.0]),
#  *2013 NEW*
extname = "gu256759184I",
e_dist12 = 5,
) # revised 2013-04-16
#===================================================================
cur20=dict(
obsdir='WR52',
obsid="00057920001",
ext=1,
anchor=[1425.8,1155.9],#array([1427.3,1158.7]),
field=[ 0.0660667482657, 0.0313906750153],
anchor1_2530=[1443.8,1143.3],
angle=35.195485,
zmxangle=35.195485,
flatangle=35.0,
dist12=637.3, # 634.6,
ratio12=0.957,
dis1=[[-297.,1814.3],[-248.7,1908],[-206.4,2000],[-167.2,2105],[-94.2,2297],\
     [-58.2,2405],[-15.8,2530],[2.1,2595],[66,2787],[99.5,2906],\
     [141.8,3066],[160.9,3137],[237.1,3409],[322.8,3718],\
     [380.5,3943],[556.4,4649],[849.1,5802]],
# subtract 5.7 pix from dis1  
dis1corr=-5.7,    
dis2=[[71.4,1721.2],[214,1908],[283,2000],[469.5,2297],[597.5,2530],[805,2906],[1060.3,3409]],
dispers1st=[ 2.76652e-10, -1.1694e-06, 1.51952e-03, 3.1814, 2600.16],
dispers2nd=[ 4.00946e-04, 1.781495, 2600.46], # array([6.834e-4,1.914,2596.4]),
#
ratio12_2000=None,
dist12_2000=275.6,#271.65,
disp2nd_2000=[ 4.00946e-04, 1.48688, 2000.061], # array([  6.83400000e-04,   1.41791877e+00,   1.99173862e+03]), 
#
#
coef0=array([ -2.02608642e-03,  -2.73943923e+00,  -9.26048055e+02]),
dlim0L=-684,
dlim0U=-594,
present0=True,
#
coef1=array([-1.80385346e-08,3.49957538e-05,-1.724964e-02, 1.03814]),
dlim1L=-365,
dlim1U=+1150,
present1=True,
#
coef2=array([  1.68807092e-05,  -3.60687976e-02,   1.81040577e+01]),
dlim2L=45,
dlim2U=1270,
present2=True,
#
coef3=array([  1.57888278e-05,  -4.45969965e-02,   3.78135497e+01]),
dlim3L=466,
dlim3U=1270,
present3=True,
#
sig0coef = None,
sig1coef=array([ -3.31889697e-04,   3.52018687e+00]),
sig2coef=array([  1.11315886e-03,   3.59534399e+00]),
sig3coef=array([  1.46780111e-03,   5.43691432e+00]),
#  *2013 NEW*
extname = "gu256764434I",
e_dist12 = 5,
) # revised 2013-04-16
#===================================================================
cur05=dict(
obsdir='WR52',
obsid="00057905001",
ext=1,
anchor=[1598.2,513.4],#array([1599.7, 512.8]),
field=[0.0942305046654, -0.0698902747508],
anchor1_2530=[1616.5,500.5],
angle=35.21699,
zmxangle=35.21699,
flatangle=35.2,
dist12=649.7,
ratio12=None,
dis1=[[-345,1721.1],[-244.5,1908],[-204.,2000],[-166.8,2105],[-91.6,2297],\
      [-57.4,2405],[-16.4,2530],[6,2596],[63,2787],[99.1,2906],\
      [145.7,3066],[163.8,3137],[241.4,3409],[328.5,3718],\
      [569.3,4649],[867.8,5802]],
# subtract 6.3 pix from dis1 
dis1corr=-6.3,    
dis2=[[830.5,2906],[1108,3409],[1280,3731]],
# overlap and problem with second order trackwidth preclude getting blue part 2nd order
dispers1st=[ 5.6895e-10, -1.09425e-06, 1.2213e-03, 3.1185, 2600.22],
dispers2nd=[ 1.323e-04, 1.7297, 2600.1],
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=None,
dlim0L=None,
dlim0U=None,
present0=False,
#
coef1=array([ -5.37545487e-09,   1.60276003e-05,  -1.18881157e-02, 7.42779975e-01]),
dlim1L=-357,
dlim1U=1150,
present1=True,
#  complete overlap
coef2=array([  2.98223703e-06,  -3.90315334e-03,  -1.57260468e-03]),
dlim2L=30, 
dlim2U=1597,
present2=True,
#   complete overlap
coef3=array([  1.80015288e-06,  -1.31820428e-03,  -1.26753187e+00]),
dlim3L=450,
dlim3U=1597,
present3=True,
#
sig0coef = None,
sig1coef=array([  3.30145068e-06,   3.95702875e-04,   3.25074312e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256741333I",
e_dist12 = 25,
) # revised 2013-04-16
#===================================================================
cur14=dict(
obsdir='WR52',
obsid="00057914001",
ext=1,
anchor=[1148.3,1093.8],#array([1147.8,1094.5]),
field=[ 0.0214946973095, 0.0219180713098],
anchor1_2530=[1165.9,1081.4],
angle=34.9232,
zmxangle=34.9232,
flatangle=35.2,
dist12=622.0,# 609.4,
ratio12=0.913,
dis1=[[-293.9,1814.3],[-244.9,1908],[-200.9,2000],[-163.7,2105],\
      [-90.3,2297],[-55.2,2405],[-16.3,2530],[5.9,2596],[63.1,2787],\
      [102.1,2906],[144,3066],[162,3137],[232.6,3409],[318,3718],\
      [374,3943],[553.3,4649],[839.3,5802]],
# subtract 7.5 pix from dis1 
dis1corr=-7.5,
dis2=[[221,1908],[264.4,2000],[461,2297],[518.5,2405],[592.5,2530],[791,2906],[1048,3409],[1204,3731]],
dispers1st=[6.1935e-10, -1.57235e-06, 1.5963e-03, 3.25267, 2600.2],
dispers2nd=[3.03537e-04, 1.798457, 2600.63],#array([6.639e-4,1.796,2566.4]),
#
ratio12_2000=None,
dist12_2000=274.3, # 256.41,
disp2nd_2000=[3.03537e-04, 1.582824,2000.1],# array([  6.63900000e-04,   1.32730293e+00,   2.01515624e+03]), 
#
coef0=array([  1.83235832e-03,   2.16425230e+00,   6.32897445e+02]),
dlim0L=-795,
dlim0U=-581,
present0=True,
#
coef1=array([ -9.88251047e-09,   2.24554324e-05,  -1.83995449e-02,4.48549229e-01]),
dlim1L=-361,
dlim1U=1137,
present1=True,
#
#coef2=array([  2.68157613e-05,  -4.04998325e-02,   9.51636712e+00]) 
coef2=array([  1.40957946e-05,  -2.59653047e-02,   8.17580547e+00]),
dlim2L=14,
dlim2U=1137,
present2=True,
#
coef3=array([ -2.33628874e-03,   6.96360246e+00]),
dlim3L=470,
dlim3U=1137,
present3=True,
#
sig0coef = None,
sig1coef=array([  3.22033222e-06,   3.93609180e-05,   3.20961327e+00]),
sig2coef=array([3.85]),
sig3coef=None,
#  *2013 NEW*
extname = "gu256752854I",
e_dist12 = 6,
) # revised 2013-04-16
#===================================================================
cur11=dict(
obsdir='WR52',
obsid="00057911001",
ext=1,
anchor=[860.1,1061.9],#array([859.9,1062.3]),
field=[ -0.0257018794074, 0.0179665293698],
anchor1_2530=[877.7,1049.4],
angle=34.645775,
zmxangle=34.645775,
flatangle=35.4,
dist12=607.5,
ratio12=None,
dis1=[[-359.3,1721.1],[-296.6,1814.3],[-247.8,1908],[-204.7,2000],\
      [-167,2105],[-93.7,2297],[-59.4,2405],[-18.4,2530],[3.2,2596],\
      [60.8,2787],[93.4,2906],[138,3066],[159,3137],[231,3409],\
      [315,3718],[369,3943],[539.5,4649],[818,5802]],
# subtract 4.5 pix from dis1 
dis1corr=-4.5,
dis2=[[441,2297],[574,2530],[767,2906]],
# none disp 2 lines are trusted - best guess
dispers1st=[ 4.75939e-10, -1.3953e-06, 1.66878e-03, 3.24645, 2600.65],
dispers2nd=[  6.02168e-04, 1.8777, 2600.484],
#
ratio12_2000=None,
dist12_2000=250.5,
disp2nd_2000=[6.02168e-04, 1.442365, 2000.38], 
#
# minus first order is present
#coefm1=array([  2.80122206e-04,   7.28032967e-01,   4.92134668e+02]),
#dlimm1L=-1325,
#dlimm1U=-1218,
#presentm1=True,
#
coef0=array([ -2.48131041e-03,  -3.25043057e+00,  -1.05325123e+03]),
dlim0L=-675,
dlim0U=-580,
present0=True,
#
#coef1=array([  1.06836198e-05,  -2.24731489e-02,   1.28028148e+00]),
coef1=array([ -1.74605e-08,   1.93936e-05, -2.06504e-02,   5.3042e-01]),
dlim1L=-371,
dlim1U=877,
present1=True,
#
#coef2=array([  1.18088202e-05,  -2.46847296e-02,   2.08077113e+00]),
coef2=array([  3.76525796e-06,  -2.08994309e-02,   5.99696499e+00]),
dlim2L=50,
dlim2U=877,
present2=True,
#
coef3=array([  1.18088202e-05,  -2.46847296e-02,   2.08077113e+00]),
#coef3=array([-0.01666105,  5.16591417]),
dlim3L=440,
dlim3U=877,
present3=True,
#
sig0coef = None,
sig1coef=array([  5.25179740e-06,  -1.07944880e-03,   2.82580754e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256747093I",
e_dist12 = 20,
) # revised 2013-04-16
#===================================================================
cur15=dict(
obsdir='WR52',
obsid="00057915001",
ext = 1,
anchor=[429.2,805.0], # array([ 427.8,  805.3]),
field=[-0.0971659263863, -0.0210218972834],
anchor1_2530=[446.7,793.3],
zmxangle= 33.800616479,  
angle = 36.1,
flatangle=36.1,
dist12=565.5,# 588,
ratio12=0.866,
dis1=[[-348,1721.1],[-294,1814.3],[-244.,1908],[-199,2000],[-164.5,2105],\
      [-92.8,2297],[-56,2405],[-18.5,2530],[4,2596],[62.,2787],[94,2906],\
      [138.4,3066],[154.,3137],[225,3409],[305,3718],[360,3943],],
# subtract 5.0 pix from dis1 to align anchor   
dis1corr=-5.0,
dis2=[[121,1814],[191,1908],[417.4,2297],[477,2405]],  
# dis2 determined pretty badly 
dispers1st=[-2.211159e-11, -1.48489e-06, 1.83258e-03, 3.32006, 2600.03],
dispers2nd=[8.598e-04, 2.141, 2600.45],#array([  4.12211766e-04,   1.90173349e+00,   2.59544237e+03]),
# 
ratio12_2000=None,dist12_2000=248.5,# 249.62,
disp2nd_2000=[8.5980e-04, 1.58745, 2000.15], # array([  4.12211766e-04,   1.62276678e+00,   1.99913585e+03]), 
#
# minus first order present
#coefm1=array([-2.7]),
#dlimm1L=-1361,
#dlimm1U=-1200,
#presentm1=True,
#
coef0=array([ -3.45705498e-04,  -4.46774871e-01,  -1.45944009e+02]),
dlim0L=-758,
dlim0U=-568,
present0=True,
#
coef1=array([  4.76744070e-09,  -1.29063460e-05,   5.74140373e-03,-1.49110]),
dlim1L=-366,
dlim1U=506,
present1=True,
#
coef2=array([  1.53482317e-04,  -2.30654062e+00]),
dlim2L=30,
dlim2U=506,
present2=True,
#
coef3=array([-4.5]),
dlim3L=470,
dlim3U=506,
present3=True,
#
sig0coef = None,
sig1coef=array([ -3.81971815e-04,   3.28397901e+00]),
sig2coef=None,
sig3coef=None,
##  *2013 NEW*
extname = "gu256753393I",
e_dist12 = 20,
) # revised 2013-04-17
#===================================================================
cur13=dict(
obsdir='WR52',
obsid="00057913001",
ext=1,
anchor=[1262.0,807.1],#array([1262.3, 810.2]),
field=[0.0408582431246, -0.0224466657854],
anchor1_2530=[1280,794.5],
angle=34.987517,
zmxangle=34.987517,
flatangle=35.1,
dist12=638.2,# 645.6,
ratio12=0.949,
dis1=[[-350,1721.1],[-295.4,1814.3],[-250.6,1908],[-206.7,2000],[-167.6,2105],\
     [-94.2,2297],[-58.0,2405],[-18.9,2530],[4.,2596],[64.0,2787],\
     [96.5,2906],[143.3,3066.3],[161.4,3137],[236.8,3409],[323,3718],\
     [379.8,3943],[560.5,4649],[851.8,5802]],
# subtract 5 pix from dis1 
dis1corr=-5.,
dis2=[[85,1721],[168,1814],[199,1908],[478,2297],[600,2530],[810,2906],[1065,3409],[1235,3731]],
# dis2 very uncertain, based below 800 on subtracting count rates with smaller and larger extraction width
dispers1st=[ 4.7954e-10, -1.32453e-06, 1.47499e-03, 3.19054, 2600.57],
dispers2nd=[ 3.0898e-04, 1.75523, 2600.20],#array([3.25197762e-04, 1.7589, 2602.1]),
#
ratio12_2000=None,dist12_2000=278, #275.49,
disp2nd_2000=[3.0898e-04, 1.5296, 2000.40], #array([  3.25197762e-04,   1.51818027e+00,   1.99565527e+03]), 
#
coef0=array([ -4.47837796e-04,  -7.08383244e-01,  -2.68632689e+02]),
dlim0L=-781,
dlim0U=-592,
present0=True,
#
coef1=array([ -4.997956e-09,   1.38095e-05,  -1.47428e-02, 2.4625]),
#array([  8.02051693e-06,  -1.42584375e-02,   2.86207786e+00]),
dlim1L=-362,
dlim1U=1150,
present1=True,
#
coef2=array([  8.02051693e-06,  -1.42584375e-02,   2.86207786e+00]),
dlim2L=30,
dlim2U=3000,
present2=True,
#
coef3=array([  -0.00802432,  3.58089]),
dlim3L=470,
dlim3U=3000,
present3=True,
#
sig0coef = None,
sig1coef=array([  2.27170072e-06,   4.11807140e-04,   3.41312313e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256748173I",
e_dist12 = 20,
) # revised 2013-04-17
#===================================================================
cur26=dict(
obsdir='WR52',
obsid="00057926001",
ext=1,
anchor=[545.5,1078.3], #array([546.9,1078.4]),
field=[-0.078008295378, 0.0222562223205],
anchor1_2530=[562.7,1066],
zmxangle=34.1,
angle=35.7,
flatangle=35.7,
dist12=None,
ratio12=None,
dis1=[[-345,1721.1],[-293,1814.3],[-245,1908],[-200.8,2000],[-162,2105],\
      [-90.7,2297],[-58,2405],[-17.3,2530],[5,2596],[61.5,2787],[94,2906],\
      [135.7,3066.3],[157.,3137],[226.8,3409],[305.8,3718],[359,3943],[527.6,4649],],
# subtract 5.2 pix from dis1 
dis1corr=-5.2,
#dis2=[[553,2530]],??? not a clear sign of 2nd order to find
dispers1st=[ 3.1629e-10, -1.3422e-06, 1.788e-03, 3.3093, 2600.13],
dispers2nd=None,
ratio12_2000=None,
dist12_2000=None,
disp2nd_2000=None, 
#
# minus first order present
#coefm1=array([0.]), 
#dlimm1L=-1440,
#dlimm1U=-1200,
#presentm1=True,
#
coef0=array([0.]),
dlim0L=-769,
dlim0U=-558,
present0=True,
#
#coef1=array([  2.77709271e-08,  -1.26743466e-05,  -6.55197135e-03, 1.73004396e+00]),
coef1=array([  9.52566075e-09,   2.30607428e-07,  -4.95328588e-03,  9.9513e-01]),
#coef1=array([  7.38534186e-07,  -4.21886366e-03,   7.53499058e-01]),
dlim1L=-357,
dlim1U=609,
present1=True,
#  complete overlap
#coef2=array([  3.54629302e-06,  -6.01360769e-03,   9.61692048e-01]),
coef2= array([  4.18695458e-06,  -3.94169892e-03,   7.86258747e-01]),
dlim2L=30,
dlim2U=609,
present2=True,
#  complete overlap
#coef3=array([ 0.00526316, -3.86315789]),
coef3=array([-0.00274994,  1.13335903]),
dlim3L=450, 
dlim3U=609,
present3=True,
#
sig0coef = None,
sig1coef=array([  2.61650312e-06,  -1.14686229e-03,   3.25534907e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu256776821I",
e_dist12 = 0,
) # revised 2013-04-17
#===================================================================
cur27=dict(
obsdir='WR52',
obsid="00056950004",
ext=1,
anchor=[1057.4,955.9],#array([1058.3,956.9]),
field=[ 0.00681607448454, 0.001390831412],
anchor1_2530=[1075,943.3],
angle= 34.81579,
zmxangle= 34.81579,
flatangle=35.6,
dist12=639.0,# 628.8,
ratio12=0.933,
dis1=[[-351,1721.1],[-301.8,1814.3],[-249,1908],[-209,2000],[-170,2105],\
      [-97,2297],[-60,2405],[-22.2,2530],[-0.2,2596],[56.2,2787],[94,2906],\
      [136,3066.3],[155.9,3137],[229.1,3409],[314.7,3718],\
      [369,3943],[549.8,4649],[832.3,5802]],
# subtract 0.8 from dis1 for anchor 
dis1corr=-0.8,
dis2 = [[138,1814.3],[263,2000],[363.5,2105],[482,2297],[582.4,2530],[1039.,3409],[1194.,3731]], 
# nearly complete overlap      
dispers1st=[ 6.8344e-10,-1.52084e-06, 1.5235e-03, 3.24217, 2600.293],
dispers2nd=[ 4.5557e-04, 1.80929, 2600.1], #array([ 4.2590e-04, 1.823,   2603.2]),
#
ratio12_2000=None,dist12_2000=274.6,#262.64,
disp2nd_2000=[ 4.556e-04, 1.4765, 2000.10],#array([  4.25900000e-04,   1.51110264e+00,   1.99278805e+03]), 
#
# used wrong anker position from ext=2
coef0=array([ -3.28803424e-03,  -4.34885525e+00,  -1.43079441e+03]),
dlim0L=-800,
dlim0U=-587,
present0=True,
#
#coef1=array([  5.86064578e-06,  -1.58924684e-02,   1.90796929e+00])
coef1=array([ -6.62555570e-09,   1.44269659e-05,  -1.87230909e-02,  1.9305]),
dlim1L=-373,
dlim1U=1150,
present1=True,
#  overlaps completely
coef2=array([  3.79215773e-06,  -1.44526461e-02,   5.188729]),
#coef2=array([  5.86064578e-06,  -1.58924684e-02,   1.90796929e+00]),
dlim2L=30,
dlim2U=1261,
present2=True,
#
#coef3=array([  5.86064578e-06,  -1.58924684e-02,   1.90796929e+00]),
coef3=array([ -0.00906866,  3.82922664]) ,
dlim3L=470,
dlim3U=1261,
present3=True,
#
sig0coef = None,
sig1coef=array([  2.60649917e-06,   2.91926922e-04,   3.39369885e+00]),
sig2coef=None,
sig3coef=None,
#  *2013 NEW*
extname = "gu230531544I",
e_dist12 = 15,
) # revised 2013-04-17
#===================================================================
cur28=dict(
obsdir='WR52',
obsid="00057927002",
# used this spectrum to make plot for first plus second orders in paper
ext=1,
anchor=[1479.3,1490.7], 
field=[ 0.07320, 0.08359],
anchor1_2530=[1497,1478],
angle= 35.6,
zmxangle= 35.6,
flatangle= 35.6,
dist12=618.2,
ratio12=1.,
dis1=[[-301.5,1814.],[-251.7,1908],[-214.5,2000],[-98.5,2297],[-63.1,2405],[-21.7,2530],\
      [0.1,2596],[56.0,2787],[91.8,2906],[136.6,3066.3],[154,3137],[230.1,3409],[305.5,3718],\
      [371.5,3943],[546.4,4649],[831.5,5802]],
dis1corr=0,
dis2 = [[141.5,1814.],[207.7,1909],[276.5,2000], [357,2105],[455.5,2297],], 
# nearly complete overlap      
dispers1st=[  6.47420302e-10, -1.56916357e-06, 1.58687572e-03, 3.243986,2600.038],
dispers2nd=[  1.19688407e-03, 2.23705,2600.13],
#
ratio12_2000=None,dist12_2000=280.6,
disp2nd_2000=[ 1.196884e-03, 1.460,2000.0], 
#
coef0=None,
dlim0L=-575,
dlim0U=-566,
present0=True,
#
coef1=array([ -3.23478066e-08,   5.83752270e-05,  -1.63860283e-02, 0.]),
dlim1L=-360,
dlim1U=950,
present1=True,
#  no overlap
coef2=array([  3.36814243e-05,  -4.56213128e-02, 26.03]),
dlim2L=30,
dlim2U=950,
present2=True,
#
coef3=array([ -2.09801e-02, 41.826 ]) ,
dlim3L=443,
dlim3U=950,
present3=True,
#
sig0coef = None,
sig1coef=array([  -1.67e-9,  3.06e-6,   3.465]),
sig2coef=None,
sig3coef=None,
#
extname = "gu380907203I",
e_dist12 = 0,
) # revised 2013-12-09
#===================================================================
cur29=dict(
obsdir='WR52',
obsid="00057927003",
ext=1,
anchor=[1463.3,1500.7], 
field=[ 0.071060, 0.085353],
anchor1_2530=[1481,1488],
angle= 35.6,
zmxangle= 35.6,
flatangle= 35.6,
dist12=625.3,
ratio12=1.,
dis1=[[-349.4,1730],[-302,1814.],[-251.6,1908],[-210.2,2000],[-98.2,2297],[-64.0,2405],[-23.6,2530],\
      [-2.4,2596],[61.1,2787],[90.4,2906],[135.8,3066.3], [156,3137],[228.6,3409],[305.2,3718],\
      [368.8,3943],[546.3,4649],[827.7,5802]],
dis1corr=0,
dis2 = [[73,1730],[108,1760],[148,1814],[216,1909],[281,2010],[341,2150],[459.5,2297],[523.5,2405],[592.5,2530],], 
#      
dispers1st=[  8.88043100e-10, -1.806399e-06, 1.616075e-03, 3.266550, 2600.01],
dispers2nd=[3.03984171e-04,   1.74711827,  2600.05],
#
ratio12_2000=None,dist12_2000=267.2,
disp2nd_2000=[ 3.03984171e-04, 1.52405468,2000.0], 
#
coef0=array([ 7.24724592e-04, 6.19816021e-01, 95.239847]),
dlim0L=-634,
dlim0U=-581,
present0=True,
#
coef1=array([ -3.87508e-08, 6.368835e-05, -1.47110e-02,0.]),
dlim1L=-367,
dlim1U=926,
present1=True,
#  
coef2=array([ 3.371039e-05, -4.060369e-02,24.949]),
dlim2L=30,
dlim2U=926,
present2=True,
#
coef3=array([ -1.1238e-02, 35.22 ]) ,
dlim3L=455,
dlim3U=926,
present3=True,
#
sig0coef = None,
sig1coef=array([  -1.67e-9,  3.06e-6,   3.465]),
sig2coef=None,
sig3coef=None,
#
extname = "gu380843846I",
e_dist12 = 0,
) # revised 2013-12-10
#===================================================================
cur30=dict(
obsdir='WR52',
obsid="00057928002",
ext=1,
anchor=[1512.3,1556.7], 
field=[ 0.0780864, 0.095343],
anchor1_2530=[1530,1544],
angle= 35.7,
zmxangle= 35.7,
flatangle= 35.7,
dist12=615.1,
ratio12=1.,
dis1=[[-252,1908],[-209,2000],[-96.3,2297],[-64,2405],[-23,2530],\
      [0.1,2596],[56.4,2787],[89.6,2906],[136,3066.3], [156,3137],[228.4,3409],[309,3718],\
      [368,3943],[545,4649],[829,5802]],
dis1corr=0,
dis2 = [[135,1814],[211,1909],[271,2010],[454,2297],[512,2405],[579,2530],], 
#      
dispers1st=[  5.87010809e-10,  -1.48595105e-06,   1.56448925e-03,  3.25276,   2600.0],
dispers2nd=[ 6.9145e-04,   1.972277,  2600.0],
#
ratio12_2000=None,dist12_2000=269,
disp2nd_2000=[ 6.91452310e-04, 1.493653, 2000.131], 
#
coef0=None,
dlim0L=-634,
dlim0U=-581,
present0=True,
#
coef1=array([ -3.47287610e-08,   6.16198013e-05,  -1.65382421e-02,0.000]),
dlim1L=-355,
dlim1U=834,
present1=True,
#  
coef2=array([  5.75084835e-05,  -5.97922150e-02,   29.32]),
dlim2L=30,
dlim2U=834,
present2=True,
#
coef3= array([ -3.46868137e-02,   53.63]),
dlim3L=463,
dlim3U=834,
present3=True,
#
sig0coef = None,
sig1coef=array([  -1.67e-9,  3.06e-6,   3.465]),
sig2coef=None,
sig3coef=None,
#
extname = "gu380617883I",
e_dist12 = 0,
) # revised 2013-12-11
#===================================================================
cur31=dict(
obsdir='WR52',
obsid="00057928004",
ext=1,
anchor=[1514.3,1539.6], 
field=[0.07929725, 0.0923272],
anchor1_2530=[1532,1526.9],
angle= 35.7,
zmxangle= 35.7,
flatangle= 35.7,
dist12=621.8,
ratio12=1.,
dis1=[ [-302.88, 1814.0], [-253.28, 1908.0],
 [-209.18, 2000.0], [-99.08, 2297.0], [-64.28, 2405.0],
 [-23.18, 2530.0], [1.52, 2596.0], [55.82, 2787.0],
 [89.92, 2906.0], [134.82, 3066.3], [156.52, 3137.0],
 [228.52, 3409.0], [302.22, 3718.0], [369.52, 3943.0],
 [543.02, 4649.0], [549.22, 4686.0], [827.22, 5802.0]],
dis1corr=0,
dis2 = [[136.4,1814],[208.7,1909],[277.4,2010],[452.6,2297],[512.3,2405],[583.7,2530],[778.8,2906],], 
#      
dispers1st=[  5.27585487e-10,  -1.51996615e-06,   1.64244182e-03,  3.2538646,   2600.0],
dispers2nd=[  5.10535506e-04, 1.876597, 2600.1],
#
ratio12_2000=None,dist12_2000=268,
disp2nd_2000=[ 5.105355e-04,1.51534, 2000.1], 
#
coef0=None,
dlim0L=-546,
dlim0U=-547,
present0=True,
#
coef1=array([ -4.39058816e-08,   7.05367413e-05,  -1.65758745e-02, 0.000]),
dlim1L=-375,
dlim1U=872,
present1=True,
#  
coef2=array([  5.75084835e-05,  -5.97922150e-02,   2.62786597e+01]),
dlim2L=30,
dlim2U=872,
present2=True,
#
coef3=  array([ -2.79635399e-02,   4.83838426e+01]), 
dlim3L=428,
dlim3U=872,
present3=True,
#
sig0coef = None,
sig1coef=array([  4.1157e-06,-7.1275e-04,3.1745]),
sig2coef=None,
sig3coef=None,
#
extname = "gu380907203I",
e_dist12 = 0,
) # revised 2013-12-11
#===================================================================
cur40=dict(              # GD153 offset observation for curvature 
obsdir='fluxcal2/GD153',
obsid="00055505002",
ext=1,
anchor=[1679.41, 1593.32],
field=[0,0],
angle= 35.42913,
zmxangle= 35.42913,
flatangle= 35.42913,
dist12=None,
ratio12=None,
dispers1st=None,
dispers2nd=None,
##
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
 
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -4.04692580e-08,   6.90820779e-05,  -2.38237710e-02, -1.24142574e+01]),
dlim1L=-366,
dlim1U=726,
present1=True,
#  
coef2=array([  4.74108897e-05,  -5.91917601e-02,   1.86938200e+01]),
dlim2L=21,
dlim2U=780,
present2=True,
#
coef3=array([ -2.69879350e-02,   3.60380074e+01]),
dlim3L=420,
dlim3U=790,
present3=True,
#
sig0coef = None,
sig1coef=3.0,
sig2coef=3.0,
sig3coef=3.0,
#  *2013 NEW*
extname = "",
)
#===================================================================
cur41=dict(              # WD0320-539 offset observation 
obsdir='fluxcal2/WD0320-539',
obsid="00054255002",
ext=1,
anchor=[1655.31, 1854.49],
field=[0,0],
angle= 35.40849,
zmxangle= 35.40849,
flatangle= 35.40849,
dist12=None,
ratio12=None,
dispers1st=None,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
# 
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -5.14588374e-08,   8.16284036e-05,  -2.42311638e-02, -1.33457211e+01]),
dlim1L=-333,
dlim1U=653,
present1=True,
#  
coef2=array([  6.24919617e-05,  -6.78884829e-02,   1.89543562e+01]),
dlim2L=33,
dlim2U=680,
present2=True,
#
coef3=array([ -3.29236735e-02,   4.03443396e+01]),
dlim3L=422,
dlim3U=705,
present3=True,
#
sig0coef = None,
sig1coef=array([ -2.421351e-09,   2.28656e-06,  -7.39636e-04,  3.5695  ]),
sig2coef=array([ -2.992549e-06,   3.20185e-03,   3.5738  ]),
sig3coef=array([  8.889489e-04,   4.5621  ]),
#  *2013 NEW*
extname = "",
)
#===================================================================
cur42=dict(              # WD1657+343 offset observation 
obsdir='fluxcal2/WD1657+343',
obsid='00055901002',
ext=1,
anchor=[ 1726.57, 1841.37],
field=[0,0],
angle=  35.47475,
zmxangle=  35.474753,
flatangle= 35.474753,
dist12=None,
ratio12=None,
dispers1st=None,
dispers2nd=None,
#
ratio12_2000=None,
dist12_2000=None,
disp2nd_2000=None, 
#
# 
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -5.58921545e-08,   7.58938879e-05,  -2.19389426e-02, -1.65166129e+01]),
dlim1L=-315,
dlim1U=663,
present1=True,
#  
coef2=array([  5.30732489e-05,  -6.26823020e-02,   1.57046197e+01]),
dlim2L=31,
dlim2U=666,
present2=True,
#
coef3=array([ -2.07649181e-02,   2.94420452e+01]),
dlim3L=422,
dlim3U=720,
present3=True,
#
sig0coef = None,
sig1coef=3.0,
sig2coef=3.0,
sig3coef=3.0,
#  *2013 NEW*
extname = "",
)
#===================================================================
cur43=dict(              # T Pyx + HST#2 offset observation 
obsdir='TPyx',
obsid='00032043001',
ext=1,
anchor=[1575,1530], # from measured line positions
field=[ 0.0878773305733, 0.0911833765739],
angle=  35.324,
zmxangle=  35.324,
flatangle= 35.324,
dist12=626,
#>>dist12=641.25,
ratio12=0.979,
#dis1 =[[1750.,-333.2], [1908,-251.] , [2143,-153.2],  [2323,-89.2],   [2800,60.7],  \
#       [4105,410.5],  [4341,475.6],  [4641,547.6], [4860,602.2],  [5018,638.2]], 
dis1 =[[-333.2,1750.], [-251.,1908.] , [-153.2,2143],  [-89.2,2323],   [60.7,2800],  \
       [410.5,4105],  [475.6,4105],  [547.6,4641], [602.2,4860],  [638.2,5018]], 
dis1corr=0.,       
dispers1st=[ 1.11860e-09, -1.65723e-06, 1.4347e-03,3.2515, 2600.56],
         #[4.512e-10,-1.3446e-6,1.5571e-3,3.2185,2600.18],
dispers2nd=[5.544e-4,1.869,2595.0], # array([5.6462e-4,1.889,2605.53]),
#
ratio12_2000=None,dist12_2000=267.64,
disp2nd_2000=[  5.54400000e-04,   1.47164999e+00,   1.99642167e+03], 
#
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#  from model fit
coef1=array([ -4.07105612e-08,   6.46616829e-05,  -2.06722917e-02,   0.00000000e+00]),
dlim1L=-315,
dlim1U=663,
present1=True,
#  from model fit
coef2=array([  4.18107903e-05,  -5.36217724e-02,   2.84265586e+01]),
dlim2L=31,
dlim2U=666,
present2=True,
#  from model fit
coef3=array([ -2.36556482e-02,   4.40569412e+01]),
dlim3L=422,
dlim3U=720,
present3=True,
#
sig0coef = None,
sig1coef=3.0,
sig2coef=3.0,
sig3coef=3.0,
##  *2013 NEW*
extname = "",
e_dist12 = 0,
) # updated anchor, disper1st, dis1 2013-04-26 NPMK
#===================================================================
cur44=dict(              # T Pyx + HST#2 offset observation 
obsdir='TPyx',
obsid='00032043002',
ext=1,
anchor=[ 1357.0,1328.4], #[1360,1314]
field=[ 0.0545008336015, 0.0584643134457],
angle=  35.160410,
zmxangle= 35.160410,
flatangle= 35.160410,
dist12=621.0,
ratio12=0.9264,
#dis1 = [[1750.,-337.], [1908,-255], [2143,-155], [2323,-88], [2800,60],\
#        [4105,407.8], [4340,474], [4641,544.2], [4861,600.7], [5018,637.3]],
dis1 = [[-337.,1750], [-255,1908], [-155,2143], [-88,2323], [60,2800],\
        [407.8,4105], [474,4340], [544.2,4641], [600.7,4861], [637.3,5018]],
dis1corr=0.,    
dispers1st=[8.13638258e-10,-1.577453e-06, 1.53182e-03, 3.2424, 2600.57],
   #[0, -7.04797622e-07,   1.47995951e-03,   3.136932,2599.78],
#dispers2nd=[9.58522412e-04,   2.14955132,   2.59917059e+03], unadjusted  
dispers2nd=None, 
# overlap second order starts to dominate at 2200 A (leaking first order line at ~4127, at d=~410pix) 
#
ratio12_2000=None,
#dist12_2000=280.6,
#disp2nd_2000=[  9.58522412e-04,   1.49698260e+00,   1.97852418e+03], # measured 2013
dist12_2000=None, 
disp2nd_2000=None,
#
#
coef0= array([ -2.84509955e-04,  -5.19721481e-01,  -2.27904947e+02]),
#array([-0.098,-67.9]),
dlim0L=-760,
dlim0U=-570,
present0=True,
#  from model+opt fit
coef1=array([ -2.05179598e-08,   4.62997128e-05,  -2.28646686e-02, -9.44946458]),
dlim1L=-386,
dlim1U=1228,
present1=True,
#  from model+opt fit
coef2=array([  2.47520764e-05,  -4.24969260e-02,   1.17423217e+01]),
dlim2L=16,
dlim2U=1228,
present2=True,
#  from model fit
coef3=array([ -2.03037328e-02,   2.16535184e+01]),
dlim3L=465,
dlim3U=1228,
present3=True,
#
#  fourth order present 
#
sig0coef=array([6.8]),
sig1coef=array([  4.61533055e-10 , -5.72189011e-07  , 2.16948729e-04 ,  4.12189]),
sig2coef=array([ -1.42555694e-06 ,  1.26413475e-03  , 5.027154]),
sig3coef=array([  4.43747900e-04 ,  6.101932]), 
##  *2013 NEW*
extname = "",
e_dist12 = 0,
) # updated anchor, disper1st, dis1 2013-04-26 NPMK
#===================================================================
cur45=dict(              # V5668 Sgr offset observation 2015-08-01
obsdir='/data/novae/NovaSgr2015b',
obsid='00033875007',
ext=1,
anchor=[1139.80,1676.48], # offset = [0,5.9]arcsec
field=[ 0.0009, 0.1041], 
angle=  35.5,
zmxangle= 35.5,
flatangle= 35.5,
dist12=603.55,
dist13=1300.,
ratio12=None,
dis1 = [[ -343.32,  1724.65], [ -324.22,  1750. ], [ -264.62,  1862.3 ],
       [ -243.52,  1908.  ],  [ -173.92,  2043. ], [ -145.42,  2151.15],
       [  -79.42,  2332.1 ],  [  -30.82,  2471. ], [   49.28,  2802.  ],
       [  352.48,  3890.2 ],  [  373.58,  3971.2], [  404.68,  4102.9 ],
       [  462.48,  4341.7 ],  [  585.38,  4862.7], [  617.48,  5007.  ]],
dis1corr=0.,    
dispers1st=[2.85851e-09, -2.63314e-06, 1.323960e-03, 3.415721, 2600.047],
dis2 = [[   58.08,    81.98,   168.18,   199.28,   300.38,   348.38,
          453.88,   531.98,   706.98],
       [ 1724.5 ,  1750.  ,  1862.3 ,  1909.  ,  2043.  ,  2151.15,
         2332.1 ,  2471.  ,  2802.  ]],
dispers2nd=[  6.07916740e-04,   1.94945337e+00,   2.60006977e+03], 
dis3 = [[  442.98,   469.68,   583.58,   621.58],
       [ 1724.5 ,  1750.  ,  1863.3 ,  1909.  ]],
dispers3rd= [ 1.0230485,   2600.00],
#
ratio12_2000=None,
dist12_2000=None, 
disp2nd_2000=None,
#
#   curve measurements
#
coef0data=array([[-794.00917632, -772.91426678, -754.6700207 , -735.85564192,
        -649.19547301, -630.38109423, -613.84724621, -592.75233667,
        -573.36782521, -562.53530409, -557.40410988, -557.40410988],
       [ 110.61603757,  108.3355068 ,  106.62510873,  102.6341799 ,
          94.08218955,   92.94192417,   91.2315261 ,   87.81072996,
          85.5301992 ,   83.24966843,   81.53927036,   81.53927036]]),
coef0= array([ -1.19569090e-04,  -2.78288237e-01,  -1.35571886e+02]),
dlim0L=-760,
dlim0U=-570,
present0=True,
coef1data=array([[-361.59744764, -347.32875108, -325.92570625, -267.95912648,
        -244.77249458, -203.74999197, -174.32080532, -146.67520574,
         -79.79069063,  -47.68612338,   -5.77182724,   35.25067536,
          67.35524262,  103.91877755,  167.23611852,  209.15041466,
         264.44161382,  317.94922591,  376.80759921,  406.23678586,
         460.63619148,  498.09151995,  539.11402255,  576.56935101,
         627.4015825 ,  650.58821441],
       [ 112.39496099,  112.39496099,  109.71958039,  107.93599332,
         107.04419978,  105.26061271,  104.36881918,  102.58523211,
         101.69343858,  100.80164504,   99.90985151,   99.01805797,
          99.01805797,   99.01805797,   99.90985151,   99.90985151,
          99.01805797,   99.90985151,   99.90985151,  100.80164504,
         101.69343858,  102.58523211,  102.58523211,  103.47702564,
         100.80164504,  100.80164504]]),
coef1=array([-4.74272660e-08,5.18575309e-05,-1.08896162e-02,0.00]),
dlim1L=-381,
dlim1U=721,
present1=True,
coef2data=array([[  13.29871304,   23.56110146,   38.38455141,   52.63786866,
          72.02238013,   82.28476855,  110.79140306,  138.72790488,
         164.38387594,  197.45157198,  227.66860456,  250.47391217,
         345.68607144,  451.16061913,  525.84800155,  704.29953359,
         706.58006435],
       [ 120.3082933 ,  119.73816061,  119.16802792,  118.59789523,
         116.88749716,  117.45762985,  115.17709909,  115.17709909,
         114.6069664 ,  113.46670102,  113.46670102,  112.32643564,
         110.61603757,  110.61603757,  110.04590488,  110.61603757,
         110.04590488]]),
coef2=array([  3.75185480e-05,  -4.02339359e-02,   2.02976244e+01]),
dlim2L=7,
dlim2U=717,
present2=True,
coef3data=array([[ 383.88496168,  410.68119812,  435.1969038 ,  450.02035375,
         465.41393638,  485.36858054,  495.06083628,  508.74402084,
         522.9973381 ,  537.82078804,  554.35463606,  576.01967829,
         595.40418975,  619.34976274,  635.88361076,  643.86546842],
       [ 126.57975289,  125.43948751,  126.0096202 ,  124.86935482,
         124.86935482,  124.29922213,  124.86935482,  124.86935482,
         125.43948751,  123.72908944,  122.58882406,  122.58882406,
         123.15895675,  122.58882406,  123.15895675,  123.15895675]]),
coef3=array([ -1.41226518e-02,   3.16074263e+01]),
dlim3L=381,
dlim3U=648,
present3=True,
#
#  fourth order present 
#
sig0coef=array([6.8]),
sig1coef=array([  4.61533055e-10 , -5.72189011e-07  , 2.16948729e-04 ,  4.12189]),
sig2coef=array([ -1.42555694e-06 ,  1.26413475e-03  , 5.027154]),
sig3coef=array([  4.43747900e-04 ,  6.101932]), 
##  *2013 NEW*
extname = "",
e_dist12 = 0,
) # updated anchor, disper1st, dis1 2013-04-26 NPMK
#===================================================================2013-04-10
#p041clines=[[2606],[2743],[2796],[2848],[3006],[3092],[3156],[3230],\
#[3578],[3744],[3833],[3934],[3966],[4060],[4095],[4300],[4868],[5183],[6566]],
#p041cbiglines=[[2796],[2848],[3009],[3230],[3578],[3744],[3833],[3934],[3966],[4060],[4095],[4300],]
##set55=dict(obsid='00056800006',wheelpos=160, ext=1, anchor=[1109.,1025.],dateobs='2006-06-17',exposure=1572,roll= None,groups=['aA'],name='P041C',nproc=0,wlshift=3,exclude=[[1600,3000],[5200,7000]],notes=''),
#set55=dict(
#obsdir='P041C',  
#obsid='00056800006',  
#ext=1, anchor=[1005.8,947.1],
#angle=  34.8, 
#zmxangle= 34.8,
#flatangle= 34.8, 
#extname = "gu172274004I",
#dis1=[[7,2606],[42.9,2743],[58,2796],[75,2848],[119,3006],[141,3092],\
#      [180,3230],[272,3578],[315,3744],[337.5,3833],[364,3934],[373,3966],\
#      [396,4060],[407,4095],[457,4300],[600,4868],[1002,6566]],
# subtract 3.75 pix to line up anchor
#dispers1st=[ 3.0651e-10, -6.47671e-07, 8.53302e-04, 3.4615, 2600.1],
#  Note P041C does not add any useful extension to the wavecal, since the S/N is too low above 4500A even at 1500s exposure
#  For G63-26 with long exposure times, the number of counts at wavelengths > 6000 is about 50, which is too low to see H-alpha.
#===================================================================
#set172=dict(obsid='00055610002',wheelpos=160, ext=1, anchor=[962.7,1377.9],dateobs='2012-11-23',exposure=969.,roll=138.,name='G63-26',nproc=0,wlshift=-70,exclude=[[1600,2200],[3700,7000]],notes='no lenticular'),            
#===================================================================
#set173=dict(obsid='00055611002',wheelpos=160, ext=1, anchor=[742.1,1258.0],dateobs='2012-11-24',exposure=1322.6,roll=138.,name='G63-26',nproc=0,wlshift=-22,exclude=[[1600,2200],[3700,7000]],notes=''),               
#===================================================================
#set168=dict(obsid='00055606002',wheelpos=160, ext=1, anchor=[1575.1,1564.5],dateobs='2012-11-19',exposure=1384.6,roll=138.0,name='G63-26',nproc=0,wlshift=-37,exclude=[[1600,2600],[3700,7000]],notes='contam. below spectrum?'),
#===================================================================
#set169=dict(obsid='00055607002',wheelpos=160, ext=1, anchor=[1221.5,1287.6],dateobs='2012-11-20',exposure=1277.3,roll=138.0,name='G63-26',nproc=0,wlshift=-22,exclude=[[1600,2600],[3700,7000]],notes=''),     
#dis1=[]  
#===================================================================
#set170=dict(obsid='00055608002',wheelpos=160, ext=1, anchor=[1559.2,1300.2],dateobs='2012-11-21',exposure=1224.2,roll=138.,name='G63-26',nproc=0,wlshift=-41,exclude=[[1600,2740],[3700,7000]],notes=''),              
#===================================================================
#set171=dict(obsid='00055609002',wheelpos=160, ext=1, anchor=[1546.8,727.9],dateobs='2012-11-22',exposure=1384.6,roll=138.,name='G63-26',nproc=0,wlshift=-30,exclude=[[1600,2730],[3700,7000]],notes=''),               
#===================================================================

#==============================================================================================
#
# WR121 uv nominal 00057500003/uvot/image/
#
#
#
#
#
#
#
#===================================================================
#
#  Next: the UV nominal grism 
#
#===================================================================
cur100=dict(
obsdir='WR52',
obsid="00056950007",
ext=1,
anchor=array([ 907.0,996.1]),
angle= 28.718,
zmxangle= 28.718,
flatangle=0.00,
dist12=683.47,
ratio12=0.96405,
dispers1st=[ 4.638e-10, -1.3245e-06, 1.362e-03, 3.203, 2600.67],
dis1=[(-341, 1725), (-294.6, 1814), (-248.7, 1910),
 (-207.1, 2008), (-98.8, 2296.0), (-63.4, 2405.0),
 (-23.5, 2530), (91.9, 2906), (160.9, 3137),
 (233.6, 3409), (314.4, 3718), (379, 3943),
 (415, 4083), (481, 4329), (563, 4649),
 (867.4, 5802)],
dis1corr=0,
dispers2nd= array([  3.17887263e-04,   1.74378991e+00,   2.60017520e+03]),
#
ratio12_2000=None,dist12_2000=313.1,
disp2nd_2000=array([  3.17887263e-04,   1.50829450e+00,   1.99787760e+03]), 
#
#
coef0=array([ -0.07306251, -49.52916596]),
dlim0L=-826,
dlim0U=-636,
present0=True,
#
coef1=array([ -1.60251884e-08,   2.87370471e-05,  -1.42229141e-02,0.400]),
dlim1L=-362,
dlim1U=1016,
present1=True,
#  
coef2=array([  1.28291304e-05,  -1.98587783e-02,   1.07024223e+01]),
dlim2L=60,
dlim2U=1016,
present2=True,
#
coef3= array([ -9.43449138e-03,   1.59246322e+01]),
dlim3L=470,
dlim3U=1016,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)

#===================================================================
cur101=dict(
obsdir='WR52',
obsid="00056950007",
ext=2,
anchor=array([ 917.1,992.0]), 
angle= 28.718,
zmxangle= 28.718,
flatangle=0.00,
dist12=683.90,
ratio12=0.9660,
dispers1st=[ 6.994e-10, -1.498e-06, 1.332e-03, 3.22096, 2600.1],
dis1=[[ -338.9,  1725. ], [ -295.6,  1814. ],  [ -247.6,  1910. ],
 [ -206. ,  2008. ], [  -97.4,  2296. ], [  -63.4,  2405. ],
 [  -21.8,  2530. ], [   92.6,  2906. ], [  158.7,  3137. ],
 [  235.1,  3409. ], [  311.5,  3718. ], [  378.7,  3943. ],
 [  419.7,  4083. ], [  482.7,  4329. ], [  563.7,  4649. ],
 [  864.1,  5802. ]],
dis1corr=0,
dispers2nd= array([ -3.29918856e-04,   1.40741237e+00,   2.60110212e+03]),
#
ratio12_2000=None,dist12_2000=291.1,
disp2nd_2000=array([ -3.29918856e-04,   1.66660687e+00,   1.99734088e+03]), 
#
#
coef0= array([ -0.09048493, -58.51779507]),
dlim0L=-810,
dlim0U=-620,
present0=True,
#
coef1=array([-1.70020236e-08,3.08486605e-05,-1.53504943e-02,2.8714]),
dlim1L=-350,
dlim1U=1030,
present1=True,
#  
coef2=array([  9.10254631e-06,  -1.55312617e-02,   1.22260816e+01]),
dlim2L=60,
dlim2U=1030,
present2=True,
#
coef3=array([ -1.40060658e-02,   2.18080962e+01]),
dlim3L=480,
dlim3U=1030,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur102=dict(
obsdir='WR1',
obsid="00037905003",
ext=1,
anchor=array([ 877.4,1009.0]),
angle= 28.71612,
zmxangle= 28.71612,
flatangle=0.00,
dist12= 682.449,
ratio12=0.9310,
dispers1st=None,
dis1=None,
dis1corr=None,
dispers2nd=array([ -2.93498611e-04,   1.44749290e+00,   2.60083419e+03]),
#
ratio12_2000=None,dist12_2000=296.2,
disp2nd_2000=array([ -2.93498611e-04,   1.67424400e+00,   1.99788648e+03]), 
#
#
coef0=array([ -0.0758324 , -51.17768958]),
dlim0L=-820,
dlim0U=-621,
present0=True,
#
coef1=array([-1.54322966e-08,2.82196185e-05,-1.40599624e-02,0.6542798]),
dlim1L=-380,
dlim1U= 987,
present1=True,
#  
coef2=array([  2.41121296e-06,  -9.52517583e-03,   8.1165]),
dlim2L=65,
dlim2U=987,
present2=True,
#  not measurable
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur103=dict(
obsdir='WR86',
obsid="00057011002",
ext=1,
anchor=array([ 947.2,1021.2]),#[]950.5,1020.9
angle= 28.75908,
zmxangle= 28.75908,
flatangle=0.00,
dist12= 684.34,
ratio12=0.970,
dispers1st=None,
dis1=[ (-251.1, 1910), (-97.6, 2296),
 (-62.3, 2405), (-24.1, 2530), (92.7, 2906),
 (231.5, 3409), (315.2, 3705), (416.3, 4064),
 (566.0, 4649), (834.5, 5696), (867.0, 5812.0)],
dis1corr=0,
dispers2nd=array([ -3.36804522e-04,   1.45301146e+00,   2.59872327e+03]),
#
ratio12_2000=None,dist12_2000=305.65,
disp2nd_2000=array([ -3.36804522e-04,   1.70810105e+00,   2.00018107e+03]), 
#
#
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -1.31495451e-08,2.63870178e-05,-1.35623775e-02,0.8676]),
dlim1L=-352,
dlim1U=1057,
present1=True,
# 
coef2=array([  2.33616973e-05,  -3.32248488e-02,   12.837]),
dlim2L=60,
dlim2U=1057,
present2=True,
# can't be measured but is faintly present
coef3=array([ -8.24292151e-03,   1.63906784e+01]),
dlim3L=510,
dlim3U=680,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur104=dict(
obsdir='WR86',
obsid="00057012002",
ext=1,
anchor=array([ 801.9,1073.9]),#[800.7, 1074.4]
angle= 28.73,
zmxangle= 28.73,
flatangle=0.00,
dist12=678.53,
ratio12=0.9550,
dispers1st=[  8.94805227e-10,  -1.70278899e-06,   1.40037682e-03,
         3.23317359e+00,   2.60064090e+03],
dis1=[ (-251.9, 1910), (-98.6, 2296),
 (-61.6, 2405), (-22.9, 2530), (93.5, 2906),
 (230.2, 3409), (311.4, 3705), (410.6, 4064),
 (562.3, 4649), (828.5, 5696), (857.7, 5812.0)],
dis1corr=0,
dispers2nd=array([ -3.36804522e-04,   1.45301146e+00,   2.59872327e+03]), # assumed 
# only 1810 and 1908 lines in second order found
#
ratio12_2000=None,dist12_2000=299.75,
disp2nd_2000=array([ -3.36804522e-04,   1.70816070e+00,   2.00002979e+03]), 
#
#
coef0=array([ -0.08171562, -53.6187171 ]),
dlim0L=-800,
dlim0U=-587,
present0=True,
#
coef1=array([ -2.54714322e-08,   3.61853924e-05,  -1.47199954e-02,.3535472]),
dlim1L=-360,
dlim1U= 894,
present1=True,
#  
coef2= array([  2.01775058e-05,  -2.67120487e-02,   1.13147188e+01]),
dlim2L= 55,
dlim2U=894,
present2=True,
#
coef3=array([ -1.26294728e-02,   1.65708462e+01]),
dlim3L=478,
dlim3U=894,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur105=dict(
obsdir='WR86',
obsid="00057013002",
ext=1,
anchor=array([ 895.,1118.]),
angle= 28.8168,
zmxangle= 28.8168 ,
flatangle=0.00,
dist12=680.93,
ratio12=0.948,
dispers1st=[ 5.9154e-10, -1.2915e-06, 1.3616e-03, 3.154, 2600.011],
dis1=[ (-251.46, 1910), (-97.86, 2296),
 (-65.06, 2405), (-23.16, 2530), (91.74, 2906),
 (232.14, 3409), (321.04, 3705), (419.24, 4064),
 (563.04, 4649), (827.04, 5696), (860.04, 5812.0)],
dis1corr=0,
dispers2nd=array([  4.99630241e-04,   1.84473513e+00,   2.59667009e+03]),
#
ratio12_2000=None,dist12_2000=321.1,
disp2nd_2000=array([  4.99630241e-04,   1.48520430e+00,   1.99761914e+03]), 
#
#
coef0=array([ -0.09905335, -69.40479469]),
dlim0L=-806,
dlim0U=-622,
present0=True,
#
coef1=array([ -1.61092e-08,3.09456e-05,-1.4648e-02,-0.0705]),
dlim1L=-360,
dlim1U=1000,
present1=True,
#  
coef2=array([  2.42082807e-05,  -3.49419540e-02,   14.48]),
dlim2L=53,
dlim2U=1000,
present2=True,
#
coef3=array([ -1.86511027e-02,   2.55813599e+01]),
dlim3L=482,
dlim3U=1002,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur106=dict(
obsdir='WR86',
obsid="00057014002",
ext=1,
anchor=[880.8,979.5],#array([ 879.5,978.9]),#[880.8,979.5]
angle=  28.6798,
zmxangle=  28.6798 ,
flatangle=0.00,
dist12= 682.8,
ratio12=0.941,
dispers1st=[8.650e-10, -1.753e-06,1.418e-03, 3.241,2600.23],
dis1=[ (-249.6, 1910), (-98.7, 2296),
 (-62.8, 2405), (-21.9, 2530), (92.5, 2906),
 (232.6, 3409), (309.5, 3705), (409.0, 4064),
 (564.0, 4649), (833.5, 5696), (865.7, 5812.0)],
dis1corr=0,
dispers2nd=array([  7.27061348e-04,   1.97392644e+00,   2.59249042e+03]),
##
ratio12_2000=None,dist12_2000=335.63,
disp2nd_2000= array([  7.27061348e-04,   1.46910865e+00,   1.99484300e+03]), 
#

coef0= array([ -0.10377563, -69.33138089]),
dlim0L=-808,
dlim0U=-628,
present0=True,
#
coef1=array([-1.064e-08,2.59454e-05,-1.649e-02,0.868]),
dlim1L=-368,
dlim1U= 984,
present1=True,
#  
coef2=array([  2.09205660e-05,  -2.70030974e-02,   11.07]),
dlim2L= 52,
dlim2U=985,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur107=dict(
obsdir='WR86',
obsid="00057018001",
ext=1,
anchor=[1049.8,767.4],#array([ 1046.5,767.9]),
angle= 28.584,
zmxangle= 28.584 ,
flatangle=0.00,
dist12=619.1,
ratio12=0.931,
dispers1st=[6.8925e-10,-1.4609e-06, 1.2573e-03, 3.2320, 2600.217],
dis1=[(-284.2, 1814), (-248.8, 1910), (-97.3, 2296),
 (-61.9, 2405), (-21.1, 2530), (92.8, 2906),
 (231.4, 3409), (313.4, 3705), (415.5, 4064),
 (569.0, 4649), (841.5, 5696), (876.1, 5812.0)],
dis1corr=0,
dispers2nd=array([  1.34283579e-03,   2.29146809e+00,   2.54063370e+03]),
#
ratio12_2000=None,dist12_2000=335.27,
disp2nd_2000= array([  1.34283579e-03,   1.52919680e+00,   1.99842609e+03]), 
#
#
coef0=array([ -0.08170072, -55.25480409]),
dlim0L=-811,
dlim0U=-634,
present0=True,
#
coef1=array([-1.22334299e-08,2.10907432e-05,-1.33601354e-02,-0.58]),
dlim1L=-360,
dlim1U=1150,
present1=True,
#  
coef2=array([  8.05150238e-06,  -1.51192810e-02,   6.98203344e+00]),
dlim2L=70,
dlim2U=1169,
present2=True,
# cant see the third order 
coef3=None,
dlim3L=470,
dlim3U=1169,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur108=dict(
obsdir='WR86',
obsid="00057018002",
ext=1,
anchor=[ 1068.75,770.58],#array([1068.0,769.1]),
angle= 28.61,
zmxangle= 28.61,
flatangle=0.00,
dist12= 692.0,
ratio12=0.987,
dispers1st=[ -2.923e-09, 1.907e-06, 1.2699e-03, 2.9235, 2600.07],
dis1=[(-287.85, 1814), (-248.35, 1910), (-98.05, 2296),
 (-62.05, 2405), (-22.55, 2530), (92.85, 2906), (234.15, 3409),
 (315.55, 3705), (414.65, 4064), (569.35, 4649), (941.35, 5696),
 (875.85, 5812.0)],
dis1corr=0,
dispers2nd=array([  6.66476323e-05,   1.55900045e+00,   2.60018963e+03]),
#
ratio12_2000=None,dist12_2000=300.29,
disp2nd_2000=array([  6.66476323e-05,   1.50678768e+00,   1.99974333e+03]), 
#
#
coef0=array([ -0.07182329, -47.24764628]),
dlim0L=-810,
dlim0U=-629,
present0=True,
#
coef1=array([ -8.57087748e-09,   1.97647275e-05,  -1.43065523e-02,0.621]),
dlim1L=-360,
dlim1U=1150,
present1=True,
#  overlaps mostly
coef2=array([  1.22622959e-05,  -1.98040713e-02,   7.062]),
dlim2L=70,
dlim2U=1230,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur109=dict(
obsdir='WR86',
obsid="00057015001",
ext=1,
anchor=array([1135.3,1167.0]), #[ 1138.8, 1166.5]
angle= 28.9626,
zmxangle= 28.9626,
flatangle=0.00,
dist12=687.42,
ratio12=0.926,
dispers1st=[8.348e-10,-1.623e-06, 1.349e-03, 3.2054, 2600.1],
dis1=[ (-250.5, 1910), (-97.9, 2296),
 (-61.6, 2405), (-23.8, 2530), (92.8, 2906),
 (231.7, 3409), (314.2, 3705), (419.2, 4064),
 (567.1, 4649), (836.5, 5696), (870.7, 5812.0)],
dis1corr=0,
dispers2nd=array([  2.28561533e-04,   1.76339446e+00,   2.60055232e+03]),
#
ratio12_2000=None,dist12_2000=329.43,
disp2nd_2000= array([  2.28561533e-04,   1.59974855e+00,   1.99856497e+03]), 
#
#
coef0=array([ -0.11772986, -85.66368875]),
dlim0L=-809,
dlim0U=-635,
present0=True,
#
coef1=array([-1.6792e-08,3.714349e-05,-1.750639e-02,0.031]),
dlim1L=-332,
dlim1U=1135,
present1=True,
#  
coef2= array([  1.53829569e-05,  -2.92911681e-02,   1.78743182e+01]),
dlim2L=70,
dlim2U=1260,
present2=True,
#  cannot see it = low s/n
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur110=dict(
obsdir='WR86',
obsid="00057015002",
ext=1,
anchor=[1132.8,1189.4],#array([1130.1,1190.8]),
angle=  28.997,
zmxangle=  28.997 ,
flatangle=0.00,
dist12=687.23,
ratio12=0.969,
dispers1st=[1.0672e-09, -1.950e-06, 1.4085e-03, 3.244, 2600.14],
dis1=[ (-250.05, 1910), (-97.955, 2296),
 (-61.855, 2405), (-22.95, 2530), (93.05, 2906),
 (229.655, 3409), (311.355, 3705), (414.155, 4064),
 (567.65, 4649), (836.35, 5696), (872.25, 5812.0)],
dis1corr=0,
dispers2nd=array([  4.63125225e-04,   1.81315862e+00,   2.60138221e+03]),
#
ratio12_2000=None,dist12_2000=317.73,
disp2nd_2000=array([  4.63125225e-04,   1.47090512e+00,   1.99464442e+03]), 
#
#
coef0=array([ -0.12781686, -92.68666524]),
dlim0L=-818,
dlim0U=-631,
present0=True,
#
coef1=array([-2.168e-08,4.24157e-05,-1.77538e-02,6.754e-01]),
dlim1L=-340,
dlim1U=1140,
present1=True,
#  no overlap
coef2= array([  1.64996098e-05,  -2.72432984e-02,   1.90428330e+01]),
dlim2L=70,
dlim2U=2670,
present2=True,
#
coef3=array([ -1.05214837e-02,   2.60471982e+01]),
dlim3L=489,
dlim3U=1272,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur111=dict(
obsdir='WR86',
obsid="00057016001",
ext=1,
anchor=[808.6,1155.6],#array([ 809.3, 1155.4]),
angle=  28.7913,
zmxangle=  28.7913,
flatangle=0.00,
dist12=677.05,
ratio12=0.965,
dispers1st=[  9.1137e-10,-1.6501e-06, 1.383e-03, 3.2055, 2600.40],
dis1=[(-252.7, 1910), (-99.7, 2296), (-61.7, 2405), (-22.4, 2530), 
(92.5, 2906), (231.3, 3409),(314.3, 3705), (416.5, 4064), 
(562.3, 4649), (829.3, 5696), (855.3, 5812.0)],
dis1corr=0,
dispers2nd=array([ -1.06576047e-03,   1.06189673e+00,   2.58966032e+03]),
#
ratio12_2000=None,dist12_2000=299.34,
disp2nd_2000= array([ -1.06576047e-03,   1.86699883e+00,   2.03652009e+03]), 
#
#
coef0=array([ -0.08383281, -56.58863706]),
dlim0L=-811,
dlim0U=-626,
present0=True,
#
coef1=array([-1.59405e-08,2.99045966e-05,-1.48291e-02,7.7195e-01]),
dlim1L=-358,
dlim1U= 895,
present1=True,
#  
coef2=array([2.42770704e-05,-3.30655320e-02,13.0574]),
dlim2L=61,
dlim2U=895,
present2=True,
#  bad S/N expected but not seen --
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur112=dict(
obsdir='WR86',
obsid="00057016002",
ext=1,
anchor=[727.8,1261.5],#array([724.3,1262.4]),
angle= 28.87,
zmxangle= 28.87,
flatangle=0.00,
dist12=672.2,
ratio12=0.957,
dispers1st=[8.6915e-10,-1.8549e-06, 1.5194e-03, 3.2500, 2600.44],
dis1=[(-252.35, 1910), (-99.15, 2296), (-61.95, 2405), (-22.65, 2530),
 (92.045, 2906), (230.35, 3409), (307.45, 3705), (405.95, 4064),
 (557.75, 4649)],
dis1corr=0,
dispers2nd=array([  1.54545716e-04,   1.71548779e+00,   2.60009220e+03]),
#
ratio12_2000=None,dist12_2000=310.26,
disp2nd_2000=array([  1.54545716e-04,   1.60361547e+00,   1.99943533e+03]), 
#
#
coef0=array([ -0.08622075, -60.4019745 ]),
dlim0L=-808,
dlim0U=-623,
present0=True,
#
coef1=array([-1.325227e-08,2.98579085e-05,-1.4901498e-02,0.0553]),
dlim1L=-334,
dlim1U= 815,
present1=True,
#  
coef2=array([ 1.6733e-05, -2.5043e-02, 13.837]),
dlim2L=50,
dlim2U= 815,
present2=True,
#
coef3=array([ -0.04469911,  40.8]),
dlim3L=490,
dlim3U=815,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur113=dict(
obsdir='WR86',
obsid="00057017001",
ext=1,
anchor=[666.9,781.4],#array([ 667.4, 779.9]),
angle= 28.31,
zmxangle= 28.31,
flatangle=0.00,
dist12= None, 
#dist12 = 678.4, One line ?
ratio12=None,
dispers1st=[  9.05909935e-10,  -1.81949560e-06,   1.46386905e-03,
         3.24617184e+00,   2.60009398e+03],
dis1=[(-250.7, 1910), (-99.2, 2296), (-61.7, 2405), (-22.7, 2530),
 (93.0, 2906), (230.2, 3409), (309.8, 3705), (407.8, 4064),
 (560.4, 4649)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0= array([ -0.04279914, -21.65044914]),
dlim0L=-800,
dlim0U=-620,
present0=True,
#
coef1=array([-8.19577e-09,1.01338e-05,-1.31486e-02,-0.0230]),
dlim1L=-360,
dlim1U= 730,
present1=True,
#  overlaps completely
coef2=None,
dlim2L=30,
dlim2U=2520,
present2=False,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur114=dict(
obsdir='WR86',
obsid="00057017002",
ext=1,
anchor=[646.0,852.3],#array([646, 851.1]),
angle=  28.3936,
zmxangle=  28.3936 ,
flatangle=0.00,
dist12=None,
#dist12=677.03,
ratio12=0.992,
dispers1st=[ 1.5412e-09, -2.295e-06, 1.4588e-03, 3.2864, 2600.],
dis1=[(-250.5, 1910), (-98.5, 2296), (-60.8, 2405),
 (-22.2, 2530), (92.2, 2906), (227.9, 3409),
 (309.7, 3705), (407.2, 4064), (560.4, 4649)],
dis1corr=0,
dispers2nd=None,
#  array([ -3.99545077e-04,   1.40238063e+00,   2.59924353e+03]),
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([ -0.04231248, -22.76222705]),
dlim0L=-806,
dlim0U=-626,
present0=True,
#
coef1=array([ -1.43531792e-08,1.70169756e-05,-1.42386930e-02,0.77]),
dlim1L=-360,
dlim1U= 728,
present1=True,
#  overlaps mostly - very uncertain fit
coef2=None,
#     array([  7.32099540e-06,  -1.39297392e-02,   5.37407070e+00]),
dlim2L=30,
dlim2U=730,
present2=False,
# nothing to be seen
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur115=dict(
obsdir='WR86',
obsid="00057019001",
ext=1,
anchor=[1165.5,683.6],#array([ 1162.1,685]),
angle= 28.60,
zmxangle= 28.60,
flatangle=0.00,
dist12=695.32,
ratio12=0.914,
dispers1st=[9.4373e-10,-1.7488e-06, 1.32544e-03,3.23957, 2601.43],
dis1=[(-248.3, 1910), (-98.1, 2296), (-62.1, 2405),
 (-23.0, 2530), (91.9, 2906), (233.1, 3409), (310.6, 3705),
 (414.5, 4064), (568.8, 4649), (838.5, 5696), (872.1, 5812.0)],
dis1corr=0,
dispers2nd=array([  6.87353700e-04,   1.96692731e+00,   2.59661313e+03]),
#
ratio12_2000=None,dist12_2000=346.14,
disp2nd_2000=array([  6.87353700e-04,   1.48690527e+00,   1.99360635e+03]), 
#
#
coef0=array([ -0.05769482, -37.81254432]),
dlim0L=-810,
dlim0U=-643,
present0=True,
#
coef1=array([-8.436364e-09,1.9837e-05,-1.4217e-02,0.608]),
dlim1L=-372,
dlim1U=1150,
present1=True,
#  
coef2=array([  2.76618447e-06,  -7.31446744e-03,   5.868]),
dlim2L=60,
dlim2U=1295,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur116=dict(
obsdir='WR86',
obsid="00057019002",
ext=1,
anchor=[1494.4,668.2],#array([ 1491.7, 668.5]),
angle= 28.74,
zmxangle= 28.74,
flatangle=0.00,
dist12=704.0,
ratio12=0.975,
dispers1st=[1.1728e-09, -2.097e-06, 1.3767e-03,3.27043, 2600.63],
dis1=[(-247.7, 1910),
 (-98.7, 2296),
 (-62.3, 2405),
 (-22.0, 2530),
 (92.5, 2906),
 (231.7, 3409),
 (312.3, 3705),
 (408.8, 4064),
 (572.0, 4649),
 (848.0, 5696),
 (879.0, 5812.0)],
dis1corr=0,
dispers2nd=array([  2.11531467e-04,   1.66246867e+00,   2.60041379e+03]),
#
ratio12_2000=None,dist12_2000=323.80,
disp2nd_2000=array([  2.11531467e-04,   1.50161971e+00,   1.99891899e+03]), 
#
#
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -1.89637297e-08,   3.14346716e-05,  -1.30395831e-02,-1.2730]),
dlim1L=-364,
dlim1U=1150,
present1=True,
#  
coef2=array([  1.32109771e-05,  -2.11183912e-02,   1.14628620e+01]),
dlim2L=70,
dlim2U=1678,
present2=True,
#  can't find it now
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur117=dict(
obsdir='WR86',
obsid="00057020001",
ext=1,
anchor=[1105.8,1081.9],#array([ 1102.3, 1082.3]),
angle= 28.91,
zmxangle= 28.91,
flatangle=0.00,
dist12=688.15,
ratio12=0.978,
dis=[(-251.3, 1910), (-98.7, 2296), (-62.2, 2405), (-22.7, 2530),
 (91.8, 2906), (233.3, 3409), (316.1, 3705), (415.5, 4064),
 (566.1, 4649), (835.3, 5696), (863.3, 5812.0)],
#dis1=[[-364,1642],[-303.6,1764],[-254.4,1909],[-101.1,2295.],[-66.3,2400],
#[-25.1,2525],[25.3,2693],[40.9,2729],[89.5,2900.],[111.6,2980],[174.1,3200],
#[230.4,3406],[313.1,3686],[565.5,4649.],[835.2,5687],[865.3,5801],[877.4,5860]],
#dispers1st=[ 2.1448e-10,-8.3315e-07, 1.1555e-03,3.1782, 2600.13],
dispers1st=[ 8.2192e-10,-1.5701e-06, 1.3429e-03,3.2019, 2600.40],
#dis1corr=2.65,
dis1corr=0,
# correct dis1 by adding 2.65 pix
dispers2nd=array([  2.72162520e-04,   1.68311008e+00,   2.60051070e+03]),
#
ratio12_2000=None,dist12_2000=306.75,
disp2nd_2000=array([  2.72162520e-04,   1.47550562e+00,   1.99816590e+03]), 
#
#
coef0=array([ -0.09679915, -69.82832267]),
dlim0L=-820,
dlim0U=-633,
present0=True,
#
coef1=array([ -1.460887e-08,3.219637e-05,-1.63888e-02,-0.0916]),
dlim1L=-356,
dlim1U=1132,
present1=True,
#  
coef2=array([  1.66648965e-05,  -2.94705731e-02, 16.97]),
dlim2L=60,
dlim2U=1240,
present2=True,
#  cannot find sufficient S/N for 3rd order
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur118=dict(
obsdir='WR86',
obsid="00057020002",
ext=1,
anchor=[1405.9,1216.0],#array([ 1403.5, 1217.1]),
angle= 29.13,
zmxangle= 29.13,
flatangle=0.00,
dist12=693.93,
ratio12=0.979,
dispers1st=[ 4.1032e-10,-1.124e-06, 1.257e-03, 3.1683, 2600.24],
dis1=[(-249.2, 1910), (-99.3, 2296), (-62.0, 2405), (-22.2, 2530),
 (92.3, 2906), (233.6, 3409), (320.0, 3705), (418.0, 4064),
 (569.0, 4649), (840.0, 5696), (874.5, 5812.0)],
dis1corr=0,
dispers2nd=array([  3.29618669e-04,   1.71389070e+00,   2.60158591e+03]),
#
ratio12_2000=None,dist12_2000=313.08,
disp2nd_2000= array([  3.29618669e-04,   1.46281620e+00,   1.99665196e+03]), 
#
#  just partly present
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([-1.87800016e-08,4.28986229e-05,-1.85800299e-02,7.55e-01]),
dlim1L=-386,
dlim1U=1200,
present1=True,
#  overlap?
coef2=array([  1.36663560e-05,  -2.40441881e-02,   2.16053336e+01]),
dlim2L=71,
dlim2U=1546,
present2=True,
#
coef3=array([  2.68774420e-05,  -6.47525302e-02,   5.41937586e+01]),
dlim3L=450,
dlim3U=1546,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur119=dict(
obsdir='WR86',
obsid="00057021001",
ext=1,
anchor=[769.3,1593.5],#array([ 758.5, 1598]),
angle= 29.27,
zmxangle= 29.27,
flatangle=0.00,
dist12= 665.04,
ratio12=0.953,
dispers1st=[1.0607e-09,-1.9456e-06,1.551e-03,3.2621, 2600.632],
dis1=[(-253.55, 1910), (-97.75, 2296), (-62.05, 2405),
 (-22.85, 2530), (91.85, 2906), (227.95, 3409),
 (306.35, 3705), (403.05, 4064), (554.35, 4649),
 (809.45, 5696), (843.15, 5812.0)],
dis1corr=0,
dispers2nd=array([  4.60966761e-04,   1.86117780e+00,   2.60172631e+03]),
#
ratio12_2000=None,dist12_2000=306.61,
disp2nd_2000=array([  4.60966761e-04,   1.53073130e+00,   1.99384923e+03]), 
#
#
coef0=array([ -0.12521468, -95.07272605]),
dlim0L=-816,
dlim0U=-624,
present0=True,
#
coef1=array([ -2.15703e-08,4.40807736e-05,-1.3277843e-02,-3.763248e-01]),
dlim1L=-378,
dlim1U=846,
present1=True,
#  no overlap
coef2= array([4.39254336e-05,-4.6140696e-02,2.2569e+01]),
dlim2L=21,
dlim2U=846,
present2=True,
#
coef3=array([ -1.16913635e-02,   3.26654376e+01]),
dlim3L=461,
dlim3U=846,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur120=dict(
obsdir='WR86',
obsid="00057021002",
ext=1,
anchor=[1028.1,1653.1],#array([ 1028.7, 1652.4]),
angle=  29.4797,
zmxangle=  29.4797,
flatangle=0.00,
dist12= 673.42,
ratio12=0.972,
dispers1st=[-3.039e-10,  -1.281e-06,   1.6161e-03, 3.2193, 2600.3],
dis1=[(-253.0, 1910), (-97.1, 2296), (-62.1, 2405), (-25.4, 2530),
 (90.0, 2906), (234.0, 3409), (304.9, 3705), (402.1, 4064),
 (558.4, 4649)],
dis1corr=0,
dispers2nd=array([  3.36997792e-04,   1.71503104e+00,   2.60040678e+03]),
#
ratio12_2000=None,dist12_2000=293.38,
disp2nd_2000=array([  3.36997792e-04,   1.45888300e+00,   1.99729313e+03]), 
#
#
coef0=array([  -0.15153889, -115.14109003]),
dlim0L=-808,
dlim0U=-615,
present0=True,
#
coef1=array([-3.08437e-08,5.744974e-05,-1.342668e-02,-1.261482e-01]),
dlim1L=-321,
dlim1U= 795,
present1=True,
# no overlap
coef2=array([  5.02717942e-05,  -5.65528751e-02,   2.90827624e+01]),
dlim2L=37,
dlim2U=795,
present2=True,
#
coef3=array([ -1.46530025e-02,   4.11871383e+01]),
dlim3L=471,
dlim3U=795,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur121=dict(
obsdir='WR86',
obsid="00057022001",
ext=1,
anchor=[559.1,1531.5],#array([ 540.4, 1538.9])
angle= 29.0921,
zmxangle= 29.0921 ,
flatangle=0.00,
dist12= 636.6,
ratio12=0.962,
dispers1st=[ 1.860e-09, -2.462e-06, 1.5570e-03,3.3120,2600.433],
dis1=[(-252.7, 1910), (-97.3, 2296), (-61.18, 2405), (-21.5, 2530),
 (90.0, 2906), (225.68, 3409), (306.88, 3705), (399.88, 4064),
 (548.8, 4649)],
dis1corr=0,
dispers2nd=array([  2.59483805e-04,   1.66355317e+00,   2.59698921e+03]),
# does not look right: array([ -3.88527569e-04,   1.42760141e+00,   2.59614350e+03]),
#
ratio12_2000=None,dist12_2000=255.76,
disp2nd_2000=array([  2.59483805e-04,   1.46591153e+00,   2.00108251e+03]), 
#
#
coef0=array([ -0.10199124, -77.32292979]),
dlim0L=-823,
dlim0U=-622,
present0=True,
#
coef1=array([ -2.49173692e-08,3.95730e-05,-1.2395e-02,-2.42610]),
dlim1L=-361,
dlim1U= 609,
present1=True,
#  
coef2=array([  3.27071114e-05,  -3.73800056e-02,   1.67320600e+01]),
dlim2L=14,
dlim2U= 609,
present2=True,
#
coef3=array([ -1.98735531e-02,   2.87561597e+01]),
dlim3L=439,
dlim3U=609,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur122=dict(
obsdir='WR86',
obsid="00057022002",
ext=1,
anchor=[704.,1495.6],#array([ 691.7,1499.9]),
angle= 29.10,
zmxangle= 29.10,
flatangle=0.00,
dist12= 665.4,
#ratio12=0.929,
#dispers1st=None,
#dispers2nd=array([  4.21764658e-04,   1.81527397e+00,   2.59657851e+03]),
ratio12=0.9602,
dispers1st=[ -1.351e-09, -5.815e-07,1.686e-03, 3.182, 2600.0],
dis1=[(-252.1, 1910), (-99.1, 2296), (-63.7, 2405), (-23.2, 2530),
 (90.9, 2906), (228.1, 3409), (311.8, 3705), (394.7, 4064), (553.1, 4649)],
dis1corr=0,
dispers2nd=array([  8.51862539e-05,   1.65569909e+00,   2.59994660e+03]),
#
ratio12_2000=None,dist12_2000=295.97,
disp2nd_2000=array([  8.51862539e-05,   1.59275761e+00,   1.99990062e+03]), 
#
#
coef0=array([ -0.10948465, -83.23821302]),
dlim0L=-820,
dlim0U=-615,
present0=True,
#
coef1=array([-2.9329657e-08,4.534893e-05,-1.2135395e-02,-2.3298778]),
dlim1L=-376,
dlim1U=774,
present1=True,
#  
coef2=array([  5.13926684e-05,  -5.30172549e-02,   1.87624734e+01]),
dlim2L=34,
dlim2U=774,
present2=True,
#
coef3=array([ -8.44161439e-06,  -2.10608438e-02,   3.30963121e+01]),
dlim3L=447,
dlim3U=774,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur123=dict(
obsdir='WR86',
obsid="00057023001",
ext=1,
anchor=[224.3,1679.2],#array([ 237.3,1672.7]),
angle= 29.04245,
zmxangle= 29.04245,
flatangle=0.00,
dist12= 580., #from matching plot spectral orders
ratio12=0.952,
dispers1st=[1.4537e-08,-1.292e-06,9.0488e-04, 3.301,2600.09],
dis1=[(-249.1, 1910), (-95.1, 2296), (-60.7, 2405),
 (-21.1, 2530), (90.4, 2906), (224.5, 3409)],
dis1corr=0,
dispers2nd=array([  8.95180891e-04,   2.20906502e+00,   2.60075366e+03]),
#
ratio12_2000=None,dist12_2000=255.89,
disp2nd_2000=array([  8.95180891e-04,   1.62878291e+00,   1.97880269e+03]), 
#
# 
coef0=array([ -0.10523993, -70.49457016]),
dlim0L=-800,
dlim0U=-587,
present0=True,
#
coef1=array([  5.84047579e-05,  -1.57798253e-02,   9.64470605e-01]),
dlim1L=-350,
dlim1U= 260,
present1=True,
#  overlaps completely
coef2=array([ -1.69429175e-02,   1.69815356e+01]),
dlim2L=40,
dlim2U=260,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur124=dict(
obsdir='WR86',
obsid="00057023002",
ext=1,
anchor=[497.8,1783.5],#array([ 489.2,1786.6]),
angle=  29.322,
zmxangle=  29.322,
flatangle=0.00,
dist12= 586.5 ,
#ratio12=0.9898,
ratio12=1.0366,
dispers1st=[-3.7627e-09, -8.1237e-07, 2.059e-03,3.237, 2600.01],
dis1=[(-253.5, 1910), (-98.4, 2296), (-62.8, 2405),
 (-24.9, 2530), (91.1, 2906), (223.8, 3409),
 (299.8, 3705), (396.8, 4064)],
dis1corr=0,
#dispers2nd=array([  9.48901939e-04,  -2.16327310e+00,   2.60273428e+03]),
dispers2nd=array([  7.64135773e-04,   2.02961161e+00,   2.60061565e+03]),
#
ratio12_2000=None,dist12_2000=237.37,
disp2nd_2000= array([  7.64135773e-04,   1.49604596e+00,   1.98515900e+03]), 
#
#
coef0=array([ -0.11967311, -91.31921996]),
dlim0L=-810,
dlim0U=-600,
present0=True,
#
coef1=array([  4.81961501e-05,  -1.21126895e-02,  -4.90553501e-01]),
dlim1L=-370,
dlim1U= 492,
present1=True,
#  no overlap Line 1703A @ 58pix , line 1900A @ 195pix from anker1 
coef2=array([  9.09705647e-05,  -5.89389487e-02,   2.28211781e+01]),
dlim2L=26,
dlim2U=492,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur125=dict(
obsdir='WR86',
obsid="00057024001",
ext=1,
anchor=array([ 62.1,1228]),
angle= 28.44837,
zmxangle= 28.44837,
flatangle=0.00,
#dist12=429.4,
dist12=None,
ratio12=None,
dispers1st=None,
dis1=None,
dis1corr=None,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([  1.53251991e-04,   1.67731837e-01,   5.12086411e+01]),
dlim0L=-805,
dlim0U=-595,
present0=True,
# 
coef1=array([ -1.121e-07,  -1.207e-05,  -1.40e-02, 0.290]),
dlim1L=-365,
dlim1U=  60,
present1=True,
#  overlaps completely
coef2=None,
dlim2L=30,
dlim2U=2520,
present2=False,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur126=dict(
obsdir='WR86',
obsid="00057024002",
ext=1,
anchor=[318.0,1345.3],#array([ 307.7,1349.6]),
angle= 28.7140,
zmxangle= 28.7140,
flatangle=0.00,
#dist12=619.64,   # vanilla zemax
dist12=None,   
ratio12=None,
# cannot find a good solution for r12 and dispersion
# dx = 74,125,209.4,284   wav = 1725,1810,1908,2010
dispers1st=[  3.554e-09,-2.044e-06,   1.311e-03,3.2687, 2600.057],
dis1=[(-250.4, 1910), (-97.1, 2296), (-60.9, 2405), (-22.9, 2530),
 (91.4, 2906), (230.4, 3409), (308.5, 3705)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([ -0.07527406, -51.64816618]),
dlim0L=-820,
dlim0U=-610,
present0=True,
#
coef1=array([  2.87428032e-05,  -1.28414187e-02,  -2.09569412e+00]),
dlim1L=-370,
dlim1U= 340,
present1=True,
#  overlaps completely
coef2=array([  2.39589824e-07,  -1.68094815e-02,   6.19478244e+00]),
dlim2L=30,
dlim2U=340,
present2=True,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur127=dict(
obsdir='WR86',
obsid="00057025001",
ext=1,
anchor=[84.3,705.7],#array([ 86,703.3]),
angle= 27.87841,
zmxangle= 27.87841,
flatangle=0.00,
#dist12= 655.31,
dist12= None,
ratio12=None,
dispers1st=[ 1.117e-08, -7.197e-06,-1.59e-03,2.966,2600.17],
dis1=[(-253.2, 1910), (-100.0, 2296), (-64.3, 2405), (-23.4, 2530)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([  1.28294182e-04,   1.80795469e-01,   7.96602336e+01]),
dlim0L=-800,
dlim0U=-590,
present0=True,
#
coef1= array([ -2.78006319e-05,  -2.19120639e-02,   3.54486204e-01]),
dlim1L=-373,
dlim1U=  78,
present1=True,
#  
coef2=None,
dlim2L=30,
dlim2U=2520,
present2=False,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur128=dict(
obsdir='WR86',
obsid="00057025002",
ext=1,
anchor=[386.2,819.7],#array([ 373.5,824.1]),
angle= 28.144,
zmxangle= 28.144,
flatangle=0.00,
#dist12= 667.13,
dist12= None,
ratio12=None,
dispers1st=[3.3668e-10,-1.9496e-06, 1.5599e-03,3.2814,2600.89],
dis1=[(-250.0, 1910), (-97.0, 2296), (-62.4, 2405), (-22.6, 2530),
 (90.5, 2906), (227.6, 3409), (308.2, 3705), (404.6, 4064)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([-0.02186922, -2.0]),
dlim0L=-820,
dlim0U=-624,
present0=True,
#
coef1=array([-0.02031744, -0.70771821]),
dlim1L=-370,
dlim1U= 410,
present1=True,
#  overlaps completely
coef2=None,
dlim2L=30,
dlim2U=2520,
present2=False,
#
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur129=dict(
obsdir='WR86',
obsid="00057026001",
ext=1,
anchor=[443.7,371.1],#array([ 444.2,368.5]),
angle= 27.724006068,
zmxangle= 27.724006068,
flatangle=0.00,
#dist12= 675.2,
dist12= None,
ratio12=None,
dispers1st=[-8.284e-10,-1.2147e-06,1.6213e-03,3.2226,2600.71],
dis1=[(-251.0, 1910), (-99.8, 2296), (-61.6, 2405), (-24.2, 2530),
 (92.1, 2906), (228.6, 3409), (308.8, 3705), (403.6, 4064)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([-0.02438574,  0.04223806]),
dlim0L=-687,
dlim0U=-613,
present0=True,
#  uv curvature spectrum at 1900: y offset -4.8 pix
coef1=array([ -1.78507014e-05,  -1.31612066e-02,  -9.89595759e-01]),
dlim1L=-360,
dlim1U= 487,
present1=True,
#  overlaps maybe  at pixel anker+ 97, 115, 127, 143, 230
coef2=None,
dlim2L=30,
dlim2U=487,
present2=False,
#
coef3=None,
dlim3L=470,
dlim3U=487,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
cur130=dict(
obsdir='WR86',
obsid="00057026002",
ext=1,
anchor=[705.0,409.4],#array([ 702.4,409.5]),
angle=  28.0260339,
zmxangle=  28.0260339,
flatangle=0.00,
#dist12=684.6,
dist12=None,
ratio12=None,
#dispers1st=array([5.35475166e-10,-1.40087781e-06,1.36460638e-03,3.20593,2600.23]),
dispers1st=[ 1.64386e-09,-2.173e-06,1.352e-03, 3.272, 2600.41],
dis1=[(-248.7, 1910), (-96.5, 2296), (-63.0, 2405),
 (-22.0, 2530), (91.7, 2906), (231.5, 3409),
 (311.8, 3705), (410.3, 4064), (563.1, 4649)],
dis1corr=0,
dispers2nd=None,
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([-0.01388531,  0.93784238]),
dlim0L=-810,
dlim0U=-632,
present0=True,
#  uv curvature estimate at 1900: spectrum [y]~2 pixels below main line zeroth-first 
coef1=array([-0.01315736, -0.7312849 ]),
dlim1L=-355,
dlim1U= 785,
present1=True,
#  overlaps completely
coef2=None,
dlim2L=30,
dlim2U=785,
present2=False,
#  overlaps 
coef3=None,
dlim3L=470,
dlim3U=785,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
#===================================================================

#===================================================================
cur216=dict(
obsdir='WR86',
obsid="00057019002",
ext=1,
anchor=array([ 1592., 400.]),
angle= 28.74,
zmxangle= 28.74,
flatangle=0.00,
dist12=704.0,
ratio12=0.975,
dispers1st=None,
dis1=None,
dis1corr=None,
dispers2nd=array([  2.11531467e-04,   1.66246867e+00,   2.60041379e+03]),
#
ratio12_2000=None,dist12_2000=232.8,
disp2nd_2000=array([  2.11531467e-04,   1.50161971e+00,   1.99891899e+03]), 
#
#
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -1.89637297e-08,   3.14346716e-05,  -1.30395831e-02,-1.2730]),
dlim1L=-364,
dlim1U=1150,
present1=True,
#  
coef2=array([  1.32109771e-05,  -2.11183912e-02,   1.14628620e+01]),
dlim2L=70,
dlim2U=1678,
present2=True,
#  can't find it now
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================

#===================================================================
cur218=dict(
obsdir='WR86',
obsid="00057020002",
ext=1,
anchor=array([ 1404., 1200.]),
angle= 29.13,
zmxangle= 29.13,
flatangle=0.00,
dist12=693.93,
ratio12=0.979,
dispers1st=None,
dis1=None,
dis1corr=None,
dispers2nd=array([  3.29618669e-04,   1.71389070e+00,   2.60158591e+03]),
#
ratio12_2000=None,dist12_2000=313.07,
disp2nd_2000=array([  3.29618669e-04,   1.46281620e+00,   1.99665196e+03]), 
#
#  just partly present
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([-1.87800016e-08,4.28986229e-05,-1.85800299e-02,7.55e-01]),
dlim1L=-386,
dlim1U=1200,
present1=True,
#  overlap?
coef2=array([  1.36663560e-05,  -2.40441881e-02,   2.16053336e+01]),
dlim2L=71,
dlim2U=1546,
present2=True,
#
coef3=array([  2.68774420e-05,  -6.47525302e-02,   5.41937586e+01]),
dlim3L=450,
dlim3U=1546,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================

#===================================================================
cur316=dict(
obsdir='WR86',
obsid="00057019002",
ext=1,
anchor=array([ 1492., 868.]),
angle= 28.74,
zmxangle= 28.74,
flatangle=0.00,
dist12=704.0,
ratio12=0.975,
dispers1st=None,
dis1=None,
dis1corr=None,
dispers2nd=array([  2.11531467e-04,   1.66246867e+00,   2.60041379e+03]),
#
ratio12_2000=None,dist12_2000=323.8,
disp2nd_2000= array([  2.11531467e-04,   1.50161971e+00,   1.99891899e+03]), 
#
#
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([ -1.89637297e-08,   3.14346716e-05,  -1.30395831e-02,-1.2730]),
dlim1L=-364,
dlim1U=1150,
present1=True,
#  
coef2=array([  1.32109771e-05,  -2.11183912e-02,   1.14628620e+01]),
dlim2L=70,
dlim2U=1678,
present2=True,
#  can't find it now
coef3=None,
dlim3L=470,
dlim3U=2520,
present3=False,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================

#===================================================================
cur318=dict(
obsdir='WR86',
obsid="00057020002",
ext=1,
anchor=array([ 1404., 1417.]),
angle= 29.13,
zmxangle= 29.13,
flatangle=0.00,
dist12=693.93,
ratio12=0.979,
dispers1st=None,
dis1corr=None,
dispers2nd=array([  3.29618669e-04,   1.71389070e+00,   2.60158591e+03]),
#
ratio12_2000=None,dist12_2000=313.07,
disp2nd_2000=array([  3.29618669e-04,   1.46281620e+00,   1.99665196e+03]), 
#
#  just partly present
coef0=None,
dlim0L=-800,
dlim0U=-587,
present0=False,
#
coef1=array([-1.87800016e-08,4.28986229e-05,-1.85800299e-02,7.55e-01]),
dlim1L=-386,
dlim1U=1200,
present1=True,
#  overlap?
coef2=array([  1.36663560e-05,  -2.40441881e-02,   2.16053336e+01]),
dlim2L=71,
dlim2U=1546,
present2=True,
#
coef3=array([  2.68774420e-05,  -6.47525302e-02,   5.41937586e+01]),
dlim3L=450,
dlim3U=1546,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================
#===================================================================
cur220=dict(
obsdir='WR86',
obsid="00057021002",
ext=1,
anchor=array([ 1228., 1452.]),
angle=  29.4797,
zmxangle=  29.4797,
flatangle=0.00,
dist12= 673.42,
ratio12=0.972,
dispers1st=None,
dis1corr=None,
dispers2nd=array([  3.36997792e-04,   1.71503104e+00,   2.60040678e+03]),
#
ratio12_2000=None,dist12_2000=None,
disp2nd_2000=None, 
#
#
coef0=array([  -0.15153889, -115.14109003]),
dlim0L=-808,
dlim0U=-615,
present0=True,
#
coef1=array([-3.08437e-08,5.744974e-05,-1.342668e-02,-1.261482e-01]),
dlim1L=-321,
dlim1U= 795,
present1=True,
# no overlap
coef2=array([  5.02717942e-05,  -5.65528751e-02,   2.90827624e+01]),
dlim2L=37,
dlim2U=795,
present2=True,
#
coef3=array([ -1.46530025e-02,   4.11871383e+01]),
dlim3L=471,
dlim3U=795,
present3=True,
#
sig0coef = None,
sig1coef=None,
sig2coef=None,
sig3coef=None,
# array([ -1.59887823e-03,   6.84074871e+00])
#  *2013 NEW*
extname = "",
e_dist12 = 0,
obslines = [],
v2_ratio12_200 = 0,
)
#===================================================================


#===================================================================
# visual nominal
#===================================================================
#vis1000=[
# wheelpos=1000
WR121_001=dict(obsdir='WR121',obsid='00057500001',ext=1,
anchor=[904.8,1074.6],angle=31.9,
dis1=[[ -235.7,  2907.0],
[ -167.2,  3260.0],[ -155.4,  3325.0],
[ -104.4,  3609.0],[  -84.4,  3717.0],[  -54.5,  3889.0],
[  -48.3,  3921.0],[  -23.6,  4070.0],[   11.3,  4267.0],
[   21.1,  4316.0],[   76.7,  4649.0],[   98.1,  4785.0],
[  110.1,  4859.0],[  137.6,  5016.0],[  157.6,  5130.0],
[  177.0,  5251.0],[  248.1,  5696.0],[  268.4,  5808.0],
[  277.2,  5880.0],[  390.0,  6580.0],[  416.4,  6740.0],],
dispers1st=[-8.18184e-07, 1.07544e-03,    5.7991, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_002=dict(obsdir='WR121',obsid='00057501002',ext=1,
#anchor=[49.0,1710.9],angle=31.87,
anchor=[32.0,1738.5],angle=31.87,
dis1=[[ -171.3,  3230.0],[ -117.9,  3511.0],[  -99.3,  3609.0],
[  -22.2,  4070.0],[   -9.6,  4131.0],[   10.2,  4267.0],],
dispers1st=[-5.63103e-07, 3.51588e-03,    6.2798, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_003 = dict(obsdir='WR121',obsid='00057501002',ext=2,
#anchor=[48.88,1710.92],angle=31.86,
anchor=[23,1750],angle=31.86,
dis1=[[ -244.4,  2840.0],[ -212.3,  2998.0],
[ -167.4,  3230.0],[ -158.2,  3268.0],[ -129.5,  3416.0],
[ -115.0,  3511.0],[  -98.6,  3609.0],[  -80.3,  3717.0],
[  -71.7,  3760.0],[  -21.2,  4070.0],[  -10.8,  4131.0],],
dispers1st=[-5.24064e-06, 1.52771e-03,    6.2410, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_004 = dict(obsdir='WR121',obsid='00057501006',ext=1,
anchor=[93.43,1737.17],angle=31.88,
dis1=[[ -246.6,  2840.0],[  -99.2,  3609.0],[  -51.9,  3889.0],
[  -44.8,  3921.0],[  -23.4,  4070.0],[  -11.5,  4131.0],
[   26.9,  4367.0],[   56.0,  4516.0],[   72.8,  4649.0],],
dispers1st=[-8.49041e-06, 2.13321e-05,    6.0368, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_005 = dict(obsdir='WR121',obsid='00057502006',ext=1,
anchor=[95.13,1252.13],angle=31.87,
dis1=[[ -102.3,  3609.0],
[  -20.4,  4070.0],[   10.5,  4267.0],[   74.7,  4649.0],],
dispers1st=[-6.09423e-05,-3.74739e-04,    6.3758, 4200.00],
dis1corr=   -0.00,
)

# wheelpos=1000
WR121_006 = dict(obsdir='WR121',obsid='00057503002',ext=1,
anchor=[440.36,1841.41],angle=32.0,
dis1=[
[ -101.4,  3609.0],[  -48.8,  3921.0],[  -21.3,  4070.0],
[   10.8,  4267.0],[   46.1,  4472.0],[   73.9,  4649.0],
[  133.4,  5016.0],[  150.3,  5130.0],[  237.9,  5696.0],
[  266.6,  5880.0],],
dispers1st=[-5.00669e-07, 1.47251e-03,    5.9536, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_007 = dict(obsdir='WR121',obsid='00057503002',ext=2,
anchor=[439.36,1841.41],angle=32.0,
dis1=[[ -102.6,  3609.0],[  -21.3,  4070.0],[   11.1,  4267.0],
[   74.6,  4649.0],[  151.2,  5130.0],[  185.0,  5321.0],
[  237.9,  5696.0],[  266.5,  5880.0],],
dispers1st=[ 1.79378e-06, 1.16011e-03,    5.8751, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_008 = dict(obsdir='WR121',obsid='00057503006',ext=1,
anchor=[471.27,1830.44],angle=32.02,
dis1=[[ -246.6,  2840.0],[ -131.7,  3441.0],[ -101.0,  3609.0],
[  -52.4,  3889.0],[  -47.2,  3921.0],[   11.7,  4267.0],
[   42.7,  4472.0],[   75.3,  4649.0],[  106.1,  4859.0],
[  132.7,  5016.0],[  150.7,  5130.0],[  170.6,  5251.0],
[  239.5,  5696.0],[  258.5,  5808.0],[  267.0,  5880.0],],
dispers1st=[-1.82927e-06, 1.47888e-03,    5.9944, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_009 = dict(obsdir='WR121',obsid='00057504002',ext=1,
anchor=[111.05,735.28],angle=31.73,
dis1=[[ -283.4,  2596.0],[ -236.0,  2840.0],[ -102.5,  3609.0],
[  -75.6,  3760.0],[  -52.6,  3889.0],[  -46.6,  3921.0],
[   11.4,  4267.0],[   20.2,  4316.0],[   45.0,  4472.0],
[   54.1,  4516.0],[   75.8,  4649.0],[   98.4,  4785.0],],
dispers1st=[-1.46413e-06, 3.88862e-04,    5.9032, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_010 = dict(obsdir='WR121',obsid='00057505002',ext=1,
anchor=[431.85,1463.7],angle=31.94,
dis1=[[ -102.3,  3609.0],[  -22.3,  4070.0],
[   11.0,  4267.0],[   74.9,  4649.0],[  153.5,  5130.0],
[  241.9,  5696.0],[  271.4,  5880.0],[  380.3,  6580.0],],
dispers1st=[-1.02729e-06, 1.30877e-03,    5.9108, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_011 = dict(obsdir='WR121',obsid='00057506002',ext=1,
anchor=[467.40,1079.07],angle=31.8,
dis1=[[ -104.5,  3609.0],[  -53.5,  3889.0],
[   10.6,  4267.0],[   75.4,  4649.0],[  153.8,  5130.0],
[  243.5,  5696.0],[  272.4,  5880.0],[  384.5,  6580.0],],
dispers1st=[-2.13097e-06, 1.68816e-03,    5.8570, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_012 = dict(obsdir='WR121',obsid='00057507002',ext=1,
anchor=[243.35,265.19],angle=31.35,
dis1=[[ -101.1,  3609.0],
[  -82.5,  3717.0],[  -74.0,  3760.0],[  -54.0,  3889.0],
[  -48.0,  3921.0],[  -21.5,  4070.0],[   11.5,  4267.0],
[   28.0,  4367.0],[   45.7,  4472.0],[   76.7,  4649.0],
[   99.0,  4785.0],[  108.4,  4859.0],[  121.5,  4924.0],
[  137.0,  5016.0],[  154.5,  5130.0],[  244.2,  5696.0],],
dispers1st=[ 2.04914e-06, 4.84006e-04,    5.8848, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_013 = dict(obsdir='WR121',obsid='00057508002',ext=1,
anchor=[950.52,1927.5],angle=32.19,
dis1=[[ -277.0,  2740.0],[ -103.0,  3609.0],[   11.6,  4267.0],
[   75.6,  4649.0],[  119.3,  4924.0],[  186.8,  5321.0],],
dispers1st=[-4.08851e-06, 1.22887e-03,    5.9243, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_014 = dict(obsdir='WR121',obsid='00057508002',ext=2,
anchor=[929.89,1933.33],angle=32.18,
dis1=[[ -102.1,  3609.0],
[  -53.6,  3889.0],[  -21.0,  4070.0],[   10.7,  4267.0],
[   74.7,  4649.0],[  135.2,  5016.0],[  153.7,  5130.0],],
dispers1st=[-3.42490e-06, 1.20780e-03,    5.9411, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_015 = dict(obsdir='WR121',obsid='00057508006',ext=1,
anchor=[983.86,1919.14],angle=32.21,
dis1=[[ -301.0,  2596.0],[ -252.1,  2840.0],
[ -248.0,  2860.0],[ -203.8,  3070.0],[ -182.2,  3169.0],
[ -170.2,  3230.0],[ -102.5,  3609.0],[  -84.9,  3717.0],
[  -48.5,  3921.0],[  -20.9,  4070.0],[   10.3,  4267.0],
[   20.3,  4316.0],[   28.4,  4367.0],[   75.1,  4649.0],
[  133.5,  5016.0],[  154.4,  5130.0],[  173.7,  5251.0],],
dispers1st=[-2.63859e-06, 1.26326e-03,    5.9205, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_016 = dict(obsdir='WR121',obsid='00057509002',ext=1,
anchor=[885.31,1519.64],angle=32.11,
dis1=[[ -270.7,  2740.0],
[ -202.2,  3070.0],[ -119.3,  3511.0],[ -104.1,  3609.0],
[  -54.7,  3889.0],[  -48.5,  3921.0],[  -21.8,  4070.0],
[   10.8,  4267.0],[   45.1,  4472.0],[   75.7,  4649.0],
[  135.8,  5016.0],[  154.1,  5130.0],[  245.3,  5696.0],
[  264.2,  5808.0],[  274.2,  5880.0],[  385.9,  6580.0],],
dispers1st=[-1.37180e-06, 1.29449e-03,    5.8660, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_017 = dict(obsdir='WR121',obsid='00057510006',ext=1,
anchor=[749.08,1197.76],angle=31.93,
dis1=[[ -235.0,  2907.0],[ -103.0,  3609.0],
[  -86.5,  3717.0],[  -22.8,  4070.0],[   20.6,  4316.0],
[   47.0,  4472.0],[   76.5,  4649.0],[  137.1,  5016.0],
[  155.3,  5130.0],[  174.1,  5251.0],[  246.2,  5696.0],
[  265.6,  5808.0],[  276.1,  5880.0],[  388.8,  6580.0],],
dispers1st=[-1.10259e-06, 1.18184e-03,    5.8327, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_018=dict(obsdir='WR121',obsid='00057511002',ext=1,
anchor=[530.84,586.50],angle=31.54,
dis1=[[ -103.4,  3609.0],[  -22.8,  4070.0],
[   11.6,  4267.0],[   76.0,  4649.0],[  185.7,  5321.0],
[  246.2,  5696.0],[  275.8,  5880.0],[  387.4,  6580.0],],
dispers1st=[-1.12792e-06, 1.20396e-03,    5.8468, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_019 = dict(obsdir='WR121',obsid='00057512006',ext=1,
anchor=[780.85,914.81],angle=31.79,
dis1=[[ -248.1,  2840.0],[ -220.2,  2984.0],[ -104.6,  3609.0],
[  -53.4,  3889.0],[  -47.3,  3921.0],[  -22.5,  4070.0],
[  -13.2,  4131.0],[   11.0,  4267.0],[   20.4,  4316.0],
[   46.5,  4472.0],[   76.9,  4649.0],[   98.5,  4785.0],
[  110.8,  4859.0],[  121.3,  4924.0],[  138.7,  5016.0],
[  156.3,  5130.0],[  187.1,  5321.0],[  247.5,  5696.0],
[  267.2,  5808.0],[  278.0,  5880.0],[  390.0,  6580.0],],
dispers1st=[-8.81492e-07, 1.09606e-03,    5.8078, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_020 = dict(obsdir='WR121',obsid='00057513002',ext=1,
anchor=[665.05,310.32],angle=31.38,
dis1=[[ -104.8,  3609.0],[  -22.6,  4070.0],
[   12.3,  4267.0],[   75.4,  4649.0],[  153.3,  5130.0],
[  246.2,  5696.0],[  276.0,  5880.0],[  390.6,  6580.0],],
dispers1st=[-2.41796e-06, 1.58595e-03,    5.8404, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_021 = dict(obsdir='WR121',obsid='00057514006',ext=1,
anchor=[1009.19,1238.02],angle=32.04,
dis1=[[ -104.3,  3609.0],[  -77.0,  3760.0],
[  -53.7,  3889.0],[  -48.0,  3921.0],[   10.6,  4267.0],
[   19.8,  4316.0],[   46.4,  4472.0],[   77.0,  4649.0],
[  121.1,  4924.0],[  137.4,  5016.0],[  156.6,  5130.0],
[  247.7,  5696.0],[  267.5,  5808.0],[  277.9,  5880.0],],
dispers1st=[-1.18523e-06, 1.13782e-03,    5.8132, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_022 = dict(obsdir='WR121',obsid='00057515002',ext=1,
anchor=[1265.39,1562.28],angle=32.31,
dis1=[[ -104.3,  3609.0],[  -22.6,  4070.0],
[   10.3,  4267.0],[   76.8,  4649.0],[  156.0,  5130.0],
[  246.6,  5696.0],[  275.6,  5880.0],[  388.8,  6580.0],],
dispers1st=[-1.51204e-06, 1.40338e-03,    5.8085, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_023 = dict(obsdir='WR121',obsid='00057516002',ext=1,
anchor=[967.57,628.86],angle=31.67,
dis1=[[ -105.0,  3609.0],[  -22.9,  4070.0],
[   11.8,  4267.0],[   75.7,  4649.0],[  156.1,  5130.0],
[  247.2,  5696.0],[  276.8,  5880.0],[  391.4,  6580.0],],
dispers1st=[-1.96397e-06, 1.49236e-03,    5.7995, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_024 = dict(obsdir='WR121',obsid='00057517002',ext=1,
anchor=[1119.95,914.57],angle=31.90,
dis1=[[ -103.8,  3609.0],[  -23.0,  4070.0],[   76.5,  4649.0],
[  247.8,  5696.0],[  278.0,  5880.0],[  392.5,  6580.0],],
dispers1st=[-1.51860e-06, 1.20347e-03,    5.8264, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_025 = dict(obsdir='WR121',obsid='00057518002',ext=1,
anchor=[1431.50,1974.46],angle=32.38,
dis1=[[ -102.8,  3609.0],[  -46.9,  3921.0],
[  -21.1,  4070.0],[   10.8,  4267.0],[   75.6,  4649.0],],
dispers1st=[-2.63405e-05, 3.73130e-04,    6.0640, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_026 = dict(obsdir='WR121',obsid='00057518002',ext=2,
anchor=[1431.52,1972.47],angle=32.37,
dis1=[[ -102.6,  3609.0],[  -53.9,  3889.0],
[  -22.8,  4070.0],[   11.6,  4267.0],[   74.8,  4649.0],],
dispers1st=[ 1.23431e-05, 1.67077e-03,    5.8051, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_027 = dict(obsdir='WR121',obsid='00057518008',ext=1,
anchor=[1455,1987],angle=32.6, # estimated anchor=[1462.9,1982.4]
dis1=[[ -248.1,  2840.0],[ -146.1,  3441.0],[ -104.1,  3609.0],
[  -86.1,  3717.0],[  -54.1,  3889.0],[  -49.1,  3921.0],
[  -24.1,  4070.0],[    7.9,  4267.0],[   73.9,  4649.0],],
dispers1st=[ 1.25894e-05, 4.38775e-03,    5.7790, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_028 = dict(obsdir='WR121',obsid='00057519002',ext=1,
anchor=[1381.28,1107.32],angle=32.13,
dis1=[[ -104.3,  3609.0],[  -23.0,  4070.0],
[   10.6,  4267.0],[   77.3,  4649.0],[  156.1,  5130.0],
[  247.8,  5696.0],[  277.6,  5880.0],[  392.5,  6580.0],],
dispers1st=[-1.70044e-06, 1.33733e-03,    5.8036, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_029 = dict(obsdir='WR121',obsid='00057520002',ext=1,
anchor=[1518.46,1416.12],angle=32.38,
dis1=[[ -270.7,  2699.0],[ -214.6,  2984.0],[ -103.6,  3609.0],
[  -85.4,  3717.0],[  -75.9,  3760.0],[  -54.3,  3889.0],
[  -49.0,  3921.0],[  -21.7,  4070.0],[   11.8,  4267.0],
[   19.7,  4316.0],[   29.5,  4367.0],[   47.1,  4472.0],
[   76.6,  4649.0],[  110.5,  4859.0],[  123.0,  4924.0],
[  155.9,  5130.0],[  186.9,  5321.0],[  248.0,  5696.0],
[  268.5,  5808.0],[  277.6,  5880.0],[  391.5,  6580.0],],
dispers1st=[-5.84327e-07, 8.61951e-04,    5.8347, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_030 = dict(obsdir='WR121',obsid='00057520008',ext=1,
anchor=[1556,1510],angle=32.4,  # anchor estimate 1560.7, 1506.1
dis1=[[ -104.1,  3616.0],[  -53.1,  3913.0],[  -23.1,  4080.0],
[   19.9,  4278.0],[   74.9,  4668.0],[  156.9,  5147.0],
[  246.9,  5709.0],[  266.9,  5834.0],[  277.9,  5891.0],],
dispers1st=[-4.35482e-06, 2.22561e-03,    5.8205, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_031 = dict(obsdir='WR121',obsid='00057521002',ext=1,
anchor=[1208.14,332.02],angle=31.54,
dis1=[[ -103.0,  3609.0],
[  -53.0,  3889.0],[  -23.9,  4070.0],[   10.4,  4267.0],
[   20.7,  4316.0],[   76.7,  4649.0],[  155.7,  5130.0],
[  247.6,  5696.0],[  277.8,  5880.0],[  392.4,  6580.0],],
dispers1st=[-1.29353e-06, 1.06426e-03,    5.8486, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_032 = dict(obsdir='WR121',obsid='00057522002',ext=1,
anchor=[1367.21,666.63],angle=31.84,
dis1=[[ -104.9,  3609.0],[  -54.7,  3889.0],[  -22.0,  4070.0],
[   10.7,  4267.0],[   76.7,  4649.0],[  156.3,  5130.0],
[  248.4,  5696.0],[  277.3,  5880.0],[  391.2,  6580.0],],
dispers1st=[-1.49020e-06, 1.34857e-03,    5.7863, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_033 = dict(obsdir='WR121',obsid='00057523002',ext=1,
anchor=[1627.74,900.36],angle=32.12,
dis1=[[ -101.6,  3609.0],[  -21.3,  4070.0],
[    9.9,  4267.0],[   75.8,  4649.0],[  155.7,  5130.0],
[  247.8,  5696.0],[  278.1,  5880.0],[  393.7,  6580.0],],
dispers1st=[-1.03970e-06, 7.76298e-04,    5.9020, 4200.00],
dis1corr=    0.00,
)

# wheelpos=1000
WR121_034 = dict(obsdir='WR121',obsid='00057524002',ext=1,
anchor=[1682.85,419.79],angle=31.8,
dis1=[[  -83.6,  3717.0],
[  -48.5,  3921.0],[  -22.0,  4070.0],[    9.8,  4267.0],
[   77.2,  4649.0],[  153.4,  5130.0],[  187.2,  5321.0],
[  245.2,  5696.0],[  276.0,  5880.0],[  388.9,  6580.0],],
dispers1st=[-1.42324e-06, 1.23159e-03,    5.8584, 4200.00],
dis1corr=    0.00,
)

#] # end vis1000
#===================================================================
#visual clocked
#===================================================================
#vis0955=[
# wheelpos=955
WR121_101 = dict(obsdir='WR121',obsid='00057500002',ext=1,
anchor=[802.02,977.72],angle=38.81,
dis1=[[ -100.9,  3609.0],
[  -84.4,  3717.0],[  -51.6,  3889.0],[  -46.9,  3921.0],
[  -22.2,  4070.0],[   10.4,  4267.0],[   20.2,  4316.0],
[   44.1,  4472.0],[   74.4,  4649.0],[  132.8,  5016.0],
[  150.6,  5130.0],[  197.3,  5395.0],[  237.5,  5696.0],
[  256.6,  5808.0],[  265.3,  5880.0],[  371.2,  6580.0],],
dispers1st=[-1.21918e-07, 1.31102e-03, 5.9483, 4200.00],
dis1corr=0,
)

# wheelpos=955
WR121_102 = dict(obsdir='WR121',obsid='00057501004',ext=1,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_103 = dict(obsdir='WR121',obsid='00057502002',ext=1,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_104 = dict(obsdir='WR121',obsid='00057502002',ext=2,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_105 = dict(obsdir='WR121',obsid='00057502004',ext=1,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_106 = dict(obsdir='WR121',obsid='00057502004',ext=2,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_107 = dict(obsdir='WR121',obsid='00057503004',ext=1,
anchor=[0,0],angle=0, # dark region
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR121_108 = dict(obsdir='WR121',obsid='00057504004',ext=1,
anchor=[269.4,652.15],angle=38.85,
dis1=[[ -182.6,  3070.0],[  -95.9,  3609.0],[  -50.6,  3889.0],
[  -21.9,  4070.0],[   72.4,  4649.0],[  145.9,  5130.0],],
dispers1st=[ 4.85737e-06, 7.25520e-04,    6.1604, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_109 = dict(obsdir='WR121',obsid='00057505004',ext=1,
anchor=[568.7,1377],angle=39.2,
dis1=[[ -206.8,  2998.0],[ -122.1,  3511.0],
[  -98.7,  3609.0],[  -20.7,  4070.0],[   10.5,  4267.0],
[   19.7,  4316.0],[   70.8,  4649.0],[  142.6,  5130.0],],
dispers1st=[ 7.55072e-06, 2.60166e-03,    6.0127, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_110 = dict(obsdir='WR121',obsid='00057506004',ext=1,
anchor=[593.44,1004.],angle=38.93,
dis1=[[ -299.3,  2511.0],[ -293.6,  2530.0],
[ -206.9,  2998.0],[ -193.7,  3070.0],[  -99.6,  3609.0],
[  -82.8,  3717.0],[  -72.4,  3760.0],[  -52.0,  3889.0],
[  -45.9,  3921.0],[  -21.3,  4070.0],[  -11.3,  4131.0],
[   10.5,  4267.0],[   20.8,  4316.0],[   72.8,  4649.0],
[  129.4,  5016.0],[  147.2,  5130.0],[  166.5,  5251.0],
[  232.5,  5696.0],[  252.1,  5808.0],[  260.0,  5880.0],],
dispers1st=[-2.66554e-07, 1.38090e-03,    6.0989, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_111 = dict(obsdir='WR121',obsid='00057507004',ext=1,
anchor=[382.85,183.7],angle=38.32,
dis1=[[  -99.5,  3609.0],[  -51.6,  3889.0],[  -47.1,  3921.0],
[   12.1,  4267.0],[   73.7,  4649.0],[  132.8,  5016.0],
[  152.0,  5130.0],[  195.1,  5395.0],[  236.6,  5696.0],],
dispers1st=[ 4.19263e-06, 4.06087e-04,    5.9626, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_112 = dict(obsdir='WR121',obsid='00057508004',ext=1,
anchor=[1114.4,1801],angle=39.22,
dis1=[[ -176.1,  3169.0],[  -99.0,  3609.0],
[  -80.0,  3717.0],[  -50.5,  3889.0],[  -45.3,  3921.0],
[   10.2,  4267.0],[   19.0,  4316.0],[   50.2,  4516.0],
[   71.1,  4649.0],[  128.5,  5016.0],[  145.2,  5130.0],
[  228.9,  5696.0],[  246.3,  5808.0],[  256.3,  5880.0],],
dispers1st=[-8.74030e-07, 1.68946e-03,    6.1781, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_113 = dict(obsdir='WR121',obsid='00057509004',ext=1,
anchor=[993.8,1436.5],angle=39.06,
dis1=[[ -237.6,  2840.0],[ -101.0,  3609.0],
[  -51.6,  3889.0],[  -46.1,  3921.0],[  -21.5,  4070.0],
[  -11.9,  4131.0],[   11.2,  4267.0],[   73.2,  4649.0],
[  132.1,  5016.0],[  148.4,  5130.0],[  234.2,  5696.0],
[  252.1,  5808.0],[  261.7,  5880.0],[  366.2,  6580.0],],
dispers1st=[-4.97955e-07, 1.37628e-03,    6.0699, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_114 = dict(obsdir='WR121',obsid='00057510004',ext=1,
anchor=[879.47,1153.1],angle=38.91,
dis1=[[ -167.7,  3230.0],[ -134.7,  3416.0],[ -117.2,  3511.0],
[ -100.5,  3609.0],[  -51.6,  3889.0],[  -13.5,  4131.0],
[   11.1,  4267.0],[   19.9,  4316.0],[   73.5,  4649.0],
[  133.1,  5016.0],[  149.3,  5130.0],[  236.2,  5696.0],
[  254.8,  5808.0],[  264.6,  5880.0],[  369.7,  6580.0],],
dispers1st=[-5.64767e-07, 1.33970e-03,    6.0225, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_115 = dict(obsdir='WR121',obsid='00057511004',ext=1,
anchor=[675.8,500.6],angle=38.47,
dis1=[[ -217.1,  2998.0],[ -182.7,  3169.0],[ -165.1,  3268.0],
[ -135.8,  3416.0],[ -102.8,  3609.0],[  -84.5,  3717.0],
[  -52.9,  3889.0],[  -47.2,  3921.0],[  -22.0,  4070.0],
[  -12.4,  4131.0],[   10.5,  4267.0],[   20.1,  4316.0],
[   27.8,  4367.0],[   44.9,  4472.0],[   74.8,  4649.0],
[  107.2,  4859.0],[  119.8,  4924.0],[  134.5,  5016.0],
[  151.8,  5130.0],[  239.7,  5696.0],[  268.3,  5880.0],],
dispers1st=[-8.75875e-07, 1.49677e-03,    5.9235, 4200.01],
dis1corr=    0.01,
)

# wheelpos=955
WR121_116 = dict(obsdir='WR121',obsid='00057512004',ext=1,
anchor=[934.45,847.62],angle=38.67,
dis1=[[ -266.4,  2699.0],
[ -169.4,  3230.0],[ -136.5,  3416.0],[ -102.7,  3609.0],
[  -74.5,  3760.0],[  -21.9,  4070.0],[   11.2,  4267.0],
[   53.2,  4516.0],[   74.5,  4649.0],[  119.7,  4924.0],
[  134.0,  5016.0],[  154.4,  5130.0],[  205.9,  5462.0],
[  241.0,  5696.0],[  259.9,  5808.0],[  270.3,  5880.0],],
dispers1st=[ 2.72772e-07, 1.08384e-03,    5.9019, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_117 = dict(obsdir='WR121',obsid='00057513004',ext=1,
anchor=[802.8,225.2],angle=38.17,
dis1=[[ -221.2,  2998.0],[ -167.2,  3268.0],[ -103.9,  3609.0],
[  -54.5,  3889.0],[  -48.9,  3921.0],[  -21.8,  4070.0],
[   10.9,  4267.0],[   20.1,  4316.0],[   29.0,  4367.0],
[   75.0,  4649.0],[  152.7,  5130.0],[  242.4,  5696.0],],
dispers1st=[-1.42218e-06, 1.61000e-03,    5.8658, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_118 = dict(obsdir='WR121',obsid='00057514004',ext=1,
anchor=[1147.4,1179.3],angle=38.86,
dis1=[[ -240.4,  2840.0],
[ -197.4,  3070.0],[ -117.5,  3511.0],[ -101.8,  3609.0],
[  -53.2,  3889.0],[  -47.3,  3921.0],[  -21.6,  4070.0],
[   11.5,  4267.0],[   51.8,  4516.0],[   74.6,  4649.0],
[  109.3,  4859.0],[  151.4,  5130.0],[  240.2,  5696.0],
[  259.2,  5808.0],[  269.0,  5880.0],[  377.9,  6580.0],],
dispers1st=[-7.44313e-07, 1.14568e-03,    5.9751, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_119 = dict(obsdir='WR121',obsid='00057515004',ext=1,
anchor=[1330.3,1475.1],angle=39.01,
dis1=[[ -247.3,  2840.0],
[ -102.5,  3609.0],[  -83.3,  3760.0],[  -51.9,  3889.0],
[  -46.7,  3921.0],[  -20.5,  4070.0],[  -12.2,  4131.0],
[   11.5,  4267.0],[   21.6,  4316.0],[   45.5,  4472.0],
[   75.6,  4649.0],[  134.8,  5016.0],[  152.4,  5130.0],
[  170.4,  5251.0],[  193.1,  5395.0],[  239.2,  5696.0],
[  258.7,  5808.0],[  267.2,  5880.0],[  374.8,  6580.0],],
dispers1st=[-7.55451e-07, 1.51958e-03,    5.9002, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_120 = dict(obsdir='WR121',obsid='00057516004',ext=1,
anchor=[1053.9,576.18],angle=38.43,
dis1=[[ -250.0,  2840.0],[ -185.2,  3169.0],
[ -136.9,  3416.0],[ -103.7,  3609.0],[  -53.6,  3889.0],
[  -48.8,  3921.0],[   11.2,  4267.0],[   20.4,  4316.0],
[   52.3,  4516.0],[   75.0,  4649.0],[  110.3,  4859.0],
[  135.4,  5016.0],[  153.5,  5130.0],[  173.2,  5251.0],
[  243.5,  5696.0],[  263.7,  5808.0],[  273.2,  5880.0],],
dispers1st=[-1.62681e-06, 1.36000e-03,    5.8845, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_121 = dict(obsdir='WR121',obsid='00057517004',ext=1,
anchor=[1167,860],angle=38.6,
dis1=[[ -216.3,  2992.0],
[ -101.3,  3616.0],[  -48.3,  3913.0],[  -21.3,  4080.0],
[   13.3,  4278.0],[   77.7,  4668.0],[  155.7,  5147.0],
[  244.7,  5709.0],[  263.7,  5834.0],[  273.7,  5891.0],],
dispers1st=[-1.06827e-06, 1.28964e-03,    5.9151, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_122 = dict(obsdir='WR121',obsid='00057518006',ext=1,
anchor=[1579.5,1894.6],angle=39.22,
dis1=[[ -260.1,  2840.0],
[ -214.4,  3070.0],[ -100.7,  3609.0],[  -52.6,  3889.0],
[  -47.3,  3921.0],[  -23.0,  4070.0],[   18.8,  4316.0],
[   72.3,  4649.0],[  129.1,  5016.0],[  147.1,  5130.0],],
dispers1st=[-2.09903e-06, 2.59345e-03,    6.0075, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_123 = dict(obsdir='WR121',obsid='00057519004',ext=1,
anchor=[1431.9,1050.1],angle=38.73,
dis1=[[ -282.6,  2699.0],[ -185.6,  3169.0],[ -139.9,  3416.0],
[ -103.8,  3609.0],[  -54.3,  3889.0],[  -49.3,  3921.0],
[   10.2,  4267.0],[   76.6,  4649.0],[  154.9,  5130.0],
[  244.0,  5696.0],[  264.5,  5808.0],[  384.5,  6580.0],],
dispers1st=[-1.25895e-06, 1.42253e-03,    5.8259, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_122 = dict(obsdir='WR121',obsid='00057520006',ext=1,
anchor=[1703.7,1427.9],angle=38.96,
dis1=[[ -237.9,  2840.0],[ -165.8,  3230.0],
[ -102.6,  3609.0],[  -82.6,  3717.0],[  -52.0,  3889.0],
[  -46.9,  3921.0],[  -22.8,  4070.0],[  -13.2,  4131.0],
[   11.0,  4267.0],[   22.0,  4316.0],[   75.8,  4649.0],
[  118.6,  4924.0],[  153.7,  5130.0],[  242.3,  5696.0],
[  263.1,  5808.0],[  272.2,  5880.0],[  382.6,  6580.0],],
dispers1st=[-3.79513e-07, 8.58600e-04,    5.9496, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_124 = dict(obsdir='WR121',obsid='00057521004',ext=1,
anchor=[1310.4,305.8],angle=38.15,
dis1=[[ -103.7,  3609.0],[  -88.0,  3717.0],[  -53.3,  3889.0],
[  -48.3,  3921.0],[  -22.8,  4070.0],[   11.1,  4267.0],
[   21.5,  4316.0],[   76.1,  4649.0],[  108.7,  4859.0],
[  136.8,  5016.0],[  158.1,  5130.0],[  176.4,  5251.0],
[  247.5,  5696.0],[  267.4,  5808.0],[  277.2,  5880.0],],
dispers1st=[-1.56360e-06, 1.35591e-03,    5.7873, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_125 = dict(obsdir='WR121',obsid='00057522004',ext=1,
anchor=[1453.6,604.92],angle=38.4,
dis1=[[ -274.4,  2740.0],[ -188.0,  3169.0],
[ -104.5,  3609.0],[  -86.5,  3717.0],[  -53.2,  3889.0],
[  -49.0,  3921.0],[  -23.4,  4070.0],[   10.6,  4267.0],
[   20.6,  4316.0],[   47.4,  4472.0],[   75.8,  4649.0],
[   99.8,  4785.0],[  138.3,  5016.0],[  155.8,  5130.0],
[  246.2,  5696.0],[  269.0,  5808.0],[  390.0,  6580.0],],
dispers1st=[-1.17574e-06, 1.30437e-03,    5.7687, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_126 = dict(obsdir='WR121',obsid='00057523004',ext=1,
anchor=[1766.8,850.2],angle=38.58,
dis1=[[ -209.6,  3070.0],[ -122.1,  3511.0],[ -104.3,  3609.0],
[  -53.6,  3889.0],[  -23.7,  4070.0],[   10.8,  4267.0],
[   53.6,  4516.0],[   76.2,  4649.0],[  154.8,  5130.0],
[  247.0,  5696.0],[  267.0,  5808.0],[  275.7,  5880.0],],
dispers1st=[-2.53514e-06, 1.52637e-03,    5.8330, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR121_127 = dict(obsdir='WR121',obsid='00057524004',ext=1,
anchor=[1798.6,385.52],angle=38.21,
dis1=[[ -122.0,  3511.0],[ -104.2,  3609.0],
[  -54.5,  3889.0],[  -48.3,  3921.0],[  -23.4,  4070.0],
[   11.4,  4267.0],[   19.8,  4316.0],[   46.4,  4472.0],
[   77.3,  4649.0],[  137.8,  5016.0],[  156.2,  5130.0],
[  247.9,  5696.0],[  271.4,  5808.0],[  280.7,  5880.0],],
dispers1st=[-2.22248e-06, 1.19877e-03,    5.8155, 4200.00],
dis1corr=    0.00,
)

# wheelpos=955
WR4_100 = dict(obsdir='WR4',obsid='00056900005',ext=1,
anchor=[1153.6,956.46],angle=38.71, # flux~ 1.5e-12 - too much coi: 4-10 times
dispers1st=[ -5.04956823e-03, 7.700,   4200.02], #bad
dis1=[[-149,2901], [65,4649], [253,5803]], 
dis1corr=0,
)

# wheelpos=955
WR4_101 = dict(obsdir='WR4',obsid='00056900006',ext=1,
anchor=[1041.6,954.46],angle=38.71, # flux~ 1.5e-12 - too much coi: 4-10 times
dispers1st=[ -4.85355811e-03,   7.74141311e+00,   4.20071786e+03], #bad
dis1=[[ -154.16 , 2901 ],[ -95.16 , 3409 ],[ -57.16 , 3737 ],
[ 58.84 , 4649 ],[ 244.84 , 5803 ],],
dis1corr=0,
)

# wheelpos=955
WR52_100 = dict(obsdir='WR52',obsid='00056950013',ext=1,
anchor=[1109.4,951.46],angle=38.71,  # flux~0.5e-12 - too much coi: 2-10 times
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR52_101 = dict(obsdir='WR52',obsid='00056950014',ext=1,
anchor=[1056.6,955.46],angle=38.71, 
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR52_102 = dict(obsdir='WR52',obsid='00056950019',ext=1,
anchor=[1104.6,957.46],angle=38.71, 
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR52_103 = dict(obsdir='WR52',obsid='00056950020',ext=1,
anchor=[1046.6,959.46],angle=38.71, 
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR92_100 = dict(obsdir='WR92',obsid='00037906005',ext=1,
anchor=[1088.6,944.46],angle=38.71, 
dispers1st=None,
dis1=[],
dis1corr=None,
)

# wheelpos=955
WR81_100 = dict(obsdir='WR81',obsid='00057550002',ext=1,
anchor=[767.57,814.46],angle=38.81, 
dispers1st=None,
dis1=[],
dis1corr=None,
)

#]  # end vis0955
#===================================================================
#===================================================================

#  Some different views on the data

clockeduv = (cur25,cur19,cur09,cur24,cur02,cur04,cur23,cur21,cur16,\
             cur06,cur22,cur08,cur01,cur12,cur18,cur20,cur05,cur14,\
             cur11,cur15,cur13,cur26,cur27,\
             cur40,cur41,cur42,cur43,cur44)
             # cur40,cur41,cur42 WD for curvature
             # cur 43,cur44 T Pyx+HST 

# for new anchor calibration (mid 2013) 
uv0160 = (cur01,cur02,cur04,cur05,cur06,cur08,cur09,\
  cur11,cur12,cur13,cur14,cur15,cur16,cur18,cur19,\
  cur20,cur21,cur22,cur23,cur24,cur25,cur26,cur27, 
  cur28, cur29, cur30, cur31, cur43, cur44)
             
nominaluv = (cur100,cur101,cur102,cur103,cur104,cur105,cur106,cur107,\
             cur108,cur109,cur110,cur111,cur112,cur113,cur114,cur115,\
             cur116,cur117,cur118,cur119,cur120,cur121,cur122,cur123,\
             cur124,cur125,cur126,cur127,cur128,cur129,cur130 , \
             cur216,cur218,cur316,cur318,)           
uv0200 = (cur100,cur101,      cur103,cur104,cur105,cur106,cur107,\
             cur108,cur109,cur110,cur111,cur112,cur113,cur114,cur115,\
             cur116,cur117,cur118,cur119,cur120,cur121,cur122,cur123,\
             cur124,cur125,cur126,cur127,cur128,cur129,cur130 , )            

vis1000 = (WR121_001,WR121_002,WR121_003,WR121_004,WR121_005,
           WR121_006,WR121_007,WR121_008,WR121_009,WR121_010,
           WR121_011,WR121_012,WR121_013,WR121_014,WR121_015,
           WR121_016,WR121_017,WR121_018,WR121_019,WR121_020,
           WR121_021,WR121_022,WR121_023,WR121_024,WR121_025,
           WR121_026,WR121_027,WR121_028,WR121_029,WR121_030,
           WR121_031,WR121_032,WR121_033,WR121_034)

vis0955 = (WR121_101,WR121_108,WR121_109,WR121_110,
           WR121_111,WR121_112,WR121_113,WR121_114,WR121_115,
           WR121_116,WR121_117,WR121_118,WR121_119,WR121_120,
           WR121_121,WR121_122,WR121_123,WR121_124,WR121_125,
           WR121_126,WR121_127,)
             
def get_curvedata_asarray(arraykey, dim=None,curvedata=clockeduv):
   '''reorganize the data from arrays. 
      use dim to force the dimension'''
   if not arraykey in list(curvedata[0].keys()): 
      print("Key error in get_curvedata_asarray : "+arraykey) 
      return   
   N = len(curvedata)
   if dim == None: 
      M = len(atleast_1d(asarray(curvedata[0].get(arraykey) )))
   else: M = dim   
   item=zeros(M*N).reshape(M,N)
   for i in range(N):
      ank = atleast_1d(asarray(curvedata[i].get(arraykey)))
      for j in range(len(ank)):
        item[j,i] = ank[j]           
   return item

def get_curvedata_aslist(key,curvedata=clockeduv):
   '''reorganize the data from arrays. 
      use dim to force the dimension'''
   if not key in list(curvedata[0].keys()): 
      print("Key error in get_curvedata_aslist : "+key)
      return   
   N = len(curvedata)
   out=list()
   for i in range(N):
      item = curvedata[i].get(key)
      out.append(item)       
   return out

def get_angle(curvedata=clockeduv):
   '''angle used for the measured curvature '''   
   N = len(curvedata)
   angle=zeros(N)
   for i in range(N):
      z = curvedata[i].get('angle')
      angle[i] = z           
   return angle

def get_present(curvedata=clockeduv):
   '''True or False arrays for zeroth first second third orders '''
   N = len(curvedata)
   present0 = array(zeros(N),dtype=bool)
   present1 = array(zeros(N),dtype=bool)
   present2 = array(zeros(N),dtype=bool)
   present3 = array(zeros(N),dtype=bool)
   for i in range(N):
      present0[i] = curvedata[i].get('present0')
      present1[i] = curvedata[i].get('present1')
      present2[i] = curvedata[i].get('present2')
      present3[i] = curvedata[i].get('present3')
   return present0,present1,present2,present3
   
def analyse1(dothis,curvedata=clockeduv):   
   '''analyse the curved tracks '''
   # get the data in arrays, lists
   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   anchor = get_curvedata_asarray('anchor',curvedata=curvedata)
   dlim0L = get_curvedata_asarray('dlim0L',curvedata=curvedata)
   dlim0U = get_curvedata_asarray('dlim0U',curvedata=curvedata)
   coef0  = get_curvedata_aslist('coef0',curvedata=curvedata)
   dlim1L = get_curvedata_asarray('dlim1L',curvedata=curvedata)
   dlim1U = get_curvedata_asarray('dlim1U',curvedata=curvedata)
   coef1  = get_curvedata_aslist('coef1',curvedata=curvedata)
   dlim2L = get_curvedata_asarray('dlim2L',curvedata=curvedata)
   dlim2U = get_curvedata_asarray('dlim2U',curvedata=curvedata)
   coef2  = get_curvedata_aslist('coef2',curvedata=curvedata)
   dlim3L = get_curvedata_asarray('dlim3L',curvedata=curvedata)
   dlim3U = get_curvedata_asarray('dlim3U',curvedata=curvedata)
   coef3  = get_curvedata_aslist('coef3',curvedata=curvedata)
   angle  = get_angle(curvedata=curvedata)
   flatangle = get_curvedata_asarray('flatangle',curvedata=curvedata)
   zmxangle  = get_curvedata_asarray('zmxangle',curvedata=curvedata)
   obsid     = get_curvedata_asarray('obsid',curvedata=curvedata)
   ext       = get_curvedata_asarray('ext',curvedata=curvedata)
   sig0coef  = get_curvedata_aslist('sig0coef',curvedata=curvedata)
   sig1coef  = get_curvedata_aslist('sig1coef',curvedata=curvedata)
   sig2coef  = get_curvedata_aslist('sig2coef',curvedata=curvedata)
   sig3coef  = get_curvedata_aslist('sig3coef',curvedata=curvedata)
   #
   Nobs = len(dlim1L.flatten())
   #   first order always present
   #   select fixed points
   X1 = arange(-375,1150,25)
   #   results list 
   Y1 = list() 
   q1 = list()
   for i in range(Nobs):
      # pixel coordinate relative to anchor point
      #  X = arange(dlim1L[i],dlim1U[i])
      # the remaining offset at the anchor point 
      Yref = coef1[i][-1]   
      Y1.append( (polyval(coef1[i],X1) - Yref) )   
      q1.append( ( (X1 > dlim1L.flatten()[i]) & (X1 < dlim1U.flatten()[i]) ) )
   Y1 = array(Y1)
   q1 = array(q1)
   
   X2 = arange(25,1150,25)
   Y2 = list() 
   q2 = list()
   for i in range(Nobs):
     if present2[i]: 
       Yref = coef1[i][-1]   
       Y2.append( (polyval(coef2[i],X2) - Yref) )   
       q2.append( ( (X2 > dlim2L.flatten()[i]) & (X2 < dlim2U.flatten()[i]) ) )
   Y2 = array(Y2)
   q2 = array(q2)   
   
   X3 = arange(450,1150,25)   
   Y3 = list() 
   q3 = list()
   qank = list()
   for i in range(Nobs):
     if present3[i]: 
       Yref = coef1[i][-1]   
       Y3.append( (polyval(coef3[i],X3) - Yref) )   
       q3.append( ( (X3 > dlim3L.flatten()[i]) & (X3 < dlim3U.flatten()[i]) ) )
       qank.append( i )
   Y3 = array(Y3)
   q3 = array(q3)   
   #   now return the result 
   if dothis == None:
      return anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank
   if dothis == 0:
      X0 = arange(-800,-550,25)
      #   results list 
      Y0 = list() 
      q0 = list()
      for i in range(Nobs):
       if present0[i]:
         # pixel coordinate relative to anchor point
         #  X = arange(dlim1L[i],dlim1U[i])
         # the remaining offset at the anchor point 
         Yref = coef1[i][-1]   
         Y0.append( (polyval(coef0[i],X0) - Yref) )   
         q0.append( ( (X0 > dlim0L.flatten()[i]) & (X0 < dlim0U.flatten()[i]) ) )
      Y0 = array(Y0)
      q0 = array(q0)
      return X0,Y0,q0      
   
def get_1stOrderFit(xin=None,yin=None,curvedata=clockeduv,kxky=1):
   '''Make a fit to the coefficients of the first order. 
      option to rotate the measured fit  first (to zemax angle; 
      or to original positions on det image; then make  
      a bilinear fit.
   '''   
   from uvotmisc import uvotrotvec
   from scipy import interpolate
   print("starting cal3.get_1stOrderFit")
   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank = analyse1(None,curvedata=curvedata)
   angle = get_angle(curvedata=curvedata)
   flatangle = (get_curvedata_asarray('flatangle',curvedata=curvedata)).flatten()
   zmxangle  = (get_curvedata_asarray('zmxangle',curvedata=curvedata)).flatten()
   qa = angle != zmxangle
   qi = where(qa)[0]
   dang = angle-zmxangle
   for i in qi:
      print("adjusting angle data point index i= ",i,' with delta angle= ',dang[i])
      x1_prime, y1_prime = uvotrotvec(X1,Y1[i,:], dang[i] )
      coef = polyfit(x1_prime,y1_prime,3)
      Y1[i,:] = polyval(coef,X1) 
      angle[i] = zmxangle[i]  
   # now all data are on the zemax angle (actually on the angle provided 
   # from the calibration file, interpolated, but based on the zemax angle.
   #
   # first compute a third order polynomial for each measurement
   Xank = anchor[0][present1]
   Yank = anchor[1][present1]
   M = len(Xank)
   c0 = zeros(M,dtype=float)
   c1 = zeros(M,dtype=float)
   c2 = zeros(M,dtype=float)
   c3 = zeros(M,dtype=float)
   for i in range(M):
      coef = polyfit( X1[q1[i]],Y1[i,q1[i]],3) 
      c3[i], c2[i], c1[i], c0[i] = coef
   # then compute a bilinear fit to the 4 coefficients of the polynomials.
   # c0 is essentially zero 
   tck_c1 = interpolate.bisplrep(Xank,Yank,c1,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   tck_c2 = interpolate.bisplrep(Xank,Yank,c2,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   tck_c3 = interpolate.bisplrep(Xank,Yank,c3,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   #     Z = interpolate.bisplev(xx, yy, tck)
   if (xin == None) | (yin == None) :
      XX, YY = linspace(1,2047),linspace(1,2047)
      cc0 = zeros(2500).reshape(50,50)
      cc1 = interpolate.bisplev(XX,YY,tck_c1)
      cc2 = interpolate.bisplev(XX,YY,tck_c2)
      cc3 = interpolate.bisplev(XX,YY,tck_c3)
      XX, YY = meshgrid(XX,YY)
      return XX,YY,(cc3,cc2,cc1,cc0),(tck_c3,tck_c2,tck_c1), (c3,c2,c1,c0)
   else:
      coef = array([interpolate.bisplev(xin,yin,tck_c3),interpolate.bisplev(xin,yin,tck_c2),\
      interpolate.bisplev(xin,yin,tck_c1), 0.]) 
      # return the coefficients of the fitting polynomial Y1(X1) for anchor point (xin,yin)
      return coef  
   
def get_2ndOrderFit(xin=None,yin=None,curvedata=clockeduv,kxky=1):
   '''Make a fit to the coefficients of the second order. 
      option to rotate the measured fit  first (to zemax angle; 
      or to original positions on det image; then make  
      a bilinear fit.
   '''   
   from uvotmisc import uvotrotvec
   from scipy import interpolate
   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank = analyse1(None,curvedata=curvedata)
   angle = get_angle(curvedata=curvedata)
   flatangle = (get_curvedata_asarray('flatangle',curvedata=curvedata)).flatten()
   zmxangle  = (get_curvedata_asarray('zmxangle',curvedata=curvedata)).flatten()
   qa = angle != zmxangle
   qi = where(qa)[0]
   dang = angle-zmxangle
   for i in qi: 
      x1_prime, y1_prime = uvotrotvec(X1,Y1[i,:], dang[i] )
      coef = polyfit(x1_prime,y1_prime,3)
      Y1[i,:] = polyval(coef,X1) 
      x2_prime, y2_prime = uvotrotvec(X2,Y2[i,:], dang[i] )
      coef = polyfit(x2_prime,y2_prime,3)
      Y2[i,:] = polyval(coef,X2) 
      angle[i] = zmxangle[i]  
   # now all data are on the zemax angle (actually on the angle provided 
   # from the calibration file, interpolated, but based on the zemax angle.
   #
   # first compute a second order polynomial for each measurement
   Xank = anchor[0][present2]
   Yank = anchor[1][present2]
   M = len(q2[:,0])
   c0 = zeros(M,dtype=float)
   c1 = zeros(M,dtype=float)
   c2 = zeros(M,dtype=float)
   for i in range(M):
      coef = polyfit( X2[q2[i,:]],Y2[i,q2[i,:]],2) 
      c2[i], c1[i], c0[i] = coef
   # then compute a bilinear fit to the 4 coefficients of the polynomials.
   tck_c0 = interpolate.bisplrep(Xank,Yank,c0,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   tck_c1 = interpolate.bisplrep(Xank,Yank,c1,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   tck_c2 = interpolate.bisplrep(Xank,Yank,c2,kx=kxky,ky=kxky,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   #     Z = interpolate.bisplev(xx, yy, tck)
   if (xin == None) | (yin == None) :
      XX, YY = linspace(1,2047),linspace(1,2047)
      cc0 = interpolate.bisplev(XX,YY,tck_c0)
      cc1 = interpolate.bisplev(XX,YY,tck_c1)
      cc2 = interpolate.bisplev(XX,YY,tck_c2)
      XX, YY = meshgrid(XX,YY)
      return XX,YY,(cc2,cc1,cc0),(tck_c2,tck_c1,tck_c0), (c2,c1,c0)
   else:    
      coef = array([interpolate.bisplev(xin,yin,tck_c2),interpolate.bisplev(xin,yin,tck_c1),\
      interpolate.bisplev(xin,yin,tck_c0)]) 
      # return the coefficients of the fitting polynomial Y1(X1) for anchor point (xin,yin)
      return coef  
     
def get_3rdOrderFit(xin=None,yin=None,curvedata=clockeduv):
   '''Make a fit to the coefficients of the third order. 
      option to rotate the measured fit  first (to zemax angle; 
      or to original positions on det image; then make  
      a bilinear fit.
   '''   
   from uvotmisc import uvotrotvec
   from scipy import interpolate

   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank = analyse1(None,curvedata=curvedata)
   angle = get_angle(curvedata=curvedata)
   flatangle = (get_curvedata_asarray('flatangle',curvedata=curvedata)).flatten()
   zmxangle  = (get_curvedata_asarray('zmxangle',curvedata=curvedata)).flatten()
   qa = abs(angle - zmxangle) > 1.0e-5
   qi = where(qa)[0]
   qp = where(present3)[0]
   NN = len(qp)
   print('NN',NN)
   print('qi',qi)
   print('qp',qp)
   dang = angle-zmxangle
   for i in range(NN):
      k = qp[i]
      x3_prime, y3_prime = uvotrotvec(X3,Y3[i,:], dang[k] )
      coef = polyfit(x3_prime,y3_prime,2)
      Y3[i,:] = polyval(coef,X3) 
      angle[k] = zmxangle[k]  
   # now all data are on the zemax angle (actually on the angle provided 
   # from the calibration file, interpolated, but based on the zemax angle.
   #
   # first compute a third order polynomial for each measurement
   Xank = anchor[0][present3]
   Yank = anchor[1][present3]
   M = len(q3[:,0])
   c0 = zeros(M,dtype=float)
   c1 = zeros(M,dtype=float)
   print('M',M)
   #c2 = zeros(M,dtype=float)
   for i in range(M):
      coef = polyfit( X3[q3[i,:]],Y3[i,q3[i,:]],1) 
      c1[i], c0[i] = coef
   # then compute a bilinear fit to the 3 coefficients of the polynomials.
   #  
   tck_c0 = None
   tck_c1 = None
   print('Xank',Xank)
   print('Yank',Yank)
   print('c0',c0)
   tck_c0 = interpolate.bisplrep(Xank,Yank,c0,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048) 
   print('tck_c0',tck_c0)
   print('c1',c1)    
   tck_c1 = interpolate.bisplrep(Xank,Yank,c1,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048)   
   print('tck_c1',tck_c1)  
   #tck_c2 = interpolate.bisplrep(Xank[qank],Yank[qank],c2,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   #     Z = interpolate.bisplev(xx, yy, tck)
   if (xin == None) | (yin == None) :
      XX, YY = linspace(0,2047),linspace(0,2047)
      cc0 = interpolate.bisplev(XX,YY,tck_c0)
      cc1 = interpolate.bisplev(XX,YY,tck_c1)
      #cc2 = interpolate.bisplev(XX,YY,tck_c2)
      XX, YY = meshgrid(XX,YY)
      return XX,YY,(cc1,cc0),(tck_c1,tck_c0), (c1,c0)
   else:    
      coef = array([interpolate.bisplev(xin,yin,tck_c1),interpolate.bisplev(xin,yin,tck_c0)]) 
      # return the coefficients of the fitting polynomial Y1(X1) for anchor point (xin,yin)
      return coef  
   
def get_0thOrderFit(xin=None,yin=None,curvedata=clockeduv):
   '''Make a fit to the coefficients of the zeroth order. 
      option to rotate the measured fit  first (to zemax angle; 
      or to original positions on det image; then make  
      a bilinear fit.
   '''   
   from uvotmisc import uvotrotvec
   from scipy import interpolate
   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank = analyse1(None,curvedata=curvedata)
   angle = get_angle(curvedata=curvedata)
   flatangle = (get_curvedata_asarray('flatangle',curvedata=curvedata)).flatten()
   zmxangle  = (get_curvedata_asarray('zmxangle',curvedata=curvedata)).flatten()
   qa = angle != zmxangle
   qi = where(qa)[0]
   dang = angle-zmxangle
   for i in qi: 
      x1_prime, y1_prime = uvotrotvec(X1,Y1[i,:], dang[i] )
      coef = polyfit(x1_prime,y1_prime,3)
      Y1[i,:] = polyval(coef,X1) 
      angle[i] = zmxangle[i]  
   # now all data are on the zemax angle (actually on the angle provided 
   # from the calibration file, interpolated, but based on the zemax angle.
   #
   # first compute a second order polynomial for each measurement
   Xank = anchor[0][present0]
   Yank = anchor[1][present0]
   X0,Y0,q0 = analyse1(0)
   M = len(q0)
   c0 = zeros(M,dtype=float)
   c1 = zeros(M,dtype=float)
   c2 = zeros(M,dtype=float)
   #c3 = zeros(M,dtype=float)
   for i in range(M):
      coef = polyfit( X0[q0[i]],Y0[i,q0[i]],1) 
      c1[i], c0[i] = coef
   # then compute a bilinear fit to the 2 coefficients of the polynomials.
   # c0 is non-zero 
   tck_c0 = interpolate.bisplrep(Xank,Yank,c0,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   tck_c1 = interpolate.bisplrep(Xank,Yank,c1,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   #tck_c2 = interpolate.bisplrep(Xank[MM],Yank[MM],c2,kx=1,ky=1,s=None,xb=0,xe=2048,yb=0,ye=2048)     
   #     Z = interpolate.bisplev(xx, yy, tck)
   if (xin == None) | (yin == None) :
      XX, YY = linspace(0,2047),linspace(0,2047)
      cc0 = interpolate.bisplev(XX,YY,tck_c0)
      cc1 = interpolate.bisplev(XX,YY,tck_c1)
      #cc2 = interpolate.bisplev(XX,YY,tck_c2)
      cc2 = cc1.copy() * 0.0
      XX, YY = meshgrid(XX,YY)
      return XX,YY,(cc2,cc1,cc0),(tck_c1,tck_c0), (c2,c1,c0)
   else:    
      coef = array([interpolate.bisplev(xin,yin,tck_c1),\
                    interpolate.bisplev(xin,yin,tck_c0)]) 
      # return the coefficients of the fitting polynomial Y0(X0) for anchor point (xin,yin)
      return coef  
  
""" 
COMPARISON BETWEEN THE MEASURED CURVATURES AND 
BILINEAR FITS TO THE COEFFICIENTS FOR THE DISTANCE 
BETWEEN FIRST AND SECOND ORDER AT PIXEL POSITION 30
WHICH IS AROUND THE START OF THE SECOND ORDER.

                                   at pixel position = +30
               anchor              measured  fit   difference
               observed            2nd-1st   2nd-1st
obsid                          
00057925001 [ 1070.9, 1733.6 ]  -  24.5 -  22.6 :  1.94
00057919001 [ 1365.3, 1510.5 ]  -  22.9 -  22.6 :  0.30
00057909001 [  658.1, 1568.5 ]  -   9.5 -  11.6 : -2.15
00057924001 [ 1247.6,  464.3 ]  -   0.3 -   1.0 : -0.72
00057902001 [  631.3,  169.2 ]  - -14.3 - -14.8 :  0.53
00057904001 [ 1289.0,  145.5 ]  -  -1.0 -   0.3 : -1.30
00057923001 [ 1042.0,  204.4 ]  -   0.8 -   0.0 :  0.72
00057921001 [ 1822.5, 1179.6 ]  -  23.8 -  23.3 :  0.49
00057916001 [  917.1,  455.3 ]  -  -3.2 -  -6.1 :  2.96
00057906001 [ 1854.2,  871.4 ]  -  18.8 -  19.0 : -0.12
00057922001 [  571.9,  482.3 ]  -  -9.2 -  -9.4 :  0.22
00057908001 [  445.0, 1305.6 ]  -   2.4 -   0.3 :  2.17
00057901001 [  308.7,  466.1 ]  - -14.4 - -13.3 : -1.11
00057912001 [  901.0,  752.9 ]  -  -2.7 -   0.0 : -2.67
00057918001 [  929.7, 1349.0 ]  -  10.7 -  12.5 : -1.87
00057920001 [ 1427.3, 1158.7 ]  -  14.5 -  16.5 : -1.99
00057905001 [ 1599.7,  512.8 ]  -   0.2 -  -0.5 :  0.76
00057914001 [ 1147.8, 1094.5 ]  -   6.1 -   7.5 : -1.37
00057911001 [  859.9, 1062.3 ]  -   3.2 -   0.7 :  2.42
00057915001 [  427.8,  805.3 ]  -  -1.2 -  -1.0 : -0.19
00057913001 [ 1262.3,  810.2 ]  -   0.5 -   0.0 :  0.48
00057926001 [  546.9, 1078.4 ]  -   1.5 -   0.2 :  1.37
00056950004 [ 1058.3,  956.9 ]  -   1.2 -   0.0 :  1.17

And ~160A further down the spectrum: 

                                   at pixel position = +300
               anchor              measured  fit   difference
               observed            2nd-1st   2nd-1st
 obsid                             
00057925001 [ 1070.9, 1733.6 ]  -  17.6 -  12.9 :  4.68
00057919001 [ 1365.3, 1510.5 ]  -  14.7 -  13.9 :  0.72
00057909001 [  658.1, 1568.5 ]  -   7.3 -   6.2 :  1.13
00057924001 [ 1247.6,  464.3 ]  -   0.2 -   1.0 : -0.78
00057902001 [  631.3,  169.2 ]  -  -8.5 - -10.4 :  1.87
00057904001 [ 1289.0,  145.5 ]  -   0.8 -  -0.3 :  1.14
00057923001 [ 1042.0,  204.4 ]  -   3.7 -  -0.2 :  3.92
00057921001 [ 1822.5, 1179.6 ]  -  14.0 -  15.6 : -1.58
00057916001 [  917.1,  455.3 ]  -  -0.9 -  -6.0 :  5.08
00057906001 [ 1854.2,  871.4 ]  -  11.3 -   9.4 :  1.84
00057922001 [  571.9,  482.3 ]  -  -4.7 -  -5.5 :  0.83
00057908001 [  445.0, 1305.6 ]  -   3.6 -  -0.1 :  3.70
00057901001 [  308.7,  466.1 ]  -  -8.0 -  -8.3 :  0.34
00057912001 [  901.0,  752.9 ]  -  -2.1 -   0.0 : -2.10
00057918001 [  929.7, 1349.0 ]  -   7.4 -   6.5 :  0.93
00057920001 [ 1427.3, 1158.7 ]  -   8.2 -  10.3 : -2.08
00057905001 [ 1599.7,  512.8 ]  -  -2.7 -   0.6 : -3.27
00057914001 [ 1147.8, 1094.5 ]  -   2.6 -   5.0 : -2.38
00057911001 [  859.9, 1062.3 ]  -   2.2 -   0.2 :  1.95
00057915001 [  427.8,  805.3 ]  -   2.7 -  -1.5 :  4.17
00057913001 [ 1262.3,  810.2 ]  -  -2.1 -   0.0 : -2.06
00057926001 [  546.9, 1078.4 ]  -   3.0 -  -0.1 :  3.11
00056950004 [ 1058.3,  956.9 ]  -  -0.7 -   0.0 : -0.74

"""

def plot_orders(curvedata=clockeduv):
   from . import cal3
   from uvotmisc import uvotrotvec
   from pylab import plot
   
   anchor, X1, Y1, q1, X2, Y2, q2, X3, Y3, q3, qank = cal3.analyse1(None,curvedata=curvedata)
   X0, Y0, q0 = cal3.analyse1(0,curvedata=curvedata)
   present0,present1,present2,present3 = get_present(curvedata=curvedata)
   dlim0L = get_curvedata_asarray('dlim0L',curvedata=curvedata).flatten()
   dlim0U = get_curvedata_asarray('dlim0U',curvedata=curvedata).flatten()
   coef0  = get_curvedata_aslist('coef0',curvedata=curvedata)
   dlim1L = get_curvedata_asarray('dlim1L',curvedata=curvedata).flatten()
   dlim1U = get_curvedata_asarray('dlim1U',curvedata=curvedata).flatten()
   coef1  = get_curvedata_aslist('coef1',curvedata=curvedata)
   dlim2L = get_curvedata_asarray('dlim2L',curvedata=curvedata).flatten()
   dlim2U = get_curvedata_asarray('dlim2U',curvedata=curvedata).flatten()
   coef2  = get_curvedata_aslist('coef2',curvedata=curvedata)
   dlim3L = get_curvedata_asarray('dlim3L',curvedata=curvedata).flatten()
   dlim3U = get_curvedata_asarray('dlim3U',curvedata=curvedata).flatten()
   coef3  = get_curvedata_aslist('coef3',curvedata=curvedata)
   angle  = get_angle(curvedata=curvedata) - 180.0
   flatangle = get_curvedata_asarray('flatangle',curvedata=curvedata)
   zmxangle  = get_curvedata_asarray('zmxangle',curvedata=curvedata)
   obsid     = get_curvedata_asarray('obsid',curvedata=curvedata)
   ext       = get_curvedata_asarray('ext',curvedata=curvedata)
   #
   indix = where(present0)[0]
   n0 = len(indix)
   for j in range(n0):
      i = indix[j]
      xank = anchor[0,i]
      yank = anchor[1,i]
      q = where( (X0 > dlim0L[i]) & (X0 < dlim0U[i]) )
      xr,yr = uvotrotvec(X0[q],Y0[j,q],angle[i]) 
      xr = xank + xr.flatten()
      yr = yank + yr.flatten()
      plot(xr,yr,'y-',linewidth=1.3)
   #
   indix = where(present1)[0]
   n0 = len(indix)
   for j in range(n0):
      i = indix[j]
      xank = anchor[0,i]
      yank = anchor[1,i]
      q = where( (X1 > dlim1L[i]) & (X1 < dlim1U[i]) )
      xr,yr = uvotrotvec(X1[q],Y1[j,q],angle[i]) 
      xr = xank + xr.flatten()
      yr = yank + yr.flatten()
      print('===============order 1===',i)
      print(xr)
      print(yr)
      plot(xr,yr,'k-',linewidth=1.2)
   #
   indix = where(present2)[0]
   n0 = len(indix)
   for j in range(n0):
      i = indix[j]
      xank = anchor[0,i]
      yank = anchor[1,i]
      q = where( (X2 > dlim2L[i]) & (X2 < dlim2U[i]) )
      plot(xr,yr,'k-',linewidth=1.2)
      xr,yr = uvotrotvec(X2[q],Y2[j,q],angle[i]) 
      xr = xank+xr.flatten()
      yr = yank+yr.flatten()
      print('===============order 2===',i)
      print(xr)
      print(yr)
      plot(xr,yr,'b-',linewidth=1.3)
   #
   indix = where(present3)[0]
   n0 = len(indix)
   for j in range(n0):
      i = indix[j]
      xank = anchor[0,i]
      yank = anchor[1,i]
      q = where( (X3 > dlim3L[i]) & (X3 < dlim3U[i]) )
      xr,yr = uvotrotvec(X3[q],Y3[j,q],angle[i]) 
      xr = xank+xr.flatten()
      yr = yank+yr.flatten()
      plot(xr,yr,'r-',linewidth=1.3)
      
def make_test_dir(rootdir='/calibration/grism/test/',what= vis1000, 
    auxil=True,raw=True, det=True, sky=True, decompress=True, lent=True,
    local_archive='/Volumes/data1/wavecal.data/', chatter=0):
   ''' 
   walk through observations in list in "what", and refresh directory 
   with new data files of the type listed (auxil, raw, det, sky) by first 
   removing the files, then putting a new in place, and 
   decompressing the data
    
   '''    
   import shutil
   import os
   from pull_swift_observation import get_uvot_data_from_obsid
   
   # check root exists
   if not os.access(rootdir,os.F_OK):
      raise IOError('Parameter rootdir=%s does not exist on filesystem '%(rootdir))
   # run through list
   for obs in what:
      name = obs['obsdir']
      obsid = obs['obsid']
      if len(obsid) < 11:
         print("problem with copying the following observation, obsid looks wrong.")
         print(obs)
         raise RuntimeError("stopping") 
      # check if data directory present
      #  Yes => delete
      if os.access(rootdir+name+'/'+obsid+'/',os.F_OK):
         shutil.rmtree(rootdir+name+'/'+obsid+'/',ignore_errors=False,onerror=None)
      # check local archive has data
      if not os.access(local_archive+'/'+name+'/',os.F_OK):
         #  no => retrieve from archive at GSFC into archive, then copy selected
         get_uvot_data_from_obsid(obsid,
                      dateobs=None,
                      rootdir=local_archive+'/'+name+'/',
                      cleanup=True,
                      chatter=1)                        
      _copy_selected(obsid,
                    fromm=local_archive+name+'/'+obsid+'/',
                    to=rootdir+name+'/'+obsid+'/',
                    auxil=auxil, raw=raw, det=det, sky=sky,
                    lent=lent,decompress=decompress,)
      
      
      
def _copy_selected(obsid,fromm=None, to=None, auxil=True,raw=True, det=True, 
    sky=True, lent=True, decompress=True,) :
    import os, shutil
    
    f1 = fromm+'auxil/'
    t1 = to + 'auxil/'
    if not os.access(t1,os.F_OK):
        os.makedirs(t1,mode=0o755)     
    f2 = fromm+'uvot/image/'
    t2 = to + 'uvot/image/'
    if not os.access(t2,os.F_OK):
        os.makedirs(t2,mode=0o755)     
    if auxil:
       os.system('cp '+f1+'* '+t1)
       os.system('gunzip '+t1+'*.gz')
    if raw:
       os.system('cp '+f2+'sw*_rw.im* '+t2)
       os.system('gunzip '+t2+'sw*_rw.img.gz')
    if det:
       os.system('cp '+f2+'sw*_dt.im* '+t2)
       os.system('gunzip '+t2+'sw*_dt.img.gz')
    if sky and lent: 
       os.system('cp '+f2+'sw*_sk.im* '+t2)
       os.system('gunzip '+t2+'sw*_sk.img.gz')
    if sky and not lent: 
       os.system('cp '+f2+'sw'+obsid+'ug?_sk.im* '+t2)
       os.system('gunzip '+t2+'sw*_sk.img.gz')
           
def fixed_wave(wave):
   """
   correct wavelength
   
   Parameters
   ==========
   wave : float or float array
   
   Returns
   =======
   wave: float or array
      corrected 
   
   Notes
   ===== 
   The wavelengths for the offset spectra [1500,1500] are a bit off, 
   The correction is twofold:
   for wave < 2133A: linear polynomial with coef=[ -3.0105e-02, -6.42027e+01]
   for wave > 4018A: linear polynomial with coef=[ -2.30762e-02, 9.27202e+01])
 
   used for nova V339 Del 
   2014-05-09 NPMK    
   """
   import numpy as np
   coef1=np.array([  3.0105e-02, -6.42027e+01])
   coef2=np.array([ -2.30762e-02, 9.27202e+01])
   w = wave.copy()
   q = w < 2133.
   w[q] = w[q] - polyval(coef1,w[q])
   q = w > 4018.
   w[q] = w[q] + polyval(coef2,w[q])
   return w
