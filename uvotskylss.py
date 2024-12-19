def uvotskylss(infile=None,test=False):
   ''' 
   The routine creates the lss file in the directories 
   given a list as an infile.
   
   Make sure the attitude files have been unzipped and 
   are of the  uat kind.
   
   2011-04-08 NPMK rough version
   '''
   import os
   if infile == None: 
      print "required input missing"
      return
      
   pwd = os.getcwd()
   x = infile.split('@',1)
   if len(x) == 2: 
      # list 
      f = open(x[1])
      lines = f.readlines()
      f.close()
      for line in lines:
         x2 = line.rsplit('/',1)
	 if len(x2) == 1:
	   dir1="."
	   _makelss(pwd,dir1,x2[0],test=test)
	 elif len(x2) == 2:
	   _makelss(pwd,x2[0],x2[1],test=test)  
      
   elif len(x) == 1:
      # not a list, just a filename
      x2 = line.rsplit('/',1)
      if len(x2) == 1:
	   dir1="."
	   _makelss(pwd,dir1,x2[0],test=test)
      elif len(x2) == 2:
	   _makelss(pwd,x2[0],x2[1],test=test) 
	   
   else:	   
      print "error in skyfile" 
            
   
def _makelss(pwd, dir1, file1, test=False ):
   '''create the lss file assuming the attitude file is the uat kind '''
   import os
   obsid = ((file1.rsplit("sw",1))[1].rsplit("u"))[0]
   filestub = file1.split("img")[0]
   attfile='sw'+obsid+'uat.fits'
   outfile=filestub+'lss'
   infile = filestub+'img'
   os.chdir(dir1)
   command="uvotskylss infile="+infile+"  outfile="+outfile+ \
     " attfile=../../auxil/"+attfile+" lssfile=CALDB clobber=YES "
   
   if test:
     print dir1,'  ', file1
     print command
   else:    
     print "\n=================",dir1,' ..... ', file1
     print command
     if os.system(command) :
        print "*** >>> problem creating lss file: "+outfile
	
   os.chdir(pwd)     
