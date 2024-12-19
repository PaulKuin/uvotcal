
#  2014-02-21 Paul Kuin MSSL
#  Study of the extra coincidence loss when very bright source (T Pyx 2011 nova) 
#  or 
#  bright extended source (SN2011fe in Messier 101)
#
import pylab as plt

# questions
# =========
#
# (1) after making coi-correction to spectrum is whole spectrum scaled 
#     or does that depend on flux(lambda) ?
#
# (2) is the flux in other field sources affected as well, or just that 
#     in the bright source?
# 
# (3) is this consistent with the lower flux in WD photometry in the white band ?
#
# (4) is there any evidence of a potential drop over the image intensifier 
#     i.e., the ground-cathode or cathode-MCP1 or MCP2-MCP23 voltage?
#
# (5) can channel depletion affect the whole MCP multiplier ? 
#
# (6) is the error in the background due to the readout streak contribution the cause? 
#

# 
#  tests
# =======
# (1) scale spectrum UVOT to STIS
# (2) determine flux of field sources 
# (3) what is WD magnitude, is there a galaxy in field?
# (4) check voltages
# (5) check image intesifier studies for high illumination
#     how much charge is deposited/extracted
#     compare to result (2)
# (6) small change in background should not affect bright star that much - compute 
#     need to know the real incident c/r to correct 
#

#
# Data
# ======
#
# T PYX 2011  - selected observation (UVOT + STIS)
# ---first observation-----------------------------------------
# offset uvot obs
# 00031973016 2011-05-07 JD2455689.0 MJD55688.5 
# AASVO B = 6.9, U = 7.3 
# field stars for comparison
# 136.1229833 -32.4484056  B1= 9.65 B2= 8.64 R1= 7.05 R2= 6.98 (second brightest *)
# 136.1825458 -32.4129944  B1=11.96 B2=11.68 R1=11.08 R2=11.05
# 136.2253458 -32.4456639  B1=11.70 B2=11.47 R1=10.95 R2=10.92
# STIS 
#
#--------------------------------------------------------------
#
# SN2011fe 
# ===========
#
# 
#
#
#
#

