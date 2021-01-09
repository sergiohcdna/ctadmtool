import matplotlib.pyplot as plt
import numpy as np
from ctaAnalysis.dmspectrum.dmflux import dmflux_anna

mass    = 1.e+4
z       = 0.018
channel = 'Tau'
emin    = 30.0
emax    = 0.9 * mass
#sigmav  = 3.6e-26
sigmav  = 1.e-28
jfactor = 9.31e+17

dmflux = dmflux_anna( sigmav , jfactor , mass , emin , emax , channel , z )

dphide = dmflux.flux()

#   By default, I computing the flux in units of GeV
#   Converting to MeV
dphide = dphide * 1.e-3

fig , ax = plt.subplots( figsize=( 9 , 6 ) )

ax.plot( dmflux.energy * 1.e+3 , dphide , color=( 0.82 , 0.10 , 0.26 ) ,\
    lw=4 , label='My dmclass' )

ax.set_xlim( 1.e+4 , 1.e+8 )
ax.set_ylim( 1.e-30 , 1.e-22 )
ax.set_xscale( 'log' )
ax.set_yscale( 'log' )
ax.set_xlabel( 'Energy $(MeV)$' )
ax.set_ylabel( '$\\frac{\\mathrm{d}\\Phi}{\\mathrm{d}E}$ ($cm^{-2} MeV^{-1} s^{-1}$)' )

ax.legend( loc='best' , prop={ 'size': 10 } )
pfile = 'PerseusFlux10TeVTau.png'
plt.savefig( pfile )
