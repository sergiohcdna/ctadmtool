import matplotlib.pyplot as plt
import numpy as np
from ctaAnalysis.dmspectrum.dmspectra import dmspectrum

masses   = ( 1.e+3 , 1.e+4 , 1.e+5 )
z        = 0.018
channel  = 'Tau'

#   Create instance of dmspectrum class with partial init
dmspec = dmspectrum.dminterp2d( channel , z )

#   pyplot
fig , ax = plt.subplots( figsize=( 9 , 6 ) )
#   get colors from cmap
cmap = plt.get_cmap( 'jet' )
colors = [ cmap( i ) for i in np.linspace( 0 , 1 , len( masses ) ) ]

#   Loop over masses to compute spectrum
for i , mass in enumerate( masses ) :

    #   Set value of mass and energy array
    energies      = np.logspace( np.log10( 30 ) , np.log10( mass ) , 100 )
    dmspec.mass   = mass
    dmspec.energy = energies

    #   Compute the value spectrum

    dnde   = dmspec.spectra()

    #   Plotting
    ax.plot( energies , dnde , color=colors[ i ] , lw=2 ,\
        label='Mass: {:.2e}GeV'.format( mass ) )

#   Making things nice :)
ax.set_xlim( 10 , 1.e+5 )
ax.set_ylim( 1.e-5 , 1.0 )
ax.set_xscale( 'log' )
ax.set_yscale( 'log' )
ax.set_xlabel( 'Energy (GeV)' )
ax.set_ylabel( '$\\frac{\\mathrm{d}N}{\\mathrm{d}E}$ (GeV$^{-1}$)' )

ax.legend( loc='best' , prop={ 'size': 10 } )
pfile = 'dmspectrumPerseus.png'
plt.savefig( pfile )
