import gammalib
import numpy as np
import argparse
from tqdm import tqdm

from ctaAnalysis.dmspectrum.dmspectra import dmspectrum

if __name__ == '__main__':

    #   List of available models
    eblmodels = ['franceschini','franceschini2017','kneiske', 'finke',
        'dominguez', 'dominguez-upper', 'dominguez-lower','inuoe',
        'inuoe-low-pop3', 'inuoe-up-pop3','gilmore', 'gilmore-fixed']

    # This is the dictionary of channels available in PPPC4DMID tables
    channels = {0:'eL', 1:'eR', 2:'e', 3:'MuL', 4:'MuR',
        5:'Mu', 6:'TauL', 7:'TauR', 8:'Tau', 9:'q',
        10:'c', 11:'b', 12:'t', 13:'WL', 14:'WT',
        15:'W', 16:'ZL', 17:'ZT', 18:'Z', 19:'g',
        20:'Gamma', 21:'h', 22:'Nue', 23:'NuMu', 24:'NuTau',
        25:'Ve', 26:'VMu', 27:'VTau'}

    msg = 'Create a fits using GModelSpectralTable'
    options = argparse.ArgumentParser(description=msg)

    source  = options.add_argument_group('Source',
        'Information about the source')
    source.add_argument('--srcname',
        help='Name of the source',
        type=str, required=False,
        default='Perseus', metavar='Perseus')

    model = options.add_argument_group('Model',
        'Model parameters')
    model.add_argument('--jfactor',
        help='Astrophysical Factor (GeV**2/cm**5)',
        type=float, required=True, metavar='1.0e+19')
    model.add_argument('--sigmav',
        help='Annihilation cross-section (cm**3/s)',
        type=float, required=True, metavar='3.0e-26')
    model.add_argument('--z', help='Redshift',
        type=float, required=True, metavar='0.01')
    model.add_argument('--eblmodel', help='EBL attenuation model',
        type=str, default='dominguez',
        metavar='dominguez', choices=eblmodels)
    model.add_argument('--mmin', help='Minimum Mass (GeV)',
        type=float, default=100.0, metavar='100.0')
    model.add_argument('--mmax', help='Maximum Mass (GeV)',
        type=float, default=1.0e+5, metavar='1.0e+5')
    model.add_argument('--mpoints', help='Number of mass points',
        type=int, required=False, default=50, metavar='50')

    esetup  = options.add_argument_group('Energy setup',
        'Min, Max and number of bins')
    esetup.add_argument('--emin', help='Minimum Energy (GeV)',
        type=float, required=False, default=50.0, metavar='50.0')
    esetup.add_argument('--emax', help='Maximum Energy (GeV)',
        type=float, required=False, default=1.0e+5, metavar='1.0e+5')
    esetup.add_argument('--nebins', help='Number of energy bins',
        type=int, required=False, default=250, metavar='250')

    args = options.parse_args()

    print('****************************************')
    print('Energy setup:\n')

    # Array with definitions of energy bins
    temin = gammalib.GEnergy(args.emin, 'GeV')
    temax = gammalib.GEnergy(args.emax, 'GeV')
    ebins = gammalib.GEbounds(args.nebins, temin, temax)

    # I need a numpy array for the energies
    # when creating an instance of dmspectrum class
    energies = np.zeros((args.nebins))
    for i in range(args.nebins):
        energies[i] = ebins.elogmean(i).GeV()

    print('****************************************')
    print('Creating Model Parameters\n')
    #   Then create the GModelPar objects for mass and channel
    dmmass    = gammalib.GModelPar('Mass', 1.0e+3)
    dmmass.unit('GeV')
    dmchannel = gammalib.GModelPar('Channel', 8)

    #   Create the GSpectralTablePar objects
    k_ch   = [k for k in channels.keys()]
    w      = (np.log10(args.mmax) - np.log10(args.mmin)) / args.mpoints
    masses = [np.ceil(10**(np.log10(args.mmin)+x*w)) for x in range(args.mpoints+1)]

    par_mass    = gammalib.GModelSpectralTablePar(dmmass, masses)
    par_channel = gammalib.GModelSpectralTablePar(dmchannel, k_ch)

    #   Create the container GSpectralTablePars and append the pars
    pars = gammalib.GModelSpectralTablePars()
    pars.append(par_mass)
    pars.append(par_channel)

    print('****************************************')
    print('Filling the spectrum table\n')
    #   GNdarray to save the spectra
    spectra = gammalib.GNdarray(len(masses), len(k_ch), args.nebins)

    #   filling the spectrum
    for index, mass in tqdm(enumerate(masses)):
        for cindex, thisch in channels.items():
            #    Create an instance of dmspectrum
            dminterp = dmspectrum(mass, energies, thisch,
                args.z, eblmod=args.eblmodel)
            spec     = dminterp.spectra()
            for eindex in range(args.nebins):
                spectra[index, cindex, eindex] = spec[eindex]
                
            del dminterp

    print('****************************************')
    print('Saving the file\n')
    #   Saving the model
    fname = 'DMModelAnnihilation{0}.fits'.format(args.srcname)
    model = gammalib.GModelSpectralTable(ebins, pars, spectra)
    model.table_par('Mass').method(0)
    model.table_par('Channel').method(0)
    model['Mass'].fix()
    model['Channel'].fix()
    model.save(fname, True)

    print('****************************************')
    print('This is the END!\n')