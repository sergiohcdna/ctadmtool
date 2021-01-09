import ctools
from csdmatter import csdmatter

inmodel = 'XMLTemplates/Perseus_DMPS_anna_1000GeV_b.xml'
evcube  = 'eventsPerseus.fits'

dm_anna = csdmatter()

dm_anna[ 'inobs' ] =  evcube
dm_anna[ 'inmodel' ] = inmodel
dm_anna[ 'srcname' ] = 'Perseus'
dm_anna[ 'caldb' ]   = 'prod3b-v2'
dm_anna[ 'irf' ]     = 'North_z20_S_5h'
dm_anna[ 'dmass' ]   = 1.0
dm_anna[ 'outfile' ] = 'Perseus_dmresults.fits'

dm_anna.logFileOpen( clobber=False )
dm_anna.execute()

print( 'END :)' )

