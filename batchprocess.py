import logging
import sys
sys.path.append('scripts')
import qcio
import qcls
import qcutils

log = qcutils.startlog('batch','logfiles/batch.log')

# get the batch processing control file
cf_batch = qcio.load_controlfile(path='controlfiles')
if len(cf_batch)==0: sys.exit()
for i in cf_batch['Files'].keys():
    cfname = cf_batch['Files'][i]
    if 'L1' in cfname:
        # L1 processing
        log.info('Starting L1 processing with '+cfname)
        cf = qcio.get_controlfilecontents(cfname)
        qcio.xl2nc(cf,'L1')
        log.info('Finished L1 processing with '+cfname)
    elif 'L2' in cfname:
        # L2 processing
        log.info('Starting L2 processing with '+cfname)
        cf = qcio.get_controlfilecontents(cfname)
        infilename = qcio.get_infilenamefromcf(cf)
        ds1 = qcio.nc_read_series(infilename)
        ds2 = qcls.l2qc(cf,ds1)
        outfilename = qcio.get_outfilenamefromcf(cf)
        ncFile = qcio.nc_open_write(outfilename)
        qcio.nc_write_series(ncFile,ds2)
        log.info('Finished L2 processing with '+cfname)
    elif 'L3' in cfname:
        # L3 processing
        log.info('Starting L3 processing with '+cfname)
        cf = qcio.get_controlfilecontents(cfname)
        infilename = qcio.get_infilenamefromcf(cf)
        ds2 = qcio.nc_read_series(infilename)
        ds3 = qcls.l3qc(cf,ds2)
        outfilename = qcio.get_outfilenamefromcf(cf)
        outputlist = qcio.get_outputlistfromcf(cf,'nc')
        ncFile = qcio.nc_open_write(outfilename)
        qcio.nc_write_series(ncFile,ds3,outputlist=outputlist)
        log.info('Finished L3 processing with '+cfname)
    