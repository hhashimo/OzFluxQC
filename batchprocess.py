import ast
import datetime
import logging
import os
import sys
sys.path.append('scripts')
import time
import qcclim
import qccpd
import qcio
import qcls
import qcutils

#log = qcutils.startlog('batch','logfiles/batch.log')
t = time.localtime()
rundatetime = datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]).strftime("%Y%m%d%H%M")
log_filename = 'logfiles/batchprocess_'+rundatetime+'.log'
logging.basicConfig(filename=log_filename,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt = '%H:%M:%S',
                    level=logging.DEBUG)
console = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
console.setFormatter(formatter)
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

# get the batch processing control file
if len(sys.argv)==1:
    cf_batch = qcio.load_controlfile(path='controlfiles')
    if len(cf_batch)==0: sys.exit()
else:
    cfname = sys.argv[1]
    if os.path.exists(cfname):
        cf_batch = qcio.get_controlfilecontents(cfname)
    else:
        logging.error("Control file "+cfname+" does not exist")
        sys.exit()

level_list = ['L1','L2','L3','concatenate','climatology','cpd','L4','L5','L6']
if "Options" in cf_batch:
    if "levels" in cf_batch["Options"]: level_list = ast.literal_eval(cf_batch["Options"]["levels"])
for level in level_list:
    if level.lower()=="l1":
        # L1 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L1 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            qcio.xl2nc(cf,'L1')
            logging.info('Finished L1 processing with '+cfname)
            logging.info('')
    elif level.lower()=="l2":
        # L2 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L2 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            infilename = qcio.get_infilenamefromcf(cf)
            ds1 = qcio.nc_read_series(infilename)
            ds2 = qcls.l2qc(cf,ds1)
            outfilename = qcio.get_outfilenamefromcf(cf)
            ncFile = qcio.nc_open_write(outfilename)
            qcio.nc_write_series(ncFile,ds2)
            logging.info('Finished L2 processing with '+cfname)
            logging.info('')
    elif level.lower()=="l3":
        # L3 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L3 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            infilename = qcio.get_infilenamefromcf(cf)
            ds2 = qcio.nc_read_series(infilename)
            ds3 = qcls.l3qc(cf,ds2)
            outfilename = qcio.get_outfilenamefromcf(cf)
            outputlist = qcio.get_outputlistfromcf(cf,'nc')
            ncFile = qcio.nc_open_write(outfilename)
            qcio.nc_write_series(ncFile,ds3,outputlist=outputlist)
            logging.info('Finished L3 processing with '+cfname)
            logging.info('')
    elif level.lower()=="fluxnet":
        # convert netCDF files to FluxNet CSV files
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting FluxNet output with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            qcio.fn_write_csv(cf)
            logging.info('Finished FluxNet output with '+cfname)
            logging.info('')
    elif level.lower()=="concatenate":
        # concatenate netCDF files
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting concatenation with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            qcio.nc_concatenate(cf)
            logging.info('Finished concatenation with '+cfname)
            logging.info('')
    elif level.lower()=="climatology":
        # climatology
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting climatology with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            qcclim.climatology(cf)
            logging.info('Finished climatology with '+cfname)
            logging.info('')
    elif level.lower()=="cpd":
        # ustar threshold from change point detection
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting CPD with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            if "Options" not in cf: cf["Options"]={}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = False
            qccpd.cpd_main(cf)
            logging.info('Finished CPD with '+cfname)
            logging.info('')
    elif level.lower()=="l4":
        # L4 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L4 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            if "Options" not in cf: cf["Options"]={}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = False
            infilename = qcio.get_infilenamefromcf(cf)
            ds3 = qcio.nc_read_series(infilename)
            ds4 = qcls.l4qc(cf,ds3)
            outfilename = qcio.get_outfilenamefromcf(cf)
            outputlist = qcio.get_outputlistfromcf(cf,'nc')
            ncFile = qcio.nc_open_write(outfilename)
            qcio.nc_write_series(ncFile,ds4,outputlist=outputlist)
            logging.info('Finished L4 processing with '+cfname)
            logging.info('')
    elif level.lower()=="l5":
        # L5 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L5 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            if "Options" not in cf: cf["Options"]={}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = False
            infilename = qcio.get_infilenamefromcf(cf)
            ds4 = qcio.nc_read_series(infilename)
            ds5 = qcls.l5qc(cf,ds4)
            outfilename = qcio.get_outfilenamefromcf(cf)
            outputlist = qcio.get_outputlistfromcf(cf,'nc')
            ncFile = qcio.nc_open_write(outfilename)
            qcio.nc_write_series(ncFile,ds5,outputlist=outputlist)
            logging.info('Finished L5 processing with '+cfname)
            logging.info('')
    elif level.lower()=="l6":
        # L6 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            logging.info('Starting L6 processing with '+cfname)
            cf = qcio.get_controlfilecontents(cfname)
            if "Options" not in cf: cf["Options"]={}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = False
            infilename = qcio.get_infilenamefromcf(cf)
            ds5 = qcio.nc_read_series(infilename)
            ds6 = qcls.l6qc(cf,ds5)
            outfilename = qcio.get_outfilenamefromcf(cf)
            outputlist = qcio.get_outputlistfromcf(cf,'nc')
            ncFile = qcio.nc_open_write(outfilename)
            qcio.nc_write_series(ncFile,ds6,outputlist=outputlist)
            logging.info('Finished L6 processing with '+cfname)
            logging.info('')
