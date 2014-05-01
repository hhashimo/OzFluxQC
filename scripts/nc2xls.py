import logging
import qcio
import qcutils
import sys
import xlwt

log = qcutils.startlog('nc2xls','../logfiles/nc2xls.log')

ncfullname = qcio.get_filename_dialog(path='.',title='Choose a netCDF file to open')
if len(ncfullname)==0: sys.exit()
xlfullname= ncfullname.replace('.nc','.xls')
# read the netCDF file
log.info(' Opening and reading netCDF file '+ncfullname)
ds = qcio.nc_read_series(ncfullname)
# open the Excel file
log.info(' Opening and writing Excel file '+xlfullname)
xlfile = xlwt.Workbook()
# add sheets to the Excel file
xlAttrSheet = xlfile.add_sheet('Attr')
# write the global attributes
log.info(' Writing the global attributes to Excel file '+xlfullname)
xlcol = 0
xlrow = 0
xlAttrSheet.write(xlrow,xlcol,'Global attributes')
xlrow = xlrow + 1
globalattrlist = ds.globalattributes.keys()
globalattrlist.sort()
for ThisOne in [x for x in globalattrlist if 'Flag' not in x]:
    xlAttrSheet.write(xlrow,xlcol,ThisOne)
    xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
    xlrow = xlrow + 1
for ThisOne in [x for x in globalattrlist if 'Flag' in x]:
    xlAttrSheet.write(xlrow,xlcol,ThisOne)
    xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
    xlrow = xlrow + 1
# write the variable attributes
log.info(' Writing the variable attributes to Excel file '+xlfullname)
xlrow = xlrow + 1
xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
xlrow = xlrow + 1
xlcol_varname = 0
xlcol_attrname = 1
xlcol_attrvalue = 2
variablelist = ds.series.keys()
variablelist.sort()
for ThisOne in ["DateTime","DateTime_UTC"]:
    if ThisOne in variablelist: variablelist.remove(ThisOne)
for ThisOne in variablelist:
    xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
    attributelist = ds.series[ThisOne]['Attr'].keys()
    attributelist.sort()
    for Attr in attributelist:
        xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
        xlAttrSheet.write(xlrow,xlcol_attrvalue,str(ds.series[ThisOne]['Attr'][Attr]))
        xlrow = xlrow + 1
# get the number of records in the file
nRecs = qcutils.get_nrecs(ds)
# add the Data and Flag sheets
xlDataSheet = xlfile.add_sheet('Data')
xlFlagSheet = xlfile.add_sheet('Flag')
# write the Excel date/time to the data and the QC flags as the first column
if 'xlDateTime' in ds.series.keys():
    log.info(' Writing the datetime to Excel file '+xlfullname)
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    xlDataSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
        xlFlagSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
    # remove xlDateTime from the list of variables to be written to the Excel file
    variablelist.remove('xlDateTime')
# now start looping over the other variables in the xl file
xlcol = xlcol + 1
# loop over variables to be output to xl file
for ThisOne in variablelist:
    # put up a progress message
    log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
    # write the units and the variable name to the header rows in the xl file
    attrlist = ds.series[ThisOne]['Attr'].keys()
    if 'long_name' in attrlist:
        longname = ds.series[ThisOne]['Attr']['long_name']
    elif 'Description' in attrlist:
        longname = ds.series[ThisOne]['Attr']['Description']
    else:
        longname = None
    if 'units' in attrlist:
        units = ds.series[ThisOne]['Attr']['units']
    elif 'Units' in attrlist:
        units = ds.series[ThisOne]['Attr']['Units']
    else:
        units = None
    xlDataSheet.write(0,xlcol,longname)
    xlDataSheet.write(1,xlcol,units)
    xlDataSheet.write(2,xlcol,ThisOne)
    # loop over the values in the variable series (array writes don't seem to work)
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,float(ds.series[ThisOne]['Data'][j]))
    # check to see if this variable has a quality control flag
    if 'Flag' in ds.series[ThisOne].keys():
        # write the QC flag name to the xk file
        xlFlagSheet.write(2,xlcol,ThisOne)
        # specify the format of the QC flag (integer)
        d_xf = xlwt.easyxf(num_format_str='0')
        # loop over QV flag values and write to xl file
        for j in range(nRecs):
            xlFlagSheet.write(j+3,xlcol,int(ds.series[ThisOne]['Flag'][j]),d_xf)
    # increment the column pointer
    xlcol = xlcol + 1

xlfile.save(xlfullname)
