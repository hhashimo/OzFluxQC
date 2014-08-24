import ast
import constants as c
import csv
import logging
import matplotlib.pyplot as plt
import numpy
import os
import platform
import qcio
import qcutils
import subprocess
import Tkinter

log = logging.getLogger('qc.rp')

def CalculateNEE(cf,ds):
    # now get a single series of NEE using the following rules:
    #  - when Fsd > 10 W/m2 (day time)
    #    - use Fc gap filled by chosen method
    #  - when Fsd < 10 W/m2 (night time)
    #    - use observed Fc (Fc_flag=0) when u* > ustar_threshold
    #    - use modelled Reco when Fc_flag!=0 or u* < ustar_threshold
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,'Fsd')
    ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,'ustar')
    Fc_label = str(cf['Derived']['NEE']['Fc'])
    Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,Fc_label)
    Reco_label = str(cf['Derived']['NEE']['Reco'])
    Reco,Reco_flag,Reco_attr = qcutils.GetSeriesasMA(ds,Reco_label)
    Fsd_threshold = float(cf['Params']['Fsd_threshold'])
    ustar_threshold = float(cf['Params']['ustar_threshold'])
    # create a series for NEE and the NEE QC flag
    nRecs = int(ds.globalattributes['nc_nrecs'])
    NEE = numpy.ma.array([c.missing_value]*nRecs,dtype=numpy.float64)
    NEE_flag = numpy.ma.array([0]*nRecs,dtype=numpy.int32)
    # fill the day time NEE
    index = numpy.ma.where(Fsd>Fsd_threshold)[0]
    NEE[index] = Fc[index]
    # fill the night time NEE when ustar is above the threshold
    index = numpy.ma.where((Fsd<=Fsd_threshold)&(ustar>ustar_threshold))[0]
    NEE[index] = Fc[index]
    # fill the night time NEE when ustar is below the threshold
    index = numpy.ma.where((Fsd<=Fsd_threshold)&(ustar<=ustar_threshold))[0]
    NEE[index] = Reco[index]
    # create the NEE series in the data structure
    units=qcutils.GetUnitsFromds(ds,Fc_label)
    attr = qcutils.MakeAttributeDictionary(long_name='Net ecosystem exchange',units=units)
    qcutils.CreateSeries(ds,'NEE',NEE,Flag=NEE_flag,Attr=attr)

def EstimateReco(cf,ds):
    # check to see if there is an Reco section in the control file
    section = qcutils.get_cfsection(cf,series='Reco',mode='quiet')
    if len(section)==0: return
    # step through the methods given in the control file
    method_list = cf[section]['Reco'].keys()
    for ThisMethod in method_list:
        if ThisMethod == 'LloydTaylor': EstimateRecoUsingLloydTaylor(cf,ds)
        if ThisMethod == 'SOLO': EstimateRecoUsingSOLO(cf,ds)

def EstimateRecoUsingLloydTaylor(cf,ds):
    # need to check that all information required is in the control file
    section = qcutils.get_cfsection(cf,series='Reco',mode='quiet')
    if len(section)==0:
        log.error('Reco section not found in control file')
        return
    if 'LloydTaylor' not in cf[section]['Reco']:
        log.error('LloydTaylor section not found in control file')
        return
    # get the driver
    if 'drivers' not in cf[section]['Reco']['LloydTaylor']:
        log.error('drivers key not found in LloydTaylor section of control file')
        return
    else:
        driver_list = eval(cf[section]['Reco']['LloydTaylor']['drivers'])
        driver_label = driver_list[0]
    # get the monthly values for the activation energy, E0
    if 'E0' not in cf[section]['Reco']['LloydTaylor']:
        log.error('E0 key not found in LloydTaylor section of control file')
        return
    else:
        E0_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['E0']))
    # get the monthly values for the base respiration, rb
    if 'rb' not in cf[section]['Reco']['LloydTaylor']:
        log.error('rb key not found in LloydTaylor section of control file')
        return
    else:
        rb_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['rb']))
    # get the output label
    if 'output' not in cf[section]['Reco']['LloydTaylor']:
        log.error('output key not found in LloydTaylor section of control file')
        return
    else:
        out_label = cf[section]['Reco']['LloydTaylor']['output']
    # ... and make an array of values for each month
    nRecs = int(ds.globalattributes['nc_nrecs'])
    E0 = numpy.ma.ones(nRecs)
    rb = numpy.ma.zeros(nRecs)
    month = ds.series['Month']['Data']
    lwr = numpy.max(numpy.array([numpy.min(month),1]))
    upr = numpy.min(numpy.array([numpy.max(month),12]))
    for m in range(lwr,upr+1):
        index = numpy.where(month==m)[0]
        E0[index] = E0_monthly[m-1]
        rb[index] = rb_monthly[m-1]
    # get the driver data
    Ts,flag,attr = qcutils.GetSeriesasMA(ds,driver_label)
    # estimate Reco using the Lloyd-Taylor expression
    t1 = 1/(c.Tref-c.T0)
    t2 = 1/(Ts-c.T0)
    Reco = rb*numpy.exp(E0*(t1-t2))
    # put the estimated Reco into the data structure
    units=qcutils.GetUnitsFromds(ds, 'Fc')
    attr = qcutils.MakeAttributeDictionary(long_name='Reco estimated using Lloyd-Taylor',units=units)
    qcutils.CreateSeries(ds,out_label,Reco,Flag=flag,Attr=attr)
    
def PartitionNEE(cf,ds):
    # check that there is a GPP section in the control file
    section = qcutils.get_cfsection(cf,series='GPP',mode='quiet')
    if len(section)==0:
        log.info('PartitionNEE: GPP section not found in control file')
        return
    # now check that the GPP section in the control file contains the
    # information needed.
    if 'NEE' not in cf[section]['GPP']:
        log.info('PartitionNEE: NEE key not in GPP section of control file')
        return
    else:
        NEE_label = str(cf[section]['GPP']['NEE'])
    if 'Reco' not in cf[section]['GPP']:
        log.info('PartitionNEE: Reco key not in GPP section of control file')
        return
    else:
        Reco_label = str(cf[section]['GPP']['Reco'])
    # now check that the requested series are in the data structure
    if NEE_label not in ds.series.keys():
        log.info('PartitionNEE: requested series '+NEE_label+' is not in the data structure')
        return
    else:
        NEE,f,a = qcutils.GetSeriesasMA(ds,NEE_label)
    if Reco_label not in ds.series.keys():
        log.info('PartitionNEE: requested series '+Reco_label+' is not in the data structure')
        return
    else:
        Reco,f,a = qcutils.GetSeriesasMA(ds,Reco_label)
    # check the units
    NEE_units=qcutils.GetUnitsFromds(ds,NEE_label)
    Reco_units=qcutils.GetUnitsFromds(ds,Reco_label)
    if NEE_units!=Reco_units:
        log.error('Units of NEE ('+NEE_label+') and Reco ('+Reco_label+') are different')
        return
    # ... at last, we can do the partitioning
    GEP = NEE - Reco
    # create the GEP QC flag
    nRecs = int(ds.globalattributes['nc_nrecs'])
    GEP_flag = numpy.ma.array([0]*nRecs,dtype=numpy.int32)
    # now check for non-zero values of GPP at night
    attr = qcutils.MakeAttributeDictionary(long_name='Gross ecosystem productivity',units=NEE_units)
    qcutils.CreateSeries(ds,'GEP',GEP,Flag=GEP_flag,Attr=attr)

def RecoUsingLloydTaylor(cf,ds):
    log.info('Estimating Reco using Lloyd-Taylor is not implemented yet')
    pass

def RecoUsingSOLO(cf,ds):
    """ Estimate Reco using SOLO. """
    # create a dictionary to hold the results and the input information
    rpSOLO_createdict(cf,ds,"Reco")
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    rpSOLO_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                 "file_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # set up the GUI
    rpSOLO_gui = Tkinter.Toplevel()
    rpSOLO_gui.wm_title("SOLO GUI (Reco)")
    rpSOLO_gui.grid()
    # top row
    nrow = 0
    rpSOLO_gui.nodesLabel = Tkinter.Label(rpSOLO_gui,text="Nodes")
    rpSOLO_gui.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpSOLO_gui.nodesEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpSOLO_gui.nodesEntry.insert(0,"10")
    rpSOLO_gui.trainingLabel = Tkinter.Label(rpSOLO_gui,text="Training")
    rpSOLO_gui.trainingLabel.grid(row=nrow,column=2,columnspan=1,sticky="E")
    rpSOLO_gui.trainingEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.trainingEntry.grid(row=nrow,column=3,columnspan=1,sticky="W")
    rpSOLO_gui.trainingEntry.insert(0,"500")
    rpSOLO_gui.factorLabel = Tkinter.Label(rpSOLO_gui,text="Nda factor")
    rpSOLO_gui.factorLabel.grid(row=nrow,column=4,columnspan=1,sticky="E")
    rpSOLO_gui.factorEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.factorEntry.grid(row=nrow,column=5,columnspan=1,sticky="W")
    rpSOLO_gui.factorEntry.insert(0,"5")
    # second row
    nrow = nrow + 1
    rpSOLO_gui.learningrateLabel = Tkinter.Label(rpSOLO_gui,text="Learning")
    rpSOLO_gui.learningrateLabel.grid(row=nrow,column=2,columnspan=1,sticky="E")
    rpSOLO_gui.learningrateEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.learningrateEntry.grid(row=nrow,column=3,columnspan=1,sticky="W")
    rpSOLO_gui.learningrateEntry.insert(0,"0.01")
    rpSOLO_gui.iterationsLabel = Tkinter.Label(rpSOLO_gui,text="Iterations")
    rpSOLO_gui.iterationsLabel.grid(row=nrow,column=4,columnspan=1,sticky="E")
    rpSOLO_gui.iterationsEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.iterationsEntry.grid(row=nrow,column=5,columnspan=1,sticky="W")
    rpSOLO_gui.iterationsEntry.insert(0,"500")
    # third row
    nrow = nrow + 1
    rpSOLO_gui.filestartLabel = Tkinter.Label(rpSOLO_gui,text="File start date")
    rpSOLO_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.fileendLabel = Tkinter.Label(rpSOLO_gui,text="File end date")
    rpSOLO_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    rpSOLO_gui.filestartValue = Tkinter.Label(rpSOLO_gui,text=str(ldt[0]))
    rpSOLO_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.fileendValue = Tkinter.Label(rpSOLO_gui,text=str(ldt[-1]))
    rpSOLO_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # fifth row
    nrow = nrow + 1
    rpSOLO_gui.startLabel = Tkinter.Label(rpSOLO_gui, text="Start date (YYYY-MM-DD)")
    rpSOLO_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.startEntry = Tkinter.Entry(rpSOLO_gui)
    rpSOLO_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # sixth row
    nrow = nrow + 1
    rpSOLO_gui.endLabel = Tkinter.Label(rpSOLO_gui, text="End date   (YYYY-MM-DD)")
    rpSOLO_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.endEntry = Tkinter.Entry(rpSOLO_gui)
    rpSOLO_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # seventh row
    nrow = nrow + 1
    rpSOLO_gui.peropt = Tkinter.IntVar()
    rpSOLO_gui.peropt.set(1)
    rpSOLO_gui.manualperiod = Tkinter.Radiobutton(rpSOLO_gui,text="Manual",variable=rpSOLO_gui.peropt,value=1)
    rpSOLO_gui.manualperiod.grid(row=nrow,column=0,columnspan=3,sticky="W")
    rpSOLO_gui.automonthly = Tkinter.Radiobutton(rpSOLO_gui,text="Monthly",variable=rpSOLO_gui.peropt,value=2)
    rpSOLO_gui.automonthly.grid(row=nrow,column=3,columnspan=3,sticky="W")
    # eigth row
    nrow = nrow + 1
    rpSOLO_gui.daysperiod = Tkinter.Radiobutton(rpSOLO_gui,text="No. days",variable=rpSOLO_gui.peropt,value=3)
    rpSOLO_gui.daysperiod.grid(row=nrow,column=0,sticky="W")
    rpSOLO_gui.daysentry = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.daysentry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpSOLO_gui.pointsperiod = Tkinter.Radiobutton(rpSOLO_gui,text="No. pts",variable=rpSOLO_gui.peropt,value=4)
    rpSOLO_gui.pointsperiod.grid(row=nrow,column=3,sticky="W")
    rpSOLO_gui.pointsentry = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.pointsentry.grid(row=nrow,column=4,columnspan=1,sticky="W")
    # ninth row
    nrow = nrow + 1
    rpSOLO_gui.minptsLabel = Tkinter.Label(rpSOLO_gui,text="Min points")
    rpSOLO_gui.minptsLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpSOLO_gui.minpts = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.minpts.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpSOLO_gui.minpts.insert(0,"200")
    rpSOLO_gui.owopt = Tkinter.IntVar()
    rpSOLO_gui.owopt.set(1)
    rpSOLO_gui.overwrite = Tkinter.Checkbutton(rpSOLO_gui, text="Overwrite", variable=rpSOLO_gui.owopt)
    rpSOLO_gui.overwrite.grid(row=nrow,column=3,columnspan=2,sticky="w")
    # tenth row
    nrow = nrow + 1
    rpSOLO_gui.doneButton = Tkinter.Button (rpSOLO_gui, text="Done",command=lambda:rpSOLO_done(rpSOLO_gui))
    rpSOLO_gui.doneButton.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.runButton = Tkinter.Button (rpSOLO_gui, text="Run",command=lambda:rpSOLO_run(ds,rpSOLO_gui,rpSOLO_info))
    rpSOLO_gui.runButton.grid(row=nrow,column=3,columnspan=3)
    # eleventh row
    nrow = nrow + 1
    rpSOLO_gui.progress_row = nrow
    rpSOLO_gui.progress = Tkinter.Label(rpSOLO_gui, text='Waiting for input ...')
    rpSOLO_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    rpSOLO_gui.wait_window(rpSOLO_gui)

def rpSOLO_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the SOLO data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("rpSOLO_createdict: Series "+series+" not found in control file, skipping ...")
        return
    # create the solo directory in the data structure
    if "solo" not in dir(ds): ds.solo = {}
    # create the dictionary keys for this series
    ds.solo[series] = {}
    # site name
    ds.solo[series]["site_name"] = ds.globalattributes["site_name"]
    # list of drivers
    ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["RecoUsingSOLO"]["drivers"])
    # name of SOLO output series in ds
    ds.solo[series]["output"] = cf[section][series]["RecoUsingSOLO"]["output"]
    # results of best fit for plotting later on
    ds.solo[series]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                  "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                  "Avg (obs)":[],"Avg (SOLO)":[],
                                  "Var (obs)":[],"Var (SOLO)":[],"Var ratio":[],
                                  "m_ols":[],"b_ols":[]}
    # create an empty series in ds if the SOLO output series doesn't exist yet
    if ds.solo[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.solo[series]["output"])
        qcutils.CreateSeries(ds,ds.solo[series]["output"],data,Flag=flag,Attr=attr)

def rpSOLO_done(rpSOLO_gui):
    # destroy the SOLO GUI
    rpSOLO_gui.destroy()

def rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info):
    """
    This is the main routine for running SOLO, an artifical neural network for estimating Reco.
    """
    startdate = rpSOLO_info["startdate"]
    enddate = rpSOLO_info["enddate"]
    log.info(" Estimating Reco using SOLO: "+startdate+" to "+enddate)
    # read the control file again, this allows the contents of the control file to
    # be changed with the SOLO GUI still displayed
    cfname = ds.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname,mode="quiet")
    solo_series = cf["Variables"].keys()
    for series in solo_series:
        section = qcutils.get_cfsection(cf,series=series,mode="quiet")
        if len(section)==0: continue
        if series not in ds.series.keys(): continue
        ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["RecoUsingSOLO"]["drivers"])
        ds.solo[series]["output"] = cf[section][series]["RecoUsingSOLO"]["output"]
    # get some useful things
    site_name = ds.globalattributes["site_name"]
    # get the time step and a local pointer to the datetime series
    ts = ds.globalattributes["time_step"]
    ldt = ds.series["DateTime"]["Data"]
    xldt = ds.series["xlDateTime"]["Data"]
    # get the start and end datetime indices
    si = qcutils.GetDateIndex(ldt,startdate,ts=ts,default=0,match="exact")
    ei = qcutils.GetDateIndex(ldt,enddate,ts=ts,default=-1,match="exact")
    # check the start and end indices
    if si >= ei:
        log.error(" RecoUsingSOLO: end datetime index ("+str(ei)+") smaller that start ("+str(si)+")")
        return
    if si==0 and ei==-1:
        log.error(" RecoUsingSOLO: no start and end datetime specified, using all data")
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        nRecs = ei - si + 1
    # loop over the series to be gap filled using solo
    # close any open plot windows
    if len(plt.get_fignums())!=0:
        for i in plt.get_fignums(): plt.close(i)
    fig_num = 0
    for series in solo_series:
        ds.solo[series]["results"]["startdate"].append(xldt[si])
        ds.solo[series]["results"]["enddate"].append(xldt[ei])
        d,f,a = qcutils.GetSeriesasMA(ds,series,si=si,ei=ei)
        if numpy.ma.count(d)<rpSOLO_info["min_points"]:
            log.error("rpSOLO: Less than "+str(rpSOLO_info["min_points"])+" points available for series "+series+" ...")
            ds.solo[series]["results"]["No. points"].append(float(0))
            results_list = ds.solo[series]["results"].keys()
            for item in ["startdate","enddate","No. points"]:
                if item in results_list: results_list.remove(item)
            for item in results_list:
                ds.solo[series]["results"][item].append(float(c.missing_value))
            continue
        drivers = ds.solo[series]["drivers"]
        output = ds.solo[series]["output"]
        # set the number of nodes for the inf files
        nodesAuto = rpSOLO_setnodesEntry(rpSOLO_gui,drivers,default=10)
        # write the inf files for sofm, solo and seqsolo
        rpSOLO_writeinffiles(rpSOLO_gui)
        # run SOFM
        result = rpSOLO_runsofm(ds,rpSOLO_gui,drivers,series,nRecs,si=si,ei=ei)
        if result!=1: return
        # run SOLO
        result = rpSOLO_runsolo(ds,drivers,series,nRecs,si=si,ei=ei)
        if result!=1: return
        # run SEQSOLO and put the SOLO data into the data structure
        result = rpSOLO_runseqsolo(ds,drivers,series,output,nRecs,si=si,ei=ei)
        if result!=1: return
        # plot the results
        fig_num = fig_num + 1
        title = site_name+' : Comparison of tower and SOLO data for '+series
        pd = gfSOLO_initplot(site_name=site_name,label=series,fig_num=fig_num,title=title,
                             nDrivers=len(drivers))
        gfSOLO_plot(pd,dsa,dsb,drivers,series,output,solo_gui,si=si,ei=ei)
        # reset the nodesEntry in the solo_gui
        if nodesAuto: rpSOLO_resetnodesEntry(rpSOLO_gui)
    if 'RecoUsingSOLO' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RecoUsingSOLO'

def rpSOLO_progress(rpSOLO_gui,text):
    """
        Update progress message in SOLO GUI
        """
    rpSOLO_gui.progress.destroy()
    rpSOLO_gui.progress = Tkinter.Label(rpSOLO_gui, text=text)
    rpSOLO_gui.progress.grid(row=rpSOLO_gui.progress_row,column=0,columnspan=6,sticky="W")
    rpSOLO_gui.update()

def rpSOLO_resetnodesEntry(rpSOLO_gui):
    rpSOLO_gui.nodesEntry.delete(0,Tkinter.END)
    rpSOLO_gui.nodesEntry.insert(0,"10")

def rpSOLO_run(ds,rpSOLO_gui,rpSOLO_info):
    # populate the solo_info dictionary with things that will be useful
    rpSOLO_info["peropt"] = rpSOLO_gui.peropt.get()
    rpSOLO_info["min_points"] = int(rpSOLO_gui.minpts.get())
    rpSOLO_info["site_name"] = ds.globalattributes["site_name"]
    rpSOLO_info["time_step"] = int(ds.globalattributes["time_step"])
    rpSOLO_info["nperhr"] = int(float(60)/rpSOLO_info["time_step"]+0.5)
    rpSOLO_info["nperday"] = int(float(24)*rpSOLO_info["nperhr"]+0.5)
    rpSOLO_info["maxlags"] = int(float(12)*rpSOLO_info["nperhr"]+0.5)
    rpSOLO_info["tower"] = {}
    rpSOLO_info["access"] = {}
    #log.info(" Estimating Reco using SOLO")
    if rpSOLO_gui.peropt.get()==1:
        rpSOLO_progress(rpSOLO_gui,"Starting manual run ...")
        # get the start and end datetimes entered in the SOLO GUI
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        rpSOLO_info["enddate"] = rpSOLO_gui.endEntry.get()
        if len(rpSOLO_info["enddate"])==0: rpSOLO_info["enddate"] = rpSOLO_info["file_enddate"]
        rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
        rpSOLO_progress(rpSOLO_gui,"Finished manual run ...")
    elif rpSOLO_gui.peropt.get()==2:
        rpSOLO_progress(rpSOLO_gui,"Starting auto (monthly) run ...")
        # get the start datetime entered in the SOLO GUI
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        startdate = dateutil.parser.parse(rpSOLO_info["startdate"])
        file_startdate = dateutil.parser.parse(rpSOLO_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpSOLO_info["file_enddate"])
        enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
        enddate = min([file_enddate,enddate])
        rpSOLO_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
            rpSOLO_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpSOLO_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ## plot the summary statistics
        #gfSOLO_plotsummary(dsb)
        rpSOLO_progress(rpSOLO_gui,"Finished auto (monthly) run ...")
    elif rpSOLO_gui.peropt.get()==3:
        rpSOLO_progress(rpSOLO_gui,"Starting auto (days) run ...")
        # get the start datetime entered in the SOLO GUI
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        startdate = dateutil.parser.parse(rpSOLO_info["startdate"])
        file_startdate = dateutil.parser.parse(rpSOLO_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpSOLO_info["file_enddate"])
        nDays = int(rpSOLO_gui.daysentry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
        enddate = min([file_enddate,enddate])
        rpSOLO_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
            rpSOLO_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpSOLO_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ## plot the summary statistics
        #gfSOLO_plotsummary(dsb)
        rpSOLO_progress(rpSOLO_gui,"Finished auto (days) run ...")
    elif rpSOLO_gui.peropt.get()==4:
        pass

def rpSOLO_runsofm(ds,rpSOLO_gui,driverlist,targetlabel,nRecs,si=0,ei=-1):
    '''
    Run sofm, the pre-processor for SOLO.
    '''
    # get the number of drivers
    ndrivers = len(driverlist)
    # add an extra column for the target data
    sofminputdata = numpy.zeros((nRecs,ndrivers))
    # now fill the driver data array
    i = 0
    badlines = []
    for TheseOnes in driverlist:
        driver,flag,attr = qcutils.GetSeries(ds,TheseOnes,si=si,ei=ei)
        index = numpy.where(abs(driver-float(c.missing_value)<c.eps))[0]
        if len(index)!=0:
            log.error(' SOLO_runsofm: c.missing_value found in driver '+TheseOnes+' at lines '+str(index))
            badlines = badlines+index.tolist()
        sofminputdata[:,i] = driver[:]
        i = i + 1
    if len(badlines)!=0:
        nBad = len(badlines)
        goodlines = [x for x in range(0,nRecs) if x not in badlines]
        sofminputdata = sofminputdata[goodlines,:]
        log.info(' SOLO_runsofm: removed '+str(nBad)+' lines from sofm input file')
        nRecs = len(goodlines)
    # now write the drivers to the SOFM input file
    sofmfile = open('solo/input/sofm_input.csv','wb')
    wr = csv.writer(sofmfile,delimiter=',')
    for i in range(sofminputdata.shape[0]):
        wr.writerow(sofminputdata[i,0:ndrivers])
    sofmfile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/sofm_4.out'): os.remove('solo/output/sofm_4.out')
    # now run SOFM
    sofmlogfile = open('solo/log/sofm.log','wb')
    if platform.system()=="Windows":
        subprocess.call(['./solo/bin/sofm.exe','solo/inf/sofm.inf'],stdout=sofmlogfile)
    else:
        subprocess.call(['./solo/bin/sofm','solo/inf/sofm.inf'],stdout=sofmlogfile)
    sofmlogfile.close()
    # check to see if the sofm output file exists, this is used to indicate that sofm ran correctly
    if os.path.exists('solo/output/sofm_4.out'):
        return 1
    else:
        log.error(' SOLO_runsofm: SOFM did not run correctly, check the GUI and the log files')
        return 0

def rpSOLO_setnodesEntry(rpSOLO_gui,drivers,default=2):
    nodesAuto = False
    if str(rpSOLO_gui.nodesEntry.get()).lower()=="auto":
        nodesAuto = True
        rpSOLO_gui.nodesEntry.delete(0,Tkinter.END)
        num_nodes = max([len(drivers)+1,default])
        rpSOLO_gui.nodesEntry.insert(0,str(num_nodes))
    return nodesAuto

def rpSOLO_writeinffiles(rpSOLO_gui):
    # sofm inf file
    f = open('solo/inf/sofm.inf','w')
    f.write(str(rpSOLO_gui.nodesEntry.get())+'\n')
    f.write(str(rpSOLO_gui.trainingEntry.get())+'\n')
    f.write(str(20)+'\n')
    f.write(str(0.01)+'\n')
    f.write(str(1234)+'\n')
    f.write('solo/input/sofm_input.csv'+'\n')
    f.write('solo/output/sofm_1.out'+'\n')
    f.write('solo/output/sofm_2.out'+'\n')
    f.write('solo/output/sofm_3.out'+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write(str(50)+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: No. of training iterations - default is 500 (changeable via GUI if used)\n')
    f.write('Line 3: No. of iterations per screen output - default is 20\n')
    f.write('Line 4: Spacing between initial weights - default is 0.01\n')
    f.write('Line 5: Seed for random number generator - default is 1234\n')
    f.write('Line 6: input data filename with path relative to current directory\n')
    f.write('Line 7: first output filename with path relative to current directory\n')
    f.write('Line 8: second output filename with path relative to current directory\n')
    f.write('Line 9: third output filename with path relative to current directory\n')
    f.write('Line 10: fourth output filename with path relative to current directory (used by SOLO)\n')
    f.write('Line 11: No. iterations per write of weights to screen - default is 50\n')
    f.close()
    # solo inf file
    f = open('solo/inf/solo.inf','w')
    f.write(str(rpSOLO_gui.nodesEntry.get())+'\n')
    f.write(str(rpSOLO_gui.factorEntry.get())+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write('solo/input/solo_input.csv'+'\n')
    f.write('training'+'\n')
    f.write(str(5678)+'\n')
    f.write(str(0)+'\n')
    f.write('solo/output/eigenValue.out'+'\n')
    f.write('solo/output/eigenVector.out'+'\n')
    f.write('solo/output/accumErr.out'+'\n')
    f.write('solo/output/accumRR.out'+'\n')
    f.write('solo/output/trainProcess.out'+'\n')
    f.write('solo/output/freqTable.out'+'\n')
    f.write('solo/output/hidOutputWt.out'+'\n')
    f.write('solo/output/errorMap.out'+'\n')
    f.write('solo/output/finResult.out'+'\n')
    f.write('solo/output/trainWin.out'+'\n')
    f.write('solo/output/trainWout.out'+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: multiplier for minimum number of points per node (NdaFactor) - default is 5 (ie 5*(no. of drivers+1) (changeable via GUI if used)\n')
    f.write('Line 3: fourth output file from SOFM, used as input to SOLO\n')
    f.write('Line 4: input data filename with path relative to current directory\n')
    f.write('Line 5: type of run ("training" or "simulation", always "training" for SOLO)\n')
    f.write('Line 6: seed for random number generator - default is 5678\n')
    f.write('Line 7: "calThreshold", not used by SOLO\n')
    f.write('Lines 8 to 18: output files from SOLO with path relative to current directory\n')
    f.close()
    # seqsolo inf file
    f = open('solo/inf/seqsolo.inf','w')
    f.write(str(rpSOLO_gui.nodesEntry.get())+'\n')
    f.write(str(0)+'\n')
    f.write(str(rpSOLO_gui.learningrateEntry.get())+'\n')
    f.write(str(rpSOLO_gui.iterationsEntry.get())+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write('solo/input/seqsolo_input.csv'+'\n')
    f.write('simulation'+'\n')
    f.write(str(9100)+'\n')
    f.write(str(0)+'\n')
    f.write('solo/output/eigenValue.out'+'\n')
    f.write('solo/output/eigenVector.out'+'\n')
    f.write('solo/output/trainWout.out'+'\n')
    f.write('solo/output/freqTable.out'+'\n')
    f.write('solo/output/errorMap.out'+'\n')
    f.write('solo/output/finResult.out'+'\n')
    f.write('solo/output/trainingRMSE.out'+'\n')
    f.write('solo/output/seqOut0.out'+'\n')
    f.write('solo/output/seqOut1.out'+'\n')
    f.write('solo/output/seqOut2.out'+'\n')
    f.write('solo/output/seqHidOutW.out'+'\n')
    f.write('solo/output/seqFreqMap.out'+'\n')
    f.write(str(c.missing_value)+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: NdaFactor - not used by SEQSOLO, default value is 0\n')
    f.write('Line 3: learning rate - default value 0.01 (must be between 0.0 1nd 1.0, changeable via GUI if used)\n')
    f.write('Line 4: number of iterations for sequential training, default value is 500 (changeable via GUI if used)\n')
    f.write('Line 5: fourth output file from SOFM, used as input file by SEQSOLO\n')
    f.write('Line 6: input data filename with path relative to current directory\n')
    f.write('Line 7: type of run ("training" or "simulation", always "simulation" for SEQSOLO)\n')
    f.write('Line 8: seed for random number generator - default is 9100\n')
    f.write('Line 9: "calThreshold" - minimum number of data points for SOLO node to be used in simulation, default value is 0 (use all nodes)\n')
    f.write('Lines 10 to 21: output files from SEQSOLO with path relative to current directory\n')
    f.write('Line 22: missing data value, default value is c.missing_value.0\n')
    f.close()
