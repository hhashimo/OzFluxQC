import ast
import constants as c
import csv
import datetime
import dateutil
import logging
import matplotlib.pyplot as plt
import numpy
import os
import platform
import qcio
import qcutils
import subprocess
import sys
import Tkinter

log = logging.getLogger('qc.rp')
# lets see if ffnet is installed
try:
    import ffnet
except ImportError:
    #log.error("FreUsingFFNET: Unable to import module ffnet")
    pass

def CalculateNEE(cf,ds):
    """
    Purpose:
     Calculate NEE from observed Fc and observed/modeled Fre.
     Input and output names are held in ds.nee.
    Usage:
     qcrp.CalculateNEE(cf,ds)
      where cf is a conbtrol file object
            ds is a data structure
    Side effects:
     Series to hold the NEE data are created in ds.
    Author: PRI
    Date: August 2014
    """
    if "nee" not in dir(ds): return
    # get the Fsd and ustar thresholds
    Fsd_threshold = float(qcio.get_keyvaluefromcf(cf,["Params"],"Fsd_threshold",default=10))
    # get the incoming shortwave radiation and friction velocity
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        Fsd[index] = Fsd_syn[index]
    ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,"ustar")
    for label in ds.nee.keys():
        Fc_label = ds.nee[label]["Fc"]
        Fre_label = ds.nee[label]["Fre"]
        output_label = ds.nee[label]["output"]
        Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,Fc_label)
        Fre,Fre_flag,Fre_attr = qcutils.GetSeriesasMA(ds,Fre_label)
        # put the day time Fc into the NEE series
        index = numpy.ma.where(Fsd>=Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = Fc[index]
        ds.series[output_label]["Flag"][index] = Fc_flag[index]
        # put the night time Fre into the NEE series
        index = numpy.ma.where(Fsd<Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = Fre[index]
        ds.series[output_label]["Flag"][index] = Fre_flag[index]
        # copy the attributes
        attr = ds.series[output_label]["Attr"]
        attr["units"] = Fc_attr["units"]
        attr["long_name"] = "Net Ecosystem Exchange calculated from "+Fc_label+" (Fc) "
        attr["long_name"] = attr["long_name"]+" and "+Fre_label+" (Fre)"
        attr["comment1"] = "Fsd threshold used was "+str(Fsd_threshold)
    del ds.nee

def CalculateNEP(cf,ds):
    """
    Purpose:
     Calculate NEP from NEE
    Usage:
     qcrp.CalculateNEP(cf,ds)
      where cf is a conbtrol file object
            ds is a data structure
    Side effects:
     Series to hold the NEP data are created in ds.
    Author: PRI
    Date: May 2015
    """
    for nee_name in cf["NEE"].keys():
        nep_name = nee_name.replace("NEE","NEP")
        nee,flag,attr = qcutils.GetSeriesasMA(ds,nee_name)
        nep = float(-1)*nee
        attr["long_name"] = "Net Ecosystem Productivity calculated as -1*NEE"
        qcutils.CreateSeries(ds,nep_name,nep,Flag=flag,Attr=attr)

def FreUsingFFNET(cf,ds):
    """
    Purpose:
     Estimate ecosystem respiration using the ffnet neural network.
    Usage:
     qcrp.FreUsingFFNET(cf,ds)
      where cf is a control file object
            ds is a data structure
    Author: PRI
    Date: August 2014
    """
    if "ffnet" not in dir(ds): return
    if "ffnet" not in sys.modules.keys():
        log.error("FreUsingFFNET: I don't think ffnet is installed ...")
        return
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    rpFFNET_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                 "file_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # set up the GUI
    rpFFNET_gui = Tkinter.Toplevel()
    rpFFNET_gui.wm_title("FFNET GUI (Reco)")
    rpFFNET_gui.grid()
    # top row
    nrow = 0
    rpFFNET_gui.nodesLabel = Tkinter.Label(rpFFNET_gui,text="Hidden Nodes")
    rpFFNET_gui.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpFFNET_gui.nodesEntry = Tkinter.Entry(rpFFNET_gui,width=6)
    rpFFNET_gui.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpFFNET_gui.nodesEntry.insert(0,"2")
    rpFFNET_gui.trainingLabel = Tkinter.Label(rpFFNET_gui,text="Training")
    rpFFNET_gui.trainingLabel.grid(row=nrow,column=2,columnspan=1,sticky="E")
    rpFFNET_gui.trainingEntry = Tkinter.Entry(rpFFNET_gui,width=6)
    rpFFNET_gui.trainingEntry.grid(row=nrow,column=3,columnspan=1,sticky="W")
    rpFFNET_gui.trainingEntry.insert(0,"500")
    # second row
    nrow = nrow + 1
    rpFFNET_gui.trainOptionVar = Tkinter.StringVar()
    rpFFNET_gui.trainOptionVar.set("TNC")
    choices = ["BFGS","CG","Genetic","Back","Rprop","TNC"]
    rpFFNET_gui.trainOptionLabel = Tkinter.Label(rpFFNET_gui,text="Training type")
    rpFFNET_gui.trainOptionLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpFFNET_gui.trainOption = Tkinter.OptionMenu(rpFFNET_gui,rpFFNET_gui.trainOptionVar,*choices)
    rpFFNET_gui.trainOption.grid(row=nrow,column=1,columnspan=1,sticky="E")
    rpFFNET_gui.trainTypeVar = Tkinter.IntVar()
    rpFFNET_gui.trainTypeVar.set(1)
    rpFFNET_gui.trainTypeStd = Tkinter.Radiobutton(rpFFNET_gui,text="Standard",variable=rpFFNET_gui.trainTypeVar,value=1)
    rpFFNET_gui.trainTypeStd.grid(row=nrow,column=2,columnspan=1,sticky="W")
    rpFFNET_gui.trainTypeFC = Tkinter.Radiobutton(rpFFNET_gui,text="Fully connected",variable=rpFFNET_gui.trainTypeVar,value=2)
    rpFFNET_gui.trainTypeFC.grid(row=nrow,column=3,columnspan=1,sticky="W")
    # third row
    nrow = nrow + 1
    rpFFNET_gui.filestartLabel = Tkinter.Label(rpFFNET_gui,text="File start date")
    rpFFNET_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    rpFFNET_gui.fileendLabel = Tkinter.Label(rpFFNET_gui,text="File end date")
    rpFFNET_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    rpFFNET_gui.filestartValue = Tkinter.Label(rpFFNET_gui,text=str(ldt[0]))
    rpFFNET_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    rpFFNET_gui.fileendValue = Tkinter.Label(rpFFNET_gui,text=str(ldt[-1]))
    rpFFNET_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # fifth row
    nrow = nrow + 1
    rpFFNET_gui.startLabel = Tkinter.Label(rpFFNET_gui, text="Start date (YYYY-MM-DD)")
    rpFFNET_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    rpFFNET_gui.startEntry = Tkinter.Entry(rpFFNET_gui)
    rpFFNET_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # sixth row
    nrow = nrow + 1
    rpFFNET_gui.endLabel = Tkinter.Label(rpFFNET_gui, text="End date   (YYYY-MM-DD)")
    rpFFNET_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    rpFFNET_gui.endEntry = Tkinter.Entry(rpFFNET_gui)
    rpFFNET_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # seventh row
    nrow = nrow + 1
    rpFFNET_gui.peropt = Tkinter.IntVar()
    rpFFNET_gui.peropt.set(1)
    rpFFNET_gui.manualperiod = Tkinter.Radiobutton(rpFFNET_gui,text="Manual",variable=rpFFNET_gui.peropt,value=1)
    rpFFNET_gui.manualperiod.grid(row=nrow,column=0,columnspan=1,sticky="W")
    rpFFNET_gui.daysperiod = Tkinter.Radiobutton(rpFFNET_gui,text="No. days",variable=rpFFNET_gui.peropt,value=2)
    rpFFNET_gui.daysperiod.grid(row=nrow,column=1,sticky="W")
    rpFFNET_gui.daysentry = Tkinter.Entry(rpFFNET_gui,width=5)
    rpFFNET_gui.daysentry.grid(row=nrow,column=2,sticky="W")
    rpFFNET_gui.daysentry.insert(0,"30")
    rpFFNET_gui.monthsperiod = Tkinter.Radiobutton(rpFFNET_gui,text="No. Months",variable=rpFFNET_gui.peropt,value=3)
    rpFFNET_gui.monthsperiod.grid(row=nrow,column=3,sticky="W")
    rpFFNET_gui.monthsentry = Tkinter.Entry(rpFFNET_gui,width=5)
    rpFFNET_gui.monthsentry.grid(row=nrow,column=4,sticky="W")
    rpFFNET_gui.monthsentry.insert(0,"1")
    # eigth row
    nrow = nrow + 1
    rpFFNET_gui.autoyearly = Tkinter.Radiobutton(rpFFNET_gui,text="Yearly",variable=rpFFNET_gui.peropt,value=4)
    rpFFNET_gui.autoyearly.grid(row=nrow,column=0,columnspan=1,sticky="W")
    rpFFNET_gui.pointsperiod = Tkinter.Radiobutton(rpFFNET_gui,text="No. pts",variable=rpFFNET_gui.peropt,value=5)
    rpFFNET_gui.pointsperiod.grid(row=nrow,column=1,sticky="W")
    rpFFNET_gui.pointsentry = Tkinter.Entry(rpFFNET_gui,width=5)
    rpFFNET_gui.pointsentry.grid(row=nrow,column=2,columnspan=1,sticky="W")
    # ninth row
    nrow = nrow + 1
    rpFFNET_gui.minptsLabel = Tkinter.Label(rpFFNET_gui,text="Min points")
    rpFFNET_gui.minptsLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpFFNET_gui.minpts = Tkinter.Entry(rpFFNET_gui,width=5)
    rpFFNET_gui.minpts.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpFFNET_gui.minpts.insert(0,"200")
    rpFFNET_gui.owopt = Tkinter.IntVar()
    rpFFNET_gui.owopt.set(1)
    rpFFNET_gui.overwrite = Tkinter.Checkbutton(rpFFNET_gui, text="Overwrite", variable=rpFFNET_gui.owopt)
    rpFFNET_gui.overwrite.grid(row=nrow,column=2,columnspan=2,sticky="w")
    # tenth row
    nrow = nrow + 1
    rpFFNET_gui.doneButton = Tkinter.Button (rpFFNET_gui, text="Done",command=lambda:rpFFNET_done(ds,rpFFNET_gui))
    rpFFNET_gui.doneButton.grid(row=nrow,column=0,columnspan=3)
    rpFFNET_gui.runButton = Tkinter.Button (rpFFNET_gui, text="Run",command=lambda:rpFFNET_run(ds,rpFFNET_gui,rpFFNET_info))
    rpFFNET_gui.runButton.grid(row=nrow,column=3,columnspan=3)
    # eleventh row
    nrow = nrow + 1
    rpFFNET_gui.progress_row = nrow
    rpFFNET_gui.progress = Tkinter.Label(rpFFNET_gui, text='Waiting for input ...')
    rpFFNET_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    rpFFNET_gui.wait_window(rpFFNET_gui)

def FreUsingLasslop(cf,ds):
    log.info('Estimating Fre using Lasslop et al is not implemented yet')
    pass

def FreUsingLloydTaylor(cf,ds):
    log.info('Estimating Fre using Lloyd-Taylor is not implemented yet')
    pass
    ## need to check that all information required is in the control file
    #section = qcutils.get_cfsection(cf,series='Reco',mode='quiet')
    #if len(section)==0:
        #log.error('Reco section not found in control file')
        #return
    #if 'LloydTaylor' not in cf[section]['Reco']:
        #log.error('LloydTaylor section not found in control file')
        #return
    ## get the driver
    #if 'drivers' not in cf[section]['Reco']['LloydTaylor']:
        #log.error('drivers key not found in LloydTaylor section of control file')
        #return
    #else:
        #driver_list = eval(cf[section]['Reco']['LloydTaylor']['drivers'])
        #driver_label = driver_list[0]
    ## get the monthly values for the activation energy, E0
    #if 'E0' not in cf[section]['Reco']['LloydTaylor']:
        #log.error('E0 key not found in LloydTaylor section of control file')
        #return
    #else:
        #E0_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['E0']))
    ## get the monthly values for the base respiration, rb
    #if 'rb' not in cf[section]['Reco']['LloydTaylor']:
        #log.error('rb key not found in LloydTaylor section of control file')
        #return
    #else:
        #rb_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['rb']))
    ## get the output label
    #if 'output' not in cf[section]['Reco']['LloydTaylor']:
        #log.error('output key not found in LloydTaylor section of control file')
        #return
    #else:
        #out_label = cf[section]['Reco']['LloydTaylor']['output']
    ## ... and make an array of values for each month
    #nRecs = int(ds.globalattributes['nc_nrecs'])
    #E0 = numpy.ma.ones(nRecs)
    #rb = numpy.ma.zeros(nRecs)
    #month = ds.series['Month']['Data']
    #lwr = numpy.max(numpy.array([numpy.min(month),1]))
    #upr = numpy.min(numpy.array([numpy.max(month),12]))
    #for m in range(lwr,upr+1):
        #index = numpy.where(month==m)[0]
        #E0[index] = E0_monthly[m-1]
        #rb[index] = rb_monthly[m-1]
    ## get the driver data
    #Ts,flag,attr = qcutils.GetSeriesasMA(ds,driver_label)
    ## estimate Reco using the Lloyd-Taylor expression
    #t1 = 1/(c.Tref-c.T0)
    #t2 = 1/(Ts-c.T0)
    #Reco = rb*numpy.exp(E0*(t1-t2))
    ## put the estimated Reco into the data structure
    #units=qcutils.GetUnitsFromds(ds, 'Fc')
    #attr = qcutils.MakeAttributeDictionary(long_name='Reco estimated using Lloyd-Taylor',units=units)
    #qcutils.CreateSeries(ds,out_label,Reco,Flag=flag,Attr=attr)

def FreUsingSOLO(cf,ds):
    """ Estimate Fre using SOLO. """
    if "solo" not in dir(ds): return
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    rpSOLO_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                 "file_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # set up the GUI
    rpSOLO_gui = Tkinter.Toplevel()
    rpSOLO_gui.wm_title("SOLO GUI (Fre)")
    rpSOLO_gui.grid()
    # top row
    nrow = 0
    rpSOLO_gui.nodesLabel = Tkinter.Label(rpSOLO_gui,text="Nodes")
    rpSOLO_gui.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    rpSOLO_gui.nodesEntry = Tkinter.Entry(rpSOLO_gui,width=6)
    rpSOLO_gui.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    rpSOLO_gui.nodesEntry.insert(0,"2")
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
    rpSOLO_gui.manualperiod.grid(row=nrow,column=0,columnspan=2,sticky="W")
    rpSOLO_gui.daysperiod = Tkinter.Radiobutton(rpSOLO_gui,text="No. days",variable=rpSOLO_gui.peropt,value=2)
    rpSOLO_gui.daysperiod.grid(row=nrow,column=2,sticky="W")
    rpSOLO_gui.daysentry = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.daysentry.grid(row=nrow,column=3,sticky="W")
    rpSOLO_gui.daysentry.insert(0,"30")
    rpSOLO_gui.monthsperiod = Tkinter.Radiobutton(rpSOLO_gui,text="No. months",variable=rpSOLO_gui.peropt,value=3)
    rpSOLO_gui.monthsperiod.grid(row=nrow,column=4,sticky="W")
    rpSOLO_gui.monthsentry = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.monthsentry.grid(row=nrow,column=5,sticky="W")
    rpSOLO_gui.monthsentry.insert(0,"1")
    # eigth row
    nrow = nrow + 1
    rpSOLO_gui.autoyearly = Tkinter.Radiobutton(rpSOLO_gui,text="Yearly",variable=rpSOLO_gui.peropt,value=4)
    rpSOLO_gui.autoyearly.grid(row=nrow,column=0,columnspan=1,sticky="W")
    rpSOLO_gui.pointsperiod = Tkinter.Radiobutton(rpSOLO_gui,text="No. pts",variable=rpSOLO_gui.peropt,value=5)
    rpSOLO_gui.pointsperiod.grid(row=nrow,column=2,sticky="W")
    rpSOLO_gui.pointsentry = Tkinter.Entry(rpSOLO_gui,width=5)
    rpSOLO_gui.pointsentry.grid(row=nrow,column=3,columnspan=1,sticky="W")
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
    rpSOLO_gui.doneButton = Tkinter.Button (rpSOLO_gui, text="Done",command=lambda:rpSOLO_done(ds,rpSOLO_gui))
    rpSOLO_gui.doneButton.grid(row=nrow,column=0,columnspan=3)
    rpSOLO_gui.runButton = Tkinter.Button (rpSOLO_gui, text="Run",command=lambda:rpSOLO_run(ds,rpSOLO_gui,rpSOLO_info))
    rpSOLO_gui.runButton.grid(row=nrow,column=3,columnspan=3)
    # eleventh row
    nrow = nrow + 1
    rpSOLO_gui.progress_row = nrow
    rpSOLO_gui.progress = Tkinter.Label(rpSOLO_gui, text='Waiting for input ...')
    rpSOLO_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    rpSOLO_gui.wait_window(rpSOLO_gui)

def GetFreFromFc(cf,ds):
    """
    Purpose:
     Get tehe observed ecosystem respiration from measurements of Fc by
     filtering out daytime periods and periods when ustar is less than
     a threshold value.
     The Fsd threshold for determining day time and night time and the
     ustar threshold are set in the [Params] section of the L5 control
     file.
    Usage:
     qcrp.GetFreFromFc(cf,ds)
     where cf is a control file object
           ds is a data structure
    Side effects:
     A new series called "Fre" is created in the data structure.
    Author: PRI
    Date: August 2014
    """
    ts = int(ds.globalattributes["time_step"])
    ldt = ds.series['DateTime']['Data']
    if "Params" not in cf.keys():
        log.error("GetFreFromFc: no [Params] section in control file")
        raise Exception("GetFreFromFc: no [Params] section in control file")
    if "Fsd_threshold" not in cf["Params"]:
        log.error("GetFreFromFc: no Fsd threshold in control file")
        raise Exception("GetFreFromFc: no Fsd threshold in control file")
    else:
        Fsd_threshold = float(cf["Params"]["Fsd_threshold"])
    ustar_threshold_list = []
    if "ustar_threshold" in cf["Params"]:
        ustar_threshold = float(cf["Params"]["ustar_threshold"])
        ustar_threshold_list.append([ldt[0].strftime("%Y-%m-%d %H:%M"),
                    ldt[-1].strftime("%Y-%m-%d %H:%M"),
                    str(ustar_threshold)])
    else:
        if "ustar_threshold" in cf:
            for n in cf["ustar_threshold"].keys():
                ustar_threshold_list.append(ast.literal_eval(cf["ustar_threshold"][str(n)]))
        else:
            log.error("GetFreFromFc: no ustar threshold in control file")
            raise Exception("GetFreFromFc: no ustar threshold in control file")
    # get the data
    Fsd,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        Fsd[index] = Fsd_syn[index]
    ustar,flag,attr = qcutils.GetSeriesasMA(ds,"ustar")
    Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,"Fc")
    # check for any missing data
    for item,label in zip([Fsd,ustar,Fc],["Fsd","ustar","Fc"]):
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        if len(index)!=0:
            log.error(" GetFreFromFc: missing data in series "+label)
            raise Exception("GetFreFromFc: missing data in series "+label)
    # apply the day/night filter
    Fre1 = numpy.ma.masked_where(Fsd>Fsd_threshold,Fc,copy=True)
    # get a copy of the day/night filtered data
    Fre2 = numpy.ma.array(Fre1)
    # get a copy of the Fc flag
    Fre_flag = numpy.array(Fc_flag)
    # loop over the list of ustar thresholds
    # make the attribute dictionary first so we can add the ustar thresholds to it
    attr = qcutils.MakeAttributeDictionary(long_name='Ecosystem respiration (observed)',units=Fc_attr["units"],
                                           Fsd_threshold=str(Fsd_threshold))
    for i,list_item in enumerate(ustar_threshold_list):
        attr["ustar_threshold_"+str(i)] = str(list_item)
        # get the start and end datetime indices
        si = qcutils.GetDateIndex(ldt,list_item[0],ts=ts,default=0,match='exact')
        ei = qcutils.GetDateIndex(ldt,list_item[1],ts=ts,default=len(ldt),match='exact')
        # get the ustar thrhold
        ustar_threshold = float(list_item[2])
        # filter out the low ustar conditions
        Fre2[si:ei] = numpy.ma.masked_where(ustar[si:ei]<ustar_threshold,Fre1[si:ei])
        # set the QC flag
        index = numpy.where(numpy.ma.getmaskarray(Fre2[si:ei])==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fre2[si:ei])==True)[0]
        Fre_flag[si:ei][index] = numpy.int32(9)
    # put the nocturnal, filtered Fc data into the data structure
    qcutils.CreateSeries(ds,"Fre",Fre2,Flag=Fre_flag,Attr=attr)
    return True

def ParseL6ControlFile(cf,ds):
    """ Parse the L6 control file. """
    # start with the repiration section
    if "Respiration" in cf.keys():
        for ThisOne in cf["Respiration"].keys():
            if "FreUsingSOLO" in cf["Respiration"][ThisOne].keys():
                rpSOLO_createdict(cf,ds,ThisOne)      # create the SOLO dictionary in ds
            if "FreUsingFFNET" in cf["Respiration"][ThisOne].keys():
                rpFFNET_createdict(cf,ds,ThisOne)     # create the FFNET dictionary in ds
            if "MergeSeries" in cf["Respiration"][ThisOne].keys():
                rpMerge_createdict(cf,ds,ThisOne)      # create the merge dictionary in ds
    if "NEE" in cf.keys():
        for ThisOne in cf["NEE"].keys():
            rpNEE_createdict(cf,ds,ThisOne)
    if "GPP" in cf.keys():
        for ThisOne in cf["GPP"].keys():
            rpGPP_createdict(cf,ds,ThisOne)

def PartitionNEE(cf,ds):
    """
    Purpose:
     Partition NEE into GPP and Fre.
     Input and output names are held in ds.nee.
    Usage:
     qcrp.PartitionNEE(cf,ds)
      where cf is a conbtrol file object
            ds is a data structure
    Side effects:
     Series to hold the GPP data are created in ds.
    Author: PRI
    Date: August 2014
    """
    if "gpp" not in dir(ds): return
    # get the Fsd thresholds
    Fsd_threshold = float(cf['Params']['Fsd_threshold'])
    # get the incoming shortwave radiation
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        Fsd[index] = Fsd_syn[index]
    # calculate GPP from NEE and Fre
    for label in ds.gpp.keys():
        NEE_label = ds.gpp[label]["NEE"]
        Fre_label = ds.gpp[label]["Fre"]
        output_label = ds.gpp[label]["output"]
        NEE,NEE_flag,NEE_attr = qcutils.GetSeriesasMA(ds,NEE_label)
        Fre,Fre_flag,Fre_attr = qcutils.GetSeriesasMA(ds,Fre_label)
        # calculate GPP
        # here we use the conventions from Chapin et al (2006)
        #  NEP = -1*NEE
        #  GPP = NEP + Reco ==> GPP = -1*NEE + Reco
        GPP = float(-1)*NEE + Fre
        # put the day time data into the GPP series
        index = numpy.ma.where(Fsd>=Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = GPP[index]
        ds.series[output_label]["Flag"][index] = NEE_flag[index]
        # put the night time Fre into the NEE series
        # This force nocturnal GPP to be 0!  Not sure this is the right thing to do.
        index = numpy.ma.where(Fsd<Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = numpy.float64(0)
        ds.series[output_label]["Flag"][index] = numpy.int32(1)
        # copy the attributes
        attr = ds.series[output_label]["Attr"]
        attr["units"] = NEE_attr["units"]
        attr["long_name"] = "Gross Primary Productivity calculated from "+NEE_label+" as -NEE+Fre "
        attr["long_name"] = attr["long_name"]+" and "+Fre_label+" (Fre)"

def rp_getdiurnalstats(DecHour,Data,dt):
    nInts = 24*int((60/dt)+0.5)
    Hr = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Av = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Sd = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Mx = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Mn = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    for i in range(nInts):
        Hr[i] = float(i)*dt/60.
        li = numpy.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
        if numpy.size(li)!=0:
            Av[i] = numpy.mean(Data[li])
            Sd[i] = numpy.std(Data[li])
            Mx[i] = numpy.max(Data[li])
            Mn[i] = numpy.min(Data[li])
    Av = numpy.ma.masked_equal(Av,numpy.float64(-9999))
    Sd = numpy.ma.masked_equal(Sd,numpy.float64(-9999))
    Mx = numpy.ma.masked_equal(Mx,numpy.float64(-9999))
    Mn = numpy.ma.masked_equal(Mn,numpy.float64(-9999))
    return Hr, Av, Sd, Mx, Mn

def rpFFNET_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the FFNET data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("FreUsingFFNET: Series "+series+" not found in control file, skipping ...")
        return
    # check that none of the drivers have missing data
    driver_list = ast.literal_eval(cf[section][series]["FreUsingFFNET"]["drivers"])
    target = cf[section][series]["FreUsingFFNET"]["target"]
    for label in driver_list:
        data,flag,attr = qcutils.GetSeriesasMA(ds,label)
        if numpy.ma.count_masked(data)!=0:
            log.error("FreUsingFFNET: driver "+label+" contains missing data, skipping target "+target)
            return
    # create the ffnet directory in the data structure
    if "ffnet" not in dir(ds): ds.ffnet = {}
    # create the dictionary keys for this series
    ds.ffnet[series] = {}
    # site name
    ds.ffnet[series]["site_name"] = ds.globalattributes["site_name"]
    # target series name
    ds.ffnet[series]["target"] = cf[section][series]["FreUsingFFNET"]["target"]
    # list of drivers
    ds.ffnet[series]["drivers"] = ast.literal_eval(cf[section][series]["FreUsingFFNET"]["drivers"])
    # name of ffnet output series in ds
    ds.ffnet[series]["output"] = cf[section][series]["FreUsingFFNET"]["output"]
    # results of best fit for plotting later on
    ds.ffnet[series]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                  "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                  "Avg (obs)":[],"Avg (FFNET)":[],
                                  "Var (obs)":[],"Var (FFNET)":[],"Var ratio":[],
                                  "m_ols":[],"b_ols":[]}
    # create an empty series in ds if the SOLO output series doesn't exist yet
    if ds.ffnet[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.ffnet[series]["output"])
        qcutils.CreateSeries(ds,ds.ffnet[series]["output"],data,Flag=flag,Attr=attr)

def rpFFNET_done(ds,rpFFNET_gui):
    # destroy the FFNET GUI
    rpFFNET_gui.destroy()
    if "ffnet" in dir(ds): del ds.ffnet

def rpFFNET_initplot(**kwargs):
    # set the margins, heights, widths etc
    pd = {"margin_bottom":0.075,"margin_top":0.075,"margin_left":0.05,"margin_right":0.05,
          "xy_height":0.20,"xy_width":0.20,"xyts_space":0.05,"xyts_space":0.05,
          "ts_width":0.9}
    # set the keyword arguments
    for key, value in kwargs.iteritems():
        pd[key] = value
    # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])/float(pd["nDrivers"]+1)
    return pd

def rpFFNET_main(ds,rpFFNET_gui,rpFFNET_info):
    """
    This is the main routine for running FFNET, an artifical neural network for estimating Reco.
    """
    startdate = rpFFNET_info["startdate"]
    enddate = rpFFNET_info["enddate"]
    log.info(" Estimating Reco using FFNET: "+startdate+" to "+enddate)
    # read the control file again, this allows the contents of the control file to
    # be changed with the FFNET GUI still displayed
    cfname = ds.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname,mode="quiet")
    ffnet_series = ds.ffnet.keys()
    for series in ffnet_series:
        section = qcutils.get_cfsection(cf,series=series,mode="quiet")
        if len(section)==0: continue
        if series not in ds.series.keys(): continue
        ds.ffnet[series]["target"] = cf[section][series]["FreUsingFFNET"]["target"]
        ds.ffnet[series]["drivers"] = ast.literal_eval(cf[section][series]["FreUsingFFNET"]["drivers"])
        ds.ffnet[series]["output"] = cf[section][series]["FreUsingFFNET"]["output"]
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
        log.error(" FreUsingFFNET: end datetime index ("+str(ei)+") smaller that start ("+str(si)+")")
        return
    if si==0 and ei==-1:
        log.error(" FreUsingFFNET: no start and end datetime specified, using all data")
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        nRecs = ei - si + 1
    # get the figure number
    if len(plt.get_fignums())==0:
        fig_num = 0
    else:
        #fig_nums = plt.get_fignums()
        #fig_num = fig_nums[-1]
        fig_num = plt.get_fignums()[-1]
    # loop over the series to be gap filled using ffnet
    for series in ffnet_series:
        ds.ffnet[series]["results"]["startdate"].append(xldt[si])
        ds.ffnet[series]["results"]["enddate"].append(xldt[ei])
        target = ds.ffnet[series]["target"]
        d,f,a = qcutils.GetSeriesasMA(ds,target,si=si,ei=ei)
        if numpy.ma.count(d)<rpFFNET_info["min_points"]:
            log.error("rpFFNET: Less than "+str(rpFFNET_info["min_points"])+" points available for series "+series+" ...")
            ds.ffnet[series]["results"]["No. points"].append(float(0))
            results_list = ds.ffnet[series]["results"].keys()
            for item in ["startdate","enddate","No. points"]:
                if item in results_list: results_list.remove(item)
            for item in results_list:
                ds.ffnet[series]["results"][item].append(float(c.missing_value))
            continue
        drivers = ds.ffnet[series]["drivers"]
        ndrivers = len(drivers)
        output = ds.ffnet[series]["output"]
        # prepare the input and target data for training
        Reco,f,a = qcutils.GetSeriesasMA(ds,target,si=si,ei=ei)
        mask = numpy.ma.getmask(Reco)
        for val in drivers:
            d,f,a = qcutils.GetSeriesasMA(ds,val,si=si,ei=ei)
            mask = numpy.ma.mask_or(mask,d.mask)
        Reco.mask = mask
        nRecs = numpy.ma.count(Reco)
        data_nm = numpy.empty((nRecs,len(drivers)+1))
        for idx,val in enumerate(drivers):
            d,f,a = qcutils.GetSeriesasMA(ds,val,si=si,ei=ei)
            d.mask = mask
            data_nm[:,idx] = numpy.ma.compressed(d)
        data_nm[:,idx+1] = numpy.ma.compressed(Reco)
        input_train = data_nm[:,0:idx+1]
        target_train = data_nm[:,idx+1]
        # design the network
        hidden_layers = rpFFNET_info["hidden"].split(",")
        if len(hidden_layers)==1:
            arch = (ndrivers,int(hidden_layers[0]),1)
        elif len(hidden_layers)==2:
            arch = (ndrivers,int(hidden_layers[0]),int(hidden_layers[1]),1)
        else:
            log.error("FreUsingFFNET: more than 2 hidden layers specified, using 1 ("+str(ndrivers)+")")
            arch = (ndrivers,ndrivers,1)
        conec = ffnet.mlgraph(arch,biases=True)
        net = ffnet.ffnet(conec)
        # train the network
        if rpFFNET_info["train_type"].lower()=="tnc":
            net.train_tnc(input_train,target_train)
        elif rpFFNET_info["train_type"].lower()=="bfgs":
            net.train_bfgs(input_train,target_train)
        elif rpFFNET_info["train_type"].lower()=="cg":
            net.train_cg(input_train,target_train)
        elif rpFFNET_info["train_type"].lower()=="genetic":
            net.train_genetic(input_train,target_train)
        elif rpFFNET_info["train_type"].lower()=="back":
            net.train_momentum(input_train,target_train)
        elif rpFFNET_info["train_type"].lower()=="rprop":
            net.train_rprop(input_train,target_train)
        else:
            raise Exception("rpFFNET: unrecognised FFNET training option")
        #output,regress=net.test(input_train,target_train)
        # get the predictions
        input_predict = numpy.empty((len(Reco),len(drivers)))
        for idx,val in enumerate(drivers):
            d,f,a = qcutils.GetSeries(ds,val,si=si,ei=ei)
            input_predict[:,idx] = d[:]
        output_predict = net.call(input_predict)
        if ei==-1:
            ds.series[output]['Data'][si:] = output_predict[:,0]
            ds.series[output]['Flag'][si:] = numpy.int32(30)
        else:
            ds.series[output]['Data'][si:ei+1] = output_predict[:,0]
            ds.series[output]['Flag'][si:ei+1] = numpy.int32(30)
        # set the attributes
        ds.series[output]["Attr"]["units"] = ds.series[target]["Attr"]["units"]
        if "modelled by FFNET" not in ds.series[output]["Attr"]["long_name"]:
            ds.series[output]["Attr"]["long_name"] = "Ecosystem respiration modelled by FFNET (ANN)"
            ds.series[output]["Attr"]["comment1"] = "Target was "+str(target)
            ds.series[output]["Attr"]["comment2"] = "Drivers were "+str(drivers)
        # plot the results
        fig_num = fig_num + 1
        title = site_name+" : "+series+" estimated using FFNET"
        pd = rpFFNET_initplot(site_name=site_name,label=target,fig_num=fig_num,title=title,
                             nDrivers=len(drivers),startdate=startdate,enddate=enddate)
        rpFFNET_plot(pd,ds,series,drivers,target,output,rpFFNET_info,si=si,ei=ei)
    if 'FreUsingFFNET' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', FreUsingFFNET'

def rpFFNET_plot(pd,ds,series,driverlist,targetlabel,outputlabel,rpFFNET_info,si=0,ei=-1):
    """ Plot the results of the FFNET run. """
    # get the time step
    ts = int(ds.globalattributes['time_step'])
    # get a local copy of the datetime series
    xdt = numpy.array(ds.series['DateTime']['Data'][si:ei+1])
    Hdh,f,a = qcutils.GetSeriesasMA(ds,'Hdh',si=si,ei=ei)
    # get the observed and modelled values
    obs,f,a = qcutils.GetSeriesasMA(ds,targetlabel,si=si,ei=ei)
    mod,f,a = qcutils.GetSeriesasMA(ds,outputlabel,si=si,ei=ei)
    # make the figure
    plt.ion()
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    fig.canvas.set_window_title(targetlabel+" (FFNET): "+pd["startdate"]+" to "+pd["enddate"])
    plt.figtext(0.5,0.95,pd["title"],ha='center',size=16)
    # XY plot of the diurnal variation
    rect1 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    # get the diurnal stats of the observations
    mask = numpy.ma.mask_or(obs.mask,mod.mask)
    obs_mor = numpy.ma.array(obs,mask=mask)
    Hr1,Av1,Sd1,Mx1,Mn1 = rp_getdiurnalstats(Hdh,obs_mor,ts)
    ax1.plot(Hr1,Av1,'b-',label="Obs")
    # get the diurnal stats of all FFNET predictions
    mod_mor = numpy.ma.array(mod,mask=mask)
    #Hr2,Av2,Sd2,Mx2,Mn2 = rp_getdiurnalstats(Hdh,mod_mor,ts)
    Hr2,Av2,Sd2,Mx2,Mn2 = rp_getdiurnalstats(Hdh,mod,ts)
    ax1.plot(Hr2,Av2,'r-',label="FFNET(all)")
    if numpy.ma.count_masked(obs)!=0:
        index = numpy.where(numpy.ma.getmaskarray(obs)==False)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(obs)==False)[0]
        # get the diurnal stats of FFNET predictions when observations are present
        Hr3,Av3,Sd3,Mx3,Mn3=rp_getdiurnalstats(Hdh[index],mod_mor[index],ts)
        ax1.plot(Hr3,Av3,'g-',label="FFNET(obs)")
    plt.xlim(0,24)
    plt.xticks([0,6,12,18,24])
    ax1.set_ylabel(targetlabel)
    ax1.set_xlabel('Hour')
    ax1.legend(loc='upper right',frameon=False,prop={'size':8})
    # XY plot of the 30 minute data
    rect2 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    ax2.plot(mod,obs,'b.')
    ax2.set_ylabel(targetlabel+'_obs')
    ax2.set_xlabel(targetlabel+'_FFNET')
    # plot the best fit line
    coefs = numpy.ma.polyfit(numpy.ma.copy(mod),numpy.ma.copy(obs),1)
    xfit = numpy.ma.array([numpy.ma.minimum(mod),numpy.ma.maximum(mod)])
    yfit = numpy.polyval(coefs,xfit)
    r = numpy.ma.corrcoef(mod,obs)
    ax2.plot(xfit,yfit,'r--',linewidth=3)
    eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(coefs[0],coefs[1],r[0][1])
    ax2.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=ax2.transAxes)
    # write the fit statistics to the plot
    numpoints = numpy.ma.count(obs)
    numfilled = numpy.ma.count(mod)-numpy.ma.count(obs)
    diff = mod - obs
    bias = numpy.ma.average(diff)
    ds.ffnet[series]["results"]["Bias"].append(bias)
    rmse = numpy.ma.sqrt(numpy.ma.mean((obs-mod)*(obs-mod)))
    plt.figtext(0.65,0.225,'No. points')
    plt.figtext(0.75,0.225,str(numpoints))
    ds.ffnet[series]["results"]["No. points"].append(numpoints)
    plt.figtext(0.65,0.200,'Hidden nodes')
    plt.figtext(0.75,0.200,str(rpFFNET_info["hidden"]))
    plt.figtext(0.65,0.175,'Training')
    plt.figtext(0.75,0.175,str(rpFFNET_info["iterations"]))
    plt.figtext(0.65,0.150,'Training type')
    plt.figtext(0.75,0.150,str(rpFFNET_info["train_type"]))
    #plt.figtext(0.65,0.125,'Learning rate')
    #plt.figtext(0.75,0.125,str(rpSOLO_gui.learningrateEntry.get()))
    #plt.figtext(0.65,0.100,'Iterations')
    #plt.figtext(0.75,0.100,str(rpSOLO_gui.iterationsEntry.get()))
    plt.figtext(0.815,0.225,'No. filled')
    plt.figtext(0.915,0.225,str(numfilled))
    plt.figtext(0.815,0.200,'Slope')
    plt.figtext(0.915,0.200,str(qcutils.round2sig(coefs[0],sig=4)))
    ds.ffnet[series]["results"]["m_ols"].append(coefs[0])
    plt.figtext(0.815,0.175,'Offset')
    plt.figtext(0.915,0.175,str(qcutils.round2sig(coefs[1],sig=4)))
    ds.ffnet[series]["results"]["b_ols"].append(coefs[1])
    plt.figtext(0.815,0.150,'r')
    plt.figtext(0.915,0.150,str(qcutils.round2sig(r[0][1],sig=4)))
    ds.ffnet[series]["results"]["r"].append(r[0][1])
    plt.figtext(0.815,0.125,'RMSE')
    plt.figtext(0.915,0.125,str(qcutils.round2sig(rmse,sig=4)))
    ds.ffnet[series]["results"]["RMSE"].append(rmse)
    var_obs = numpy.ma.var(obs)
    ds.ffnet[series]["results"]["Var (obs)"].append(var_obs)
    var_mod = numpy.ma.var(mod)
    ds.ffnet[series]["results"]["Var (FFNET)"].append(var_mod)
    ds.ffnet[series]["results"]["Var ratio"].append(var_obs/var_mod)
    ds.ffnet[series]["results"]["Avg (obs)"].append(numpy.ma.average(obs))
    ds.ffnet[series]["results"]["Avg (FFNET)"].append(numpy.ma.average(mod))    
    # time series of drivers and target
    ts_axes = []
    rect = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    ts_axes.append(plt.axes(rect))
    ts_axes[0].plot(xdt,obs,'b.',xdt,mod,'r-')
    ts_axes[0].set_xlim(xdt[0],xdt[-1])
    TextStr = targetlabel+'_obs ('+ds.series[targetlabel]['Attr']['units']+')'
    ts_axes[0].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[0].transAxes)
    TextStr = outputlabel+'('+ds.series[outputlabel]['Attr']['units']+')'
    ts_axes[0].text(0.85,0.85,TextStr,color='r',horizontalalignment='right',transform=ts_axes[0].transAxes)
    for ThisOne,i in zip(driverlist,range(1,pd["nDrivers"]+1)):
        this_bottom = pd["ts_bottom"] + i*pd["ts_height"]
        rect = [pd["margin_left"],this_bottom,pd["ts_width"],pd["ts_height"]]
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))
        data,flag,attr = qcutils.GetSeriesasMA(ds,ThisOne,si=si,ei=ei)
        ts_axes[i].plot(xdt,data)
        plt.setp(ts_axes[i].get_xticklabels(),visible=False)
        TextStr = ThisOne+'('+ds.series[ThisOne]['Attr']['units']+')'
        ts_axes[i].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i].transAxes)
    # save a hard copy of the plot
    sdt = xdt[0].strftime("%Y%m%d")
    edt = xdt[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_FFNET_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    # draw the plot on the screen
    plt.draw()
    # turn off interactive plotting
    plt.ioff()
    
def rpFFNET_progress(rpFFNET_gui,text):
    """
        Update progress message in FFNET GUI
        """
    rpFFNET_gui.progress.destroy()
    rpFFNET_gui.progress = Tkinter.Label(rpFFNET_gui, text=text)
    rpFFNET_gui.progress.grid(row=rpFFNET_gui.progress_row,column=0,columnspan=4,sticky="W")
    rpFFNET_gui.update()

def rpFFNET_run(ds,rpFFNET_gui,rpFFNET_info):
    # populate the rpFFNET_info dictionary with things that will be useful
    rpFFNET_info["hidden"] = rpFFNET_gui.nodesEntry.get()
    rpFFNET_info["iterations"] = rpFFNET_gui.trainingEntry.get()
    rpFFNET_info["train_type"] = str(rpFFNET_gui.trainOptionVar.get())
    rpFFNET_info["peropt"] = rpFFNET_gui.peropt.get()
    rpFFNET_info["min_points"] = int(rpFFNET_gui.minpts.get())
    rpFFNET_info["site_name"] = ds.globalattributes["site_name"]
    rpFFNET_info["time_step"] = int(ds.globalattributes["time_step"])
    rpFFNET_info["nperhr"] = int(float(60)/rpFFNET_info["time_step"]+0.5)
    rpFFNET_info["nperday"] = int(float(24)*rpFFNET_info["nperhr"]+0.5)
    rpFFNET_info["maxlags"] = int(float(12)*rpFFNET_info["nperhr"]+0.5)
    rpFFNET_info["tower"] = {}
    rpFFNET_info["access"] = {}
    #log.info(" Estimating Reco using SOLO")
    if rpFFNET_gui.peropt.get()==1:
        rpFFNET_progress(rpFFNET_gui,"Starting manual run ...")
        # get the start and end datetimes entered in the SOLO GUI
        rpFFNET_info["startdate"] = rpFFNET_gui.startEntry.get()
        if len(rpFFNET_info["startdate"])==0: rpFFNET_info["startdate"] = rpFFNET_info["file_startdate"]
        rpFFNET_info["enddate"] = rpFFNET_gui.endEntry.get()
        if len(rpFFNET_info["enddate"])==0: rpFFNET_info["enddate"] = rpFFNET_info["file_enddate"]
        rpFFNET_main(ds,rpFFNET_gui,rpFFNET_info)
        rpFFNET_progress(rpFFNET_gui,"Finished manual run ...")
    elif rpFFNET_gui.peropt.get()==2:
        rpFFNET_progress(rpFFNET_gui,"Starting auto (days) run ...")
        # get the start datetime entered in the SOLO GUI
        rpFFNET_info["startdate"] = rpFFNET_gui.startEntry.get()
        if len(rpFFNET_info["startdate"])==0: rpFFNET_info["startdate"] = rpFFNET_info["file_startdate"]
        startdate = dateutil.parser.parse(rpFFNET_info["startdate"])
        file_startdate = dateutil.parser.parse(rpFFNET_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpFFNET_info["file_enddate"])
        nDays = int(rpFFNET_gui.daysentry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
        enddate = min([file_enddate,enddate])
        rpFFNET_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpFFNET_main(ds,rpFFNET_gui,rpFFNET_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
            rpFFNET_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpFFNET_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ## plot the summary statistics
        #gfSOLO_plotsummary(ds)
        rpFFNET_progress(rpFFNET_gui,"Finished auto (days) run ...")
    elif rpFFNET_gui.peropt.get()==3:
        rpFFNET_progress(rpFFNET_gui,"Starting auto (monthly) run ...")
        # get the start datetime entered in the SOLO GUI
        rpFFNET_info["startdate"] = rpFFNET_gui.startEntry.get()
        if len(rpFFNET_info["startdate"])==0: rpFFNET_info["startdate"] = rpFFNET_info["file_startdate"]
        startdate = dateutil.parser.parse(rpFFNET_info["startdate"])
        file_startdate = dateutil.parser.parse(rpFFNET_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpFFNET_info["file_enddate"])
        nMonths = int(rpFFNET_gui.monthsentry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(months=nMonths)
        enddate = min([file_enddate,enddate])
        rpFFNET_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpFFNET_main(ds,rpFFNET_gui,rpFFNET_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=nMonths)
            rpFFNET_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpFFNET_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ## plot the summary statistics
        #gfSOLO_plotsummary(ds)
        rpFFNET_progress(rpFFNET_gui,"Finished auto (monthly) run ...")
    elif rpFFNET_gui.peropt.get()==4:
        # automatic run with yearly datetime periods
        rpFFNET_progress(rpFFNET_gui,"Starting auto (yearly) run ...")
        # get the start date
        rpFFNET_info["startdate"] = rpFFNET_gui.startEntry.get()
        if len(rpFFNET_info["startdate"])==0: rpFFNET_info["startdate"] = rpFFNET_info["file_startdate"]
        startdate = dateutil.parser.parse(rpFFNET_info["startdate"])
        # get the start year
        start_year = startdate.year
        enddate = dateutil.parser.parse(str(start_year+1)+"-01-01 00:00")
        file_enddate = dateutil.parser.parse(rpFFNET_info["file_enddate"])
        enddate = min([file_enddate,enddate])
        rpFFNET_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpFFNET_main(ds,rpFFNET_gui,rpFFNET_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(years=1)
            rpFFNET_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpFFNET_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ### plot the summary statistics
        ##gfSOLO_plotsummary(ds)
        rpFFNET_progress(rpFFNET_gui,"Finished auto (yearly) run ...")
    elif rpFFNET_gui.peropt.get()==5:
        pass

def rpGPP_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about calculating GPP."""
    # create the ffnet directory in the data structure
    if "gpp" not in dir(ds): ds.gpp = {}
    # create the dictionary keys for this series
    ds.gpp[series] = {}
    # output series name
    ds.gpp[series]["output"] = series
    # CO2 flux
    ds.gpp[series]["NEE"] = cf["GPP"][series]["NEE"]
    # ecosystem respiration
    ds.gpp[series]["Fre"] = cf["GPP"][series]["Fre"]
    # create an empty series in ds if the output series doesn't exist yet
    if ds.gpp[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.gpp[series]["output"])
        qcutils.CreateSeries(ds,ds.gpp[series]["output"],data,Flag=flag,Attr=attr)

def rpMerge_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the merging of gap filled
        and tower data."""
    merge_prereq_list = []
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # create the ffnet directory in the data structure
    if "merge" not in dir(ds): ds.merge = {}
    # check to see if this series is in the "merge first" list
    # series in the "merge first" list get merged first so they can be used with existing tower
    # data to re-calculate Fg, Fn and Fa
    merge_order = "standard"
    if series in merge_prereq_list: merge_order = "prerequisite"
    if merge_order not in ds.merge.keys(): ds.merge[merge_order] = {}
    # create the dictionary keys for this series
    ds.merge[merge_order][series] = {}
    # output series name
    ds.merge[merge_order][series]["output"] = series
    # source
    ds.merge[merge_order][series]["source"] = ast.literal_eval(cf[section][series]["MergeSeries"]["Source"])
    # create an empty series in ds if the output series doesn't exist yet
    if ds.merge[merge_order][series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.merge[merge_order][series]["output"])
        qcutils.CreateSeries(ds,ds.merge[merge_order][series]["output"],data,Flag=flag,Attr=attr)

def rpNEE_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about calculating NEE."""
    # create the ffnet directory in the data structure
    if "nee" not in dir(ds): ds.nee = {}
    # create the dictionary keys for this series
    ds.nee[series] = {}
    # output series name
    ds.nee[series]["output"] = series
    # CO2 flux
    ds.nee[series]["Fc"] = cf["NEE"][series]["Fc"]
    # ecosystem respiration
    ds.nee[series]["Fre"] = cf["NEE"][series]["Fre"]
    # create an empty series in ds if the output series doesn't exist yet
    if ds.nee[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.nee[series]["output"])
        qcutils.CreateSeries(ds,ds.nee[series]["output"],data,Flag=flag,Attr=attr)

def rpSOLO_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the SOLO data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("FreUsingSOLO: Series "+series+" not found in control file, skipping ...")
        return
    # check that none of the drivers have missing data
    driver_list = ast.literal_eval(cf[section][series]["FreUsingSOLO"]["drivers"])
    target = cf[section][series]["FreUsingSOLO"]["target"]
    for label in driver_list:
        data,flag,attr = qcutils.GetSeriesasMA(ds,label)
        if numpy.ma.count_masked(data)!=0:
            log.error("FreUsingSOLO: driver "+label+" contains missing data, skipping target "+target)
            return
    # create the solo directory in the data structure
    if "solo" not in dir(ds): ds.solo = {}
    # create the dictionary keys for this series
    ds.solo[series] = {}
    # site name
    ds.solo[series]["site_name"] = ds.globalattributes["site_name"]
    # target series name
    ds.solo[series]["target"] = cf[section][series]["FreUsingSOLO"]["target"]
    # list of drivers
    ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["FreUsingSOLO"]["drivers"])
    # name of SOLO output series in ds
    ds.solo[series]["output"] = cf[section][series]["FreUsingSOLO"]["output"]
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

def rpSOLO_done(ds,rpSOLO_gui):
    # destroy the SOLO GUI
    rpSOLO_gui.destroy()
    if "solo" in dir(ds): del ds.solo

def rpSOLO_initplot(**kwargs):
    # set the margins, heights, widths etc
    pd = {"margin_bottom":0.075,"margin_top":0.075,"margin_left":0.05,"margin_right":0.05,
          "xy_height":0.20,"xy_width":0.20,"xyts_space":0.05,"xyts_space":0.05,
          "ts_width":0.9}
    # set the keyword arguments
    for key, value in kwargs.iteritems():
        pd[key] = value
    # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])/float(pd["nDrivers"]+1)
    return pd

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
    solo_series = ds.solo.keys()
    for series in solo_series:
        section = qcutils.get_cfsection(cf,series=series,mode="quiet")
        if len(section)==0: continue
        if series not in ds.series.keys(): continue
        ds.solo[series]["target"] = cf[section][series]["FreUsingSOLO"]["target"]
        ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["FreUsingSOLO"]["drivers"])
        ds.solo[series]["output"] = cf[section][series]["FreUsingSOLO"]["output"]
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
        log.error(" FreUsingSOLO: end datetime index ("+str(ei)+") smaller that start ("+str(si)+")")
        return
    if si==0 and ei==-1:
        log.error(" FreUsingSOLO: no start and end datetime specified, using all data")
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        nRecs = ei - si + 1
    # get the figure number
    if len(plt.get_fignums())==0:
        fig_num = 0
    else:
        #fig_nums = plt.get_fignums()
        #fig_num = fig_nums[-1]
        fig_num = plt.get_fignums()[-1]
    # loop over the series to be gap filled using solo
    for series in solo_series:
        ds.solo[series]["results"]["startdate"].append(xldt[si])
        ds.solo[series]["results"]["enddate"].append(xldt[ei])
        target = ds.solo[series]["target"]
        d,f,a = qcutils.GetSeriesasMA(ds,target,si=si,ei=ei)
        if numpy.ma.count(d)<rpSOLO_info["min_points"]:
            log.error("rpSOLO: Less than "+str(rpSOLO_info["min_points"])+" points available for series "+target+" ...")
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
        result = rpSOLO_runsofm(ds,rpSOLO_gui,drivers,target,nRecs,si=si,ei=ei)
        if result!=1: return
        # run SOLO
        result = rpSOLO_runsolo(ds,drivers,target,nRecs,si=si,ei=ei)
        if result!=1: return
        # run SEQSOLO and put the SOLO data into the data structure
        result = rpSOLO_runseqsolo(ds,drivers,target,output,nRecs,si=si,ei=ei)
        if result!=1: return
        # plot the results
        fig_num = fig_num + 1
        title = site_name+" : "+series+" estimated using SOLO"
        pd = rpSOLO_initplot(site_name=site_name,label=target,fig_num=fig_num,title=title,
                             nDrivers=len(drivers),startdate=startdate,enddate=enddate)
        rpSOLO_plot(pd,ds,series,drivers,target,output,rpSOLO_gui,si=si,ei=ei)
        # reset the nodesEntry in the rpSOLO_gui
        if nodesAuto: rpSOLO_resetnodesEntry(rpSOLO_gui)
    if 'FreUsingSOLO' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', FreUsingSOLO'

def rpSOLO_plot(pd,ds,series,driverlist,targetlabel,outputlabel,rpSOLO_gui,si=0,ei=-1):
    """ Plot the results of the SOLO run. """
    # get the time step
    ts = int(ds.globalattributes['time_step'])
    # get a local copy of the datetime series
    xdt = numpy.array(ds.series['DateTime']['Data'][si:ei+1])
    Hdh,f,a = qcutils.GetSeriesasMA(ds,'Hdh',si=si,ei=ei)
    # get the observed and modelled values
    obs,f,a = qcutils.GetSeriesasMA(ds,targetlabel,si=si,ei=ei)
    mod,f,a = qcutils.GetSeriesasMA(ds,outputlabel,si=si,ei=ei)
    # make the figure
    plt.ion()
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    fig.canvas.set_window_title(targetlabel+" (SOLO): "+pd["startdate"]+" to "+pd["enddate"])
    plt.figtext(0.5,0.95,pd["title"],ha='center',size=16)
    # XY plot of the diurnal variation
    rect1 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    # get the diurnal stats of the observations
    mask = numpy.ma.mask_or(obs.mask,mod.mask)
    obs_mor = numpy.ma.array(obs,mask=mask)
    Hr1,Av1,Sd1,Mx1,Mn1 = rp_getdiurnalstats(Hdh,obs_mor,ts)
    ax1.plot(Hr1,Av1,'b-',label="Obs")
    # get the diurnal stats of all SOLO predictions
    mod_mor = numpy.ma.array(mod,mask=mask)
    #Hr2,Av2,Sd2,Mx2,Mn2 = rp_getdiurnalstats(Hdh,mod_mor,ts)
    Hr2,Av2,Sd2,Mx2,Mn2 = rp_getdiurnalstats(Hdh,mod,ts)
    ax1.plot(Hr2,Av2,'r-',label="SOLO(all)")
    if numpy.ma.count_masked(obs)!=0:
        index = numpy.where(numpy.ma.getmaskarray(obs)==False)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(obs)==False)[0]
        # get the diurnal stats of SOLO predictions when observations are present
        Hr3,Av3,Sd3,Mx3,Mn3=rp_getdiurnalstats(Hdh[index],mod_mor[index],ts)
        ax1.plot(Hr3,Av3,'g-',label="SOLO(obs)")
    plt.xlim(0,24)
    plt.xticks([0,6,12,18,24])
    ax1.set_ylabel(targetlabel)
    ax1.set_xlabel('Hour')
    ax1.legend(loc='upper right',frameon=False,prop={'size':8})
    # XY plot of the 30 minute data
    rect2 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    ax2.plot(mod,obs,'b.')
    ax2.set_ylabel(targetlabel+'_obs')
    ax2.set_xlabel(targetlabel+'_SOLO')
    # plot the best fit line
    coefs = numpy.ma.polyfit(numpy.ma.copy(mod),numpy.ma.copy(obs),1)
    xfit = numpy.ma.array([numpy.ma.minimum(mod),numpy.ma.maximum(mod)])
    yfit = numpy.polyval(coefs,xfit)
    r = numpy.ma.corrcoef(mod,obs)
    ax2.plot(xfit,yfit,'r--',linewidth=3)
    eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(coefs[0],coefs[1],r[0][1])
    ax2.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=ax2.transAxes)
    # write the fit statistics to the plot
    numpoints = numpy.ma.count(obs)
    numfilled = numpy.ma.count(mod)-numpy.ma.count(obs)
    diff = mod - obs
    bias = numpy.ma.average(diff)
    ds.solo[series]["results"]["Bias"].append(bias)
    rmse = numpy.ma.sqrt(numpy.ma.mean((obs-mod)*(obs-mod)))
    plt.figtext(0.65,0.225,'No. points')
    plt.figtext(0.75,0.225,str(numpoints))
    ds.solo[series]["results"]["No. points"].append(numpoints)
    plt.figtext(0.65,0.200,'Nodes')
    plt.figtext(0.75,0.200,str(rpSOLO_gui.nodesEntry.get()))
    plt.figtext(0.65,0.175,'Training')
    plt.figtext(0.75,0.175,str(rpSOLO_gui.trainingEntry.get()))
    plt.figtext(0.65,0.150,'Nda factor')
    plt.figtext(0.75,0.150,str(rpSOLO_gui.factorEntry.get()))
    plt.figtext(0.65,0.125,'Learning rate')
    plt.figtext(0.75,0.125,str(rpSOLO_gui.learningrateEntry.get()))
    plt.figtext(0.65,0.100,'Iterations')
    plt.figtext(0.75,0.100,str(rpSOLO_gui.iterationsEntry.get()))
    plt.figtext(0.815,0.225,'No. filled')
    plt.figtext(0.915,0.225,str(numfilled))
    plt.figtext(0.815,0.200,'Slope')
    plt.figtext(0.915,0.200,str(qcutils.round2sig(coefs[0],sig=4)))
    ds.solo[series]["results"]["m_ols"].append(coefs[0])
    plt.figtext(0.815,0.175,'Offset')
    plt.figtext(0.915,0.175,str(qcutils.round2sig(coefs[1],sig=4)))
    ds.solo[series]["results"]["b_ols"].append(coefs[1])
    plt.figtext(0.815,0.150,'r')
    plt.figtext(0.915,0.150,str(qcutils.round2sig(r[0][1],sig=4)))
    ds.solo[series]["results"]["r"].append(r[0][1])
    plt.figtext(0.815,0.125,'RMSE')
    plt.figtext(0.915,0.125,str(qcutils.round2sig(rmse,sig=4)))
    ds.solo[series]["results"]["RMSE"].append(rmse)
    var_obs = numpy.ma.var(obs)
    ds.solo[series]["results"]["Var (obs)"].append(var_obs)
    var_mod = numpy.ma.var(mod)
    ds.solo[series]["results"]["Var (SOLO)"].append(var_mod)
    ds.solo[series]["results"]["Var ratio"].append(var_obs/var_mod)
    ds.solo[series]["results"]["Avg (obs)"].append(numpy.ma.average(obs))
    ds.solo[series]["results"]["Avg (SOLO)"].append(numpy.ma.average(mod))    
    # time series of drivers and target
    ts_axes = []
    rect = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    ts_axes.append(plt.axes(rect))
    ts_axes[0].plot(xdt,obs,'b.',xdt,mod,'r-')
    ts_axes[0].set_xlim(xdt[0],xdt[-1])
    TextStr = targetlabel+'_obs ('+ds.series[targetlabel]['Attr']['units']+')'
    ts_axes[0].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[0].transAxes)
    TextStr = outputlabel+'('+ds.series[outputlabel]['Attr']['units']+')'
    ts_axes[0].text(0.85,0.85,TextStr,color='r',horizontalalignment='right',transform=ts_axes[0].transAxes)
    for ThisOne,i in zip(driverlist,range(1,pd["nDrivers"]+1)):
        this_bottom = pd["ts_bottom"] + i*pd["ts_height"]
        rect = [pd["margin_left"],this_bottom,pd["ts_width"],pd["ts_height"]]
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))
        data,flag,attr = qcutils.GetSeriesasMA(ds,ThisOne,si=si,ei=ei)
        ts_axes[i].plot(xdt,data)
        plt.setp(ts_axes[i].get_xticklabels(),visible=False)
        TextStr = ThisOne+'('+ds.series[ThisOne]['Attr']['units']+')'
        ts_axes[i].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i].transAxes)
    # save a hard copy of the plot
    sdt = xdt[0].strftime("%Y%m%d")
    edt = xdt[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_SOLO_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    # draw the plot on the screen
    plt.draw()
    # turn off interactive plotting
    plt.ioff()

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
    # populate the rpSOLO_info dictionary with things that will be useful
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
        # manual run using start and end datetime entered via GUI
        rpSOLO_progress(rpSOLO_gui,"Starting manual run ...")
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        rpSOLO_info["enddate"] = rpSOLO_gui.endEntry.get()
        if len(rpSOLO_info["enddate"])==0: rpSOLO_info["enddate"] = rpSOLO_info["file_enddate"]
        rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
        rpSOLO_progress(rpSOLO_gui,"Finished manual run ...")
    elif rpSOLO_gui.peropt.get()==2:
        # automatc run with number of days specified by user via the GUI
        rpSOLO_progress(rpSOLO_gui,"Starting auto (days) run ...")
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
        #gfSOLO_plotsummary(ds)
        rpSOLO_progress(rpSOLO_gui,"Finished auto (days) run ...")
    elif rpSOLO_gui.peropt.get()==3:
        # automatic run with monthly datetime periods
        rpSOLO_progress(rpSOLO_gui,"Starting auto (monthly) run ...")
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        startdate = dateutil.parser.parse(rpSOLO_info["startdate"])
        file_startdate = dateutil.parser.parse(rpSOLO_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpSOLO_info["file_enddate"])
        nMonths = int(rpSOLO_gui.monthsentry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(months=nMonths)
        enddate = min([file_enddate,enddate])
        rpSOLO_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=nMonths)
            rpSOLO_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpSOLO_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ## plot the summary statistics
        #gfSOLO_plotsummary(ds)
        rpSOLO_progress(rpSOLO_gui,"Finished auto (monthly) run ...")
    elif rpSOLO_gui.peropt.get()==4:
        # automatic run with yearly datetime periods
        rpSOLO_progress(rpSOLO_gui,"Starting auto (yearly) run ...")
        # get the start date
        rpSOLO_info["startdate"] = rpSOLO_gui.startEntry.get()
        if len(rpSOLO_info["startdate"])==0: rpSOLO_info["startdate"] = rpSOLO_info["file_startdate"]
        startdate = dateutil.parser.parse(rpSOLO_info["startdate"])
        # get the start year
        start_year = startdate.year
        enddate = dateutil.parser.parse(str(start_year+1)+"-01-01 00:00")
        #file_startdate = dateutil.parser.parse(rpSOLO_info["file_startdate"])
        file_enddate = dateutil.parser.parse(rpSOLO_info["file_enddate"])
        #enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
        enddate = min([file_enddate,enddate])
        rpSOLO_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            rpSOLO_main(ds,rpSOLO_gui,rpSOLO_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(years=1)
            rpSOLO_info["startdate"] = startdate.strftime("%Y-%m-%d")
            rpSOLO_info["enddate"] = enddate.strftime("%Y-%m-%d")
        ### plot the summary statistics
        ##gfSOLO_plotsummary(ds)
        rpSOLO_progress(rpSOLO_gui,"Finished auto (yearly) run ...")
    elif rpSOLO_gui.peropt.get()==5:
        pass

def rpSOLO_runseqsolo(ds,driverlist,targetlabel,outputlabel,nRecs,si=0,ei=-1):
    '''
    Run SEQSOLO.
    '''
    # get the number of drivers    
    ndrivers = len(driverlist)
    # add an extra column for the target data
    seqsoloinputdata = numpy.zeros((nRecs,ndrivers+1))
    # now fill the driver data array
    i = 0
    for TheseOnes in driverlist:
        driver,flag,attr = qcutils.GetSeries(ds,TheseOnes,si=si,ei=ei)
        seqsoloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the unmodified ds each time
    target,flag,attr = qcutils.GetSeries(ds,targetlabel,si=si,ei=ei)
    # now load the target data into the data array
    seqsoloinputdata[:,ndrivers] = target[:]
    # now strip out the bad data
    cind = numpy.zeros(nRecs)
    iind = numpy.arange(nRecs)
    # do only the drivers not the target
    for i in range(ndrivers):
        index = numpy.where(seqsoloinputdata[:,i]==c.missing_value)[0]
        if len(index!=0): cind[index] = 1
    # index of good data
    index = numpy.where(cind==0)[0]
    nRecs_good = len(index)
    gooddata = numpy.zeros((nRecs_good,ndrivers+1))
    for i in range(ndrivers+1):
        gooddata[:,i] = seqsoloinputdata[:,i][index]
    # keep track of the good data indices
    goodindex = iind[index]
    # and then write the seqsolo input file
    seqsolofile = open('solo/input/seqsolo_input.csv','wb')
    wr = csv.writer(seqsolofile,delimiter=',')
    for i in range(gooddata.shape[0]):
        wr.writerow(gooddata[i,0:ndrivers+1])
    seqsolofile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/seqOut2.out'): os.remove('solo/output/seqOut2.out')
    # now run SEQSOLO
    #log.info(' GapFillUsingSOLO: running SEQSOLO')
    seqsolologfile = open('solo/log/seqsolo.log','wb')
    if platform.system()=="Windows":
        subprocess.call(['./solo/bin/seqsolo.exe','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    else:
        subprocess.call(['./solo/bin/seqsolo','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    seqsolologfile.close()
    # check to see if the solo output file exists, this is used to indicate that solo ran correctly
    if os.path.exists('solo/output/seqOut2.out'):
        # now read in the seqsolo results, use the seqOut2 file so that the learning capability of
        # seqsolo can be used via the "learning rate" and "Iterations" GUI options
        seqdata = numpy.genfromtxt('solo/output/seqOut2.out')
        # put the SOLO modelled data back into the data series
        if ei==-1:
            ds.series[outputlabel]['Data'][si:][goodindex] = seqdata[:,1]
            ds.series[outputlabel]['Flag'][si:][goodindex] = numpy.int32(30)
        else:
            ds.series[outputlabel]['Data'][si:ei+1][goodindex] = seqdata[:,1]
            ds.series[outputlabel]['Flag'][si:ei+1][goodindex] = numpy.int32(30)
        # set the attributes
        ds.series[outputlabel]["Attr"]["units"] = ds.series[targetlabel]["Attr"]["units"]
        if "modelled by SOLO" not in ds.series[outputlabel]["Attr"]["long_name"]:
            ds.series[outputlabel]["Attr"]["long_name"] = "Ecosystem respiration modelled by SOLO (ANN)"
            ds.series[outputlabel]["Attr"]["comment1"] = "Target was "+str(targetlabel)
            ds.series[outputlabel]["Attr"]["comment2"] = "Drivers were "+str(driverlist)
        return 1
    else:
        log.error(' SOLO_runseqsolo: SEQSOLO did not run correctly, check the SOLO GUI and the log files')
        return 0

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

def rpSOLO_runsolo(ds,driverlist,targetlabel,nRecs,si=0,ei=-1):
    '''
    Run SOLO.
    '''
    ndrivers = len(driverlist)
    # add an extra column for the target data
    soloinputdata = numpy.zeros((nRecs,ndrivers+1))
    # now fill the driver data array, drivers come from the modified ds
    i = 0
    for TheseOnes in driverlist:
        driver,flag,attr = qcutils.GetSeries(ds,TheseOnes,si=si,ei=ei)
        soloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the ds each time
    target,flag,attr = qcutils.GetSeries(ds,targetlabel,si=si,ei=ei)
    # now load the target data into the data array
    soloinputdata[:,ndrivers] = target[:]
    # now strip out the bad data
    cind = numpy.zeros(nRecs)
    for i in range(ndrivers+1):
        index = numpy.where(soloinputdata[:,i]==c.missing_value)[0]
        if len(index!=0): cind[index] = 1
    index = numpy.where(cind==0)[0]
    nRecs_good = len(index)
    gooddata = numpy.zeros((nRecs_good,ndrivers+1))
    for i in range(ndrivers+1):
        gooddata[:,i] = soloinputdata[:,i][index]
    # and then write the solo input file, the name is assumed by the solo.inf control file
    solofile = open('solo/input/solo_input.csv','wb')
    wr = csv.writer(solofile,delimiter=',')
    for i in range(gooddata.shape[0]):
        wr.writerow(gooddata[i,0:ndrivers+1])
    solofile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/eigenValue.out'): os.remove('solo/output/eigenValue.out')
    # now run SOLO
    #log.info(' GapFillUsingSOLO: running SOLO')
    solologfile = open('solo/log/solo.log','wb')
    if platform.system()=="Windows":
        subprocess.call(['./solo/bin/solo.exe','solo/inf/solo.inf'],stdout=solologfile)
    else:
        subprocess.call(['./solo/bin/solo','solo/inf/solo.inf'],stdout=solologfile)
    solologfile.close()
    # check to see if the solo output file exists, this is used to indicate that solo ran correctly
    if os.path.exists('solo/output/eigenValue.out'):
        return 1
    else:
        log.error(' SOLO_runsolo: SOLO did not run correctly, check the SOLO GUI and the log files')
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
