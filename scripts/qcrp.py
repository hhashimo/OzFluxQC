import ast
import collections
import constants as c
import datetime
import dateutil
import logging
import matplotlib.pyplot as plt
import numpy
import os
import qcio
import qcrpNN
import qcts
import qcutils
import sys
import xlrd

log = logging.getLogger('qc.rp')

def CalculateET(ds):
    """
    Purpose:
     Calculate ET from Fe
    Usage:
     qcrp.CalculateET(ds)
      where ds is a data structure
    Side effects:
     Series to hold the ET data are created in ds.
    Author: PRI
    Date: June 2015
    """
    ts = int(ds.globalattributes["time_step"])
    Fe,flag,attr = qcutils.GetSeriesasMA(ds,"Fe")
    ET = Fe*ts*60/c.Lv
    attr["long_name"] = "Evapo-transpiration calculated from latent heat flux"
    attr["units"] = "mm"
    qcutils.CreateSeries(ds,"ET",ET,Flag=flag,Attr=attr)

def CalculateNEE(cf,ds):
    """
    Purpose:
     Calculate NEE from observed Fc and observed/modeled ER.
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
    Fsd_threshold = float(qcutils.get_keyvaluefromcf(cf,["Params"],"Fsd_threshold",default=10))
    # get the incoming shortwave radiation and friction velocity
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        Fsd[index] = Fsd_syn[index]
    ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,"ustar")
    for label in ds.nee.keys():
        if "Fc" not in ds.nee[label] and "ER" not in ds.nee[label]: continue
        Fc_label = ds.nee[label]["Fc"]
        ER_label = ds.nee[label]["ER"]
        output_label = ds.nee[label]["output"]
        Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,Fc_label)
        ER,ER_flag,ER_attr = qcutils.GetSeriesasMA(ds,ER_label)
        # put the day time Fc into the NEE series
        index = numpy.ma.where(Fsd>=Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = Fc[index]
        ds.series[output_label]["Flag"][index] = Fc_flag[index]
        # put the night time ER into the NEE series
        index = numpy.ma.where(Fsd<Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = ER[index]
        ds.series[output_label]["Flag"][index] = ER_flag[index]
        # copy the attributes
        attr = ds.series[output_label]["Attr"]
        attr["units"] = Fc_attr["units"]
        attr["long_name"] = "Net Ecosystem Exchange calculated from "+Fc_label+" (Fc) "
        attr["long_name"] = attr["long_name"]+" and "+ER_label+" (ER)"
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

def ERUsingFFNET(cf,ds):
    """
    Purpose:
     Estimate ecosystem respiration using the ffnet neural network.
    Usage:
     qcrp.ERUsingFFNET(cf,ds)
      where cf is a control file object
            ds is a data structure
    Author: PRI
    Date: August 2014
    """
    if "ffnet" not in dir(ds): return
    if "ffnet" not in sys.modules.keys():
        log.error("ERUsingFFNET: I don't think ffnet is installed ...")
        return
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    FFNET_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                    "file_enddate":enddate.strftime("%Y-%m-%d %H:%M"),
                    "plot_path":cf["Files"]["plot_path"]}
    # check to see if this is a batch or an interactive run
    call_mode = qcutils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive")
    FFNET_info["call_mode"]= call_mode
    if call_mode.lower()=="interactive":
        FFNET_info["show_plots"] = True
        # call the FFNET GUI
        qcrpNN.rpFFNET_gui(cf,ds,FFNET_info)
    else:
        if "GUI" in cf:
            if "FFNET" in cf["GUI"]:
                qcrpNN.rpFFNET_run_nogui(cf,ds,FFNET_info)

def ERUsingLasslop(cf,ds):
    log.info('Estimating ER using Lasslop et al is not implemented yet, but we are working on it now ...')
    pass
    # get necessary data
    #  - PAR or Fsd
    #    - convert Fsd to PAR
    #    - check PAR units are umol/m2/s
    #  - soil or air temperature
    #
    # get a list of years in the data set
    #
    # loop over entries in 
    # loop over years to get annual E0 values
    #  - if year has > minimum good points
    #    - curve fit to get E0
    #  - else
    #    - set E0 to missing value
    # replace any missing E0 with mean values
    #  E0_dict = rpLL_get_annual_E0()
    #
    # loop through data set with window_size and window_step to get
    # curve fit parameters phi, Aopt, k and rb
    #  fit_dict = rpLL_get_fit_parameters()
    #   fit_dict["middate"],fit_dict["phi"],fit_dict["Aopt"],fit_dict["k"],fit_dict["E0_short"],fit_dict["E0_long"]
    #dt = ds.series["DateTime"]["Data"]
    #radn,f,a = qcutils.GetSeriesasMA(ds,"Fsd")
    #hd_label = cf
    #hd,f,a = qcutils.GetSeriesasMA(ds,hd_label)
    #startdate = dt[0]
    #enddate = startdate + dateutil.relativedelta.relativedelta(days=window_size)
    #while startdate>dt[-1]:
        #middate = startdate+(enddate-startdate)/2
        #E0 = E0_dict[middate.year]
        #si = qcutils.GetDateIndex(dt,str(startdate),ts=ts)
        #ei = qcutils.GetDateIndex(dt,str(enddate),ts=ts)
    # plot phi, Aopt, k and rb values
    #  rpLL_plot_fit_parameters(fit_dict)
    # interpolate phi, Aopt, k and rb values to daily time step
    #  fit_dict_daily = rpLL_interpolate_fit_parameters(fit_dict)
    # calculate ER_LL
    #  rpLL_calculateER()
    # - replicate daily values of phi, Aopt, k and rb at data time step
    # - calculate ER from temperature, rb and E0
    # - put ER in data structure

def ERUsingLloydTaylor(cf,ds):
    log.info('Estimating ER using Lloyd-Taylor is not implemented yet')
    pass
    ## need to check that all information required is in the control file
    #section = qcutils.get_cfsection(cf,series='ER',mode='quiet')
    #if len(section)==0:
        #log.error('ER section not found in control file')
        #return
    #if 'LloydTaylor' not in cf[section]['ER']:
        #log.error('LloydTaylor section not found in control file')
        #return
    ## get the driver
    #if 'drivers' not in cf[section]['ER']['LloydTaylor']:
        #log.error('drivers key not found in LloydTaylor section of control file')
        #return
    #else:
        #driver_list = eval(cf[section]['ER']['LloydTaylor']['drivers'])
        #driver_label = driver_list[0]
    ## get the monthly values for the activation energy, E0
    #if 'E0' not in cf[section]['ER']['LloydTaylor']:
        #log.error('E0 key not found in LloydTaylor section of control file')
        #return
    #else:
        #E0_monthly = numpy.array(eval(cf[section]['ER']['LloydTaylor']['E0']))
    ## get the monthly values for the base respiration, rb
    #if 'rb' not in cf[section]['ER']['LloydTaylor']:
        #log.error('rb key not found in LloydTaylor section of control file')
        #return
    #else:
        #rb_monthly = numpy.array(eval(cf[section]['ER']['LloydTaylor']['rb']))
    ## get the output label
    #if 'output' not in cf[section]['ER']['LloydTaylor']:
        #log.error('output key not found in LloydTaylor section of control file')
        #return
    #else:
        #out_label = cf[section]['ER']['LloydTaylor']['output']
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
    ## estimate ER using the Lloyd-Taylor expression
    #t1 = 1/(c.Tref-c.T0)
    #t2 = 1/(Ts-c.T0)
    #ER = rb*numpy.exp(E0*(t1-t2))
    ## put the estimated ER into the data structure
    #units=qcutils.GetUnitsFromds(ds, 'Fc')
    #attr = qcutils.MakeAttributeDictionary(long_name='ER estimated using Lloyd-Taylor',units=units)
    #qcutils.CreateSeries(ds,out_label,ER,Flag=flag,Attr=attr)

def ERUsingSOLO(cf,ds):
    """ Estimate ER using SOLO. """
    if "solo" not in dir(ds): return
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    solo_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                   "file_enddate":enddate.strftime("%Y-%m-%d %H:%M"),
                   "plot_path":cf["Files"]["plot_path"]}
    # check to see if this is a batch or an interactive run
    call_mode = qcutils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive")
    solo_info["call_mode"]= call_mode
    if call_mode.lower()=="interactive": solo_info["show_plots"] = True
    if call_mode.lower()=="interactive":
        # call the ERUsingSOLO GUI
        qcrpNN.rpSOLO_gui(cf,ds,solo_info)
    else:
        if "GUI" in cf:
            if "SOLO" in cf["GUI"]:
                qcrpNN.rpSOLO_run_nogui(cf,ds,solo_info)

def GetERFromFc(cf,ds):
    """
    Purpose:
     Get the observed ecosystem respiration from measurements of Fc by
     filtering out daytime periods and periods when ustar is less than
     a threshold value.
     The Fsd threshold for determining day time and night time and the
     ustar threshold are set in the [Params] section of the L5 control
     file.
    Usage:
     qcrp.GetERFromFc(cf,ds)
     where cf is a control file object
           ds is a data structure
    Side effects:
     A new series called "ER" is created in the data structure.
    Author: PRI
    Date: August 2014
    """
    # needs a fecking good refactor
    ts = int(ds.globalattributes["time_step"])
    ldt = ds.series['DateTime']['Data']
    # get the Fsd threshold
    if "Options" in cf.keys():
        if "Fsd_threshold" in cf["Options"].keys():
            Fsd_threshold = float(cf["Options"]["Fsd_threshold"])
        else:
            log.warning(" No Fsd threshold in [Options] section of control file")
            log.warning(" ... using default value of 10 W/m2")
            Fsd_threshold = float(10)
    else:
        log.warning(" No [Options] section of control file for Fsd threshold")
        log.warning(" ... using default value of 10 W/m2")
        Fsd_threshold = float(10)
    # get the ustar thresholds
    if "cpd_filename" in cf["Files"]:
        ustar_dict = get_ustarthreshold_from_cpdresults(cf)
    else:
        msg = " CPD results filename not in control file"
        log.warning(msg)
        ustar_dict = get_ustarthreshold_from_cf(cf,ldt)
    # make sure we have an entry in ustar_dict for all years
    start_year = ldt[0].year
    end_year = ldt[-1].year
    data_years = range(start_year,end_year+1)
    ustar_years = ustar_dict.keys()
    ustar_list = ustar_dict[ustar_years[0]]
    for year in data_years:
        if str(year) not in ustar_years:
            ustar_dict[str(year)] = {}
            for item in ustar_list:
                ustar_dict[str(year)][item] = float(c.missing_value)
    # get the data
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" not in ds.series.keys(): qcts.get_synthetic_fsd(ds)
    Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
    sa,flag,attr = qcutils.GetSeriesasMA(ds,"solar_altitude")
    ustar,ustar_flag,attr = qcutils.GetSeriesasMA(ds,"ustar")
    Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,"Fc")
    # get a copy of the Fc flag
    ER_flag = numpy.array(Fc_flag)

    # only accept Fc and ustar data when both have a QC flag value of 0
    ustar = numpy.ma.masked_where((ustar_flag!=0)|(Fc_flag!=0),ustar)
    Fc = numpy.ma.masked_where((ustar_flag!=0)|(Fc_flag!=0),Fc)
    index_notok = numpy.where((ustar_flag!=0)|(Fc_flag!=0))[0]
    #ustar_flag[index_notok] = numpy.int32(61)
    ER_flag[index_notok] = numpy.int32(61)
    # check for any missing data
    for item,label in zip([Fsd,Fsd_syn],["Fsd","Fsd_syn"]):
        index = numpy.where(numpy.ma.getmaskarray(item)==True)[0]
        if len(index)!=0:
            log.error(" GetERFromFc: missing data in series "+label)
            raise Exception("GetERFromFc: missing data in series "+label)

    # apply the day/night filter
    # get the day/night filter type from the control file
    daynightfilter_type = qcutils.get_keyvaluefromcf(cf,["Options"],"DayNightFilter",default="Fsd")
    # trap any types not implemented and set to Fsd
    if daynightfilter_type not in ["Fsd","sa"]: daynightfilter_type = "Fsd"
    # make the attribute dictionary first so we can add the ustar thresholds to it
    ER_attr = qcutils.MakeAttributeDictionary(long_name='Ecosystem respiration (observed)',units=Fc_attr["units"])
    Fsd_attr["long_name"] = "Incoming shortwave radiation, filtered"
    Fsd_attr["units"] = "W/m2"

    # apply the day/night filter
    if daynightfilter_type=="Fsd":
        # we are using Fsd and possibly Fsd_syn to define day/night
        ER_attr["Fsd_threshold"] = str(Fsd_threshold)
        #Fsd_attr["Fsd_threshold"] = str(Fsd_threshold)
        if qcutils.cfoptionskeylogical(cf,Key='UseFsdsyn_threshold',default=False):
            # we are using Fsd and Fsd_syn
            ER1 = numpy.ma.masked_where((Fsd>Fsd_threshold)|(Fsd_syn>Fsd_threshold),Fc,copy=True)
            Fsd1 = numpy.ma.masked_where((Fsd>Fsd_threshold)|(Fsd_syn>Fsd_threshold),Fsd,copy=True)
            index_daynight = numpy.ma.where((Fsd>Fsd_threshold)|(Fsd_syn>Fsd_threshold))[0]
            ER_flag[index_daynight] = numpy.int32(62)
        else:
            # we are only using Fsd
            ER1 = numpy.ma.masked_where(Fsd>Fsd_threshold,Fc,copy=True)
            Fsd1 = numpy.ma.masked_where(Fsd>Fsd_threshold,Fsd,copy=True)
            index_daynight = numpy.ma.where(Fsd>Fsd_threshold)[0]
            ER_flag[index_daynight] = numpy.int32(62)
    else:
        sa_threshold = int(qcutils.get_keyvaluefromcf(cf,["Options"],"sa_threshold",default="-5"))
        ER_attr["sa_threshold"] = str(sa_threshold)
        Fsd_attr["sa_threshold"] = str(sa_threshold)
        ER1 = numpy.ma.masked_where(sa>sa_threshold,Fc,copy=True)
        Fsd1 = numpy.ma.masked_where(sa>sa_threshold,Fsd,copy=True)
        index_daynight = numpy.ma.where(sa>sa_threshold)[0]
        ER_flag[index_daynight] = numpy.int32(63)
    # get a copy of the day/night filtered data
    ER2 = numpy.ma.array(ER1)
    Fsd2 = numpy.ma.array(Fsd1)

    # loop over the list of ustar thresholds
    year_list = ustar_dict.keys()
    year_list.sort()
    # get the average of good ustar threshold values
    good_values = []
    for year in year_list:
        ustar_threshold = float(ustar_dict[year]["ustar_mean"])
        if ustar_threshold!=float(c.missing_value):
            good_values.append(ustar_threshold)
    ustar_threshold_mean = numpy.sum(numpy.array(good_values))/len(good_values)
    # now loop over the years in the data to apply the ustar threshold
    for year in year_list:
        start_date = str(year)+"-01-01 00:30"
        if ts==60: start_date = str(year)+"-01-01 01:00"
        end_date = str(int(year)+1)+"-01-01 00:00"
        # get the ustar threshold
        ustar_threshold = float(ustar_dict[year]["ustar_mean"])
        if ustar_threshold==float(c.missing_value): ustar_threshold = ustar_threshold_mean
        ER_attr["ustar_threshold_"+str(year)] = str(ustar_threshold)
        # get the start and end datetime indices
        si = qcutils.GetDateIndex(ldt,start_date,ts=ts,default=0,match='exact')
        ei = qcutils.GetDateIndex(ldt,end_date,ts=ts,default=len(ldt),match='exact')
        # filter out the low ustar conditions
        ER2[si:ei] = numpy.ma.masked_where(ustar[si:ei]<ustar_threshold,ER1[si:ei])
        Fsd2[si:ei] = numpy.ma.masked_where(ustar[si:ei]<ustar_threshold,Fsd1[si:ei])
        # set the QC flag
        index_lowustar = numpy.ma.where(ustar[si:ei]<ustar_threshold)[0]
        ER_flag[si:ei][index_lowustar] = numpy.int32(64)

    # apply quantile filter
    if qcutils.cfoptionskeylogical(cf,Key='UseQuantileFilter',default=False):
        ER_attr["long_name"] = ER_attr["long_name"]+", quantile filter not used"
        Fsd_attr["long_name"] = Fsd_attr["long_name"]+", quantile filter not used"
        qcutils.CreateSeries(ds,"ER_nqf",ER2,Flag=ER_flag,Attr=ER_attr)
        qcutils.CreateSeries(ds,"Fsd_nqf",Fsd2,Flag=Fsd_flag,Attr=Fsd_attr)
        quantile_lower = float(qcutils.get_keyvaluefromcf(cf,["Options"],"QuantileValue",default="2.5"))
        quantile_upper = float(100) - quantile_lower
        q = numpy.percentile(numpy.ma.compressed(ER2),[quantile_lower,quantile_upper])
        ER2 = numpy.ma.masked_where((ER2<q[0])|(ER2>q[1]),ER2)
        Fsd2 = numpy.ma.masked_where((ER2<q[0])|(ER2>q[1]),Fsd2)
        index_qf = numpy.ma.where((ER2<q[0])|(ER2>q[1]))[0]
        ER_flag[index_qf] = numpy.int32(65)
        Fsd_flag[index_qf] = numpy.int32(65)
        ER_attr["long_name"].replace(", quantile filter not used",", quantile filter used")
        ER_attr["ER_quantile"] = str(quantile_lower)+","+str(quantile_upper)
        Fsd_attr["long_name"].replace(", quantile filter not used",", quantile filter used")
        Fsd_attr["ER_quantile"] = str(quantile_lower)+","+str(quantile_upper)

    # put the nocturnal, filtered Fc data into the data structure
    qcutils.CreateSeries(ds,"ER",ER2,Flag=ER_flag,Attr=ER_attr)
    qcutils.CreateSeries(ds,"Fsd_filtered",Fsd2,Flag=Fsd_flag,Attr=Fsd_attr)
    return

def GetERIndicator(cf,ds):
    """
    Purpose:
     Indicator values are:
      - 1 OK to use Fc as ER
      - 0 not OK to use Fc as ER
    Useage:
    Author: PRI
    Date: August 2015
    """
    # get the day/night indicator
    # 1 ==> night, 0 ==> day
    daynight_indicator = get_daynight_indicator(cf,ds)
    # get the evening indicator
    # 1 ==> within "evening" definition, 0 ==> outside "evening"
    evening_indicator = get_evening_indicator(cf,ds)
    # get the turbulent/non-turbulent indicator
    # 1 ==> turbulent, 0 ==> non-turbulent
    turbulence_indicator = get_turbulence_indicator(cf,ds)
    # get ER indicator
    # 1 ==> OK to use Fc as ER, 0 ==> not OK to use Fc as ER
    ER_indicator = daynight_indicator*evening_indicator*turbulence_indicator
    # put the indicator series in the data structure
    nRecs = len(ER_indicator)
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    attr = qcutils.MakeAttributeDictionary(long_name="ER indicator series")
    qcutils.CreateSeries(ds,"ER_indicator",ER_indicator,Flag=flag,Attr=attr)

def get_daynight_indicator(cf,ds):
    # get the day/night indicator
    nRecs = int(ds.globalattributes["nc_nrecs"])
    daynight_indicator = numpy.ones(nRecs,dtype=numpy.int32)
    # get the filter type
    filter_type = qcutils.get_keyvaluefromcf(cf,["Options"],"DayNightFilter",default="Fsd")
    # get the indicator series
    if filter_type.lower()=="fsd":
        # get the data
        Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        # get the Fsd threshold
        Fsd_threshold = int(qcutils.get_keyvaluefromcf(cf,["Options"],"Fsd_threshold",default=10))
        # we are using Fsd and Fsd_syn to define day/night
        index = numpy.ma.where((Fsd>Fsd_threshold)|(Fsd_syn>Fsd_threshold))[0]
        daynight_indicator[index] = numpy.int32(0)
    elif filter_type.lower()=="sa":
        # get the data
        sa,flag,attr = qcutils.GetSeriesasMA(ds,"solar_altitude")
        # get the solar altitude threshold
        sa_threshold = int(qcutils.get_keyvaluefromcf(cf,["Options"],"sa_threshold",default="-5"))
        # we are using solar altitude to define day/night
        index = numpy.ma.where(sa>sa_threshold)[0]
        daynight_indicator[index] = numpy.int32(0)
    else:
        msg = "Unrecognised DayNightFilter option in L6 control file, no filter applied ..."
        log.error(msg)
    return daynight_indicator

def get_evening_indicator(cf,ds):
    # make sure we have the synthetic downwelling shortwave and the solar altitude
    if ("Fsd_syn" not in ds.series.keys() or
        "solar_altitude" not in ds.series.keys()): qcts.get_synthetic_fsd(ds)
    
def get_turbulence_indicator(cf,ds):
    pass

def get_ustarthreshold_from_cf(cf,ldt):
    """
    Purpose:
     Returns a dictionary containing ustar thresholds for each year read from
     the control file.  If no [ustar_threshold] section is found then a
     default value of 0.25 is used.
    Usage:
     ustar_dict = qcrp.get_ustarthreshold_from_cf(cf,ldt)
     where cf is the control file object
           ldt is the Python datetime series from the data structure
    Author: PRI
    Date: July 2015
    """
    ustar_dict = collections.OrderedDict()
    ustar_threshold_list = []
    if "ustar_threshold" in cf.keys():
        msg = " Using values from ustar_threshold section"
        log.info(msg)
        for n in cf["ustar_threshold"].keys():
            ustar_threshold_list.append(ast.literal_eval(cf["ustar_threshold"][str(n)]))
        for item in ustar_threshold_list:
            startdate = dateutil.parser.parse(item[0])
            year = startdate.year
            ustar_dict[str(year)] = {}
            ustar_dict[str(year)]["ustar_mean"] = float(item[2])
    else:
        log.error(" No [ustar_threshold] section in control file")
        log.error(" ... using default value of 0.25 m/s")
        startyear = ldt[0].year
        endyear = ldt[-1].year
        years = range(startyear,endyear+1)
        for year in years:
            ustar_dict[str(year)] = {}
            ustar_dict[str(year)]["ustar_mean"] = float(0.25)
    return ustar_dict

def get_ustarthreshold_from_cpdresults(cf):
    # do some stuff
    cpd_path = cf["Files"]["file_path"]
    cpd_name = cpd_path+cf["Files"]["cpd_filename"]
    cpd_wb = xlrd.open_workbook(cpd_name)
    annual_ws = cpd_wb.sheet_by_name("Annual")
    header_list = [x for x in annual_ws.row_values(0)]
    year_list = [str(int(x)) for x in annual_ws.col_values(0)[1:]]
    ustar_dict = collections.OrderedDict()
    for i,year in enumerate(year_list):
        ustar_dict[year] = collections.OrderedDict()
        for item in header_list:
            xlcol = header_list.index(item)
            val = annual_ws.col_values(xlcol)[i+1]
            typ = annual_ws.col_types(xlcol)[i+1]
            if typ==2:
                ustar_dict[year][item] = float(val)
            else:
                ustar_dict[year][item] = float(c.missing_value)
    return ustar_dict

def L6_summary(cf,ds):
    """
    Purpose:
     Produce summaries of L6 data, write them to an Excel spreadsheet and plot them.
    Usage:
    Author: PRI
    Date: June 2015
    """
    log.info(" Doing the L6 summary")
    # set up a dictionary of lists
    series_dict = L6_summary_createseriesdict(cf,ds)
    # open the Excel workbook
    nc_name = qcio.get_outfilenamefromcf(cf)
    xl_name = nc_name.replace(".nc","_Summary.xls")
    xl_file = qcio.xl_open_write(xl_name)
    if xl_file=='':
        log.error("L6_summary: error opening Excel file "+xl_name)
        return
    # daily averages and totals
    daily_dict = L6_summary_daily(ds,series_dict)
    L6_summary_write_xlfile(xl_file,"Daily",daily_dict)
    # monthly averages and totals
    monthly_dict = L6_summary_monthly(ds,series_dict)
    L6_summary_write_xlfile(xl_file,"Monthly",monthly_dict)
    # annual averages and totals
    annual_dict = L6_summary_annual(ds,series_dict)
    L6_summary_write_xlfile(xl_file,"Annual",annual_dict)
    # cumulative totals
    cumulative_dict = L6_summary_cumulative(ds,series_dict)
    for year in cumulative_dict.keys():
        L6_summary_write_xlfile(xl_file,"Cummulative("+str(year)+")",cumulative_dict[str(year)])
    # close the Excel workbook
    xl_file.save(xl_name)
    # plot the daily averages and sums
    L6_summary_plotdaily(cf,ds,daily_dict)
    # plot the cumulative sums
    L6_summary_plotcumulative(cf,ds,cumulative_dict)

def L6_summary_plotdaily(cf,ds,daily_dict):
    """
    Purpose:
     Plot the daily averages or sums with a 30 day filter.
    Usage:
     L6_summary_plotdaily(daily_dict)
     where daily_dict is the dictionary of results returned by L6_summary_daily
    Author: PRI
    Date: June 2015
    """
    type_list = []
    for item in daily_dict.keys():
        if item[0:2]=="ER": type_list.append(item[2:])
    for item in type_list:
        if "NEE"+item not in daily_dict or "GPP"+item not in daily_dict:
            type_list.remove(item)
    # plot time series of NEE, GPP and ER
    sdate = daily_dict["DateTime"]["data"][0].strftime("%d-%m-%Y")
    edate = daily_dict["DateTime"]["data"][-1].strftime("%d-%m-%Y")
    site_name = ds.globalattributes["site_name"]
    title_str = site_name+": "+sdate+" to "+edate
    for item in type_list:
        if cf["Options"]["call_mode"].lower()=="interactive":
            plt.ion()
        else:
            plt.ioff()
        fig = plt.figure(figsize=(16,4))
        fig.canvas.set_window_title("Carbon Budget: "+item.replace("_",""))
        plt.figtext(0.5,0.95,title_str,horizontalalignment='center')
        plt.plot(daily_dict["DateTime"]["data"],daily_dict["NEE"+item]["data"],'b-',alpha=0.3)
        plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["NEE"+item]["data"],window_len=30),
                 'b-',linewidth=2,label="NEE"+item+" (30 day filter)")
        plt.plot(daily_dict["DateTime"]["data"],daily_dict["GPP"+item]["data"],'g-',alpha=0.3)
        plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["GPP"+item]["data"],window_len=30),
                 'g-',linewidth=2,label="GPP"+item+" (30 day filter)")
        plt.plot(daily_dict["DateTime"]["data"],daily_dict["ER"+item]["data"],'r-',alpha=0.3)
        plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["ER"+item]["data"],window_len=30),
                 'r-',linewidth=2,label="ER"+item+" (30 day filter)")
        plt.axhline(0)
        plt.xlabel("Date")
        plt.ylabel(daily_dict["NEE"+item]["units"])
        plt.legend(loc='upper left',prop={'size':8})
        plt.tight_layout()
        sdt = daily_dict["DateTime"]["data"][0].strftime("%Y%m%d")
        edt = daily_dict["DateTime"]["data"][-1].strftime("%Y%m%d")
        plot_path = cf["Files"]["plot_path"]+"L6/"
        if not os.path.exists(plot_path): os.makedirs(plot_path)
        figname = plot_path+site_name.replace(" ","")+"_CarbonBudget"+item
        figname = figname+"_"+sdt+"_"+edt+'.png'
        fig.savefig(figname,format='png')
        if cf["Options"]["call_mode"].lower=="interactive":
            plt.draw()
            plt.ioff()
        else:
            plt.ion()
    # plot time series of Fn,Fg,Fh,Fe
    if cf["Options"]["call_mode"].lower()=="interactive":
        plt.ion()
    else:
        plt.ioff()
    fig = plt.figure(figsize=(16,4))
    fig.canvas.set_window_title("Surface Energy Budget")
    plt.figtext(0.5,0.95,title_str,horizontalalignment='center')
    plt.plot(daily_dict["DateTime"]["data"],daily_dict["Fn"]["data"],'k-',alpha=0.3)
    plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["Fn"]["data"],window_len=30),
             'k-',linewidth=2,label="Fn (30 day filter)")
    plt.plot(daily_dict["DateTime"]["data"],daily_dict["Fg"]["data"],'g-',alpha=0.3)
    plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["Fg"]["data"],window_len=30),
             'g-',linewidth=2,label="Fg (30 day filter)")
    plt.plot(daily_dict["DateTime"]["data"],daily_dict["Fh"]["data"],'r-',alpha=0.3)
    plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["Fh"]["data"],window_len=30),
             'r-',linewidth=2,label="Fh (30 day filter)")
    plt.plot(daily_dict["DateTime"]["data"],daily_dict["Fe"]["data"],'b-',alpha=0.3)
    plt.plot(daily_dict["DateTime"]["data"],qcts.smooth(daily_dict["Fe"]["data"],window_len=30),
             'b-',linewidth=2,label="Fe (30 day filter)")
    plt.xlabel("Date")
    plt.ylabel(daily_dict["Fn"]["units"])
    plt.legend(loc='upper left',prop={'size':8})
    plt.tight_layout()
    sdt = daily_dict["DateTime"]["data"][0].strftime("%Y%m%d")
    edt = daily_dict["DateTime"]["data"][-1].strftime("%Y%m%d")
    plot_path = cf["Files"]["plot_path"]+"L6/"
    if not os.path.exists(plot_path): os.makedirs(plot_path)
    figname = plot_path+site_name.replace(" ","")+"_SEB"
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    if cf["Options"]["call_mode"].lower=="interactive":
        plt.draw()
        plt.ioff()
    else:
        plt.ion()

def L6_summary_plotcumulative(cf,ds,cumulative_dict):
    ts = int(ds.globalattributes["time_step"])
    # cumulative plots
    color_list = ["blue","red","green","yellow","magenta","black","cyan","brown"]
    year_list = cumulative_dict.keys()
    year_list.sort()
    type_list = []
    for item in cumulative_dict[year_list[0]].keys():
        if item[0:2]=="ER": type_list.append(item[2:])
    for item in type_list:
        if "NEE"+item not in cumulative_dict[year_list[0]] or "GPP"+item not in cumulative_dict[year_list[0]]:
            type_list.remove(item)
    # do the plots
    site_name = ds.globalattributes["site_name"]
    title_str = site_name+": "+year_list[0]+" to "+year_list[-1]
    for item in type_list:
        if cf["Options"]["call_mode"].lower()=="interactive":
            plt.ion()
        else:
            plt.ioff()
        fig = plt.figure(figsize=(12,12))
        fig.canvas.set_window_title("Cumulative plots: "+item.replace("_",""))
        plt.suptitle(title_str)
        plt.subplot(221)
        plt.title("NEE: "+item.replace("_",""),fontsize=12)
        for n,year in enumerate(year_list):
            x = numpy.arange(0,len(cumulative_dict[year]["NEE"+item]["data"]))*ts/float(60)
            plt.plot(x,cumulative_dict[year]["NEE"+item]["data"],color=color_list[numpy.mod(n,8)],
                     label=str(year))
            plt.xlabel("Hour of Year")
            plt.ylabel(cumulative_dict[year]["NEE"+item]["units"])
            plt.legend(loc='lower left',prop={'size':8})
        plt.subplot(222)
        plt.title("GPP: "+item.replace("_",""),fontsize=12)
        for n,year in enumerate(year_list):
            x = numpy.arange(0,len(cumulative_dict[year]["GPP"+item]["data"]))*ts/float(60)
            plt.plot(x,cumulative_dict[year]["GPP"+item]["data"],color=color_list[numpy.mod(n,8)],
                     label=str(year))
            plt.xlabel("Hour of Year")
            plt.ylabel(cumulative_dict[year]["GPP"+item]["units"])
            plt.legend(loc='lower right',prop={'size':8})
        plt.subplot(223)
        plt.title("ER: "+item.replace("_",""),fontsize=12)
        for n,year in enumerate(year_list):
            x = numpy.arange(0,len(cumulative_dict[year]["ER"+item]["data"]))*ts/float(60)
            plt.plot(x,cumulative_dict[year]["ER"+item]["data"],color=color_list[numpy.mod(n,8)],
                     label=str(year))
            plt.xlabel("Hour of Year")
            plt.ylabel(cumulative_dict[year]["ER"+item]["units"])
            plt.legend(loc='lower right',prop={'size':8})
        plt.subplot(224)
        plt.title("ET & Precip",fontsize=12)
        for n,year in enumerate(year_list):
            x = numpy.arange(0,len(cumulative_dict[year]["ET"]["data"]))*ts/float(60)
            plt.plot(x,cumulative_dict[year]["ET"]["data"],color=color_list[numpy.mod(n,8)],
                     label=str(year))
            plt.plot(x,cumulative_dict[year]["Precip"]["data"],color=color_list[numpy.mod(n,8)],
                     linestyle='--')
            plt.xlabel("Hour of Year")
            plt.ylabel(cumulative_dict[year]["ET"]["units"])
            plt.legend(loc='upper left',prop={'size':8})
        plt.tight_layout(rect=[0, 0, 1, 0.98])
        # save a hard copy of the plot
        sdt = year_list[0]
        edt = year_list[-1]
        plot_path = cf["Files"]["plot_path"]+"L6/"
        if not os.path.exists(plot_path): os.makedirs(plot_path)
        figname = plot_path+site_name.replace(" ","")+"_Cumulative_"+item.replace("_","")
        figname = figname+"_"+sdt+"_"+edt+'.png'
        fig.savefig(figname,format='png')
        if cf["Options"]["call_mode"].lower=="interactive":
            plt.draw()
            plt.ioff()
        else:
            plt.ion()

def L6_summary_createseriesdict(cf,ds):
    """
    Purpose:
     Create a dictionary containing lists of variables, operators and formats
    for use by the daily, annual and cumulative routines.
    Usage:
     series_dict = L6_summary_createseriesdict(cf,ds)
     where cf is a control file object
           ds is an OzFluxQC data structure
           series_dict is a dictionary of various variable lists
    Author: PRI
    Date: June 2015
    """
    ts = int(ds.globalattributes["time_step"])
    series_dict = {"daily":{},"annual":{},"cumulative":{},"lists":{}}
    # adjust units of NEE, NEP, GPP and ER
    sdl = series_dict["lists"]
    sdl["nee"] = [item for item in cf["NEE"].keys() if "NEE" in item and item in ds.series.keys()]
    sdl["gpp"] = [item for item in cf["GPP"].keys() if "GPP" in item and item in ds.series.keys()]
    sdl["fre"] = [item for item in cf["ER"].keys() if "ER" in item and item in ds.series.keys()]
    sdl["nep"] = [item.replace("NEE","NEP") for item in sdl["nee"]]
    sdl["co2"] = sdl["nee"]+sdl["nep"]+sdl["gpp"]+sdl["fre"]
    for item in sdl["co2"]:
        series_dict["daily"][item] = {}
        series_dict["cumulative"][item] = {}
        series_dict["daily"][item]["operator"] = "sum"
        series_dict["daily"][item]["format"] = "0.00"
        series_dict["cumulative"][item]["operator"] = "sum"
        series_dict["cumulative"][item]["format"] = "0.00"
    series_dict["daily"]["Ah"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["Cc"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fc"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["Fe"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fh"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fg"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fn"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fsd"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fsu"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Fld"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["Flu"] = {"operator":"average","format":"0.0"}
    series_dict["daily"]["ps"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["q"] = {"operator":"average","format":"0.0000"}
    series_dict["daily"]["RH"] = {"operator":"average","format":"0"}
    series_dict["daily"]["Sws"] = {"operator":"average","format":"0.000"}
    series_dict["daily"]["Ta"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["Ts"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["ustar"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["Ws"] = {"operator":"average","format":"0.00"}
    series_dict["daily"]["ET"] = {"operator":"sum","format":"0.0"}
    series_dict["daily"]["Precip"] = {"operator":"sum","format":"0.0"}
    series_dict["cumulative"]["ET"] = series_dict["daily"]["ET"]
    series_dict["cumulative"]["Precip"] = series_dict["daily"]["Precip"]
    series_dict["annual"] = series_dict["daily"]
    series_dict["monthly"] = series_dict["daily"]
    return series_dict

def L6_summary_daily(ds,series_dict,xl_file=None):
    """
    Purpose:
     Calculate the daily averages or sums of various quantities and write
     them to a worksheet in an Excel workbook.
    Usage:
     L6_summary_daily(xl_file,ds,series_dict)
     where xl_file is an Excel file object
           ds is an OzFluxQC data structure
           series_dict is a dictionary of various variable lists
    Author: PRI
    Date: June 2015
    """
    log.info(" Doing the daily summary at L6")
    dt = ds.series["DateTime"]["Data"]
    ts = int(ds.globalattributes["time_step"])
    si = qcutils.GetDateIndex(dt,str(dt[0]),ts=ts,default=0,match="startnextday")
    ei = qcutils.GetDateIndex(dt,str(dt[-1]),ts=ts,default=len(dt)-1,match="endpreviousday")
    ldt = dt[si:ei+1]
    ntsInDay = int(24.0*60.0/float(ts))
    nDays = int(len(ldt))/ntsInDay
    ldt_daily = [ldt[0]+datetime.timedelta(days=i) for i in range(0,nDays)]
    daily_dict = {}
    daily_dict["DateTime"] = {"data":ldt_daily,"units":"Days","format":"dd/mm/yyyy"}
    series_list = series_dict["daily"].keys()
    series_list.sort()
    for item in series_list:
        if item not in ds.series.keys(): continue
        daily_dict[item] = {}
        data_1d,flag,attr = qcutils.GetSeriesasMA(ds,item,si=si,ei=ei)
        if item in series_dict["lists"]["co2"]:
            data_1d = qcutils.convertunits(data_1d,attr["units"],"gC/m2",ts)
            daily_dict[item]["units"] = "gC/m2"
        else:
            daily_dict[item]["units"] = attr["units"]
        data_2d = data_1d.reshape(nDays,ntsInDay)
        if series_dict["daily"][item]["operator"].lower()=="average":
            daily_dict[item]["data"] = numpy.ma.average(data_2d,axis=1)
        elif series_dict["daily"][item]["operator"].lower()=="sum":
            daily_dict[item]["data"] = numpy.ma.sum(data_2d,axis=1)
            daily_dict[item]["units"] = daily_dict[item]["units"]+"/day"
        else:
            print "unrecognised series: ",series
        daily_dict[item]["format"] = series_dict["daily"][item]["format"]
    return daily_dict

def L6_summary_write_xlfile(xl_file,sheet_name,data_dict):
    # add the daily worksheet to the summary Excel file
    xl_sheet = xl_file.add_sheet(sheet_name)
    qcio.xl_write_data(xl_sheet,data_dict)

def L6_summary_monthly(ds,series_dict):
    """
    Purpose:
     Calculate the monthly averages or sums of various quantities and write
     them to a worksheet in an Excel workbook.
    Usage:
     L6_summary_monthly(xl_file,ds,series_dict)
     where xl_file is an Excel file object
           ds is an OzFluxQC data structure
           series_dict is a dictionary of various variable lists
    Author: PRI
    Date: July 2015
    """
    log.info(" Doing the monthly summaries at L6")
    dt = ds.series["DateTime"]["Data"]
    ts = int(ds.globalattributes["time_step"])
    si = qcutils.GetDateIndex(dt,str(dt[0]),ts=ts,default=0,match="startnextmonth")
    ldt = dt[si:]
    monthly_dict = {}
    monthly_dict["DateTime"] = {"data":[],"units":"Months","format":"dd/mm/yyyy"}
    # create arrays in monthly_dict
    series_list = series_dict["monthly"].keys()
    series_list.sort()
    # create the data arrays
    for item in series_list:
        monthly_dict[item] = {"data":numpy.array([])}
    # loop over the months in the data file
    start_date = ldt[0]
    end_date = start_date+dateutil.relativedelta.relativedelta(months=1)
    end_date = end_date-dateutil.relativedelta.relativedelta(minutes=ts)
    last_date = ldt[-1]
    while start_date<=last_date:
        si = qcutils.GetDateIndex(ldt,str(start_date),ts=ts,default=0)
        ei = qcutils.GetDateIndex(ldt,str(end_date),ts=ts,default=len(ldt)-1)
        monthly_dict["DateTime"]["data"].append(ldt[si])
        for item in series_list:
            if item not in ds.series.keys(): continue
            data_1d,flag,attr = qcutils.GetSeriesasMA(ds,item,si=si,ei=ei)
            if item in series_dict["lists"]["co2"]:
                data_1d = qcutils.convertunits(data_1d,attr["units"],"gC/m2",ts)
                monthly_dict[item]["units"] = "gC/m2"
            else:
                monthly_dict[item]["units"] = attr["units"]
            if series_dict["monthly"][item]["operator"].lower()=="average":
                monthly_dict[item]["data"] = numpy.append(monthly_dict[item]["data"],
                                                          numpy.ma.average(data_1d))
            elif series_dict["monthly"][item]["operator"].lower()=="sum":
                monthly_dict[item]["data"] = numpy.append(monthly_dict[item]["data"],
                                                          numpy.ma.sum(data_1d))
                monthly_dict[item]["units"] = monthly_dict[item]["units"]+"/month"
            else:
                print "unrecognised operator"
            monthly_dict[item]["format"] = series_dict["monthly"][item]["format"]
        start_date = end_date+dateutil.relativedelta.relativedelta(minutes=ts)
        end_date = start_date+dateutil.relativedelta.relativedelta(months=1)
        end_date = end_date-dateutil.relativedelta.relativedelta(minutes=ts)
    return monthly_dict

def L6_summary_annual(ds,series_dict):
    """
    Purpose:
     Calculate the annual averages or sums of various quantities and write
     them to a worksheet in an Excel workbook.
    Usage:
     L6_summary_annual(xl_file,ds,series_dict)
     where xl_file is an Excel file object
           ds is an OzFluxQC data structure
           series_dict is a dictionary of various variable lists
    Author: PRI
    Date: June 2015
    """
    log.info(" Doing the annual summaries at L6")
    dt = ds.series["DateTime"]["Data"]
    ts = int(ds.globalattributes["time_step"])
    nperDay = int(24/(float(ts)/60.0)+0.5)
    si = qcutils.GetDateIndex(dt,str(dt[0]),ts=ts,default=0,match="startnextday")
    ei = qcutils.GetDateIndex(dt,str(dt[-1]),ts=ts,default=len(dt)-1,match="endpreviousday")
    ldt = dt[si:ei+1]
    start_year = ldt[0].year
    end_year = ldt[-1].year
    year_list = range(start_year,end_year+1,1)
    annual_dict = {}
    annual_dict["DateTime"] = {"data":[datetime.datetime(yr,1,1) for yr in year_list],
                               "units":"Years","format":"dd/mm/yyyy"}
    annual_dict["nDays"] = {"data":numpy.array([float(-9999)]*len(year_list)),
                            "units":"","format":"0"}
    # create arrays in annual_dict
    series_list = series_dict["annual"].keys()
    series_list.sort()
    for item in series_list:
        annual_dict[item] = {"data":numpy.array([float(-9999)]*len(year_list))}
    for i,year in enumerate(year_list):
        if ts==30:
            start_date = str(year)+"-01-01 00:30"
        elif ts==60:
            start_date = str(year)+"-01-01 01:00"
        end_date = str(year+1)+"-01-01 00:00"
        si = qcutils.GetDateIndex(dt,start_date,ts=ts,default=0)
        ei = qcutils.GetDateIndex(dt,end_date,ts=ts,default=len(dt)-1)
        nDays = int((ei-si+1)/nperDay+0.5)
        annual_dict["nDays"]["data"][i] = nDays
        for item in series_list:
            if item not in ds.series.keys(): continue
            data_1d,flag,attr = qcutils.GetSeriesasMA(ds,item,si=si,ei=ei)
            if item in series_dict["lists"]["co2"]:
                data_1d = qcutils.convertunits(data_1d,attr["units"],"gC/m2",ts)
                annual_dict[item]["units"] = "gC/m2"
            else:
                annual_dict[item]["units"] = attr["units"]
            if series_dict["annual"][item]["operator"].lower()=="average":
                annual_dict[item]["data"][i] = numpy.ma.average(data_1d)
            elif series_dict["annual"][item]["operator"].lower()=="sum":
                annual_dict[item]["data"][i] = numpy.ma.sum(data_1d)
                annual_dict[item]["units"] = annual_dict[item]["units"]+"/year"
            else:
                print "unrecognised operator"
            annual_dict[item]["format"] = series_dict["annual"][item]["format"]
    return annual_dict

def L6_summary_cumulative(ds,series_dict):
    """
    Purpose:
     Calculate the cumulative sums of various quantities and write
     them to a worksheet in an Excel workbook.
    Usage:
     L6_summary_cumulative(xl_file,ds,series_dict)
     where xl_file is an Excel file object
           ds is an OzFluxQC data structure
           series_dict is a dictionary of various variable lists
    Author: PRI
    Date: June 2015
    """
    log.info(" Doing the cumulative summaries at L6")
    dt = ds.series["DateTime"]["Data"]
    ts = int(ds.globalattributes["time_step"])
    si = qcutils.GetDateIndex(dt,str(dt[0]),ts=ts,default=0,match="startnextday")
    ei = qcutils.GetDateIndex(dt,str(dt[-1]),ts=ts,default=len(dt)-1,match="endpreviousday")
    ldt = dt[si:ei+1]
    start_year = ldt[0].year
    end_year = ldt[-1].year
    year_list = range(start_year,end_year+1,1)
    series_list = series_dict["cumulative"].keys()
    cumulative_dict = {}
    for i,year in enumerate(year_list):
        cumulative_dict[str(year)] = {}
        if ts==30:
            start_date = str(year)+"-01-01 00:30"
        elif ts==60:
            start_date = str(year)+"-01-01 01:00"
        end_date = str(year+1)+"-01-01 00:00"
        si = qcutils.GetDateIndex(dt,start_date,ts=ts,default=0)
        ei = qcutils.GetDateIndex(dt,end_date,ts=ts,default=len(dt)-1)
        ldt = dt[si:ei+1]
        cumulative_dict[str(year)]["DateTime"] = {"data":ldt,"units":"Year",
                                                  "format":"dd/mm/yyyy HH:MM"}
        for item in series_list:
            cumulative_dict[str(year)][item] = {}
            data,flag,attr = qcutils.GetSeriesasMA(ds,item,si=si,ei=ei)
            if item in series_dict["lists"]["co2"]:
                data = qcutils.convertunits(data,attr["units"],"gC/m2",ts)
                cumulative_dict[str(year)][item]["units"] = "gC/m2"
            else:
                cumulative_dict[str(year)][item]["units"] = attr["units"]
            cumulative_dict[str(year)][item]["data"] = numpy.ma.cumsum(data)
            cumulative_dict[str(year)][item]["format"] = series_dict["cumulative"][item]["format"]
            cumulative_dict[str(year)][item]["units"] = cumulative_dict[str(year)][item]["units"]+"/year"
    return cumulative_dict

def ParseL6ControlFile(cf,ds):
    """ Parse the L6 control file. """
    # start with the repiration section
    if "Respiration" in cf.keys() and "ER" not in cf.keys(): cf["ER"] = cf["Respiration"]
    if "ER" in cf.keys():
        for ThisOne in cf["ER"].keys():
            if "ERUsingSOLO" in cf["ER"][ThisOne].keys():
                qcrpNN.rpSOLO_createdict(cf,ds,ThisOne)      # create the SOLO dictionary in ds
            if "ERUsingFFNET" in cf["ER"][ThisOne].keys():
                qcrpNN.rpFFNET_createdict(cf,ds,ThisOne)     # create the FFNET dictionary in ds
    if "NEE" in cf.keys():
        for ThisOne in cf["NEE"].keys():
            rpNEE_createdict(cf,ds,ThisOne)
    if "GPP" in cf.keys():
        for ThisOne in cf["GPP"].keys():
            rpGPP_createdict(cf,ds,ThisOne)

def PartitionNEE(cf,ds):
    """
    Purpose:
     Partition NEE into GPP and ER.
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
    Fsd_threshold = float(cf['Options']['Fsd_threshold'])
    # get the incoming shortwave radiation
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,"Fsd")
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd_syn")
        index = numpy.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Fsd)==True)[0]
        Fsd[index] = Fsd_syn[index]
    # calculate GPP from NEE and ER
    for label in ds.gpp.keys():
        if "NEE" not in ds.gpp[label] and "ER" not in ds.gpp[label]: continue
        NEE_label = ds.gpp[label]["NEE"]
        ER_label = ds.gpp[label]["ER"]
        output_label = ds.gpp[label]["output"]
        NEE,NEE_flag,NEE_attr = qcutils.GetSeriesasMA(ds,NEE_label)
        ER,ER_flag,ER_attr = qcutils.GetSeriesasMA(ds,ER_label)
        # calculate GPP
        # here we use the conventions from Chapin et al (2006)
        #  NEP = -1*NEE
        #  GPP = NEP + ER ==> GPP = -1*NEE + ER
        GPP = float(-1)*NEE + ER
        # put the day time data into the GPP series
        index = numpy.ma.where(Fsd>=Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = GPP[index]
        ds.series[output_label]["Flag"][index] = NEE_flag[index]
        # put the night time ER into the NEE series
        # This force nocturnal GPP to be 0!  Not sure this is the right thing to do.
        index = numpy.ma.where(Fsd<Fsd_threshold)[0]
        ds.series[output_label]["Data"][index] = numpy.float64(0)
        ds.series[output_label]["Flag"][index] = numpy.int32(1)
        # copy the attributes
        attr = ds.series[output_label]["Attr"]
        attr["units"] = NEE_attr["units"]
        attr["long_name"] = "Gross Primary Productivity calculated from "+NEE_label+" as -NEE+ER "
        attr["long_name"] = attr["long_name"]+" and "+ER_label+" (ER)"

def rpGPP_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about calculating GPP."""
    # create the ffnet directory in the data structure
    if "gpp" not in dir(ds): ds.gpp = {}
    # create the dictionary keys for this series
    ds.gpp[series] = {}
    # output series name
    ds.gpp[series]["output"] = series
    # CO2 flux
    if "NEE" in cf["GPP"][series].keys():
        ds.gpp[series]["NEE"] = cf["GPP"][series]["NEE"]
    # ecosystem respiration
    if "ER" in cf["GPP"][series].keys():
        ds.gpp[series]["ER"] = cf["GPP"][series]["ER"]
    # create an empty series in ds if the output series doesn't exist yet
    if ds.gpp[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.gpp[series]["output"])
        qcutils.CreateSeries(ds,ds.gpp[series]["output"],data,Flag=flag,Attr=attr)

def rpNEE_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about calculating NEE."""
    # create the ffnet directory in the data structure
    if "nee" not in dir(ds): ds.nee = {}
    # create the dictionary keys for this series
    ds.nee[series] = {}
    # output series name
    ds.nee[series]["output"] = series
    # CO2 flux
    if "Fc" in cf["NEE"][series].keys():
        ds.nee[series]["Fc"] = cf["NEE"][series]["Fc"]
    # ecosystem respiration
    if "ER" in cf["NEE"][series].keys():
        ds.nee[series]["ER"] = cf["NEE"][series]["ER"]
    # create an empty series in ds if the output series doesn't exist yet
    if ds.nee[series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.nee[series]["output"])
        qcutils.CreateSeries(ds,ds.nee[series]["output"],data,Flag=flag,Attr=attr)

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
