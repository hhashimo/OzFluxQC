import ast
import copy
import constants as c
import datetime
import numpy
import time
import qcrp
import qcts
import qcutils
import logging

log = logging.getLogger('qc.ck')

def ApplyQCChecks(variable):
    """
    Purpose:
     Apply the QC checks speified in the control file object to a single variable
    Usage:
     qcck.ApplyQCChecks(variable)
     where variable is a variable dictionary as returned by qcutils.GetVariableAsDict()
    Author: PRI
    Date: September 2016
    """
    # do the range check
    ApplyRangeCheckToVariable(variable)
    # do the diurnal check
    #do_diurnalcheck_variable(cf,variable)
    # do exclude dates
    #do_excludedates_variable(cf,variable)
    # do exclude hours
    #do_excludehours_variable(cf,variable)
    return

def ApplyRangeCheckToVariable(variable):
    """
    Purpose:
    Usage:
    Author: PRI
    Date: September 2016
    """
    dt = variable["DateTime"]
    # Check to see if a lower limit has been specified
    if "rangecheck_lower" in variable["Attr"]:
        attr = variable["Attr"]["rangecheck_lower"]
        lower = numpy.array(eval(attr))
        valid_lower = str(numpy.min(lower))
        month = numpy.array([dt[i].month for i in range(0,len(dt))])
        lower_series = lower[month-1]
        index = numpy.ma.where(variable["Data"]<lower_series)[0]
        variable["Data"][index] = numpy.ma.masked
        variable["Flag"][index] = numpy.int32(2)
        valid_range = variable["Attr"]["valid_range"]
        old_lower = valid_range.split(",")[0]
        valid_range = valid_range.replace(old_lower,valid_lower)
        variable["Attr"]["valid_range"] = valid_range
    if "rangecheck_upper" in variable["Attr"]:
        attr = variable["Attr"]["rangecheck_upper"]
        upper = numpy.array(eval(attr))
        valid_upper = str(numpy.min(upper))
        month = numpy.array([dt[i].month for i in range(0,len(dt))])
        upper_series = upper[month-1]
        index = numpy.ma.where(variable["Data"]>upper_series)[0]
        variable["Data"][index] = numpy.ma.masked
        variable["Flag"][index] = numpy.int32(2)
        valid_range = variable["Attr"]["valid_range"]
        old_upper = valid_range.split(",")[1]
        valid_range = valid_range.replace(old_upper,valid_upper)
        variable["Attr"]["valid_range"] = valid_range
    return

def ApplyTurbulenceFilter(cf,ds,ustar_threshold=None):
    """
    Purpose:
    Usage:
    Author:
    Date:
    """
    opt = ApplyTurbulenceFilter_checks(cf,ds)
    if not opt["OK"]: return
    # local point to datetime series
    ldt = ds.series["DateTime"]["Data"]
    # time step
    ts = int(ds.globalattributes["time_step"])
    # dictionary of utar thresold values
    if ustar_threshold==None:
        ustar_dict = qcrp.get_ustar_thresholds(cf,ldt)
    else:
        ustar_dict = qcrp.get_ustar_thresholds_annual(ldt,ustar_threshold)
    # initialise a dictionary for the indicator series
    indicators = {}
    # get data for the indicator series
    ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,"ustar")
    Fsd,f,a = qcutils.GetSeriesasMA(ds,"Fsd")
    if "solar_altitude" not in ds.series.keys(): qcts.get_synthetic_fsd(ds)
    Fsd_syn,f,a = qcutils.GetSeriesasMA(ds,"Fsd_syn")
    sa,f,a = qcutils.GetSeriesasMA(ds,"solar_altitude")
    # get the day/night indicator series
    # indicators["day"] = 1 ==> day time, indicators["day"] = 0 ==> night time
    indicators["day"] = qcrp.get_day_indicator(cf,Fsd,Fsd_syn,sa)
    ind_day = indicators["day"]["values"]
    # get the turbulence indicator series
    if opt["turbulence_filter"].lower()=="ustar":
        # indicators["turbulence"] = 1 ==> turbulent, indicators["turbulence"] = 0 ==> not turbulent
        indicators["turbulence"] = qcrp.get_turbulence_indicator_ustar(ldt,ustar,ustar_dict,ts)
    elif opt["turbulence_filter"].lower()=="ustar_evg":
        # ustar >= threshold ==> ind_ustar = 1, ustar < threshold == ind_ustar = 0
        indicators["ustar"] = qcrp.get_turbulence_indicator_ustar(ldt,ustar,ustar_dict,ts)
        ind_ustar = indicators["ustar"]["values"]
        # ustar >= threshold during day AND ustar has been >= threshold since sunset ==> indicators["turbulence"] = 1
        # indicators["turbulence"] = 0 during night once ustar has dropped below threshold even if it
        # increases above the threshold later in the night
        indicators["turbulence"] = qcrp.get_turbulence_indicator_ustar_evg(ldt,ind_day,ind_ustar,ustar,ustar_dict,ts)
    elif opt["turbulence_filter"].lower()=="l":
        #indicators["turbulence] = get_turbulence_indicator_l(ldt,L,z,d,zmdonL_threshold)
        indicators["turbulence"] = numpy.ones(len(ldt))
        msg = " Use of L as turbulence indicator not implemented, no filter applied"
        log.warning(msg)
    else:
        msg = " Unrecognised turbulence filter option ("
        msg = msg+opt["turbulence_filter"]+"), no filter applied"
        log.error(msg)
        return
    # initialise the final indicator series as the turbulence indicator
    # subsequent filters will modify the final indicator series
    # we must use copy.deepcopy() otherwise the "values" array will only
    # be copied by reference not value.  Damn Python's default of copy by reference!
    indicators["final"] = copy.deepcopy(indicators["turbulence"])
    # check to see if the user wants to accept all day time observations
    # regardless of ustar value
    if opt["accept_day_times"].lower()=="yes":
        # if yes, then we force the final indicator to be 1
        # if ustar is below the threshold during the day.
        idx = numpy.where(indicators["day"]["values"]==1)[0]
        indicators["final"]["values"][idx] = numpy.int(1)
        indicators["final"]["attr"].update(indicators["day"]["attr"])
    # get the evening indicator series
    indicators["evening"] = qcrp.get_evening_indicator(cf,Fsd,Fsd_syn,sa,ts)
    indicators["dayevening"] = {"values":indicators["day"]["values"]+indicators["evening"]["values"]}
    indicators["dayevening"]["attr"] = indicators["day"]["attr"].copy()
    indicators["dayevening"]["attr"].update(indicators["evening"]["attr"])
    if opt["use_evening_filter"].lower()=="yes":
        idx = numpy.where(indicators["dayevening"]["values"]==0)[0]
        indicators["final"]["values"][idx] = numpy.int(0)
        indicators["final"]["attr"].update(indicators["dayevening"]["attr"])
    # save the indicator series
    ind_flag = numpy.zeros(len(ldt))
    long_name = "Turbulence indicator, 1 for turbulent, 0 for non-turbulent"
    ind_attr = qcutils.MakeAttributeDictionary(long_name=long_name,units="None")
    qcutils.CreateSeries(ds,"turbulence_indicator",indicators["turbulence"]["values"],Flag=ind_flag,Attr=ind_attr)
    long_name = "Day indicator, 1 for day time, 0 for night time"
    ind_attr = qcutils.MakeAttributeDictionary(long_name=long_name,units="None")
    qcutils.CreateSeries(ds,"day_indicator",indicators["day"]["values"],Flag=ind_flag,Attr=ind_attr)
    long_name = "Evening indicator, 1 for evening, 0 for not evening"
    ind_attr = qcutils.MakeAttributeDictionary(long_name=long_name,units="None")
    qcutils.CreateSeries(ds,"evening_indicator",indicators["evening"]["values"],Flag=ind_flag,Attr=ind_attr)
    long_name = "Day/evening indicator, 1 for day/evening, 0 for not day/evening"
    ind_attr = qcutils.MakeAttributeDictionary(long_name=long_name,units="None")
    qcutils.CreateSeries(ds,"dayevening_indicator",indicators["dayevening"]["values"],Flag=ind_flag,Attr=ind_attr)
    long_name = "Final indicator, 1 for use data, 0 for don't use data"
    ind_attr = qcutils.MakeAttributeDictionary(long_name=long_name,units="None")
    qcutils.CreateSeries(ds,"final_indicator",indicators["final"]["values"],Flag=ind_flag,Attr=ind_attr)
    # loop over the series to be filtered
    for series in opt["filter_list"]:
        msg = " Applying "+opt["turbulence_filter"]+" filter to "+series
        log.info(msg)
        # get the data
        data,flag,attr = qcutils.GetSeriesasMA(ds,series)
        # continue to next series if this series has been filtered before
        if "turbulence_filter" in attr:
            msg = " Series "+series+" has already been filtered, skipping ..."
            log.warning(msg)
            continue
        # save the non-filtered data
        qcutils.CreateSeries(ds,series+"_nofilter",data,Flag=flag,Attr=attr)
        # now apply the filter
        data_filtered = numpy.ma.masked_where(indicators["final"]["values"]==0,data,copy=True)
        flag_filtered = numpy.copy(flag)
        idx = numpy.where(indicators["final"]["values"]==0)[0]
        flag_filtered[idx] = numpy.int32(61)
        # update the series attributes
        for item in indicators["final"]["attr"].keys():
            attr[item] = indicators["final"]["attr"][item]
        # and write the filtered data to the data structure
        qcutils.CreateSeries(ds,series,data_filtered,Flag=flag_filtered,Attr=attr)
        # and write a copy of the filtered datas to the data structure so it
        # will still exist once the gap filling has been done
        qcutils.CreateSeries(ds,series+"_filtered",data_filtered,Flag=flag_filtered,Attr=attr)
    return

def ApplyTurbulenceFilter_checks(cf,ds):
    """
    Purpose:
    Usage:
    Author:
    Date:
    """
    opt = {"OK":True,"turbulence_filter":"ustar","filter_list":['Fc']}
    # return if there is no Options section in control file
    if "Options" not in cf:
        msg = " ApplyTurbulenceFilter: Options section not found in control file"
        log.warning(msg)
        opt["OK"] = False
        return opt
    # get the value of the TurbulenceFilter key in the Options section
    opt["turbulence_filter"] = qcutils.get_keyvaluefromcf(cf,["Options"],"TurbulenceFilter",default="None")
    # return if turbulence filter disabled
    if opt["turbulence_filter"].lower()=="none":
        msg = " Turbulence filter disabled in control file at "+ds.globalattributes["nc_level"]
        log.info(msg)
        opt["OK"] = False
        return opt
    # check to see if filter type can be handled
    if opt["turbulence_filter"].lower() not in ["ustar","ustar_evg","l"]:
        msg = " Unrecognised turbulence filter option ("
        msg = msg+opt["turbulence_filter"]+"), no filter applied"
        log.error(msg)
        opt["OK"] = False
        return opt
    # get the list of series to be filtered
    if "FilterList" in cf["Options"]:
        opt["filter_list"] = ast.literal_eval(cf["Options"]["FilterList"])
    # check to see if the series are in the data structure
    for item in opt["filter_list"]:
        if item not in ds.series.keys():
            msg = " Series "+item+" given in FilterList not found in data stucture"
            log.warning(msg)
            opt["filter_list"].remove(item)
    # return if the filter list is empty
    if len(opt["filter_list"])==0:
        msg = " FilterList in control file is empty, skipping turbulence filter"
        log.warning(msg)
        opt["OK"] = False
        return opt
    # get the value of the DayNightFilter key in the Options section
    opt["daynight_filter"] = qcutils.get_keyvaluefromcf(cf,["Options"],"DayNightFilter",default="None")
    # check to see if filter type can be handled
    if opt["daynight_filter"].lower() not in ["fsd","sa","none"]:
        msg = " Unrecognised day/night filter option ("
        msg = msg+opt["daynight_filter"]+"), no filter applied"
        log.error(msg)
        opt["OK"] = False
        return opt
    # check to see if all day time values are to be accepted
    opt["accept_day_times"] = qcutils.get_keyvaluefromcf(cf,["Options"],"AcceptDayTimes",default="Yes")
    opt["use_evening_filter"] = qcutils.get_keyvaluefromcf(cf,["Options"],"UseEveningFilter",default="Yes")
    
    return opt

def cliptorange(data, lower, upper):
    data = rangecheckserieslower(data,lower)
    data = rangecheckseriesupper(data,upper)
    return data

def CoordinateAh7500AndFcGaps(cf,ds,Fcvar='Fc'):
    '''Cleans up Ah_7500_Av based upon Fc gaps to for QA check on Ah_7500_Av v Ah_HMP.'''
    if not qcutils.cfoptionskeylogical(cf,Key='CoordinateAh7500&FcGaps'): return
    log.info(' Doing the Ah_7500 check')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='AhcheckFc'):
        Fclist = ast.literal_eval(cf['FunctionArgs']['AhcheckFc'])
        Fcvar = Fclist[0]
    
    # index1  Index of bad Ah_7500_Av observations
    index1 = numpy.where((ds.series['Ah_7500_Av']['Flag']!=0) & (ds.series['Ah_7500_Av']['Flag']!=10))
    
    # index2  Index of bad Fc observations
    index2 = numpy.where((ds.series[Fcvar]['Flag']!=0) & (ds.series[Fcvar]['Flag']!=10))
    
    ds.series['Ah_7500_Av']['Data'][index2] = numpy.float64(c.missing_value)
    ds.series['Ah_7500_Av']['Flag'][index2] = ds.series[Fcvar]['Flag'][index2]
    ds.series['Ah_7500_Av']['Flag'][index1] = ds.series['Ah_7500_Av']['Flag'][index1]
    if 'CoordinateAh7500AndFcGaps' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CoordinateAh7500AndFcGaps'

def CoordinateFluxGaps(cf,ds,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh'):
    if not qcutils.cfoptionskeylogical(cf,Key='CoordinateFluxGaps'): return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='gapsvars'):
        vars = ast.literal_eval(cf['FunctionArgs']['gapsvars'])
        Fc_in = vars[0]
        Fe_in = vars[1]
        Fh_in = vars[2]
    Fc,f,a = qcutils.GetSeriesasMA(ds,Fc_in)
    Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    # April 2015 PRI - changed numpy.ma.where to numpy.where
    index = numpy.where((numpy.ma.getmaskarray(Fc)==True)|
                        (numpy.ma.getmaskarray(Fe)==True)|
                        (numpy.ma.getmaskarray(Fh)==True))[0]
    #index = numpy.ma.where((numpy.ma.getmaskarray(Fc)==True)|
                           #(numpy.ma.getmaskarray(Fe)==True)|
                           #(numpy.ma.getmaskarray(Fh)==True))[0]
    # the following for ... in loop is not necessary
    for i in range(len(index)):
        j = index[i]
        if Fc.mask[j]==False:
            Fc.mask[j]=True
            Fc[j] = numpy.float64(c.missing_value)
            ds.series[Fc_in]['Flag'][j] = numpy.int32(19)
        if Fe.mask[j]==False:
            Fe.mask[j]=True
            Fe[j] = numpy.float64(c.missing_value)
            ds.series[Fe_in]['Flag'][j] = numpy.int32(19)           
        if Fh.mask[j]==False:
            Fh.mask[j]=True
            Fh[j] = numpy.float64(c.missing_value)
            ds.series[Fh_in]['Flag'][j] = numpy.int32(19)
    ds.series[Fc_in]['Data']=numpy.ma.filled(Fc,float(c.missing_value))
    ds.series[Fe_in]['Data']=numpy.ma.filled(Fe,float(c.missing_value))
    ds.series[Fh_in]['Data']=numpy.ma.filled(Fh,float(c.missing_value))
    log.info(' Finished gap co-ordination')

def CreateNewSeries(cf,ds):
    '''Create a new series using the MergeSeries or AverageSeries instructions.'''
    log.info(' Checking for new series to create')
    for ThisOne in cf['Variables'].keys():
        if 'MergeSeries' in cf['Variables'][ThisOne].keys():
            qcts.MergeSeries(cf,ds,ThisOne,[0,10,20,30,40,50])
        if 'AverageSeries' in cf['Variables'][ThisOne].keys():
            qcts.AverageSeriesByElements(cf,ds,ThisOne)

def do_CSATcheck(cf,ds):
    '''Rejects data values for series specified in CSATList for times when the Diag_CSAT
       flag is non-zero.  If the Diag_CSAT flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       CSATList.'''
    if "Diag_CSAT" not in ds.series.keys():
        msg = " Diag_CSAT not found in data, skipping CSAT checks ..."
        log.warning(msg)
        return
    log.info(' Doing the CSAT check')
    CSAT_all = ['Ux','Uy','Uz',
                'Ws_CSAT','Wd_CSAT','Wd_CSAT_Compass','Tv_CSAT',
                'UzT','UxT','UyT','UzA','UxA','UyA','UzC','UxC','UyC',
                'UxUz','UyUz','UxUy','UxUx','UyUy','UzUz']
    CSAT_list = []
    for item in CSAT_all:
        if item in ds.series.keys(): CSAT_list.append(item)
    index = numpy.where(ds.series['Diag_CSAT']['Flag']!=0)
    log.info('  CSATCheck: Diag_CSAT ' + str(numpy.size(index)))
    for ThisOne in CSAT_list:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(3)
        else:
            log.error('  qcck.do_CSATcheck: series '+str(ThisOne)+' in CSATList not found in ds.series')
    if 'CSATCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CSATCheck'

def do_dependencycheck(cf,ds,section='',series='',code=23,mode="quiet"):
    if len(section)==0 and len(series)==0: return
    if len(section)==0: section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if "DependencyCheck" not in cf[section][series].keys(): return
    if "Source" not in cf[section][series]["DependencyCheck"]:
        msg = " DependencyCheck: keyword Source not found for series "+series+", skipping ..."
        log.error(msg)
        return
    if mode=="verbose":
        msg = " Doing DependencyCheck for "+series
        log.info(msg)
    # get the precursor source list from the control file
    source_list = ast.literal_eval(cf[section][series]["DependencyCheck"]["Source"])
    # get the data
    dependent_data,dependent_flag,dependent_attr = qcutils.GetSeriesasMA(ds,series)
    # loop over the precursor source list
    for item in source_list:
        # check the precursor is in the data structure
        if item not in ds.series.keys():
            msg = " DependencyCheck: "+series+" precursor series "+item+" not found, skipping ..."
            continue
        # get the precursor data
        precursor_data,precursor_flag,precursor_attr = qcutils.GetSeriesasMA(ds,item)
        # mask the dependent data where the precurso is masked
        dependent_data = numpy.ma.masked_where(numpy.ma.getmaskarray(precursor_data)==True,dependent_data)
        # get an index of masked precursor data
        index = numpy.ma.where(numpy.ma.getmaskarray(precursor_data)==True)[0]
        # set the dependent QC flag
        dependent_flag[index] = numpy.int32(code)
    # put the data back into the data structure
    dependent_attr["DependencyCheck_source"] = str(source_list)
    qcutils.CreateSeries(ds,series,dependent_data,Flag=dependent_flag,Attr=dependent_attr)
    if 'do_dependencychecks' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_dependencychecks'

def do_diurnalcheck(cf,ds,section='',series='',code=5):
    if 'DiurnalCheck' not in cf[section][series].keys(): return
    if 'NumSd' not in cf[section][series]['DiurnalCheck'].keys(): return
    dt = float(ds.globalattributes['time_step'])
    n = int((60./dt) + 0.5)             #Number of timesteps per hour
    nInts = int((1440.0/dt)+0.5)        #Number of timesteps per day
    Av = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Sd = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    NSd = numpy.array(eval(cf[section][series]['DiurnalCheck']['NumSd']),dtype=float)
    for m in range(1,13):
        mindex = numpy.where(ds.series['Month']['Data']==m)[0]
        if len(mindex)!=0:
            lHdh = ds.series['Hdh']['Data'][mindex]
            l2ds = ds.series[series]['Data'][mindex]
            for i in range(nInts):
                li = numpy.where(abs(lHdh-(float(i)/float(n))<c.eps)&(l2ds!=float(c.missing_value)))
                if numpy.size(li)!=0:
                    Av[i] = numpy.mean(l2ds[li])
                    Sd[i] = numpy.std(l2ds[li])
                else:
                    Av[i] = float(c.missing_value)
                    Sd[i] = float(c.missing_value)
            Lwr = Av - NSd[m-1]*Sd
            Upr = Av + NSd[m-1]*Sd
            hindex = numpy.array(n*lHdh,int)
            index = numpy.where(((l2ds!=float(c.missing_value))&(l2ds<Lwr[hindex]))|
                                ((l2ds!=float(c.missing_value))&(l2ds>Upr[hindex])))[0] + mindex[0]
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['diurnalcheck_numsd'] = cf[section][series]['DiurnalCheck']['NumSd']
    if 'DiurnalCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',DiurnalCheck'

def do_EC155check(cf,ds):
    """
    Purpose:
    Usage:
    Author: PRI
    Date: September 2015
    """
    # check to see if we have a Diag_IRGA series to work with
    if "Diag_IRGA" not in ds.series.keys():
        msg = " Diag_IRGA not found in data, skipping IRGA checks ..."
        log.warning(msg)
        return
    # seems OK to continue
    log.info(' Doing the EC155 check')
    # list of series that depend on IRGA data quality
    EC155_list = ['H2O_IRGA_Av','CO2_IRGA_Av','H2O_IRGA_Sd','CO2_IRGA_Sd','H2O_IRGA_Vr','CO2_IRGA_Vr',
                 'UzA','UxA','UyA','UzH','UxH','UyH','UzC','UxC','UyC']
    index = numpy.where(ds.series['Diag_IRGA']['Flag']!=0)
    log.info('  EC155Check: Diag_IRGA rejects ' + str(numpy.size(index)))
    EC155_dependents = []
    for item in ['Signal_H2O','Signal_CO2','H2O_IRGA_Sd','CO2_IRGA_Sd']:
        if item in ds.series.keys(): EC155_dependents.append(item)
    for item in EC155_dependents:
        index = numpy.where(ds.series[item]['Flag']!=0)
        log.info('  EC155Check: '+item+' rejected '+str(numpy.size(index))+' points')
        ds.series['Diag_IRGA']['Flag'] = ds.series['Diag_IRGA']['Flag'] + ds.series[item]['Flag']
    index = numpy.where((ds.series['Diag_IRGA']['Flag']!=0))
    log.info('  EC155Check: Total rejected ' + str(numpy.size(index)))
    for ThisOne in EC155_list:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.warning(' do_EC155check: series '+str(ThisOne)+' in EC155 list not found in data structure')
    if 'EC155Check' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',EC155Check'

def do_excludedates(cf,ds,section='',series='',code=6):
    if 'ExcludeDates' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeDates'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeDateList = ast.literal_eval(cf[section][series]['ExcludeDates'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeDateList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeDateList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        ds.series[series]['Data'][si:ei] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][si:ei] = numpy.int32(code)
        ds.series[series]['Attr']['ExcludeDates_'+str(i)] = cf[section][series]['ExcludeDates'][str(i)]
    if 'ExcludeDates' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeDates'

def do_excludehours(cf,ds,section='',series='',code=7):
    if 'ExcludeHours' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeHours'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeHourList = ast.literal_eval(cf[section][series]['ExcludeHours'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeHourList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeHourList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        for j in range(len(ExcludeHourList[2])):
            ExHr = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').hour
            ExMn = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').minute
            index = numpy.where((ds.series['Hour']['Data'][si:ei]==ExHr)&
                                (ds.series['Minute']['Data'][si:ei]==ExMn))[0] + si
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['ExcludeHours_'+str(i)] = cf[section][series]['ExcludeHours'][str(i)]
    if 'ExcludeHours' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeHours'

def do_IRGAcheck(cf,ds):
    """
    Purpose:
     Decide which IRGA check routine to use depending on the setting
     of the "irga_type" key in the [Options] section of the control
     file.  The default is Li7500.
    Usage:
    Author: PRI
    Date: September 2015
    """
    irga_list = ["li7500","li7500a","li7500rs","ec150","ec155","irgason"]
    # get the IRGA type from the control file
    irga_type = qcutils.get_keyvaluefromcf(cf,["Options"],"irga_type", default="li7500")
    # remove any hyphens or spaces
    for item in ["-"," "]:
        if item in irga_type: irga_type = irga_type.replace(item,"")
    # check the IRGA type against the list of suppprted devices
    if irga_type.lower() not in irga_list:
        msg = " Unrecognised IRGA type "+irga_type+" given in control file, IRGA checks skipped ..."
        log.error(msg)
        return
    # do the IRGA checks
    if irga_type.lower()=="li7500":
        ds.globalattributes["irga_type"] = irga_type
        do_li7500check(cf,ds)
    elif irga_type.lower() in ["li7500a","irgason"]:
        ds.globalattributes["irga_type"] = irga_type
        do_li7500acheck(cf,ds)
    elif irga_type.lower() in ["ec155","ec150","irgason"]:
        ds.globalattributes["irga_type"] = irga_type
        do_EC155check(cf,ds)
    else:
        msg = " Unsupported IRGA type "+irga_type+", contact the devloper ..."
        log.error(msg)
        return

def do_li7500check(cf,ds):
    '''Rejects data values for series specified in LI75List for times when the Diag_7500
       flag is non-zero.  If the Diag_7500 flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       LI75Lisat.  Additional checks are done for AGC_7500 (the LI-7500 AGC value),
       Ah_7500_Sd (standard deviation of absolute humidity) and Cc_7500_Sd (standard
       deviation of CO2 concentration).'''
    if "Diag_7500" not in ds.series.keys():
        msg = " Diag_7500 not found in data, skipping 7500 checks ..."
        log.warning(msg)
        return
    log.info(' Doing the 7500 check')
    LI75List = ['Ah_7500_Av','Cc_7500_Av','Ah_7500_Sd','Cc_7500_Sd',
                'UzA','UxA','UyA','UzC','UxC','UyC']
    index = numpy.where(ds.series['Diag_7500']['Flag']!=0)
    log.info('  7500Check: Diag_7500 ' + str(numpy.size(index)))
    LI75_dependents = []
    for item in ['AGC_7500','Ah_7500_Sd','Cc_7500_Sd','AhAh','CcCc']:
        if item in ds.series.keys(): LI75_dependents.append(item)
    if "Ah_7500_Sd" and "AhAh" in LI75_dependents: LI75_dependents.remove("AhAh")
    if "Cc_7500_Sd" and "CcCc" in LI75_dependents: LI75_dependents.remove("CcCc")
    for item in LI75_dependents:
        if item in ds.series.keys():
            index = numpy.where(ds.series[item]['Flag']!=0)
            log.info('  7500Check: '+item+' rejected '+str(numpy.size(index))+' points')
            ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[item]['Flag']
    index = numpy.where((ds.series['Diag_7500']['Flag']!=0))
    log.info('  7500Check: Total ' + str(numpy.size(index)))
    for ThisOne in LI75List:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.error('  qcck.do_7500check: series '+str(ThisOne)+' in LI75List not found in ds.series')
    if '7500Check' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',7500Check'

def do_li7500acheck(cf,ds):
    #msg = " Li-7500A check not implemented yet, contact the developer ..."
    #log.warning(msg)
    #return
    '''Rejects data values for series specified in LI75List for times when the Diag_7500
       flag is non-zero.  If the Diag_IRGA flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       LI75Lisat.  Additional checks are done for AGC_7500 (the LI-7500 AGC value),
       Ah_7500_Sd (standard deviation of absolute humidity) and Cc_7500_Sd (standard
       deviation of CO2 concentration).'''
    if "Diag_IRGA" not in ds.series.keys():
        msg = " Diag_IRGA not found in data, skipping IRGA checks ..."
        log.warning(msg)
        return
    log.info(' Doing the 7500A check')
    LI75List = ['H2O_IRGA_Av','CO2_IRGA_Av','H2O_IRGA_Sd','CO2_IRGA_Sd','H2O_IRGA_Vr','CO2_IRGA_Vr',
                'UzA','UxA','UyA','UzC','UxC','UyC']
    index = numpy.where(ds.series['Diag_IRGA']['Flag']!=0)
    log.info('  7500ACheck: Diag_IRGA ' + str(numpy.size(index)))
    LI75_dependents = []
    for item in ['Signal_H2O','Signal_CO2','H2O_IRGA_Sd','CO2_IRGA_Sd','H2O_IRGA_Vr','CO2_IRGA_Vr']:
        if item in ds.series.keys(): LI75_dependents.append(item)
    if "H2O_IRGA_Sd" and "H2O_IRGA_Vr" in LI75_dependents: LI75_dependents.remove("H2O_IRGA_Vr")
    if "CO2_IRGA_Sd" and "CO2_IRGA_Vr" in LI75_dependents: LI75_dependents.remove("CO2_IRGA_Vr")
    for item in LI75_dependents:
        if item in ds.series.keys():
            index = numpy.where(ds.series[item]['Flag']!=0)
            log.info('  7500ACheck: '+item+' rejected '+str(numpy.size(index))+' points')
            ds.series['Diag_IRGA']['Flag'] = ds.series['Diag_IRGA']['Flag'] + ds.series[item]['Flag']
    index = numpy.where((ds.series['Diag_IRGA']['Flag']!=0))
    log.info('  7500ACheck: Total ' + str(numpy.size(index)))
    for ThisOne in LI75List:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            #log.warning('  qcck.do_7500acheck: series '+str(ThisOne)+' in LI75List not found in ds.series')
            pass
    if '7500ACheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',7500ACheck'

def do_linear(cf,ds):
    level = ds.globalattributes['nc_level']
    for ThisOne in cf['Variables'].keys():
        if qcutils.haskey(cf,ThisOne,'Linear'):
            qcts.ApplyLinear(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'Drift'):
            qcts.ApplyLinearDrift(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'LocalDrift'):
            qcts.ApplyLinearDriftLocal(cf,ds,ThisOne)
    if 'do_linear' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_linear'

def do_rangecheck(cf,ds,section='',series='',code=2):
    '''Applies a range check to data series listed in the control file.  Data values that
       are less than the lower limit or greater than the upper limit are replaced with
       c.missing_value and the corresponding QC flag element is set to 2.'''
    if 'RangeCheck' not in cf[section][series].keys(): return
    if 'Lower' in cf[section][series]['RangeCheck'].keys():
        lwr = numpy.array(eval(cf[section][series]['RangeCheck']['Lower']))
        valid_lower = numpy.min(lwr)
        lwr = lwr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']<lwr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_lower'] = cf[section][series]['RangeCheck']['Lower']
    if 'Upper' in cf[section][series]['RangeCheck'].keys():
        upr = numpy.array(eval(cf[section][series]['RangeCheck']['Upper']))
        valid_upper = numpy.min(upr)
        upr = upr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']>upr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_upper'] = cf[section][series]['RangeCheck']['Upper']
        ds.series[series]['Attr']['valid_range'] = str(valid_lower)+','+str(valid_upper)
    if 'RangeCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',RangeCheck'

def do_qcchecks(cf,ds,mode="verbose"):
    if "nc_level" in ds.globalattributes:
        level = str(ds.globalattributes["nc_level"])
        if mode!="quiet": log.info(" Doing the QC checks at level "+str(level))
    else:
        if mode!="quiet": log.info(" Doing the QC checks")
    # get the series list from the control file
    series_list = []
    for item in ["Variables","Drivers","Fluxes"]:
        if item in cf:
            section = item
            series_list = cf[item].keys()
    if len(series_list)==0:
        msg = " do_qcchecks: Variables, Drivers or Fluxes section not found in control file, skipping QC checks ..."
        log.warning(msg)
        return
    # loop over the series specified in the control file
    # first time for general QC checks
    for series in series_list:
        # check the series is in the data structure
        if series not in ds.series.keys():
            if mode!="quiet":
                msg = " do_qcchecks: series "+series+" not found in data structure, skipping ..."
                log.warning(msg)
            continue
        # if so, do the QC checks
        do_qcchecks_oneseries(cf,ds,section=section,series=series)
    # loop over the series in the control file
    # second time for dependencies
    for series in series_list:
        # check the series is in the data structure
        if series not in ds.series.keys():
            if mode!="quiet":
                msg = " do_qcchecks: series "+series+" not found in data structure, skipping ..."
                log.warning(msg)
            continue
        # if so, do dependency check
        do_dependencycheck(cf,ds,section=section,series=series,code=23,mode="quiet")

def do_qcchecks_oneseries(cf,ds,section='',series=''):
    if len(section)==0:
        section = qcutils.get_cfsection(cf,series=series,mode='quiet')
        if len(section)==0: return
    # do the range check
    do_rangecheck(cf,ds,section=section,series=series,code=2)
    # do the diurnal check
    do_diurnalcheck(cf,ds,section=section,series=series,code=5)
    # do exclude dates
    do_excludedates(cf,ds,section=section,series=series,code=6)
    # do exclude hours
    do_excludehours(cf,ds,section=section,series=series,code=7)
    # do wind direction corrections
    do_winddirectioncorrection(cf,ds,section=section,series=series)
    if 'do_qcchecks' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_qcchecks'

def do_winddirectioncorrection(cf,ds,section='',series=''):
    if 'CorrectWindDirection' not in cf[section][series].keys(): return
    qcts.CorrectWindDirection(cf,ds,series)

def rangecheckserieslower(data,lower):
    if lower is None:
        log.info(' rangecheckserieslower: no lower bound set')
        return data
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data<lower,data)
    else:
        index = numpy.where((abs(data-numpy.float64(c.missing_value))>c.eps)&(data<lower))[0]
        data[index] = numpy.float64(c.missing_value)
    return data

def rangecheckseriesupper(data,upper):
    if upper is None:
        log.info(' rangecheckserieslower: no upper bound set')
        return data
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data>upper,data)
    else:
        index = numpy.where((abs(data-numpy.float64(c.missing_value))>c.eps)&(data>upper))[0]
        data[index] = numpy.float64(c.missing_value)
    return data

def UpdateVariableAttributes_QC(cf, variable):
    """
    Purpose:
    Usage:
    Side effects:
    Author: PRI
    Date: November 2016
    """
    label = variable["Label"]
    section = qcutils.get_cfsection(cf,series=label,mode='quiet')
    if label not in cf[section]:
        return
    if "RangeCheck" not in cf[section][label]:
        return
    if "Lower" in cf[section][label]["RangeCheck"]:
        variable["Attr"]["rangecheck_lower"] = cf[section][label]["RangeCheck"]["Lower"]
    if "Upper" in cf[section][label]["RangeCheck"]:
        variable["Attr"]["rangecheck_upper"] = cf[section][label]["RangeCheck"]["Upper"]
    return