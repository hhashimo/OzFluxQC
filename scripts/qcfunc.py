import meteorologicalfunctions as mf
import qcutils

def AhfromRH(ds,Ah_out,RH_in,Ta_in):
    """
    Purpose:
     Function to calculate absolute humidity given relative humidity and
     air temperature.  Absolute humidity is not calculated if any of the
     input series are missing or if the specified output series already
     exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     qcfunc.AhfromRH(ds,"Ah_HMP_2m","RH_HMP_2m","Ta_HMP_2m")
    Author: PRI
    Date: September 2015
    """
    for item in [RH_in,Ta_in]:
        if item not in ds.series.keys():
            msg = " AhfromRH: Requested series "+item+" not found, "+Ah_out+" not calculated"
            log.error(msg)
            return 0
    if Ah_out in ds.series.keys():
        msg = " AhfromRH: Output series "+Ah_out+" already exists, skipping ..."
        log.error(msg)
        return 0
    RH_data,RH_flag,RH_attr = qcutils.GetSeriesasMA(ds,RH_in)
    Ta_data,Ta_flag,Ta_attr = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah_data = mf.absolutehumidityfromRH(Ta_data,RH_data)
    Ah_attr = qcutils.MakeAttributeDictionary(long_name="Absolute humidity calculated from "+RH_in+" and "+Ta_in,
                                              height=RH_attr["height"],
                                              units="g/m3")
    qcutils.CreateSeries(ds,Ah_out,Ah_data,FList=[RH_in,Ta_in],Attr=Ah_attr)
    return 1

def AhfromMR(ds,Ah_out,MR_in,Ta_in,ps_in):
    """
    Purpose:
     Function to calculate absolute humidity given the water vapour mixing
     ratio, air temperature and pressure.  Absolute humidity is not calculated
     if any of the input series are missing or if the specified output series
     already exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     qcfunc.AhfromMR(ds,"Ah_IRGA_Av","H2O_IRGA_Av","Ta_HMP_2m","ps")
    Author: PRI
    Date: September 2015
    """
    for item in [MR_in,Ta_in,ps_in]:
        if item not in ds.series.keys():
            msg = " AhfromMR: Requested series "+item+" not found, "+Ah_out+" not calculated"
            log.error(msg)
            return 0
    if Ah_out in ds.series.keys():
        msg = " AhfromMR: Output series "+Ah_out+" already exists, skipping ..."
        log.error(msg)
        return 0
    MR_data,MR_flag,MR_attr = qcutils.GetSeriesasMA(ds,MR_in)
    Ta_data,Ta_flag,Ta_attr = qcutils.GetSeriesasMA(ds,Ta_in)
    ps_data,ps_flag,ps_attr = qcutils.GetSeriesasMA(ds,ps_in)
    Ah_data = mf.h2o_gpm3frommmolpmol(MR_data,Ta_data,ps_data)
    long_name = "Absolute humidity calculated from "+MR_in+", "+Ta_in+" and "+ps_in
    Ah_attr = qcutils.MakeAttributeDictionary(long_name=long_name,
                                              height=MR_attr["height"],
                                              units="g/m3")
    qcutils.CreateSeries(ds,Ah_out,Ah_data,FList=[MR_in,Ta_in,ps_in],Attr=Ah_attr)
    return 1

def test(arg1,arg2):
    print "got args:",arg1,arg2
    return "that worked"

    