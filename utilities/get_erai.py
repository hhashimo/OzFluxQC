#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import sys

erai_info = {}
erai_info["class"] = "ei"
erai_info["dataset"] = "interim"
erai_info["expver"] = "1"
erai_info["grid"] = "0.75/0.75"
erai_info["levtype"] = "sfc"
erai_info["param"] = "39.128/40.128/41.128/42.128/134.128/139.128/146.128/147.128/159.128/165.128/166.128/167.128/168.128/169.128/170.128/175.128/176.128/177.128/183.128/228.128/236.128"
erai_info["step"] = "3/6/9/12"
erai_info["stream"] = "oper"
erai_info["time"] = "00/12"
erai_info["type"] = "fc"
erai_info["format"] = "netcdf"

if len(sys.argv)==1:
    print "Command line syntax is:"
    print " python get_erai.py <country>"
    print "where <country> can be;"
    print " Australia"
    print " USA"
    sys.exit

if sys.argv[1].lower()=="australia":
    erai_info["area"] = "-10/110/-45/155"
    target_directory = "/mnt/OzFlux/ERAI/"
    start_year = 2014
    end_year = 2016
elif sys.argv[1].lower()=="usa":
    erai_info["area"] = "70/229.5/30/300"
    target_directory = "/mnt/AmeriFlux/ERAI/"
    start_year = 1991
    end_year = 2016
else:
    print "Unrecognised country option entered on command line"
    print "Valid country options are:"
    print " australia"
    print " usa"
    sys.exit()

server = ECMWFDataServer()
year_list = range(start_year,end_year,1)
for year in year_list:
    print " Processing year: ",str(year)
    erai_info["date"] = str(year)+"-01-01/to/"+str(year+1)+"-01-01"
    erai_info["target"] = target_directory+"ERAI_"+str(year)+".nc"
    server.retrieve(erai_info)
