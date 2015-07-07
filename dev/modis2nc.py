import netCDF4

dap_url = "http://www.auscover.org.au/thredds/dodsC/"
dap_folder = "auscover/lpdaac-aggregates/c5/v2-nc4/aust/MOD13Q1.005/"
dap_name = "MOD13Q1.aggregated.aust.005.enhanced_vegetation_index.ncml"
dap_file = dap_url+dap_folder+dap_name

nc_file = netCDF4.Dataset(dap_file,"r")

lat_resolution = getattr(nc_file,"geospatial_lat_resolution")
lon_resolution = getattr(nc_file,"geospatial_lon_resolution")

# get the site information and the AWS stations to use
xlname = "../../BoM/Locations/AWS_Locations.xls"
wb = xlrd.open_workbook(xlname)
sheet = wb.sheet_by_name("OzFlux")
xl_row = 10
xl_col = 0
bom_sites_info = {}
for n in range(xl_row,sheet.nrows):
    xlrow = sheet.row_values(n)
    bom_sites_info[str(xlrow[0])] = {}
    bom_sites_info[xlrow[0]]["latitude"] = xlrow[1]
    bom_sites_info[xlrow[0]]["longitude"] = xlrow[2]
    bom_sites_info[xlrow[0]]["elevation"] = xlrow[3]
    for i in [4,10,16,22]:
        if xlrow[i]!="":
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))] = {}
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["site_name"] = xlrow[i]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["latitude"] = xlrow[i+2]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["longitude"] = xlrow[i+3]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["elevation"] = xlrow[i+4]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["distance"] = xlrow[i+5]

site_list = bom_sites_info.keys()
