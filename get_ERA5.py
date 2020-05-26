# Guillaume Dueymes: version 29-10-2019
# Code python permettant d extraire de ERA5
#
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
# Besoin d'avoir une clef API d'ECMWF
# 

from ecmwfapi import ECMWFDataServer
 
server = ECMWFDataServer()

yearStart = 2000
yearEnd = 2001
monthStart = 1
monthEnd = 12 

for year in list(range(yearStart, yearEnd + 1)):
    requestMonthList = []
    for month in list(range(monthStart, monthEnd + 1)):
            requestMonthList.append('%04d-%02d-01' % (year, month))
    requestMonths = "/".join(requestMonthList)
        # we submit a data request for the current year
    target = "era5_edna_ea_"+str(year)+"_sfc.nc"  
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "stream": "enda",
        "expver": "1",
        "type": "an",
        "levtype": "sfc",
        "date": requestMonths,
        "param": "167.128",                      #identification de la variable 
        "number": "0/1/2/3/4/5/6/7/8/9",
        "target": target,
        "grid": "0.25/0.25",
        "area": "75/-100/40/-45",  # North/West/South/East in Geographic lat/long degrees.
        "time": "00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00",
        "format": "netcdf",      
        })