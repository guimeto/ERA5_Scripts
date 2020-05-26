# Guillaume Dueymes: version 29-10-2019
# Code python permettant d extraire de ERA5
#
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
# Besoin d'avoir une clef API d'ECMWF
# 
import cdsapi
c = cdsapi.Client()

yearStart = 2007
yearEnd = 2008
monthStart = 1
monthEnd = 12

for year in list(range(yearStart, yearEnd + 1)):
    
    requestMonthList = []
    for month in list(range(monthStart, monthEnd + 1)):
        # we submit a data request for the current year
        target = "era5_edna_ea_"+str(year)+f"{month:02d}"+"_sfc.nc"  
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type':'reanalysis',
                'variable': [
                '10m_u_component_of_wind', '10m_v_component_of_wind',
            ],
                'year':year,
                'month':month,
                'day':[
                    '01','02','03',
                    '04','05','06',
                    '07','08','09',
                    '10','11','12',
                    '13','14','15',
                    '16','17','18',
                    '19','20','21',
                    '22','23','24',
                    '25','26','27',
                    '28','29','30',
                    '31'
                ],
                'time':[
                    '00:00','01:00','02:00',
                    '03:00','04:00','05:00',
                    '06:00','07:00','08:00',
                    '09:00','10:00','11:00',
                    '12:00','13:00','14:00',
                    '15:00','16:00','17:00',
                    '18:00','19:00','20:00',
                    '21:00','22:00','23:00'
                ],
                'format':'netcdf'
            },target)