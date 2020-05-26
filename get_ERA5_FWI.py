# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:51:46 2020

@author: guillaume
"""

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'cems-fire-historical',
    {
        'format': 'netcdf',
        'variable': [
            'fire_weather_index', 
        ],
        'product_type': 'reanalysis',
        'year': [
            '2018', 
        ],
        'month': [
            '04', '05', '06',
            '07', '08',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'version': '3.1',
        'dataset': 'Consolidated dataset',
    },
    'download2.nc')