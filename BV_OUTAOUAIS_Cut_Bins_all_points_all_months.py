# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:06:25 2019

@author: guillaume
"""
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import warnings; warnings.filterwarnings(action='once')
import seaborn as sns
import warnings; warnings.filterwarnings(action='ignore')
import numpy as np

titre='Outaouais Bassin Versant: Daily mean precipitation'


legends = ['2019', '2020']

# Lecture du bassin versant de la rivière des Outaouais interpolé sur le doaime CORDEX-NAM44 
MASK = xr.open_dataset('Outaouais_ERA5_Grid.nc')
lat1d=MASK.variables['latitude'][:]
lon1d=MASK.variables['longitude'][:]
lon2d, lat2d = np.meshgrid(lon1d, lat1d)

mask = MASK['tp'][ :, :].values


# POUR 2020
year = 2020
file = 'J:/REANALYSES/ERA5/PR_1h_Outaouais/era5_pr_1h_'
multi_file = [f'{file}{year-1}{month}_sfc.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}_sfc.nc' for month in ['01','02','03']])
pr_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
pr_all = pr_all.assign_coords(longitude=(((pr_all.longitude + 180) % 360) - 180)).sortby('longitude') 
pr_all = pr_all *1000
pr_all =pr_all.resample(time = '1D').sum() 

file = 'J:/REANALYSES/ERA5/T2m_1h_Outaouais/era5_t2m_1h_'
multi_file = [f'{file}{year-1}{month}_sfc.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}_sfc.nc' for month in ['01','02','03']])
tt_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
tt_all = tt_all.assign_coords(longitude=(((tt_all.longitude + 180) % 360) - 180)).sortby('longitude') 
tt_all = tt_all - 273.15
tt_all = tt_all.resample(time = '1D').mean() 

tt_all.time[-1].values
x = pd.to_datetime(tt_all.time[-1].values)
last_day = str(x.date().day)
result_2020 = []
tt_2020 = []
pr_2020 = []                
# On va détecter les centiles des P>1mm et relever la température associée
PR_w_precip = pr_all.tp.where( (pr_all.tp >= 1) )
# On va détecter les journées avec précipitation et relever la température associée           
TT_w_precip = tt_all.t2m.where( (pr_all.tp >= 1) )

# On filtre les points de grille au-dessus du bassin versant          
PR_w_precip_BV = PR_w_precip.where(MASK.tp >= 0 )
TT_w_precip_BV = TT_w_precip.where(MASK.tp >= 0 )
                
# On effectue une moyenne spatiale par jour       
PR_masked =  PR_w_precip_BV
TT_masked =  TT_w_precip_BV
    
# On crée un masque 1D des points de grille sans nan
PR_masked = PR_w_precip_BV.stack(dim=['time','latitude','longitude']).notnull()
TT_masked = TT_w_precip_BV.stack(dim=['time','latitude','longitude']).notnull()

# On applique ce asque sur les données 1D  
PR_stacked = PR_w_precip_BV.stack(dim=['time','latitude','longitude'])[PR_masked].values
TT_stacked = TT_w_precip_BV.stack(dim=['time','latitude','longitude'])[TT_masked].values
    
tt_2020.append(pd.DataFrame(TT_stacked))          
pr_2020.append(pd.DataFrame(PR_stacked)) 
      
                            
data_pr_2020= pd.concat(pr_2020) 
data_tt_2020 = pd.concat(tt_2020) 
         
result_2020 = pd.concat([ data_pr_2020,data_tt_2020],axis=1)
result_2020.columns = ['Precipitation', 'Temperature']

# POUR 2019
year = 2019
file = 'J:/REANALYSES/ERA5/PR_1h_Outaouais/era5_pr_1h_'
multi_file = [f'{file}{year-1}{month}_sfc.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}_sfc.nc' for month in ['01','02','03']])
pr_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
pr_all = pr_all.assign_coords(longitude=(((pr_all.longitude + 180) % 360) - 180)).sortby('longitude') 
pr_all = pr_all.sel(time=slice(str(year-1)+'-11-01', str(year)+'-03-'+last_day))*1000
pr_all = pr_all.resample(time = '1D').sum() 

file = 'J:/REANALYSES/ERA5/T2m_1h_Outaouais/era5_t2m_1h_'
multi_file = [f'{file}{year-1}{month}_sfc.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}_sfc.nc' for month in ['01','02','03']])
tt_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
tt_all = tt_all.assign_coords(longitude=(((tt_all.longitude + 180) % 360) - 180)).sortby('longitude') 
tt_all = tt_all.sel(time=slice(str(year-1)+'-11-01', str(year)+'-03-'+last_day)) - 273.15
tt_all = tt_all.resample(time = '1D').mean() 

result_2019 = []
tt_2019 = []
pr_2019 = []                
# On va détecter les centiles des P>1mm et relever la température associée
PR_w_precip = pr_all.tp.where( (pr_all.tp >= 1) )
# On va détecter les journées avec précipitation et relever la température associée           
TT_w_precip = tt_all.t2m.where( (pr_all.tp >= 1) )

# On filtre les points de grille au-dessus du bassin versant          
PR_w_precip_BV = PR_w_precip.where(MASK.tp >= 0 )
TT_w_precip_BV = TT_w_precip.where(MASK.tp >= 0 )
                
# On effectue une moyenne spatiale par jour       
PR_masked =  PR_w_precip_BV
TT_masked =  TT_w_precip_BV
    
# On crée un masque 1D des points de grille sans nan
PR_masked = PR_w_precip_BV.stack(dim=['time','latitude','longitude']).notnull()
TT_masked = TT_w_precip_BV.stack(dim=['time','latitude','longitude']).notnull()

# On applique ce asque sur les données 1D  
PR_stacked = PR_w_precip_BV.stack(dim=['time','latitude','longitude'])[PR_masked].values
TT_stacked = TT_w_precip_BV.stack(dim=['time','latitude','longitude'])[TT_masked].values
    
tt_2019.append(pd.DataFrame(TT_stacked))          
pr_2019.append(pd.DataFrame(PR_stacked)) 
      
                            
data_pr_2019= pd.concat(pr_2019) 
data_tt_2019 = pd.concat(tt_2019) 
         
result_2019 = pd.concat([ data_pr_2019, data_tt_2019],axis=1)
result_2019.columns = ['Precipitation', 'Temperature']

bins_T = list(range(-20,10,2))
lst = list(range(-20,12,2))
lst2 = list(range(-18,12,2))
labels = [format(x, '02d') for x in lst]
labels2 = [format(x, '02d') for x in lst2]
label_bins = [x+':'+y for x,y in zip(labels,labels2)] 
 
result_2019['temp_bins'] = pd.cut(x=result_2019['Temperature'], bins=lst, labels=label_bins)
result_2020['temp_bins'] = pd.cut(x=result_2020['Temperature'], bins=lst, labels=label_bins)



list_tmp=[]
for label in label_bins:
    tmp = pd.DataFrame(result_2019['Precipitation'].loc[result_2019['temp_bins'] == label].values)
    tmp.columns = [label]
    list_tmp.append(tmp)
final_2019 = pd.concat(list_tmp, sort=False)

list_tmp=[]
for label in label_bins:
    tmp = pd.DataFrame(result_2020['Precipitation'].loc[result_2020['temp_bins'] == label].values)
    tmp.columns = [label]
    list_tmp.append(tmp)
final_2020 = pd.concat(list_tmp, sort=False)

final_2019 = final_2019.assign(Location=1)
final_2020 = final_2020.assign(Location=2)


cdf = pd.concat([final_2019, final_2020])

mdf = pd.melt(cdf, id_vars=['Location'], var_name=['temp_bins'])
fig = plt.figure(figsize=(28,16))
ax = sns.boxplot(x="temp_bins", 
                 y="value",
                 hue="Location", 
                 data=mdf, 
                 showfliers=False,
                 palette=[sns.xkcd_rgb["medium green"], 
                          sns.xkcd_rgb["medium blue"],
                          sns.xkcd_rgb["pale red"]],
                 )    # https://xkcd.com/color/rgb/


   


#ax.set(ylim=(0, 50))
#plt.legend(title='Smoker', loc='upper left', labels=['RCMs historical', 'RCMs rcp45 scenario', 'RCMs rcp85 scenario'])
handles, _ = ax.get_legend_handles_labels()
ax.legend(handles, legends ,prop={'size':30})
ax.set_title('Bassin versant des Outaouais: 1er Novembre au '+str(last_day)+' Mars' , fontdict={'fontsize': 30, 'fontweight': 'bold'})
ax.set_ylabel('Précipitation quotidienne [mm]', fontdict={'fontsize': 30, 'fontweight': 'bold'})
ax.set_xlabel('Température moyenne journalière [Celcius]', fontdict={'fontsize': 30, 'fontweight': 'bold'})
for item in ax.get_xticklabels():
    item.set_rotation(45)
    
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)    
    

fig = plt.gcf()
fig.subplots_adjust(top=0.9)     
plt.savefig(('./figures_update/Precipitation_Daily_by_TMEAM_Bins_All_points_November_to_March.png'), dpi=300, bbox_inches='tight')



