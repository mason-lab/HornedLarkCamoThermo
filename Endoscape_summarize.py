#-------------------------#
#------INSTRUCTIONS-------#
#-------------------------#
'''The following Python script was used to summarize the data from the biophysical simulations on Horned Larks.
The script calculates the cooling benefits between birds with average dorsal absorptance and the observed dorsal absorptance.'''

#---------------------#
#------LIBRARIES------#
#---------------------#
import glob as glob
import pandas as pandas
from collections import OrderedDict
import numpy
from math import log as log

#---------------------#
#-------SCRIPT--------#
#---------------------#
#Summarize the data for simulations that incorporated average dorsal absorptance
files = glob.glob('Output/avg_absorptance_*.csv')
hourly_results = pandas.DataFrame(columns = ['type','species','activity_pattern','site','version','climate_scenario','mass_change','seek_temperature','wind','posture','fiber_density','shade','feather_depths','orientation','shape','physiology_known','mass','surface_area','julian_day','hour','Srad','Rabs','Remi','Tes','Tref','Te_dorsal','Hi_d','Ri_d','Ks','Kf_d','Ksfi_d','Ke_d','Ke_empirical','Ke_theoretical','Tb','Qgen','EHL','MHP','EHL/MHP','energy_balance','water_loss_mass','dehydration_mass','excess_proportion_of_mass_lost','excess_water','metabolic_heat_required','water_loss_required','proportion_water_loss_required','excess_metabolic_heat','Qgen_water_avg','Qgen_water_sum','Qgen_heat_avg','Qgen_heat_sum','activity_thermal_stress','activity_water_stress','water_stress','heat_stress','hr_restriction'])
for i in range(len(files)):
    temporary = pandas.read_csv(files[i])
    hourly_results = hourly_results.append(temporary)

hourly_results.to_csv('Output/lark_avg_concatenated.csv',columns=['type','species','activity_pattern','site','version','climate_scenario','mass_change','seek_temperature','wind','posture','fiber_density','shade','feather_depths','orientation','shape','physiology_known','mass','surface_area','julian_day','hour','Srad','Rabs','Remi','Tes','Tref','Te_dorsal','Hi_d','Ri_d','Ks','Kf_d','Ksfi_d','Ke_d','Ke_empirical','Ke_theoretical','Tb','Qgen','EHL','MHP','EHL/MHP','energy_balance','water_loss_mass','dehydration_mass','excess_proportion_of_mass_lost','excess_water','metabolic_heat_required','water_loss_required','proportion_water_loss_required','excess_metabolic_heat','Qgen_water_avg','Qgen_water_sum','Qgen_heat_avg','Qgen_heat_sum','activity_thermal_stress','activity_water_stress','water_stress','heat_stress','hr_restriction'],index = False)

daily_results = hourly_results.groupby(['species','site','climate_scenario'],as_index=False)
daily_df1 = daily_results.agg(OrderedDict([
                 ('Tes' , numpy.mean),
                 ('Qgen_water_sum' , sum),
                 ('Qgen_heat_sum' , sum),
                 ('activity_water_stress' , sum),
                 ('hr_restriction' , sum)
                ]))
daily_df1['Cooling_costs_avg'] = (daily_df1['Qgen_water_sum']/12.0)/-1000.0
daily_df1['Heating_costs_avg'] = (daily_df1['Qgen_heat_sum']/12.0)/1000.0
daily_df1.to_csv('Output/lark_avg_abs.csv',index = False)

#Summarize the data for simulations that incorporated observed dorsal absorptance
files = glob.glob('Output/dorsal_absorptance*.csv')
hourly_results = pandas.DataFrame(columns = ['type','species','activity_pattern','site','version','climate_scenario','mass_change','seek_temperature','wind','posture','fiber_density','shade','feather_depths','orientation','shape','physiology_known','mass','surface_area','julian_day','hour','Srad','Rabs','Remi','Tes','Tref','Te_dorsal','Hi_d','Ri_d','Ks','Kf_d','Ksfi_d','Ke_d','Ke_empirical','Ke_theoretical','Tb','Qgen','EHL','MHP','EHL/MHP','energy_balance','water_loss_mass','dehydration_mass','excess_proportion_of_mass_lost','excess_water','metabolic_heat_required','water_loss_required','proportion_water_loss_required','excess_metabolic_heat','Qgen_water_avg','Qgen_water_sum','Qgen_heat_avg','Qgen_heat_sum','activity_thermal_stress','activity_water_stress','water_stress','heat_stress','hr_restriction'])
for i in range(len(files)):
    temporary = pandas.read_csv(files[i])
    hourly_results = hourly_results.append(temporary)

hourly_results.to_csv('Output/lark_dorsal_concatenated.csv',columns=['type','species','activity_pattern','site','version','climate_scenario','mass_change','seek_temperature','wind','posture','fiber_density','shade','feather_depths','orientation','shape','physiology_known','mass','surface_area','julian_day','hour','Srad','Rabs','Remi','Tes','Tref','Te_dorsal','Hi_d','Ri_d','Ks','Kf_d','Ksfi_d','Ke_d','Ke_empirical','Ke_theoretical','Tb','Qgen','EHL','MHP','EHL/MHP','energy_balance','water_loss_mass','dehydration_mass','excess_proportion_of_mass_lost','excess_water','metabolic_heat_required','water_loss_required','proportion_water_loss_required','excess_metabolic_heat','Qgen_water_avg','Qgen_water_sum','Qgen_heat_avg','Qgen_heat_sum','activity_thermal_stress','activity_water_stress','water_stress','heat_stress','hr_restriction'],index = False)

daily_results = hourly_results.groupby(['species','site','climate_scenario'],as_index=False)
daily_df2 = daily_results.agg(OrderedDict([
                 ('Tes' , numpy.mean),
                 ('Qgen_water_sum' , sum),
                 ('Qgen_heat_sum' , sum),
                 ('activity_water_stress' , sum),
                 ('hr_restriction' , sum)
                ]))
daily_df2['Cooling_costs_dorsal'] = (daily_df2['Qgen_water_sum']/12.0)/-1000.0
daily_df2['Heating_costs_dorsal'] = (daily_df2['Qgen_heat_sum']/12.0)/1000.0
daily_df2.to_csv('Output/lark_dorsal_abs.csv',index = False)

#Now create the final dataframe that contains the cooling benefits (Delta cooling_avg_dorsal)
newDF = pandas.DataFrame()
newDF['site'] = daily_df3['site']
newDF['Cooling_costs_avg'] = daily_df1['Cooling_costs_avg']
newDF['Heating_costs_avg'] = daily_df1['Heating_costs_avg']
newDF['Cooling_costs_dorsal'] = daily_df2['Cooling_costs_dorsal']
newDF['Heating_costs_dorsal'] = daily_df2['Heating_costs_dorsal']
newDF['Delta_cooling_avg_dorsal'] = newDF['Cooling_costs_avg'] - newDF['Cooling_costs_dorsal']
newDF['Delta_heating_avg_dorsal'] = newDF['Heating_costs_avg'] - newDF['Heating_costs_dorsal']


newDF.to_csv('Output/cooling_benefit.csv',index = False)

print('Summarizing is complete')
