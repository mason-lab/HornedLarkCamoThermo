#-------------------------#
#------INSTRUCTIONS-------#
#-------------------------#

'''To run the script, your root directory needs to be the folder
containing all of the .txt files containing the reflectance measurements.
You will also need to move the solar_specrum.csv into the same directory with
the .txt files.'''

#---------------------#
#------LIBRARIES------#
#---------------------#

import pandas
import numpy
import glob
from collections import OrderedDict

#---------------------#
#------FUNCTIONS------#
#---------------------#

def return_species_name_list():
    'returns a unique list of names in a directory using the first 7 letters of the filename'
    total_name_list = []
    list_of_files = glob.glob(str('*txt'))
    for each_file in range(len(list_of_files)):
        total_name_list.append(list_of_files[each_file][0:7])
    uniques = set(total_name_list)
    return list(uniques)
    
def return_number(species):
    'same as return_species_name function, but returns the number of individuals in a species'
    total_name_list = []
    list_of_files = glob.glob(str('*txt'))
    for each_file in range(len(list_of_files)):
        total_name_list.append(list_of_files[each_file][0:7])
    number_of = []
    for each_spp in range(len(species)):
        number_of.append(total_name_list.count(species[each_spp]))
    return(number_of)

def make_master_dataframe(list_of_files):
    'concatenates a large dataframe for each reflectance measurement'
    spp_dataframe = pandas.DataFrame(columns = ['spp','file_id','specimen_id','measurement_id','side','wavelength','reflectance'])
    measurement_counter = 0.0
    specimen_counter = 0.0
    for each_file in range(len(list_of_files)):#list_of_files
        if measurement_counter > 9:
            specimen_counter += 1.0
            measurement_counter = 0.0
        if measurement_counter >= 0 and measurement_counter <5:
            side_of_measurement = 'ventral'
        else:
            side_of_measurement = 'dorsal'
        data = pandas.read_table(list_of_files[each_file])
        data.columns = ['wavelength','reflectance']
        spp_dataframe = spp_dataframe.append(pandas.DataFrame({"spp":numpy.array([list_of_files[each_file][0:7]]*2151),"file_id":numpy.array([list_of_files[each_file][8:11]]*2151),"specimen_id":numpy.repeat(numpy.array(specimen_counter),2151),"measurement_id":numpy.repeat(numpy.array(measurement_counter),2151),"side":numpy.repeat(numpy.array(side_of_measurement),2151),"wavelength":numpy.arange(350,2501),"reflectance":numpy.array(data.reflectance)}))
        measurement_counter += 1
    return spp_dataframe
    
def reflectance_summary(master_db):
    'calculates the reflectance from the master dataframe'
    mean_df = master_db.groupby(['spp','side','wavelength']).mean()
    mean_df = mean_df.reset_index()
    std_df = master.groupby(['spp','side','wavelength']).std()
    std_df = std_df.reset_index()
    mean_df['std_reflectance'] = std_df['reflectance']
    mean_df['upper'] = mean_df['reflectance'] + mean_df['std_reflectance']
    mean_df['lower'] = mean_df['reflectance'] - mean_df['std_reflectance']
    #spp_dataframe = spp_dataframe.append(pandas.DataFrame({"spp":mean_df.spp,"side":mean_df.side,"wavelength":mean_df.wavelength,"mean_reflectance":mean_df.reflectance,"std_reflectance":std_df.reflectance}))
    return mean_df
    
def mean_reflectance(master_db):
    'calculates the mean reflectance for each subspecies on each side of the bird'
    mean_df = master_db.groupby(['spp','side']).mean()
    mean_df = mean_df.reset_index()
    std_df = master.groupby(['spp','side']).std()
    std_df = std_df.reset_index()
    mean_df['std_reflectance'] = std_df['reflectance']
    return mean_df
            
def create_massive_df(directory):
    'creates a large dataframe consisting of the species, measurement, reflectance, and wavelength'
    list_of_files = glob.glob(str(directory)+'*txt')
    sorted_files = sort(list_of_files)
    final_dataframe = pandas.DataFrame(columns = ['spp','value','wavelength','reflectance'])
    for each_file in range(len(sorted_files)):
       data = pandas.read_table(sorted_files[each_file])
       data.columns = ['wavelength','reflectance']
       file_df = pandas.DataFrame({"spp":sorted_files[each_file][12:19],"value":sorted_files[each_file][20:23],"wavelength":data.wavelength,"reflectance_avg":data.reflectance,"reflectance_std":data.reflectance})
       final_dataframe = final_dataframe.append(file_df)
    sample_size = len(list_of_files)/10
    id_ = []
    for i in range(sample_size):
        id_.append([i]*21510)
    id_ = [item for sublist in id_ for item in sublist]
    final_dataframe['individual'] = id_
    final_dataframe['side'] = (['ventral']*2151*5 + ['dorsal'] * 2151*5)*(sample_size)
    return final_dataframe
       
        
def calculate_individual_values(dataframe):
    'calculates reflectance for each individual'
    grouped = dataframe.groupby(['spp','individual','side','wavelength'],as_index=False)
    output = grouped.agg(OrderedDict([
                 ('reflectance_avg' , numpy.mean),
                 ('reflectance_std' , numpy.std)
                ]))
    return(output)
    
def calculate_individual_averages(dataframe):
    'calculates reflectance for each individual'
    grouped = dataframe.groupby(['spp','individual','side'],as_index=False)
    output = grouped.agg(OrderedDict([
                ('reflectance_avg' , numpy.mean),
                ('reflectance_std' , numpy.std)
                ]))
    return(output)
    
def solar_correction_group(dataframe):
    'calculates averages of solar-corrected absorptance'
    grouped = dataframe.groupby(['spp','individual','side'],as_index=False)
    output = grouped.agg(OrderedDict([
                ('solar_absorptance' , numpy.sum),
                ('reflectance_std' , numpy.mean),
                ('total_solar_watts', numpy.mean)
                ]))
    return(output)
    
    
#---------------------#
#-------SCRIPT--------#
#---------------------#

'This script calculates mean and std reflectance across the whole wavelength'
list_of_files = glob.glob(str('*txt'))
list_of_files.sort()
master = make_master_dataframe(list_of_files)
final = reflectance_summary(master)
final.to_csv('raw_reflectance.csv',index=False,header=True,columns=['spp','side','wavelength','reflectance','std_reflectance'])

'This script was made to create lists of species abbreviations and sample totals within the current working directory:'
spp = return_species_name_list()
num  = return_number(spp)
df = pandas.DataFrame({"spp":spp,"num":num})
df = df.sort_values(by = 'spp')

'This script was made to create lists of species abbrevs and sample totals within the current working directory:'
spp1 = return_species_name_list()
num1 = return_number(spp1)
df1 = pandas.DataFrame({"spp":spp1,"num":num1})
df1 = df1.sort_values(by = 'spp')

'This script makes dataframes to assess individual variation'
df = create_massive_df('horned_larks/')
df_c = calculate_individual_values(df)
df_i = calculate_individual_averages(df_c)
df_i_dorsal = df_i.drop(df_i[df_i.side == 'ventral'].index)
df_i_ventral = df_i.drop(df_i[df_i.side == 'dorsal'].index)

df_i_dorsal = df_i_dorsal.reset_index(drop=True)
df_i_ventral = df_i_ventral.reset_index(drop=True)

df_total = df_i_dorsal
df_total['ventral_reflectance_avg'] = df_i_ventral['reflectance_avg']
df_total['ventral_reflectance_std'] = df_i_ventral['reflectance_std']

df_total = df_total.drop(['side'],axis=1)
df_total.columns = ['spp', 'individual','dorsal_reflectance_avg','dorsal_reflectance_std','ventral_reflectance_avg','ventral_reflectance_std']
df_total.to_csv('individual_reflectances.csv',index = False)

'This script calculates the solar correction for absorptance'
solar_dataframe = pandas.read_csv('solar_spectrum.csv') 
df_c['reflectance_avg'] = 1- df_c['reflectance_avg']
combined_df = df_c.set_index('wavelength').combine_first(solar_dataframe.set_index('wavelength')).reset_index()
combined_df['solar_absorptance'] = combined_df['watts']*combined_df['reflectance_avg']
combined_df['total_solar_watts'] = sum(solar_dataframe.watts)
combined_df = combined_df.dropna()

solar_df = solar_correction_group(combined_df)
solar_df['absorptance_corrected'] = solar_df['solar_absorptance']/solar_df['total_solar_watts']
solar_df.to_csv('relflectance_corrected.csv',index = False)

solar_dorsal = solar_df.drop(solar_df[solar_df.side == 'ventral'].index)
solar_ventral = solar_df.drop(solar_df[solar_df.side == 'dorsal'].index)

solar_dorsal = solar_dorsal.reset_index(drop=True)
solar_ventral = solar_ventral.reset_index(drop=True)

solar_total = solar_dorsal
solar_total['ventral_absorptance_corrected_avg'] = solar_ventral['absorptance_corrected']
solar_total['ventral_absorptance_corrected_std'] = solar_ventral['reflectance_std']
solar_total = solar_total.drop(['side'],axis=1)
solar_total.columns = ['spp', 'individual','sum_solar_absorptance','dorsal_reflectance_std','total_solar_watts','dorsal_absorptance_corrected','ventral_absorptance_corrected_avg','ventral_reflectance_corrected_std']
solar_total.to_csv('relflectance_corrected.csv',index = False)