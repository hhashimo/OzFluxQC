import calendar
import datetime
import numpy
import matplotlib.pyplot as plt
import os
import scipy
import sys
import pdb

# code from dark_T_response_functions.py
def TRF(data_dict, Eo, rb):
    return rb * numpy.exp(Eo * (1 / (10 + 46.02) - 1 / (data_dict['TempC'] + 46.02)))

def optimise_rb(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    
    # Get drivers and response
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE']        
    
    try:
        params = scipy.optimize.curve_fit(lambda x, b:
                           TRF(x, params_dict['Eo_default'], b),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['rb_prior']])[0]                           
    except RuntimeError:
        params = [numpy.nan]

    # If negative rb returned, set to nan
    if params[0] < 0:
        error_state = 9
        params = [numpy.nan]
       
    return params, error_state

# code from Partition_NEE.py
def get_dates(datetime_array, configs_dict):

    # Assign configs to local vars
    window = configs_dict['window_size_days']

    # Create a series of continuous whole day dates that will be used for output
    # (parameter series will be interpolated between window centres)
    start_date = datetime_array[0].date()
    end_date = datetime_array[-1].date()
    num_days = (end_date - start_date).days + 1 # Add 1 so is inclusive of both end members
    all_dates_array = numpy.array([start_date + datetime.timedelta(i) for i in xrange(num_days)])

    # Create a shifted array
    shift_mins = 60 * configs_dict['measurement_interval']
    shift_datetime_array = datetime_array - datetime.timedelta(minutes = shift_mins)

    # Check that first and last days are complete and revise start and end dates if required
    temp_date = datetime.datetime.combine((shift_datetime_array[0] + datetime.timedelta(1)).date(), 
                                    datetime.datetime.min.time())
    num_obs = len(numpy.where(shift_datetime_array < temp_date)[0])
    if num_obs < 24 * (1 / configs_dict['measurement_interval']):
        start_date = start_date + datetime.timedelta(1)
    temp_date = datetime.datetime.combine(shift_datetime_array[-1].date(), 
                                    datetime.datetime.min.time())
    num_obs = len(numpy.where(shift_datetime_array >= temp_date)[0])
    if num_obs < 24 * (1 / configs_dict['measurement_interval']):
        end_date = end_date - datetime.timedelta(1)
    
    # Calculate the dates that represent the centre of the window for each step
    num_days = (end_date - start_date).days + 1 - window # Add 1 so is inclusive of both end members
    first_fit_day = start_date + datetime.timedelta(window / 2)
    step_days = numpy.arange(0, num_days, configs_dict['step_size_days'])
    step_dates_array = [first_fit_day + datetime.timedelta(i) for i in step_days]

    # Make an index dictionary for step dates    
    step_dates_index_dict = {}
    for date in step_dates_array:    
        date_time = (datetime.datetime.combine(date, datetime.datetime.min.time()) 
                     + datetime.timedelta(hours = 12))
        start_date = date_time - datetime.timedelta(window / 2.0)
        end_date = date_time + datetime.timedelta(window / 2.0)
        start_ind = numpy.where(datetime_array == start_date)[0].item() + 1
        end_ind = numpy.where(datetime_array == end_date)[0].item()
        step_dates_index_dict[date] = [start_ind, end_ind]
    
    # Make an index dictionary for all dates
    all_dates_index_dict = {}
    for date in all_dates_array:
        date_time = datetime.datetime.combine(date, datetime.datetime.min.time())
        if date == all_dates_array[0]:
            start_ind = 0
        else:
            start_date = date_time + datetime.timedelta(hours = configs_dict['measurement_interval'])
            start_ind = numpy.where(datetime_array == start_date)[0].item()
        if date == all_dates_array[-1]:
            end_ind = len(datetime_array)
        else:
            end_date = date_time + datetime.timedelta(1)
            end_ind = numpy.where(datetime_array == end_date)[0].item()
        all_dates_index_dict[date] = [start_ind, end_ind]
    
    # Make an index dictionary for years
    years_index_dict = {}
    year_array = numpy.array([i.year for i in shift_datetime_array])
    year_list = list(set(year_array))
    for yr in year_list:
        index = numpy.where(year_array == yr)[0]
        years_index_dict[yr] = [index[0], index[-1]]
    
    return step_dates_index_dict, all_dates_index_dict, years_index_dict

def make_initial_guess_dict(data_d):

    # Calculate the parameter values that are intialised from data
    index = numpy.where(data_d['PAR'] < 10)[0]
    daytime_NEE_mean = numpy.nanmean(data_d['NEE'][index])
    daytime_NEE_range = (numpy.nanpercentile(data_d['NEE'][index], 5) - 
                         numpy.nanpercentile(data_d['NEE'][index], 95))

    params_dict = {'Eo_prior': 100,
                   'k_prior': 0,
                   'alpha_prior': -0.01,
                   'rb_prior': daytime_NEE_mean,
                   'beta_prior': daytime_NEE_range,
                   'alpha_default': 0,
                   'beta_default': 0,
                   'k_default': 0 }
    
    return params_dict

def optimise_annual_Eo(data_dict, params_dict, configs_dict, year_index_dict):
    
    # Initialise local variables with configurations
    min_pct = configs_dict['minimum_pct_annual']
    msmt_int = configs_dict['measurement_interval']
    
    # Get Eo for each year and compile dictionary
    yearsEo_dict = {}
    yearsQC_dict = {}
    Eo_pass_keys = []
    Eo_range_fail_keys = []
    Eo_nan_fail_keys = []
    year_list = year_index_dict.keys()
    print 'Eo optimised using whole year is as follows:'
    for yr in year_list:

        # Calculate number of recs for year
        days = 366 if calendar.isleap(yr) else 365
        recs = days * (24 / msmt_int) / 2

        # Subset data        
        sub_dict = subset_window(data_dict, year_index_dict[yr])
        sub_dict = subset_nan(sub_dict)
        noct_flag = True
        sub_dict = subset_daynight(sub_dict, noct_flag)
        
        # Calculate percent of potential annual data that the subset contains
        pct = round(float(len(sub_dict['NEE'])) / recs * 100)

        # Fit L&T parameters if minimum data criterion satisfied, otherwise nan
        if pct > min_pct:
            params, error_code = dark.optimise_all(sub_dict, params_dict)
        else:
            params, error_code = [numpy.nan, numpy.nan], 10                                         

        # Assign year to pass, range_fail or nan_fail list for subsequent QC and fill
        Eo = params[0]
        yearsEo_dict[yr] = Eo
        yearsQC_dict[yr] = error_code
        if numpy.isnan(Eo):
            Eo_nan_fail_keys.append(yr)
        elif ((Eo < 50) | (Eo > 400)):
            Eo_range_fail_keys.append(yr)
        else:
            Eo_pass_keys.append(yr)

        print '    - ' + str(yr) + ': ' + str(round(params[0], 1))
    
    # Do QC on Eo
    if len(Eo_pass_keys) != len(yearsEo_dict):
        if len(Eo_nan_fail_keys) == len(yearsEo_dict):
            print 'Could not find any values of Eo for any years! Exiting...'
            sys.exit()
        elif len(Eo_pass_keys) != 0:
            Eo_mean = numpy.array([yearsEo_dict[i] for i in Eo_pass_keys]).mean()
            all_fail_keys = Eo_range_fail_keys + Eo_nan_fail_keys
            for i in (all_fail_keys):
                yearsEo_dict[i] = Eo_mean
            all_fail_keys = [str(key) for key in all_fail_keys]
            if len(all_fail_keys) > 1:
                all_fail_str = ', '.join(all_fail_keys)
            else:
                all_fail_str = all_fail_keys[0]
            print 'Eo optimisation failed for the following years: ' + all_fail_str
            print 'Eo for these years estimated from the mean of all other years'
        else:
            for i in Eo_range_fail_keys:
                yearsEo_dict[i] == 50 if yearsEo_dict[i] < 50 else 400
            Eo_mean = [yearsEo_dict[i] for i in Eo_range_fail_keys].mean()
            for i in Eo_nan_fail_keys:
                yearsEo_dict[i] = Eo_mean
            print 'Warning! Eo estimates were out of range for all years'
            print 'Low estimates have been set to lower limit (50);'
            print 'High estimates have been set to upper limit (400);'
            print 'Parameter estimates are unlikely to be robust!'
    else:
        print 'Eo estimates passed QC for all years'
        
    return yearsEo_dict, yearsQC_dict

def subset_window(data_dict, index_list):

    # Subset the arrays on basis of index list
    sub_dict = {}
    for i in data_dict.keys():
        sub_dict[i] = data_dict[i][index_list[0]: index_list[1] + 1]
    
    return sub_dict

def subset_daynight(data_dict, noct_flag):
    
    # Turn dictionary into an array
    temp_array = numpy.empty([len(data_dict['NEE']), len(data_dict)])
    for i, var in enumerate(data_dict.keys()):
        temp_array[:, i] = data_dict[var]

    # Create night / day subsetting index and subset data
    if noct_flag:
        daynight_index = numpy.where(data_dict['PAR'] < 10)[0]
    else:
        daynight_index = numpy.where(data_dict['PAR'] > 10)[0]
    temp_array = temp_array[daynight_index]
    
    sub_dict = {var: temp_array[:, i] for i, var in enumerate(data_dict.keys())}

    return sub_dict

def subset_nan(data_dict):
    
    # Turn dictionary into an array
    temp_array = numpy.empty([len(data_dict['NEE']), len(data_dict)])
    for i, var in enumerate(data_dict.keys()):
        temp_array[:, i] = data_dict[var]

    # Create nan subsetting index and subset data and count
    QCdata_index = numpy.where(numpy.all(~numpy.isnan(temp_array), axis=1))    
    temp_array = temp_array[QCdata_index]

    sub_dict = {var: temp_array[:, i] for i, var in enumerate(data_dict.keys())}

    return sub_dict

def estimate_Re_GPP(sub_d, params_d, GPP = False):

    return_dict = {}
    if GPP:
        GPP, Re = light.LRF_part(sub_d, params_d['Eo'], params_d['rb'],
                                 params_d['alpha'], params_d['beta'], 
                                 params_d['k'])
        return_dict['Re'] = Re
        return_dict['GPP'] = GPP
    else:
        Re = dark.TRF(sub_d, params_d['Eo'], params_d['rb'])
        return_dict['Re'] = Re
    return return_dict

def plot_windows(data_dict, configs_dict, date, noct_flag):

    # Set parameters from dicts
    path = configs_dict['window_plot_output_path']
    window = configs_dict['window_size_days']
    
    for i in range(2):
        sub_d = subset_daynight(data_dict, noct_flag)
        if noct_flag:
            daynight_ind = 'noct'
            x_lab = r'Temperature ($^{o}C$)'
            x_var = sub_d['TempC']
            y_var1 = sub_d['NEE']
            y_var2 = sub_d['Re']
        else:            
            daynight_ind = 'day'
            x_lab = r'PAR ($\mu mol\/photons\/m^{-2}s^{-1}$)'
            x_var = sub_d['PAR']
            y_var1 = sub_d['NEE']
            y_var2 = sub_d['NEE_est']
              
        # Plot
        date_str = datetime.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        plt.plot(x_var, y_var1, 'bo' , label = 'NEE_obs')
        plt.plot(x_var, y_var2, 'ro', label = 'NEE_est')
        plt.title('Fit for ' + str(window) + ' day window centred on ' + 
                  date_str + '\n', fontsize = 22)
        plt.xlabel(x_lab, fontsize = 16)
        plt.ylabel(r'NEE ($\mu mol C\/m^{-2} s^{-1}$)', fontsize = 16)
        plt.axhline(y = 0, color = 'black')
        plot_out_name = daynight_ind + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(path, plot_out_name))
        plt.close(fig)
        
    return

def interp_params(param_rslt_array):

    def do_interp(array_1D):
        xp = numpy.arange(len(arr))
        fp = array_1D[:]
        nan_index = numpy.isnan(fp)
        fp[nan_index] = numpy.interp(xp[nan_index], xp[~nan_index], fp[~nan_index])
        return fp   
    
    arr = param_rslt_array.copy()    
    num_vars = numpy.shape(arr)
    if len(num_vars) == 1:
        arr = do_interp(arr)
    else:
        num_vars = num_vars[1]
        for i in range(num_vars):
            arr[:, i] = do_interp(arr[:, i])

    return arr            
