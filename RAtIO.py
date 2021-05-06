#################
# Rapid         #
# AuTomatic     #
# Isotope       #
# Optimisation  #
#################

from my_imports import *

def velocity_correction(wavelength, rv):
    ''' 
    Transforming the wavelengths in velocity space

    Lets hope my algebra is right. Should take the velocity and scale it depending on the rv_correction array that should be defined at the top of the notebook.

    Parameters
    ----------
    wavelength : array
        The array you want to scale
    rv : int
        The radial velocity correction

    Returns
    -------
    arr
        The wavelength in the rest frame
    '''

    return wavelength * (1+(-rv/300000))

def call_pymoogi(filename):
    os.system('echo q | pymoogi ' + filename)

def get_region(r):
    if r == 0:
        lw = 5134.42
        uw = 5140.46
    elif r == 1:
        lw = 5134.42
        uw = 5134.85
    elif r == 2:
        lw = 5138.59
        uw = 5138.95
    elif r == 3:
        lw = 5140.04
        uw = 5140.46
    else:
        print('wavelength region error')
        lw = 0
        uw = 1
    return lw, uw

def get_lines(r):
    if r == 0:
        wl  = [5134.569, 5134.653, 5134.734,5138.710, 5138.768, 5138.785, 
        5138.823, 5138.860, 5140.202, 5140.251, 5140.286, 5140.302, 5140.358]
        iso = [24, 25, 26, 24, 25, 25, 26, 26, 24, 25, 25, 26, 26]
    elif r == 1:
        wl  = [5134.569, 5134.653, 5134.734]
        iso = [24, 25, 26]
    elif r == 2:
        wl  = [5138.710, 5138.768, 5138.785, 5138.823, 5138.860]
        iso = [24, 25, 25, 26, 26]
    elif r == 3:
        wl  = [5140.202, 5140.251, 5140.286, 5140.302, 5140.358]
        iso = [24, 25, 25, 26, 26]
    else:
        print('line region error')
        wl = [0]
        iso = [0]
    return wl, iso

def calc_chi(raw, r):
    
    # hard coded wavelength bounds
    lw, uw = get_region(r)

    raw_flux = raw[(raw.wavelength > lw) & (raw.wavelength < uw)]
    
    return chisquare(raw_flux.flux, raw_flux.model_flux)[0]

def interp_smooth(raw, smooth):

    # perform a cubic spline on the data to make the wavelengths line up with eachother
    tck = interpolate.splrep(smooth.wavelength, smooth.flux, s=0)

    # evaluate the value of the spline at the wavelength points of the original spectra
    new_flux = interpolate.splev(raw.wavelength, tck, der=0)
    
    # add the model flux to the dataframe
    raw['model_flux'] = pd.Series(new_flux, index=raw.index)
    
    # return the dataframe with the new interpolated column in it 
    return raw
    
def make_model_plots(raw, smooth, output_filename, region, rv):

    fig = plt.figure(constrained_layout=True, figsize = (8, 3))
    gs = fig.add_gridspec(2, 1, height_ratios = [1, 0.3])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # shift things backwards so it looks correct whe you plot things
    wavelength_shifted = velocity_correction(raw.wavelength, -rv)

    # Getting the plot bounds
    lw, uw = get_region(region)
    cropped_flux = raw[(raw.wavelength > lw) & (raw.wavelength < uw)].flux
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()

    #getting the positions to plot the lines
    wl, iso = get_lines(region)
    text_24 = False
    text_25 = False
    text_26 = False

    # square surrounding the fitting region
    ax1.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)

    col24 = '#73454E'
    col25 = '#995C68'
    col26 = '#BA8C95'

    # for plotting the vertical lines of the isotopes
    for i in range(len(wl)):
        if iso[i] == 24:
            ax1.plot([wl[i], wl[i]], [min_flux - 0.06, max_flux + 0.05], color = col24, linestyle = '--', dashes=(5, 1))
            if not text_24:
                ax1.text(wl[i] - 0.0, min_flux - 0.09, r'$^{24}\rm{MgH}$', color = col24, fontsize=10)
                text_24 = True
        if iso[i] == 25:
            ax1.plot([wl[i], wl[i]], [min_flux - 0.04, max_flux + 0.05], color = col25,linestyle = '--', dashes=(3, 1))
            if not text_25:
                ax1.text(wl[i] - 0.0, min_flux - 0.07, r'$^{25}\rm{MgH}$', color = col25, fontsize=10)
                text_25 = True
        if iso[i] == 26:
            ax1.plot([wl[i], wl[i]], [min_flux - 0.02, max_flux + 0.05], color = col26,linestyle = '--', dashes=(1, 1))
            if not text_26:
                ax1.text(wl[i] - 0.0, min_flux - 0.05, r'$^{26}\rm{MgH}$', color = col26, fontsize=10)
                text_26 = True

    ax1.plot(wavelength_shifted, raw.model_flux, color ='#E1A4A7', label = 'model')
    #ax1.plot(smooth.wavelength, smooth.flux, color ='#eca1a6', label = 'model')
    ax1.plot(wavelength_shifted, raw.flux, '.', color ='#2a9d8f', label = 'spectra', marker = '+')
    #ax1.plot(smooth.wavelength, smooth.flux, '.', color ='#d6cbd3', label = 'model')

    ax2.plot(wavelength_shifted, raw.flux - raw.model_flux, color ='#ada397')

    ax1.legend(frameon=False)
    ax1.set_xlim(lw - 0.2, uw + 0.2)
    ax2.set_xlabel(r'Wavelength ($\AA$)')

    ax1.set_ylim(min_flux - 0.1, max_flux + 0.03)
    ax1.set_ylabel('Norm. flux')

    # Limits for the residual plot
    ax2.set_ylim(-0.024,0.024)
    

    ax1.ticklabel_format(useOffset=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.ticklabel_format(useOffset=False)

    ax1.tick_params(direction='in', axis='both', which='both', bottom=True,top=True, left=True, right=True)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(direction='in', axis='both', which='both', bottom=True,top=True, left=True, right=True)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    gs.update(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    

    plt.savefig('./plots/'+ output_filename, facecolor='white', bbox_inches='tight', dpi = 300)
    plt.close()

def make_temp_file(filename):
    # will need these files in your directory - wont make them apparently...
    f  = open(filename, "a+") 
    f.write('')
    f.close() 

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par):
    # I doubt I'll ever want to change these so initialise them here
    standard_out = 'out1'
    summary_out  = 'out2'

    # will need these files in your directory - wont make them apparently...
    make_temp_file(standard_out)
    make_temp_file(summary_out)
    make_temp_file(out_filename)

  
    par_string = "synth\n" +\
    "standard_out   '" + standard_out +"'\n"                    + \
    "summary_out    '" + summary_out +"'\n"                     + \
    "smoothed_out   '" + out_filename +"'\n"                    + \
    "model_in       't4070g040m18.newmod'\n"                    + \
    "lines_in       'linelist.MgH'\n"                           + \
    "observed_in    '" + raw_spec_filename +"'\n"               + \
    "atmosphere    1\n"                                         + \
    "molecules     2\n"                                         + \
    "lines         2\n"                                         + \
    "flux/int      0\n"                                         + \
    "plotpars      1\n"                                         + \
    wavelength_region + " 0.15 1.05\n"                          + \
    str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
    "d          0.047 0.0 0.0 "+ str(par['s']) +" 0.0\n"        + \
    "abundances   4    1\n"                                     + \
    "6            0.0500000\n"                                  + \
    "12           " + str(par['mg']) + "\n"                     + \
    "24           0.50000\n"                                    + \
    "26           0.20000\n"                                    + \
    "isotopes      4    1\n"                                    + \
    "606.01212     5.0\n"                                       + \
    "112.00124     "+ str(par['i_24']) +"\n"                    + \
    "112.00125     "+ str(par['i_25']) +"\n"                    + \
    "112.00126     "+ str(par['i_26']) +"\n"                    + \
    "obspectrum 5\n"                                            + \
    "synlimits\n"                                               + \
    wavelength_region + " 0.01 5.0\n"                           + \
    "plot 2\n"                                                  + \
    "damping 2\n"


    # writing that string to a file 
    par_file  = open(in_filename, "w+") 
    par_file.write(par_string)
    par_file.close() 
    return in_filename, out_filename

def change_s(d, increase=True):
    change_by = 0.2
    ll = 4 # lower limit 
    ul = 9 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['s'] + change_by <= ul:
        d['s'] += change_by
        return d
    elif not increase and d['s'] - change_by >= ll:
        d['s'] -= change_by
        return d
    else:
        return None

def change_mg(d, increase=True):
    change_by = 0.02
    ll = -1 # lower limit 
    ul = 1.5 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['mg'] + change_by <= ul:
        d['mg'] += change_by
        d['mg'] = round(d['mg'], 2)
        return d
    elif not increase and d['mg'] - change_by >= ll:
        d['mg'] -= change_by
        d['mg'] = round(d['mg'], 2)
        return d
    else:
        return None

def change_24(d, increase=True):
    change_by = 0.1
    ll = 0.1 # lower limit 
    ul = 8 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_24'] + change_by <= ul:
        d['i_24'] += change_by
        return d
    elif not increase and d['i_24'] - change_by >= ll:
        d['i_24'] -= change_by
        return d
    else:
        return None

def change_25(d, increase=True):
    change_by = 0.5
    ll = 1 # lower limit 
    ul = 15 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_25'] + change_by <= ul:
        d['i_25'] += change_by
        d['i_25'] = round(d['i_25'],3)
        return d
    elif not increase and d['i_25'] - change_by >= ll:
        d['i_25'] -= change_by
        d['i_25'] = round(d['i_25'],3)
        return d
    else:
        return None

def change_26(d, increase=True):
    change_by = 0.5
    ll = 1 # lower limit 
    ul = 15 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_26'] + change_by <= ul:
        d['i_26'] += change_by
        d['i_26'] = round(d['i_26'],3)
        return d
    elif not increase and d['i_26'] - change_by >= ll:
        d['i_26'] -= change_by
        d['i_26'] = round(d['i_26'],3)
        return d
    else:
        return None

def change_rv(d, increase=True):
    change_by = 1
    ll =-5 # lower limit 
    ul = 5 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['rv'] + change_by <= ul:
        d['rv'] += change_by
        d['rv'] = round(d['rv'], 2)
        return d
    elif not increase and d['rv'] - change_by >= ll:
        d['rv'] -= change_by
        d['rv'] = round(d['rv'], 2)
        return d
    else:
        return None

def read_raw_spectra(filename):
    return pd.read_table(filename, sep="\s+", 
                         header=None, names = ['wavelength', 'flux'])

def read_smoothed_spectra(filename, rv):
    # different to reading raw spectra because we have to skip some headder rows
    smooth = pd.read_table(filename, sep="\s+", header=None, skiprows = [0,1],
                         names = ['wavelength', 'flux'])
    # run_interpolation for the values of the raw spectra wavelength
    smooth.wavelength = velocity_correction(smooth.wavelength, rv)
    return smooth
    
def get_chi_squared(raw, out_filename, region, guess, make_plot = True):

    # read in the smoothed data, interpolate spectra
    smooth = read_smoothed_spectra(out_filename, guess['rv'])
    raw = interp_smooth(raw, smooth)
    
    # make a plot of the model
    if make_plot:
        make_model_plots(raw, smooth, out_filename, region, guess['rv'])

    # return the chi quared value over the line
    return calc_chi(raw, region)

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv

def get_wavelength_region(raw_wavelength):
    lower_wavelength = raw_wavelength[0]
    upper_wavelength = raw_wavelength[len(raw_wavelength)-1] # -1 isnt working for some reason
    return str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) 

def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess):

    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')


    # creates a parameter string in the directory that moog can read
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess)

    # create the smoothed spectra by calling pymoogi
    call_pymoogi(in_filename)

    # read in the smoothed model spectra and calculate the chi squared value
    cs = get_chi_squared(raw_spectra, out_filename, region, guess, make_plot = True)
    
    # return a dataframe with a single row (to be added to a larger df later)
    return pd.DataFrame({'filename'   : out_filename, 
                         'chi_squared': cs, 
                         's'          : guess['s'],
                         'mg'         : guess['mg'],
                         'i_24'       : guess['i_24'],
                         'i_25'       : guess['i_25'],
                         'i_26'       : guess['i_26'],
                         'rv'         : guess['rv'],
                         'ratio'      : calc_ratio(guess['i_24'], guess['i_25'], guess['i_26'])
                         }, index=[1])

def generate_neighbours(guess, region):
    
    # return an array of dictionaries?
    # increase and decrease each value?

    # a list of dictionaries
    new_guesses = []
    # Note: must pass a copy of the dictionary!!!
    new_guesses.append(change_s(guess.copy(), True))  # increase 
    new_guesses.append(change_s(guess.copy(), False)) # decrease 

    new_guesses.append(change_mg(guess.copy(), True))  # increase 
    new_guesses.append(change_mg(guess.copy(), False)) # decrease 

    # only optimise for these for the individual regions, not the whole thing
    if not region == -1:
        new_guesses.append(change_24(guess.copy(), True))  # increase 
        new_guesses.append(change_24(guess.copy(), False)) # decrease 

        new_guesses.append(change_25(guess.copy(), True))  # increase 
        new_guesses.append(change_25(guess.copy(), False)) # decrease 

        new_guesses.append(change_26(guess.copy(), True))  # increase 
        new_guesses.append(change_26(guess.copy(), False)) # decrease 

    #new_guesses.append(change_rv(guess.copy(), True))  # increase 
    #new_guesses.append(change_rv(guess.copy(), False)) # decrease 

    # Strip the none values
    new_guesses = [i for i in new_guesses if i != None]
    
    return new_guesses

def filter_guesses(guess_arr, chi_df):
    # make sure you arent running a model you have already done
    
    # compare the strings of making the file name to what is currently
    # in the chi squared data frame

    no_duplicates_arr = []
    
    for dict in guess_arr:
        # make the string
        dict_name = make_filenames(dict, 'out')
        have_found = False # we havent found it yet

        for file in chi_df.filename:
            if dict_name == file:
                have_found = True
        
        # if you get throught the loop and you havent found it:
        if not have_found:
            no_duplicates_arr.append(dict)

    return no_duplicates_arr

def reconstruct_min_chi(min):
    # return the dictionary
    return      {'s'    : min.s, 
                 'mg'   : round(min.mg, 2), 
                 'i_24' : min.i_24, 
                 'i_25' : min.i_25, 
                 'i_26' : min.i_26, 
                 'rv'   : min.rv}

def find_minimum_neighbour(raw_spec_filename, raw_spectra, wavelength_region, region, guess, chi_df):

    # generate neighbours close to the guess (that havent already been run)
    guess_arr = generate_neighbours(guess, region)
    print('The length of the guess array before filtering is: ', len(guess_arr))
    guess_arr = filter_guesses(guess_arr, chi_df)
    print('The length of the guess array after filtering is: ', len(guess_arr))
    
    # if there are no neighbours that we can run, this must be our minimum
    if len(guess_arr) == 0:
        return chi_df

    # run optimise_model_fit on the neighbours
    for par in guess_arr:
        # add the new chi squared values to the df
        chi_of_model = optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, par)
        chi_df = chi_df.append(chi_of_model)
    
    # return chi_df with the results of the new models
    return chi_df

def model_finder():
    
    data_path = '/Users/maddie/Desktop/M22/moog_m22/scomb2 3 L/'
    region = 2
    os.chdir(data_path+'r'+str(region))
    os.system('mkdir plots')

    # initial guesses as a dictionary
    guess = initial_guess()

    raw_spec_filename = 'scomb2.M22_III-3.L.txt.mgh'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    wavelength_region = get_wavelength_region(raw_spectra.wavelength)


    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess)

    best_guess_chi = chi_df.chi_squared.iloc[0] # should only be 1 thing in the df atm

    # currently this is the best guess we have
    best_guess = guess
    while len(chi_df) < 300:
        # add the neighbours to the dataframe
        chi_df = find_minimum_neighbour(raw_spec_filename, raw_spectra, 
                        wavelength_region, region, best_guess, chi_df)
        
        # get the best chi-squared fit
        chi_df = chi_df.sort_values(by = ['chi_squared'])
        minimum_model = chi_df.iloc[0]
        print(chi_df)
        print('minimum model identified: ', minimum_model)

        new_best_guess_chi = minimum_model.chi_squared
        new_best_guess = reconstruct_min_chi(minimum_model)
        print('current best guess: ', best_guess_chi)
        print('new best guess: ', new_best_guess_chi)
        
        if best_guess_chi > new_best_guess_chi:
            best_guess = new_best_guess
            best_guess_chi = new_best_guess_chi
            print('now using ', best_guess)
        elif best_guess_chi == new_best_guess_chi:
            print('could not find a better option - exiting')
            break

    return chi_df
        
def calc_ratio(i_24, i_25, i_26):
    i24_percentage=1/(0.01*i_24)
    i25_percentage=1/(0.01*i_25)
    i26_percentage=1/(0.01*i_26)

    isotope_sum = i24_percentage + i25_percentage + i26_percentage

    i24_ratio = (i24_percentage/isotope_sum) * 100
    i25_ratio = (i25_percentage/isotope_sum) * 100
    i26_ratio = (i26_percentage/isotope_sum) * 100
    
    return str(round(i24_ratio,2)) + '_' + str(round(i25_ratio,2)) + '_' + str(round(i26_ratio,2))


def calc_moog(r_24, r_25, r_26):
    i24=1/(0.01*r_24)
    i25=1/(0.01*r_24)
    i26=1/(0.01*r_24)
    return [i24, i25, i26]

def calc_moog_string(r_24, r_25, r_26):
    i24=1/(0.01*r_24)
    i25=1/(0.01*r_24)
    i26=1/(0.01*r_24)
    return str(round(i24,2)) + '_' + str(round(i25,2)) + '_' + str(round(i26,2))

def initial_guess():

    s = 8.41
    mg = 0.55
    i_24 = 2
    i_25 = 15
    i_26 = 13
    rv = 1.8

    # return the guess as a dictionary
    return {'s'    : s, 
            'mg'   : mg, 
            'i_24' : i_24, 
            'i_25' : i_25, 
            'i_26' : i_26, 
            'rv'   : rv}

csv_out = model_finder()
print(csv_out)

csv_out.to_csv('all_fits.csv')