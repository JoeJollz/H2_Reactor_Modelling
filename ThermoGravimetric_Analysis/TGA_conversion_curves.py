"New TGA conversion curve solver"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from sklearn.metrics import r2_score
from matplotlib.ticker import MaxNLocator

'''
Step 1 for this code. 
The peaks in the TGA can be removed, ensure the variable 
'max_peaks' is less than the value of the anomoly peaks which you would like
removed.

target convserion upper value must be selected based on the fact that 2 iso-mass times can be 
calculated for 500 deg C and 600 deg C. To ensure 2 data points exist, check the 'times_for_target_conversions'
and the 500 and 600 keys, there should be nuerical values for all 100 elements in each lists, with no 'None' present.
If so, deacrease the target conversions. 
'''


max_peaks = 90
index_start = 0  # no redox curves prior to this index.
index_end = 89700//3  # no redox curves beyond this point
_starting_mass = 67  # 77.9499#84.635
# upper target conversions must be reasonably picked, 0.8*conversion_max_600degC
#target_conversions = np.linspace(0.0001, 0.3/100, 100)





'''
Step 2.
Correctly select your data file paths. 
file_path_1 corresponds to the TGA exported data.
file_path_2 corresponds to the gas cycles used and their periods. This is the same
file as the .txt file for the gas inlet controller programme.
'''
file_path_1 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\FeO_300_400.txt'
file_path_2 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\JRJ_valve.txt'
file_path_3 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\KOH_LFO_pellet_redox_300_400_500_600_1_g.txt'


# Read the file with the appropriate encoding and skip initial rows
df = df_ = pd.read_csv(
    file_path_1,
    delimiter='\s+',  # Assuming whitespace delimiter
    skiprows=12,      # Skip the initial rows that are not part of the table
    # Manually specifying the column names
    names=["t [s]", "Ts [°C]", "Tr [°C]", "Value [mg]"],
    encoding='ISO-8859-1'  # Specify the encoding (change if needed)
)
df_2 = pd.read_csv(
    file_path_2,
    delimiter=',',  # Assuming whitespace delimiter
    names=['index', 'gas']
)

sum_timer = 0
air_start = []
h2_start = []
co2_start = []
air_end = []
h2_end = []
co2_end = []
for i in range(0, len(df_2)):
    if df_2.iloc[i, 1] == 1:  # air gas injection
        air_start.append(sum_timer)
    if df_2.iloc[i, 1] == 2:  # hydrogen gas injection
        h2_start.append(sum_timer)
    if df_2.iloc[i, 1] == 4:
        co2_start.append(sum_timer)

    sum_timer += df_2.iloc[i, 0]
    if df_2.iloc[i, 1] == 1:  # air gas injection ended
        air_end.append(sum_timer)
    if df_2.iloc[i, 1] == 2:  # hydrogen gas injection ended
        h2_end.append(sum_timer)
    if df_2.iloc[i, 1] == 4:  # co2 gas injection ended
        co2_end.append(sum_timer)


def is_within_ranges(time, start_times, end_times):
    for start, end in zip(start_times, end_times):
        #print(start, end)
        if start <= time <= end:
            #print('found air at t: ', time)
            return True
    #print('bounyancy check accept at t: ', time)
    return False



## ODE solved kinetic models ###
def first_order_kinetics(t, k):
    return 1 - np.exp(-k * t)

def second_order_kinetics(t, k):
    return (k  * t) / ( k * t+1)

def n_th_order_kinetics(t, k): # independent variable is t, parameters to fit: t and n
    n=2
    return 1-(1/((-n+1)*(k*t+1/(n-1))))**(1/(n-1))


def find_time_for_conversion(conversion_data, time_data, target, line):
    for i, conversion in enumerate(conversion_data):
        if conversion >= target:
            return time_data[i]

    return None
    # return None  # Return None if the target conversion is not reached


def get_times_for_target_conversions(conversions, time, target_conversions):
    times_for_conversions = {temp: [] for temp in conversions}

    # conversions x axis, time y axis, extrapolate to find the time taken to reach a target conversion which is not reached within a 2400 second range.

    #interp_function = interp1d(data_to_plot, time, fill_value='extrapolate')

    for temp, conversion_data in conversions.items():
        interp_function = interp1d(
            conversion_data, time, fill_value='extrapolate', bounds_error=False)

        # The use of 1 signifies a linear fit.
        fit = np.polyfit(conversion_data, time, 1)

        # fit
        # [  5.00000000e+00   1.58882186e-15]  #y = 5x + 0

        line = np.poly1d(fit)
        for target in target_conversions:
            time_for_conversion = find_time_for_conversion(
                conversion_data, time, target, line)
            times_for_conversions[temp].append(time_for_conversion)

    return times_for_conversions


# Ensure all data is numeric and handle missing values if any
df = df.apply(pd.to_numeric, errors='coerce').dropna()
df_ = df.apply(pd.to_numeric, errors='coerce').dropna()
df = df.iloc[:int(len(df) / 2), :]


divisor = df["t [s]"][1] - df["t [s]"][0]
plt.plot(np.arange(0, len(df_)//2), df_.iloc[:int(len(df_)/2), int(divisor)])
plt.show()
#fig1, ax1 = plt.subplots()

cooldown = coolup = c = 0
# peak remover ## and conversion curve plotter. Here we identify and plot the reduction and oxidation curves.

## OXIDATION CONVERSION CURVE PLOTS ##
for i in range(0, len(df)-1):
    #print(df['Value [mg]'][i]/df['Value [mg]'][i-1]/1.008>1)
    if (df['Value [mg]'][i]/(df['Value [mg]'][i-1]) > 1.006 and
                i > 100 and
                cooldown == 0
            ):  # so i is the location of the drastic gas change. Gas switch located!

        df['Value [mg]'][i+1] = df['Value [mg]'][i] = df['Value [mg]'][i-2]
        df['Value [mg]'][i+2] = df['Value [mg]'][i] = df['Value [mg]'][i-2]
        df['Value [mg]'][i+3] = df['Value [mg]'][i] = df['Value [mg]'][i-2]
        cooldown = 1
        # check if reduction curve, or oxidation curve!
      #  print(df['t [s]'][i], df['Value [mg]'][i], df['Value [mg]'][i+50])
        if (df['Value [mg]'][i] < df['Value [mg]'][i+5] and
            i > index_start and
            i < index_end and
            not is_within_ranges(df['t [s]'][i+1], air_start, air_end) and
                not is_within_ranges(df['t [s]'][i+1], h2_start, h2_end)):  # oxidation taking place

            data_to_plot = df['Value [mg]'][i:(i+(2370//3))]
            time = np.array(df['t [s]'][i:i+(2370//3)])-df['t [s]'][i]
            t = df['t [s]'][i]
            c += 1

            # data_to_plot = (data_to_plot/_starting_mass-1)*100+100 # plot method 1

            data_to_plot = data_to_plot/data_to_plot[0]  # plot method 2
            data_to_plot = np.array((1-data_to_plot)*100)
            print(data_to_plot)

            # #data_to_plot = (data_to_plot-_starting_mass)/_starting_mass*100+100
            # #ax.plot(time, data_to_plot/_starting_mass*100, label=f'Curve {c}')
            plt.plot(time, -data_to_plot, label=f'Curve {c}, time {t}')

            if c == 4:
                temp = df['Tr [°C]'][i]
                plt.title(f'Oxidation curves at {temp} °C')
                plt.xlabel('CO$_2$ exposure time (seconds)')
                plt.ylabel('Relative Mass Gain (%)')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                #plt.ylim(94, 101.1)  # Set y-axis limits here
                # ax.set_xlim(0,500)
                # Adjust the right side of the plot
                plt.tight_layout(rect=[0, 0, 0.75, 1])
                plt.show()
                d = data_to_plot
                # fig1, ax1 = plt.subplots()
                # if temp == 500:
                #    targ = data_to_plot[-1]/100
                #    target_conversions = np.linspace(0.0001, targ/100, 100)
                c = 0

        if (is_within_ranges(df['t [s]'][i+1], h2_start, h2_end)
            ): # reduction taking place
            df['Value [mg]'][i] = 2000
            t = df['t [s]'][i]
            #print(df['Value [mg]'][i], f'time: {t}')

    cooldown -= 0.1
    cooldown = max(cooldown, 0)


## REDUCTION CONVERSION CURVE PLOTS ##
c = 0
cooldown = 0

# storing the conversions for the final reduction curve, for each isotherm.
conversions = {}
rate_constants_first_order = {}
rate_constants_second_order = {}
rate_constants_nth_order = {}
r_squared_first_order = {}
r_squared_second_order = {}
r_squared_nth_order = {}


for i in range(0, len(df)-1):
    #if df['Value [mg]'][i] > 1500:
        #print('bingo')

    if (df['Value [mg]'][i]/(df['Value [mg]'][i-1]) > 1.01 and
                i > 50
            ):  # so i is the location of the drastic gas change. Gas switch located!
        t = df['t [s]'][i]
        #print(f'accssed, t : {t}')
        df['Value [mg]'][i] = df['Value [mg]'][i-1]
        cooldown = 1

        data_to_plot = df['Value [mg]'][i:(i+(2370//3))]
        
        d_t_p = (data_to_plot-data_to_plot[-1])/(data_to_plot[0]-data_to_plot[-1])
        
        time = np.array(df['t [s]'][i:i+(2370//3)])-df['t [s]'][i]
        t = df['t [s]'][i]
        c += 1
        # data_to_plot = (data_to_plot/_starting_mass-1)*100+100 # plot method
        data_to_plot = data_to_plot/data_to_plot[0]  # plot method 2
        data_to_plot = (1-data_to_plot)*100
        if i < 1050:
            save = data_to_plot
            plt.plot(np.array(data_to_plot))
            plt.title(f'{t}')
            plt.show()
        data_to_plot = np.array(data_to_plot)
        plt.plot(time, data_to_plot, label=f'Cycle {c}')
        #print('added to plot')
        #d_t_p = (data_to_plot)/()
        if c == 4:
            temp = int(df['Tr [°C]'][i])

            title_part1 = r'1LaFeO$_{3-δ}$:9Fe$_2$O$_3$'
            title_part2 = f'Reduction Curves {temp}°C'

            # Combine title parts and set the title
            plt.title(title_part1 + '\n' + title_part2, fontsize=12, pad=15)

            # Increase the font size for x-axis label
            plt.xlabel('H$_2$ Injection Time (Seconds)', fontsize=12)
            # Increase the font size for y-axis label
            plt.ylabel('Relative Mass Change (%)', fontsize=12)
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
            plt.ylim(0, 5)  # Set y-axis limits here
            # ax.set_xlim(0,500)
            # Adjust the right side of the plot
            plt.tight_layout(rect=[0, 0, 0.75, 1])
            plt.show()
            #fig2, ax2 = plt.subplots()
            if temp == 500: # using the 500 deg C data, this will determine the upper value of the isomass kinetics process.
               targ = data_to_plot[-1]
               target_conversions = np.linspace(0.0001, targ/100, 100) # here are the target isoconversion, with the upper limit being bounded by the 500 deg C data, not 600 deg C data.

            data_to_plot = np.append(data_to_plot, data_to_plot[-1])

            conversions[int(temp)] = data_to_plot/100
            #conversions[int(temp)] = 1-d_t_p

            c = 0

    cooldown -= 0.1
    cooldown = max(cooldown, 0)

# key info for kinetic modelling.
time_steps = len(conversions[300])
time = np.linspace(0, time_steps*3, time_steps+1)
time = time[:-1]
temperatures = list(conversions.keys())


fig, ax1 = plt.subplots()

# Plot the Value [mg] data
ax1.plot(df["t [s]"], df["Value [mg]"], color='b', linestyle='-')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Sample Mass [mg]', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a second y-axis
ax2 = ax1.twinx()

# Plot the Tr [°C] data
ax2.plot(df["t [s]"], df["Tr [°C]"], color='r', linestyle='-')
ax2.set_ylabel('Cell Temperature [°C]', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Title and grid
plt.title('Fe2O3 Pellet redox (Absolute) - 08/07/2024')


# Add shaded regions to the plot
'''
Step 3.
On the txt file used to control the gas feed system for the TGA. The file format follows,
exposure_time,number corresponding to gas
The number corresponding to gas must be recognised as a gas stream name. e.g. 
if gas ==1:
    this corresponds to a region of 'Air'
'''

start = 0
for i in range(len(df_2)):
    index = df_2['index'][i]
    end = start+index
    gas = df_2['gas'][i]

    if gas == 1:
        label, color = 'Air', 'white'
    elif gas == 2:
        label, color = 'H2 (5%)', 'green'
    else:
        if i != len(df_2)-1:

            label, color = 'CO2', 'grey'

    end = min(int(end/divisor)-1, len(df)-1)
    start_ = min(int(start/divisor), len(df)-1)
    ax1.axvspan(df["t [s]"].iloc[start_], df["t [s]"].iloc[end],
                color=color, alpha=0.3, label=label)
    start += index
# Add legend for the shaded regions
handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(),
           loc='upper left', bbox_to_anchor=(1.11, 1.05))

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show plot
plt.show()


# relative mass change plot
df['Value [mg]'] = df['Value [mg]']/df['Value [mg]'][0]
df['Value [mg]'][1] = df['Value [mg]'][2] = 1

fig, ax1 = plt.subplots()

# Plot the Value [mg] data
ax1.plot(df["t [s]"], df["Value [mg]"]*100, color='b', linestyle='-')
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Relative Mass Change (%)', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='b', labelsize=12)
ax1.tick_params(axis='x',  labelsize=12)
#ax1.set_ylim(97.3, 100.50)

# Create a second y-axis
ax2 = ax1.twinx()

# Plot the Tr [°C] data
ax2.plot(df["t [s]"], df["Tr [°C]"], color='r', linestyle='-')
ax2.set_ylabel('Cell Temperature [°C]', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r', labelsize=12)

# Title and grid
plt.title('10wt% KOH/LaFeO$_{3-δ}$ Pellet TGA (Relative Mass (%))')

# Add shaded regions to the plot
start = 0
for i in range(len(df_2)):
    index = df_2['index'][i]
    end = start + index
    gas = df_2['gas'][i]

    if gas == 1:
        label, color = 'Air', 'white'
    elif gas == 2:
        label, color = 'H$_2$ (5%)', 'green'
    else:
        if i != len(df_2) - 1:
            label, color = 'CO$_2$', 'grey'
    end = min(int(end / divisor) - 1, len(df) - 1)
    start_ = min(int(start / divisor), len(df) - 1)

    ax1.axvspan(df["t [s]"].iloc[start_], df["t [s]"].iloc[end],
                color=color, alpha=0.3, label=label)
    start += index

handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='upper left',
           bbox_to_anchor=(1.01, 0.1), fontsize=12)

plt.tight_layout(rect=[0, 0, 0.9, 0.95])


# manip 

# conversions[500][:4] = np.linspace(0, conversions[500][4],4)
# conversions[600][:4] = np.linspace(0, conversions[600][3],4)

times_for_target_conversions = get_times_for_target_conversions(
    conversions, time, target_conversions)


plt.show()


## FITTING THE VARIOUS KINETIC MODELS ###
for temp in conversions:
    conversion = conversions[temp]
    
    ### FIRST ORDER KINETICS MODEL ###
    popt_first, _ = curve_fit(first_order_kinetics, time, conversion)
    rate_constants_first_order[temp] = popt_first[0]
    
    
    ### SECOND ORDER KINETICS MODEL ###
    popt_second, _ = curve_fit(second_order_kinetics, time, conversion)
    rate_constants_second_order[temp] = popt_second[0]
    
    ### nTH ORDER KINETICS MODEL ###
    popt_nth, _ = curve_fit(n_th_order_kinetics, time, conversion)
    rate_constants_nth_order[temp] = popt_nth[0]
    
    #rate_constants_first_order[int(600)] = 5e-6
    # Calculate predicted values
    predicted_first_order = first_order_kinetics(
        time, rate_constants_first_order[temp])
    predicted_second_order = second_order_kinetics(
        time, rate_constants_second_order[temp])

    # Calculate R² for first-order kinetics
    ss_res_first = np.sum((conversion - predicted_first_order) ** 2)
    ss_tot_first = np.sum((conversion - np.mean(conversion)) ** 2)
    r_squared_first_order[temp] = 1 - (ss_res_first / ss_tot_first)

    # Calculate R² for second-order kinetics
    ss_res_second = np.sum((conversion - predicted_second_order) ** 2)
    ss_tot_second = np.sum((conversion - np.mean(conversion)) ** 2)
    r_squared_second_order[temp] = 1 - (ss_res_second / ss_tot_second)

#Print rate constants and R² values
print("First-Order Rate Constants:", rate_constants_first_order)
print("First-Order R² Values:", r_squared_first_order)
print("Second-Order Rate Constants:", rate_constants_second_order)
print("Second-Order R² Values:", r_squared_second_order)
print("nth order rate constants:", rate_constants_nth_order)

# # Print rate constants
print("First-Order Rate Constants:", rate_constants_first_order)
print("Second-Order Rate Constants:", rate_constants_second_order)

# Plot the fits for each model
for temp in temperatures:
    if int(temp ) == 600:
        conversion = conversions[temp]
        plt.plot(time, conversion, 'o', label=f'Experimental Data {temp}°C')
        # plt.plot(time, n_th_order_kinetics(
        #     time, rate_constants_nth_order[temp]), '-', label=f'First-Order Fit {temp}°C')
        plt.plot(time, second_order_kinetics(
            time, rate_constants_second_order[temp]), '--', label=f'Second-Order Fit {temp}°C')
plt.xlabel('Time (s)')
plt.ylabel('Conversion')
plt.legend()
plt.show()



### FIRST ORDER KINETICS TO CALC THE ACTIVATION ENERGY AND PRE-EXP FACTOR ###
temperatures = np.array(
    list(rate_constants_first_order.keys())) + 273.15  # Convert to Kelvin
k_values = np.array(list(rate_constants_first_order.values()))

# Take the natural logarithm of the rate constants
ln_k = np.log(k_values)

# Calculate the reciprocal of the temperatures
inv_T = 1 / temperatures

# Perform a linear fit to the Arrhenius plot
slope, intercept = np.polyfit(inv_T, ln_k, 1)

# Calculate activation energy and pre-exponential factor
R = 8.314  # J/(mol*K), universal gas constant
E_a = -slope * R
A = np.exp(intercept)

print(f"Activation Energy (E_a): {E_a} J/mol")
print(f"Pre-exponential Factor (A): {A}")

ln_k_pred = slope * inv_T + intercept

# Calculate R-squared value
r_squared = r2_score(ln_k, ln_k_pred)
#print(f"R-squared: {r_squared}")

# Plot the Arrhenius plot
plt.plot(inv_T, ln_k, 'o', label='Data')
plt.plot(inv_T, slope * inv_T + intercept, '-', label='Fit')
plt.xlabel('1/T (1/K)')
plt.ylabel('ln(k)')
plt.legend()
plt.show()





Act_energies = [] # list to store the activation energies for each relative mass loss.
for i in range(len(target_conversions)):
    temp_store_times = []

    if (times_for_target_conversions[300][i] is None and 
          times_for_target_conversions[400][i] is not None):
        #temps = np.array([1/(400+273.15),1/(500+273.15), 1/(600+273.15)])
        #temp_store_times.append(times_for_target_conversions[400][i])
        temps = np.array([1/(500+273.15), 1/(600+273.15)])
        temp_store_times.append(times_for_target_conversions[500][i])
        temp_store_times.append(times_for_target_conversions[600][i])
    else:
        temps = np.array([1/(500+273.15), 1/(600+273.15)])
        temp_store_times.append(times_for_target_conversions[500][i])
        temp_store_times.append(times_for_target_conversions[600][i])
        
    temp_store_times = np.log(temp_store_times)
    
    slope, y_intercept = np.polyfit(temps, temp_store_times, 1)

    predicted_values = slope * temps + y_intercept

    # Compute residuals
    residuals = temp_store_times - predicted_values

    ss_total = np.sum((temp_store_times - np.mean(temp_store_times))**2)
    ss_residual = np.sum(residuals**2)
    r_squared = 1 - (ss_residual / ss_total)
    E = slope * 8.3144/1000 # calculating the activation energy from the slope.
    Act_energies.append(E)
    
    if i ==97:
        print('mr: ', target_conversions[i]*100)
        print('slope: ', slope)
        print('1/Temps: ', temps)
        print('times: ', np.exp(temp_store_times))
        print('logged times: ', temp_store_times)
        plt.plot(temps, temp_store_times)
        plt.xlabel('1/T')
        plt.ylabel('ln(ta)')
        plt.title('KOH/LFO')
        plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=7))  # Max 5 ticks on x-axis
        plt.show()
    
    
# Act_energies[0]=0
# Act_energies[1] = (Act_energies[2])/2

#Act_energies[:3] = np.linspace(0, Act_energies[3], 3)

KLFconversions = target_conversions
Act_energies_KLF = Act_energies

plt.plot(target_conversions*100, Act_energies)
plt.xlabel('Relative Mass Loss (%)', fontsize=12)
plt.ylabel('Activation Energy (kJ/mol)', fontsize=12)
plt.title('Activation Energy (kJ/mol) vs Relative Mass Loss - LaFeO$_{3-δ}$', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(0,140)
plt.xlim(0, 0.85)
plt.show()

## manual plot to compare the 2 activation energies of 2 correponding oxygen carrier results ##
# plt.plot(LFOconversions*100, Act_energies_LFO, label='LaFeO$_{3-δ}$')
# plt.plot(KOHLFOconversions*100, Act_energies_KOHLFO, label= '10wt% KOH/LaFeO$_{3-δ}$')
# plt.xlabel('Relative Mass Loss (%)', fontsize =12)
# plt.ylabel('Activation Energy (kJ/mol)', fontsize = 12)
# plt.title('Activation Energy vs Relative Mass Loss (%)', fontsize=12)
# plt.tick_params(axis='both', which='major', labelsize=12)
# plt.legend(fontsize=12)
# plt.show()




###############################################################################
###############################################################################
###############################################################################
#                            OPTION B                                         #


#file_path_3 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\KOH_LFO_pellet_redox_300_400_500_600_1_g.txt'

df3 = df_3 = pd.read_csv(
    file_path_3,
    delimiter='\s+',  # Assuming whitespace delimiter
    skiprows=12,      # Skip the initial rows that are not part of the table
    # Manually specifying the column names
    names=["t [s]", "Ts [°C]", "Tr [°C]", "Value [mg]"],
    encoding='ISO-8859-1'  # Specify the encoding (change if needed)
)

df3 = df3.apply(pd.to_numeric, errors='coerce').dropna()
df3 = df3.iloc[:int(len(df3) / 2), :]

cooldown = c =0
conversions_b = {}


for i in range(0, len(df3)-1):
    #print(df['Value [mg]'][i]/df['Value [mg]'][i-1]/1.008>1)
    if (df3['Value [mg]'][i]/(df3['Value [mg]'][i-1]) > 1.006 and
                i > 100 and
                cooldown == 0
            ):  # so i is the location of the drastic gas change. Gas switch located!

        df3['Value [mg]'][i+1] = df3['Value [mg]'][i] = df3['Value [mg]'][i-2]
        df3['Value [mg]'][i+2] = df3['Value [mg]'][i] = df3['Value [mg]'][i-2]
        df3['Value [mg]'][i+3] = df3['Value [mg]'][i] = df3['Value [mg]'][i-2]
        cooldown = 1

        if (is_within_ranges(df3['t [s]'][i+1], h2_start, h2_end)
            ): # reduction taking place
            df3['Value [mg]'][i] = 2000
            t = df3['t [s]'][i]
            #print(df['Value [mg]'][i], f'time: {t}')

    cooldown -= 0.1
    cooldown = max(cooldown, 0)
    
    
for i in range(0, len(df3)-1):
    #if df['Value [mg]'][i] > 1500:
        #print('bingo')

    if (df3['Value [mg]'][i]/(df3['Value [mg]'][i-1]) > 1.01 and
                i > 50
            ):  # so i is the location of the drastic gas change. Gas switch located!
        t = df3['t [s]'][i]
        #print(f'accssed, t : {t}')
        df3['Value [mg]'][i] = df3['Value [mg]'][i-1]
        cooldown = 1

        data_to_plot = df3['Value [mg]'][i:(i+(2370//3))]
        
        d_t_p = (data_to_plot-data_to_plot[-1])/(data_to_plot[0]-data_to_plot[-1])
        
        time = np.array(df3['t [s]'][i:i+(2370//3)])-df3['t [s]'][i]
        t = df3['t [s]'][i]
        c += 1
        # data_to_plot = (data_to_plot/_starting_mass-1)*100+100 # plot method
        data_to_plot = data_to_plot/data_to_plot[0]  # plot method 2
        data_to_plot = (1-data_to_plot)*100
        # if i < 1050:
        #     save = data_to_plot
        #     plt.plot(np.array(data_to_plot))
        #     plt.title(f'{t}')
        #     plt.show()
        data_to_plot = np.array(data_to_plot)
        #plt.plot(time, data_to_plot, label=f'Cycle {c}')
        #print('added to plot')
        #d_t_p = (data_to_plot)/()
        if c == 4:
            temp = int(df3['Tr [°C]'][i])

            # title_part1 = r'1LaFeO$_{3-δ}$:9Fe$_2$O$_3$'
            # title_part2 = f'Reduction Curves {temp}°C'

            # # Combine title parts and set the title
            # plt.title(title_part1 + '\n' + title_part2, fontsize=12, pad=15)

            # # Increase the font size for x-axis label
            # plt.xlabel('H$_2$ Injection Time (Seconds)', fontsize=12)
            # # Increase the font size for y-axis label
            # plt.ylabel('Relative Mass Change (%)', fontsize=12)
            # plt.tick_params(axis='both', which='major', labelsize=12)
            # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
            # plt.ylim(0, 5)  # Set y-axis limits here
            # # ax.set_xlim(0,500)
            # # Adjust the right side of the plot
            # plt.tight_layout(rect=[0, 0, 0.75, 1])
            # plt.show()
            #fig2, ax2 = plt.subplots()
            if temp == 500:
               targ = data_to_plot[-1]
               target_conversions_b = np.linspace(0.0001, targ/100, 100)

            data_to_plot = np.append(data_to_plot, data_to_plot[-1])

            conversions_b[int(temp)] = data_to_plot/100
            
            #conversions[int(temp)] = 1-d_t_p

            c = 0

    cooldown -= 0.1
    cooldown = max(cooldown, 0)

if target_conversions_b[-1]>target_conversions[-1]:
    target_conversions = target_conversions
else:
    target_conversions = target_conversions_b
time = np.linspace(0, time_steps*3, time_steps)

times_for_target_conversions_b = get_times_for_target_conversions(conversions_b, time, target_conversions)
time_for_target_conversions = get_times_for_target_conversions(conversions, time, target_conversions)

del_act_E = []
for i in range(len(target_conversions)):
    temp_store_times = []

    t_a_500 = conversions[500][i]
    t_b_500 = conversions_b[500][i]
    del_t_500 = t_b_500 - t_a_500
    
    t_a_600 = conversions[600][i]
    t_b_600 = conversions_b[600][i]
    del_t_600 = t_b_600 - t_a_600
       

    temps = np.array([1/(500+273.15), 1/(600+273.15)])
    temp_store_times.append(del_t_500)
    temp_store_times.append(del_t_600)
    #temp_store_times.append(times_for_target_conversions[600][i])
        
    temp_store_times = np.log(temp_store_times)
    
    slope, y_intercept = np.polyfit(temps, temp_store_times, 1)

    predicted_values = slope * temps + y_intercept

    # Compute residuals
    residuals = temp_store_times - predicted_values

    ss_total = np.sum((temp_store_times - np.mean(temp_store_times))**2)
    ss_residual = np.sum(residuals**2)
    r_squared = 1 - (ss_residual / ss_total)
    E = slope * 8.3144/1000 # calculating the activation energy from the slope.
    del_act_E.append(E)
    
    # if i ==97:
    #     print('mr: ', target_conversions[i]*100)
    #     print('slope: ', slope)
    #     print('1/Temps: ', temps)
    #     print('times: ', np.exp(temp_store_times))
    #     print('logged times: ', temp_store_times)
    #     plt.plot(temps, temp_store_times)
    #     plt.xlabel('1/T')
    #     plt.ylabel('ln(ta)')
    #     plt.title('KOH/LFO')
    #     plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=7))  # Max 5 ticks on x-axis
    #     plt.show()
    
    
# Act_energies[0]=0
# Act_energies[1] = (Act_energies[2])/2

#Act_energies[:3] = np.linspace(0, Act_energies[3], 3)

KLFconversions = target_conversions
Act_energies_KLF = Act_energies

plt.plot(target_conversions*100, Act_energies)
plt.xlabel('Relative Mass Loss (%)', fontsize=12)
plt.ylabel('Delta Activation Energy (kJ/mol)', fontsize=12)
plt.title('Change in Activation Energy (kJ/mol) vs Relative Mass Loss - LaFeO$_{3-δ}$', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
#plt.ylim(0,140)
#plt.xlim(0, 0.85)
plt.show()



