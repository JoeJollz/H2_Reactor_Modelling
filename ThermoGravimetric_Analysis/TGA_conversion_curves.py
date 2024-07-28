"New TGA conversion curve solver" 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

'''
Step 1 for this code. 
The peaks in the TGA can be removed, ensure the variable 
'max_peaks' is less than the value of the anomoly peaks which you would like
removed.
'''


max_peaks = 85.5
index_start = 0 # no redox curves prior to this index.
index_end = 89700//3 #  no redox curves beyond this point
_starting_mass = 84.635#77.9499#84.635

'''
Step 2.
Correctly select your data file paths. 
file_path_1 corresponds to the TGA exported data.
file_path_2 corresponds to the gas cycles used and their periods. This is the same
file as the .txt file for the gas inlet controller programme.
'''
file_path_1 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\10KOH_1LFO_9FE2O3_PEL.txt'
file_path_2 = r'C:\Users\jrjol\OneDrive - University of Cambridge\Documents\Cambridge\Project\TGA DATA\Doped TGA 600 500 400 pellets\JRJ_valve.txt'

# Read the file with the appropriate encoding and skip initial rows
df= df_ = pd.read_csv(
    file_path_1,
    delimiter='\s+',  # Assuming whitespace delimiter
    skiprows=12,      # Skip the initial rows that are not part of the table
    names=["t [s]", "Ts [°C]", "Tr [°C]", "Value [mg]"],  # Manually specifying the column names
    encoding='ISO-8859-1'  # Specify the encoding (change if needed)
)
df_2 = pd.read_csv(
    file_path_2,
    delimiter=',',  # Assuming whitespace delimiter
    names = ['index', 'gas']
    )

sum_timer = 0
air_start = []
h2_start = []
co2_start = []
air_end = []
h2_end = []
co2_end = []
for i in range(0, len(df_2)):
    if df_2.iloc[i, 1] == 1: # air gas injection
        air_start.append(sum_timer)
    if df_2.iloc[i,1] ==2: # hydrogen gas injection
        h2_start.append(sum_timer)
    if df_2.iloc[i,1] == 4:
        co2_start.append(sum_timer)
    
    sum_timer += df_2.iloc[i,0]
    if df_2.iloc[i,1] ==1: #air gas injection ended
        air_end.append(sum_timer)
    if df_2.iloc[i,1] ==2: # hydrogen gas injection ended
        h2_end.append(sum_timer)
    if df_2.iloc[i,1] == 4: # co2 gas injection ended
        co2_end.append(sum_timer)
        
def is_within_ranges(time, start_times, end_times):
    for start, end in zip(start_times, end_times):
        #print(start, end)
        if start <= time <= end:
            #print('found air at t: ', time)
            return True
    #print('bounyancy check accept at t: ', time)
    return False

# Ensure all data is numeric and handle missing values if any
df = df.apply(pd.to_numeric, errors='coerce').dropna()
df_ = df.apply(pd.to_numeric, errors='coerce').dropna()
df = df.iloc[:int(len(df) / 2), :]


divisor = df["t [s]"][1] - df["t [s]"][0]
plt.plot(np.arange(0,len(df_)//2),df_.iloc[:int(len(df_)/2), int(divisor)])
plt.show()
fig1, ax1 = plt.subplots()

cooldown = coolup =c = 0
## peak remover ## and conversion curve plotter. Here we identify and plot the reduction and oxidation curves. 

## OXIDATION CONVERSION CURVE PLOTS ##
for i in range(0,len(df)-1):
    #print(df['Value [mg]'][i]/df['Value [mg]'][i-1]/1.008>1)
    if (df['Value [mg]'][i]/(df['Value [mg]'][i-1])>1.008 and 
        i > 100 and 
        cooldown == 0
        ):  # so i is the location of the drastic gas change. Gas switch located! 
    
        df['Value [mg]'][i+1]= df['Value [mg]'][i]=df['Value [mg]'][i-2]
        df['Value [mg]'][i+2]= df['Value [mg]'][i]=df['Value [mg]'][i-2]
        df['Value [mg]'][i+3]= df['Value [mg]'][i]=df['Value [mg]'][i-2]
        cooldown = 1
        #check if reduction curve, or oxidation curve!
        print( df['t [s]'][i], df['Value [mg]'][i], df['Value [mg]'][i+50])
        if (df['Value [mg]'][i] < df['Value [mg]'][i+5] and
            i>index_start and 
            i<index_end and 
            not is_within_ranges(df['t [s]'][i+1], air_start, air_end) and
            not is_within_ranges(df['t [s]'][i+1], h2_start, h2_end)): # oxidation taking place
            
            data_to_plot = df['Value [mg]'][i:(i+(2370//3))]
            time = np.array(df['t [s]'][i:i+(2370//3)])-df['t [s]'][i]
            t = df['t [s]'][i]
            c +=1

            #data_to_plot = (data_to_plot/_starting_mass-1)*100+100 # plot method 1
            
            
            data_to_plot = data_to_plot/data_to_plot[0]  # plot method 2
            data_to_plot = (1-data_to_plot)*100
            
            # #data_to_plot = (data_to_plot-_starting_mass)/_starting_mass*100+100
            # #ax.plot(time, data_to_plot/_starting_mass*100, label=f'Curve {c}')
            ax1.plot(time, data_to_plot, label=f'Curve {c}, time {t}')
            
            if c == 4: 
                temp = df['Tr [°C]'][i]
                ax1.set_title(f'Oxidation curves at {temp} °C')
                ax1.set_xlabel('CO$_2$ exposure time (seconds)')
                ax1.set_ylabel('Relative Mass Change (w.r.t 100%)')
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax1.set_ylim(94, 101.1)  # Set y-axis limits here
                #ax.set_xlim(0,500)
                fig1.tight_layout(rect=[0, 0, 0.75, 1])  # Adjust the right side of the plot
                plt.show()
                fig1, ax1 = plt.subplots()
                
                c = 0   
        
                
        if (is_within_ranges(df['t [s]'][i+1], h2_start, h2_end)
            ): # reduction taking place
            df['Value [mg]'][i] = 2000
            t = df['t [s]'][i]
            print(df['Value [mg]'][i], f'time: {t}')
            
    cooldown -= 0.1
    cooldown = max(cooldown, 0)



## REDUCTION CONVERSION CURVE PLOTS ##
c=0
cooldown = 0

for i in range(0,len(df)-1):
    if df['Value [mg]'][i]>1500:
        print('bingo')
    
    if (df['Value [mg]'][i]/(df['Value [mg]'][i-1])>1.01 and 
        i >50
        ):  # so i is the location of the drastic gas change. Gas switch located! 
        t = df['t [s]'][i]
        print(f'accssed, t : {t}')
        df['Value [mg]'][i]= df['Value [mg]'][i-1]
        cooldown = 1
        
        data_to_plot = df['Value [mg]'][i:(i+(2370//3))]
        time = np.array(df['t [s]'][i:i+(2370//3)])-df['t [s]'][i]
        t = df['t [s]'][i]
        c +=1
        #data_to_plot = (data_to_plot/_starting_mass-1)*100+100 # plot method 
        data_to_plot = data_to_plot/data_to_plot[0]  # plot method 2
        data_to_plot = (1-data_to_plot)*100
        if i<1050:
            save = data_to_plot
            plt.plot(np.array(data_to_plot))
            plt.title(f'{t}')
            plt.show()
        data_to_plot = np.array(data_to_plot)
        #data_to_plot = (data_to_plot-_starting_mass)/_starting_mass*100+100
        #ax.plot(time, data_to_plot/_starting_mass*100, label=f'Curve {c}')
        plt.plot(time, data_to_plot, label=f'Curve {c}, time {t}')
        print('added to plot')

        if c == 4: 
            temp = df['Tr [°C]'][i]
            plt.title(f'Reduction curves at {temp} °C')
            plt.xlabel('H$_2$ exposure time (seconds)')
            plt.ylabel('Relative Mass Change (w.r.t 100%)')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            #ax.set_ylim(94, 101.1)  # Set y-axis limits here
            #ax.set_xlim(0,500)
            plt.tight_layout(rect=[0, 0, 0.75, 1])  # Adjust the right side of the plot
            plt.show()
            #fig2, ax2 = plt.subplots()
            
            c = 0

    cooldown -= 0.1
    cooldown = max(cooldown, 0)
        
    if df['Value [mg]'][i]>max_peaks:
        df['Value [mg]'][i] = df['Value [mg]'][i+4]
        df['Value [mg]'][i] = df['Value [mg]'][i+4]

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
    
    if gas ==1:
        label, color = 'Air', 'white'
    elif gas ==2:
        label, color = 'H2 (5%)', 'green'
    else:
        if i !=len(df_2)-1:
            
            label, color = 'CO2', 'grey'
    
    end = min(int(end/divisor)-1, len(df)-1)
    start_ = min(int(start/divisor), len(df)-1)
    ax1.axvspan(df["t [s]"].iloc[start_], df["t [s]"].iloc[end], color=color, alpha=0.3, label=label)
    start += index
# Add legend for the shaded regions
handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.11, 1.05))

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show plot
plt.show()


## relative mass change plot
df['Value [mg]'] = df['Value [mg]']/df['Value [mg]'][0]
df['Value [mg]'][1]= df['Value [mg]'][2] =1

fig, ax1 = plt.subplots()

# Plot the Value [mg] data
ax1.plot(df["t [s]"], df["Value [mg]"]*100, color='b', linestyle='-')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Relative Mass Change (%)', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a second y-axis
ax2 = ax1.twinx()

# Plot the Tr [°C] data
ax2.plot(df["t [s]"], df["Tr [°C]"], color='r', linestyle='-')
ax2.set_ylabel('Cell Temperature [°C]', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Title and grid
plt.title('3wt% KOH/LaFeO$_3-$ Pellet TGA (Relative (%))')

# Add shaded regions to the plot
start = 0
for i in range(len(df_2)):
    index = df_2['index'][i]
    end = start + index
    gas = df_2['gas'][i]
    
    if gas == 1:
        label, color = 'Air', 'white'
    elif gas == 2:
        label, color = 'H2 (5%)', 'green'
    else:
        if i != len(df_2) - 1:
            label, color = 'CO2', 'grey'
    end = min(int(end / divisor) - 1, len(df) - 1)
    start_ = min(int(start / divisor), len(df) - 1)
    
    ax1.axvspan(df["t [s]"].iloc[start_], df["t [s]"].iloc[end], color=color, alpha=0.3, label=label)
    start += index

handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.11, 1.05))

plt.tight_layout(rect=[0, 0, 0.9, 0.95])

plt.show()