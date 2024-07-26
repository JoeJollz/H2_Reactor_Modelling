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
stoichometric_starting_mass = 84.635

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
air_end = []
for i in range(0, len(df_2)):
    if df_2.iloc[i, 1] == 1: # air gas injection
        air_start.append(sum_timer)
    sum_timer += df_2.iloc[i,0]
    if df_2.iloc[i,1] ==1: #air gas injection ended
        air_end.append(sum_timer)
        
def is_within_ranges(time, start_times, end_times):
    for start, end in zip(start_times, end_times):
        if start <= time <= end:
            print('found air curve at t: ', time)
            return True
    
    return False

# Ensure all data is numeric and handle missing values if any
df = df.apply(pd.to_numeric, errors='coerce').dropna()
df_ = df.apply(pd.to_numeric, errors='coerce').dropna()
df = df.iloc[:int(len(df) / 2), :]


divisor = df["t [s]"][1] - df["t [s]"][0]
plt.plot(np.arange(0,len(df_)//2),df_.iloc[:int(len(df_)/2), int(divisor)])
plt.show()

## peak remover ## and conversion curve plotter. Here we identify and plot the reduction and oxidation curves. 
for i in range(0,len(df)-1):
    if (df['Value [mg]'][i]/(df['Value [mg]'][i-1])>1.008 and i > 400):  # so i is the location of the drastic gas change. Gas switch located! 
        df['Value [mg]'][i+1]= df['Value [mg]'][i]=df['Value [mg]'][i-1]
        df['Value [mg]'][i+2]= df['Value [mg]'][i]=df['Value [mg]'][i-1]
        
        #check if reduction curve, or oxidation curve!
        
        if (df['Value [mg]'][i] < df['Value [mg]'][i+5] and
            i>index_start and 
            i<index_end and 
            not is_within_ranges(df['t [s]'][i+1], air_start, air_end)): # oxidation taking place
            data_to_plot = df['Value [mg]'][i:(i+(2370//3))]
            time = np.array(df['t [s]'][i:i+(2370//3)])-df['t [s]'][i]
            t = df['t [s]'][i]
            plt.plot(time, data_to_plot)
            plt.title(f'Start time {t} seconds')
            plt.show()
        
        
        
        
    if df['Value [mg]'][i]>max_peaks:
        df['Value [mg]'][i] = df['Value [mg]'][i+4]
        df['Value [mg]'][i] = df['Value [mg]'][i+4]

    


#df['Value [mg]'][1]= df['Value [mg]'][2] =1

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

#divisor = df["t [s]"][1] - df["t [s]"][0]


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

# fig, ax1 = plt.subplots()

# # Plot the Value [mg] data
# ax1.plot(df["t [s]"], df["Value [mg]"]*100, color='b', linestyle='-')
# ax1.set_xlabel('Time [s]')
# ax1.set_ylabel('Relative Mass Change (%)', color='b')
# ax1.tick_params(axis='y', labelcolor='b')

# # Create a second y-axis
# ax2 = ax1.twinx()

# # Plot the Tr [°C] data
# ax2.plot(df["t [s]"], df["Tr [°C]"], color='r', linestyle='-')
# ax2.set_ylabel('Cell Temperature [°C]', color='r')
# ax2.tick_params(axis='y', labelcolor='r')

# # Title and grid
# plt.title('Fe2O3 Pellet redox (Relative (%)) - 08/07/2024')

# # Add shaded regions to the plot
# start = 0
# for i in range(len(df_2)):
#     index = df_2['index'][i]
#     end = start+index
#     gas = df_2['gas'][i]
    
#     if gas ==1:
#         label, color = 'Air', 'white'
#     elif gas ==2:
#         label, color = 'H2 (5%)', 'green'
#     else:
#         if i !=len(df_2)-1:
            
#             label, color = 'CO2', 'grey'
#     end = min(int(end/divisor)-1, len(df)-1)
#     start_ = min(int(start/divisor), len(df)-1)
    
#     ax1.axvspan(df["t [s]"].iloc[start_], df["t [s]"].iloc[end], color=color, alpha=0.3, label=label)
#     start += index

# # Add legend for the shaded regions
# handles, labels = ax1.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax1.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.11, 1.05))

# # Adjust layout to make room for the legend
# plt.tight_layout(rect=[0, 0, 0.85, 1])

# # Show plot
# plt.show()

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

# Add legend for the shaded regions
# handles, labels = ax1.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax1.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.11, 1.05))

# #Create a color bar for the phases
# cmap = ListedColormap(['red','yellow'])
# norm = BoundaryNorm([95, 95.3, 103], cmap.N)

# # Create a third axis for the color bar
# cbar_ax = fig.add_axes([0.8, 0.15, 0.01, 0.72])
# cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar_ax)
# cbar.set_ticks([95.9, 103])
# cbar.set_ticklabels(['Magnetite', 'Hematite'])
# cbar.set_label('Phase')
handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.11, 1.05))

# Create an inset axis for the color bar
# inset_ax = inset_axes(ax1,
#                       width="5%",  # width of the inset axes in fraction of the main plot
#                       height="50%",  # height of the inset axes in fraction of the main plot
#                       loc='center right',
#                       bbox_to_anchor=(0.7, -0.45, 1, 1.85),
#                       bbox_transform=ax1.transAxes,
#                       borderpad=0)

# # Create the color bar plot
# inset_ax.set_xlim(0, 1)
# inset_ax.set_ylim(90, 100)

# # Plotting the color bar segments
# inset_ax.fill_between([0, 1], 90, 95.3, color='yellow')
# inset_ax.fill_between([0, 1], 93.3, 100, color='red')

# # Adding labels for the color bar
# inset_ax.text(1.5, 90, 'Wustite + Magnetite', ha='left', va='center', color='black', fontsize=8)
# inset_ax.text(1.5, 93.3, 'Magnetite', ha='left', va='center', color='black', fontsize=8)
# inset_ax.text(1.5, 100, 'Hematite', ha='left', va='center', color='black', fontsize=8)

# Removing x-axis ticks and labels
# inset_ax.set_xticks([])
# inset_ax.set_xticklabels([])

# inset_ax.set_yticks([])
# inset_ax.set_yticklabels([])
# Setting y-axis labels
# inset_ax.set_yticks([95, 95.6, 100])
# inset_ax.set_yticklabels(['95', '95.6', '100'])

# Adjust layout to make room for the legend and color bar
plt.tight_layout(rect=[0, 0, 0.9, 0.95])

# Show plot
plt.show()