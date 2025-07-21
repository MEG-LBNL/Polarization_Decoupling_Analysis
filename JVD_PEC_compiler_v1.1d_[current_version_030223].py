#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:31:26 2020

@author: pcagbo
Faradaic efficiency compiler script for all gas products.
"""
from extractor_v2 import *
import matplotlib.pyplot as mt
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
import scipy.ndimage
import pandas as pd

mt.rcParams.update({'font.size': 20, 'lines.markersize': 5, 'font.weight': 'bold'}) #for publication.

directory0 = os.getcwd() + '/input_data' #use as input for 'direc' variable in 'Analyze'for runs when not compiling multiple datasets
electrode_area = 2.5 #cm-2
Zreal = 57 #impedance value for IR drop calculation
#analyzer for Gamry CA data and SRI GC data
#example function call: Analyze(directory0, 2, 0, 4.85, 345, 11.5, 3, 'separate', 'JVD_2_Argon_ctrl_0V.csv' )
#2 line CA data header, no header for GC data, flowrate = 5sccm, initial injection at 345 sec, GC run time = 11.5 minutes,
#3 GC injections per condition (CA) tested.
def Analyze(direc, headerCA, headerGC, flowrate, t_inject_start, GC_runtime, num_injections, datatype, saveas):
    directory2 = direc + '/*.log'
    #extract potentiostatic data:
    if datatype == 'separate': #if data is read in from separate *.dta files:
        directory1 = direc + '/*.DTA'
        current, filelist = extract_auto(directory1, headerCA, '\t', 2, 4) #import nxm array of CA data.
        front_voltages = extract_auto(directory1, headerCA, '\t', 2, 3)[0][1] #import 1xm array of front-contact voltage data.
        back_voltages = extract_auto(directory1, headerCA, '\t', 2, 6)[0][1] #import 1xm array of back-contact voltage data
    else: #if data is read in as a pre-assembled array from an *.csv file:
        directory1 = direc + '/*.csv'
        current = np.loadtxt(directory1[0], skiprows = 3, delimiter = '\t')[:,4]
        filelist = read_header(directory1, headerCA, ',')
    CA_array = pull_CA(GC_runtime, t_inject_start, current)
    #extract gas chromatographic data:
    filenames = np.sort(glob.glob(directory2)); counter = 0; GC_peaks = np.zeros((num_injections*len(filelist), 5))
    for every in filenames: #!!! Note: MAY NEED TO CHANGE DELIMITERS BETWEEN ',' AND '\t'!!!
        if 'TCD' in every:#store hydrogen peak data:
            GC_peaks[:,0] = extract_auto(every, headerGC, '\t', 0, 1)[0][:,0] #import 2xm array of SRI GC data; col 6 = integral, 7 = peak height.
        elif 'FID' in every: #assume other log file is FID data; store CO, CH4, C2H4 and C2H6 integrals:
            GC_peaks[:,1] = extract_auto(every, headerGC, '\t', 2, 3)[0][:,0] #CO integrals: #old cols: 6,7
            GC_peaks[:,2] = extract_auto(every, headerGC, '\t', 4, 5)[0][:,0] #CH4 integrals: old cols: 13,14
            GC_peaks[:,3] = extract_auto(every, headerGC, '\t', 6, 7)[0][:,0] #C2H4 integrals: old cols: 20,21
            GC_peaks[:,4] = extract_auto(every, headerGC, '\t', 8, 9)[0][:,0] #C2H6 integrals: old cols: 27,28
        else:
            pass
    #scrub out zero vectors from GC_peaks (removing empty lines from SRI GC file):
    if np.zeros((1,2)) in GC_peaks:
        counter = 0
        for every in GC_peaks:
            if np.array_equal(every, np.zeros((1, 2))[0]):
                GC_peaks = np.delete(GC_peaks, counter, 0)
            else:
                counter += 1
    else:
        pass
    #organize data, calculate faradaic efficiencies.
    f=flowrate*(7.44e-7)*96485*(0.000001) #flow rate x n x...etc. factor for H2,CO FE calculation.
    calibrations = calibration(GC_peaks[:,0], GC_peaks[:,1], GC_peaks[:,2], GC_peaks[:,3], GC_peaks[:,4]) #***may need to enter each vector input explicitly.***
    assembled = np.zeros((len(CA_array[:,0])*len(filelist), 18))
    counter = 0; counter2 = 1
    while counter+num_injections <= len(assembled):
        assembled[counter:counter+num_injections,0] = CA_array[:,0] #Tinject
        assembled[counter:counter+num_injections,1] = np.transpose(CA_array)[counter2,:] #Current at time = Tinject
        counter += num_injections; counter2 += 1
    assembled[:,2:7] = GC_peaks #H2 and CO raw integrals
    assembled[:,7] = calibrations[0] #[H2]/ppm
    assembled[:,8] = calibrations[1] #[CO]/ppm
    assembled[:,9] = calibrations[2] #[CH4]/ppm
    assembled[:,10] = calibrations[3] #[C2H4]/ppm
    assembled[:,11] = calibrations[4] #[C2H6]/ppm
    assembled[:,12] = assembled[:,7]*(f*2)/abs(assembled[:,1]) #H2 faradaic efficiency
    assembled[:,13] = assembled[:,8]*(f*2)/abs(assembled[:,1]) #CO faradaic efficiency
    assembled[:,14] = assembled[:,9]*(f*8)/abs(assembled[:,1]) #CH4 faradaic efficiency
    assembled[:,15] = assembled[:,10]*(f*12)/abs(assembled[:,1]) #C2H4 faradaic efficiency
    assembled[:,16] = assembled[:,11]*(f*14)/abs(assembled[:,1]) #C2H6 faradaic efficiency

    #clean out cells with negative F.E. values, and set them to zero:
    row,col = 0,12
    for eachrow in assembled:
        for cell in eachrow[12:17]:
            if cell < 0:
                assembled[row,col] = 0; col += 1
            else:
                col += 1
        row += 1; col = 12
    
    #write final column of total F.E.s, after generating null values.
    assembled[:,17] = assembled[:,12] + assembled[:,13] + assembled[:,14] + assembled[:,15] + assembled[:,16] #Sum faradaic efficiency

    #store in data frame, save data to file
    outputs = ['Injection Time / s', 'Current', 'H2', 'CO', 'CH4', 'C2H4', 'C2H6', \
               '[H2]/ppm', '[CO]/ppm', '[CH4]/ppm', '[C2H4]/ppm', '[C2H6]/ppm', \
               'FE(H2)', 'FE(CO)', 'FE(CH4)', 'FE(C2H4)', 'FE(C2H6)', 'Total FE']
    dataframe = pd.DataFrame(assembled, columns = outputs)  
    dataframe.to_csv(os.getcwd() + '/output/' + saveas + '_Compiled_injection_timepoints+FEs.csv')
    dataframe2 = pd.DataFrame(current, columns = ['Time / s'] + filelist)
    dataframe2.to_csv(os.getcwd() + '/output/' + saveas +'_Compiled_raw_CA_data.csv')
    counter = 1
    #mt.style.use('bmh')
    figure = mt.gcf(); figure.set_size_inches(9.6, 5.4)
    #mt.figure(figsize = (4,2))
    while counter < len(current[0,:]):
        mt.plot(current[:,0]/60, current[:,counter], 'b-', linewidth = 1.5)
        counter += 1
    mt.title(saveas)
    mt.xlabel('Time  [min]'); mt.ylabel('J [mA cm$^-$$^2$]')
    #mt.show()
    mt.savefig(os.getcwd()+ '/output/' + saveas + '_compiled_CAs.png', dpi = 180)
    mt.close()
    print (dataframe); return assembled, front_voltages, back_voltages; print (filelist)

#example function call: pull_CA(11.5, 345, array)
def pull_CA(GC_runtime, t_inject_start, array): #GC_runtime is total time between injections in minutes (typically 13 minutes total)
    counter = t_inject_start; short_array = []
    for every in array[0:,]:
        if counter % (GC_runtime*60) == t_inject_start:
            short_array.append(array[counter-1,0:])
            counter += 1
        else:
            counter += 1
    return np.array(short_array)

def calibration(H2_integral, CO_integral, CH4_integral, C2H4_integral, C2H6_integral):
    H2_ppm = 18.594*H2_integral - 23.39695; CO_ppm = 3.002*CO_integral - 55.65248 #MEG GC calibration
    CH4_ppm = 3.8582*CH4_integral - 23.6510; C2H4_ppm = 1.8499*C2H4_integral - 10.9095; C2H6_ppm = 2.0390*C2H6_integral - 152.9782

    #David's GC calibration, 02.18.2020
    #H2_ppm = 12.2196*H2_integral - 84.5549; CO_ppm = 3.7254*CO_integral - 57.4404 
    #CH4_ppm = 3.8582*CH4_integral - 23.6510; C2H4_ppm = 1.8499*C2H4_integral - 10.9095; C2H6_ppm = 2.0390*C2H6_integral - 152.9782

    return H2_ppm, CO_ppm, CH4_ppm, C2H4_ppm, C2H6_ppm

#for use with multiple input folders
#call ex: compile_arrays(2, 0, 4.8, 345, 11.5, 3, 'separate') - when input data are as individual raw dta files
#call ex: compile_arrays(2, 0, 4.8, 345, 11.5, 3, '') - when input data in each folder come as a pre-compiled csv table.
#call ex: compile_arrays(58, 0, 4.8, 390, 13, 2, 'separate') #for secondary GC settings
#compile_arrays(58, None, 4.8, 390, 13, 2, 'separate') #using new fileread procedure, header in log files = None.
def compile_arrays(headerCA, headerGC, flowrate, t_inject_start, GC_runtime, num_injections, datatype):
    Voc = 0 #keep at 0V for dark device.
    directory = glob.glob(os.getcwd() + '/input_data/*')
    #Vector of Potentials tested in volts; V0 = file0, V1 = file1, etc.
    Vall = np.arange(-3, -5, -0.5) #forward series
    #Vall = np.arange(-4.0, -1.75, 0.25) #reverse series
    #initialize lists
    A1_lst = []; B1_lst = []
    A1_avg_lst = []; B1_avg_lst = []
    A1_std_lst = []; B1_std_lst = []
    for folder in directory:
        expt_name = extract_filename(folder)[1]
        A1sub, front_voltages, back_voltages = Analyze(folder, headerCA, headerGC, flowrate, t_inject_start, GC_runtime, num_injections, datatype, expt_name)
        front_voltage = round(front_voltages[3], 1)
        back_voltage = round(back_voltages[3], 1)
        filenames = glob.glob(folder + '/*.DTA')
        A1 = np.zeros((np.shape(A1sub)[0], 9))
        #V1 = Vall[-len(filenames):] #for runs where PV biases (current too low) were skipped, truncate voltage vector accordingly
        V1 = filenames #use with procedures where voltage is automatically extracted from .DTA files.
        #A1 columns: col0 = voltages, col1 = current, col2 = FE(H2), col3 = FE(CO), col4 = total FE
        counter = 0 #add the appropriate voltages to first column of array A1.
        for every in V1:
            A1[counter:counter+num_injections,0] = front_voltage #Column of applied bias.
            A1[counter:counter+num_injections,8] = back_voltage #floating voltage.
            counter += num_injections
        A1[:,1] = A1sub[:,1]; A1[:,2:8] = A1sub[:,12:]
        
        #Clear A1 matrix of bad injection runs:
        count = 0
        for every in A1:
            #exclude runs for H2 FEs exceeding 100% or under 0% or total FEs exceeding 110% or below 50%:
            if every[7] < 0.7 or every[7] > 1.14:
                A1 = np.delete(A1, count, 0)
            else:
                count += 1
        
        B1 = np.zeros((np.shape(A1)[0], 2))
        
        #NOTE: vector (A1[:,0] - Voc) gives the absolute cell voltage = (reverse bias voltage - Voc); Voc is open circuit voltage of PV.
        B1[:,0] = front_voltage * (A1[:,1]*1000/electrode_area) #col0 = Voltage*Current density
        B1[:,1] = (A1[:,7] - A1[:,2]) / A1[:,2] # col1 = CO2RR/H2 F.E. ratio.
        A1_avg = np.zeros((np.int(len(A1)/num_injections), 9))
        A1_std = np.zeros((np.int(len(A1)/num_injections), 9))
        B1_avg = np.zeros((np.int(len(B1)/num_injections), 2)); B1_std = np.zeros((np.int(len(B1)/num_injections), 2))

        A1_lst.append(A1); B1_lst.append(B1) #Append arrays to a list of arrays.
        #Average A1 array and calculate std deviations:
        counter = 0; counter2 = 0
        for every in A1_avg[:,0]:#for each row in A1_avg:
            A1_avg[counter,0] = A1[counter2,0] #populate voltages (same voltage for n injection points)
            A1_avg[counter,1] = np.average(A1[counter2:counter2+num_injections,1])#add averages of current, H2 FE and CO FE and total FE to A1_avg array
            A1_avg[counter,2] = np.average(A1[counter2:counter2+num_injections,2]) #H2
            A1_avg[counter,3] = np.average(A1[counter2:counter2+num_injections,3]) #CO
            A1_avg[counter,4] = np.average(A1[counter2:counter2+num_injections,4]) #CH4
            A1_avg[counter,5] = np.average(A1[counter2:counter2+num_injections,5]) #C2H4
            A1_avg[counter,6] = np.average(A1[counter2:counter2+num_injections,6]) #C2H6
            A1_avg[counter,7] = np.average(A1[counter2:counter2+num_injections,7]) #Sum
            A1_avg[counter,8] = A1[counter2,8]
            #add std.deviations of current, H2 FE and CO FE and total FE to A1_avg array
            A1_std[counter,0] = A1[counter2,0]
            A1_std[counter,1] = np.std(A1[counter2:counter2+num_injections,1])
            A1_std[counter,2] = np.std(A1[counter2:counter2+num_injections,2])
            A1_std[counter,3] = np.std(A1[counter2:counter2+num_injections,3])
            A1_std[counter,4] = np.std(A1[counter2:counter2+num_injections,4])
            A1_std[counter,5] = np.std(A1[counter2:counter2+num_injections,5])
            A1_std[counter,6] = np.std(A1[counter2:counter2+num_injections,6])
            A1_std[counter,7] = np.std(A1[counter2:counter2+num_injections,7])
            A1_std[counter,8] = A1[counter2,8]
            
            A1_avg[:,1] = np.round(A1_avg[:,1], 5) #format array values to only 2 decimal places
            A1[:,1] = np.round(A1[:,1], 5)
            
            A1_avg_lst.append(A1_avg); A1_std_lst.append(A1_std)
            counter += 1; counter2 += num_injections
        #Average B1 array and calculate std deviations.
        counter = 0; counter2 = 0
        for every in B1_avg[:,0]: #for each row in B1_avg:
            B1_avg[counter,0] = np.average(B1[counter2:counter2+num_injections,0])
            B1_avg[counter,1] = np.average(B1[counter2:counter2+num_injections,1])
            B1_std[counter,0] = np.std(B1[counter2:counter2+num_injections,0])
            B1_std[counter,1] = np.std(B1[counter2:counter2+num_injections,1])
            B1_avg_lst.append(B1_avg); B1_std_lst.append(B1_std)
            counter += 1; counter2 += num_injections
            
        count = 0
        for every in A1_avg:
            #exclude runs for H2 FEs exceeding 100% or under 0% or total FEs exceeding 110% or below 50%:
            if every[7] < 0.7:
                A1_avg = np.delete(A1_avg, count, 0)
                A1_std = np.delete(A1_std, count, 0)
            else:
                count += 1
        
        #Save arrays generated from each folder as dataframes:
        x_axis = (list(A1_avg[:,1]*1000)); x_axis2 = (list(A1[:,1]*1000)) #originally, x_axis2 = (list(A1[:8]) ) for voltage labeling.
        index = pd.Index(x_axis, name = 'Current [mA]')
        index2 = pd.Index(x_axis2, name = 'Current [mA]')
        
        A1_cols = ['Cell Potential [V]', 'Current / A', 'FE(H2)', 'FE(CO)', 'FE(CH4)', 'FE(C2H4)', \
                   'FE(C2H6)', 'FE(Total)', 'Cell Potential [V]']
        B1_cols = ['Power Density, mW/cm2', 'CO2RR/H2 Faradaic Yield ratio']
        df_A1 = pd.DataFrame(A1[:,2:7], columns = A1_cols[2:7], index=index2)
        df_A1.to_csv(os.getcwd() + '/output/' + expt_name + '_A1_raw.csv')
        df_A1_avg = pd.DataFrame(A1_avg[:,2:7], columns = A1_cols[2:7], index=index)
        df_A1_avg.to_csv(os.getcwd() + '/output/' + expt_name + '_A1_avg.csv')
        df_A1_std = pd.DataFrame(A1_std[:,2:7], columns = A1_cols[2:7])
        df_A1_std.to_csv(os.getcwd() + '/output/' + expt_name + '_A1_std.csv')
        df_B1 = pd.DataFrame(B1, columns = B1_cols)
        df_B1.to_csv(os.getcwd() + '/output/' + expt_name + '_B1_raw.csv')
        df_B1_avg = pd.DataFrame(B1_avg, columns = B1_cols)
        df_B1_avg.to_csv(os.getcwd() + '/output/' + expt_name + '_B1_avg.csv')
        df_B1_std = pd.DataFrame(B1_std, columns = B1_cols)
        df_B1_std.to_csv(os.getcwd() + '/output/' + expt_name + '_B1_std.csv')
        
        #Plot averaged injections as bar graphs.
        figure = mt.gcf(); figure.set_size_inches(10.5, 4)
        #mt.style.use('bmh')
        colors = ['#B1B1B1', '#FF0000', '#103596', '#A7FFF8', '#000000' ]

        #df_A1_avg = pd.DataFrame(A1_avg, index=index)
        ax = df_A1_avg.plot(kind='bar', stacked=True, color = colors)
        #ax = df_A1_avg.plot(kind='bar', stacked=True, color = colors, yerr=A1_std[:,2:7], capsize = 4, rot=0)
        ax.set_ylabel('Faradaic Yield')
        mt.legend(title='', bbox_to_anchor=(1.0, 1), loc='upper left')
        mt.ylim(0, 1.1)
        mt.title('Averaged Injections (' + str(num_injections) + ' / illumination level)')
        mt.savefig(os.getcwd()+ '/output/' + expt_name + '_(averaged)volt_dependence_FE.png', dpi = 720, bbox_inches='tight')
        mt.grid(False)
        mt.show()
        mt.close()
        
        #Plot individual injections as bar graphs.
        figure = mt.gcf(); figure.set_size_inches(10.5, 4)
        #mt.style.use('bmh')
        colors = ['#B1B1B1', '#FF0000', '#103596', '#A7FFF8', '#000000']

        #df_A1_avg = pd.DataFrame(A1_avg, index=index)
        ax = df_A1.plot(kind='bar', stacked=True, color = colors)
        #ax = df_A1_avg.plot(kind='bar', stacked=True, color = colors, yerr=A1_std[:,2:7], capsize = 4, rot=0)
        ax.set_ylabel('Faradaic Yield')#; ax.set_xlabel('Cell Voltage / V')
        mt.legend(title='', bbox_to_anchor=(1.0, 1), loc='upper left')
        mt.ylim(0, 1.1)
        #mt.title('Individual Injections (' + str(num_injections) + '/ illumination level)')
        mt.savefig(os.getcwd()+ '/output/' + expt_name + '_volt_dependence_FE.png', dpi = 720, bbox_inches='tight')
        mt.grid(False)
        mt.show()
        mt.close()
        
        size = 0
        for every in A1_lst:
            size += len(every)
        A1_array = np.zeros((size, 9)) #super array contains all values from all powers measured for map generation
        B1_array = np.zeros((size, 2))
        #Add all values from various power tests to A1_array and B1_array:
        counter = 0
        for every in A1_lst:
            A1_array[counter:counter+len(every),:] = every
            counter += len(every)
        counter = 0
        for every in B1_lst:
            B1_array[counter:counter+len(every),:] = every
            counter += len(every)
        #purge super arrays of runs where total faradaic efficiency is over 1.2 or negative:
        """count = 0
        for every in A1_array[:,]:
            #exclude runs for H2 FEs exceeding 100% or under 0% or total FEs exceeding 110% or below 90%:
            if every[2] > 1.1 or every[2] < 0.0 or every[7] > 1.1 or every[7] < 0.8:
                A1_array = np.delete(A1_array, count, 0)
                B1_array = np.delete(B1_array, count, 0)
            else:
                count += 1
        """
        #mt.style.use('bmh')
        figure = mt.gcf(); figure.set_size_inches (10.5, 5)
        ax = mt.subplot(1,2,1) #plot CO/H2 ratios as fn of Power Density
        mt.grid(False)
        #mt.plot(B1_array[:,0], B1_array[:,1], 'k-', linewidth = 0.5)
        mt.plot(B1_array[:,0], B1_array[:,1], 'co')
        #mt.ylim(-0.5, 10)
        #mt.xlim(-0.5, 36)
        mt.xlabel('Power Density, JxV [mW cm$^-$$^2$]')
        mt.ylabel('CO$_2$RR/H$_2$ Faradaic Yield Ratio')
        #ax.set_xticks([0,10,20,30])
    
        ax = mt.subplot(1,2,2) #plot CO/H2 ratios as fn of Voltage
        mt.grid(False)
        #mt.plot(A1_array[:,0], B1_array[:,1], 'k-', linewidth = 0.5)
        mt.plot(A1_array[:,0], B1_array[:,1], 'co')
        #mt.ylim(-0.5, 10)
        #mt.xlim(-5, -1.6)
        mt.xlabel('Applied Potential [V]')
        mt.ylabel('CO$_2$RR/H$_2$ Faradaic Yield Ratio')
        #ax.set_xticks([-1.8,-2.8,-3.8,-4.8])
    
        mt.tight_layout(); mt.savefig(os.getcwd() + '/output/product_ratio.png', dpi=720, bbox_inches = 'tight')
    
        #sort arrays before printing to file
        A1_array = A1_array[A1_array[:,5].argsort()]#sort array ascending monotonically by cell potential vs Ewe
        B1_array = B1_array[B1_array[:,0].argsort()]#sort array monotonically in ascending order by v/j ratio
    
        dataframeA1_array = pd.DataFrame(A1_array, columns = A1_cols)
        dataframeA1_array.to_csv(os.getcwd() + '/output/A1_super_array(all_data).csv')
        dataframeB1_array = pd.DataFrame(B1_array, columns = B1_cols)
        dataframeB1_array.to_csv(os.getcwd() + '/output/B1_super_array(all_data).csv')
    
        print (np.array(filenames))
        print (A1_array)
        print front_voltages, front_voltage