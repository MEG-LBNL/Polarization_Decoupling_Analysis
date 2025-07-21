#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import matplotlib.pyplot as mt
import numpy as np; import pandas as pd
import glob; import os; import csv
import scipy.signal as sig
"""
Created on Thu Feb 20 17:13:07 2020

@author: pcagbo
"""
#General data extraction from a group of text files (.csv, etc.)
#General Directions:
#When entering strings, be sure to surround them in quotes.
#Filenames must be entered with the file extension i.e. "file.csv
def extract(directory, header, delim, col0, col1):  
    #initialize arrays and vectors (assumes same file structure and x axis for all data as file 0):
    init_array = fileread(np.sort((glob.glob(directory)))[0], header, delim, col0, col1)[1]
    filelist = []; filenum = len(glob.glob(directory))
    xdata = init_array[:,0]
    All_data = np.zeros((np.shape(init_array)[0], filenum+1)); All_data[:,0] = xdata
    
    responses = ['yes', 'Yes', 'y', 'Y', 'YES']
    base = raw_input("Do you need to baseline adjust? (y/n) ")
    xcolumn = raw_input("What is the name of the x column?")
    if base in responses:
        low = input("Enter lower bound for averaging: ")
        high = input("Enter upper bound for averaging: ")
        low_ind = approx_index(low, xdata)[0][0]
        high_ind = approx_index(high, xdata)[0][0]
        saveto = os.getcwd() + '/Extracted&compiled+Baseline_Adjusted.csv'
    else:
        low = ''; high = ''
        saveto = os.getcwd() + '/Extracted&compiled.csv'
    
    counter = 1
    for filename in np.sort(glob.glob(directory)):
        data = fileread(filename, header, delim, col0, col1)[1]; yarray = data[:,1]
        if type(low) == int or type(low) == float:
            adjusted_array = baseline(yarray, low_ind, high_ind)[1]
            All_data[:,counter] = adjusted_array
            filelist.append(extract_filename(filename)[1])
            counter += 1
        else:
            filelist.append(extract_filename(filename)[1])
            All_data[:,counter] = yarray
            counter += 1
    
    #plot data?
    plotdata = raw_input("Do you want to plot the compiled data?")
    if plotdata in responses:
        counter = 0
        for every in All_data[:,1:]:
            mt.plot(All_data[:,0], every, label = filelist[counter])
            counter += 1
        mt.legend(loc = 'best'); mt.show()
        
    #writing adjusted data to file:
    dataframe = pd.DataFrame(All_data, columns = [xcolumn] + filelist)
    dataframe.to_csv(saveto)
    return All_data, filelist#, init_array, yarray

def extract_auto(directory, header, delim, col0, col1): #no user input options for batch operations   
    #initialize arrays and vectors (assumes same file structure and x axis for all data as file 0):
    init_array = fileread(glob.glob(directory)[0], header, delim, col0, col1)[1]
    filelist = []; filenum = len(glob.glob(directory))
    xdata = init_array[:,0]
    All_data = np.zeros((np.shape(init_array)[0], filenum+1)); All_data[:,0] = xdata
    
    counter = 1
    #for filename in sorted(glob.glob(directory), key=lambda name: int(name[110:-4])): 
    for filename in np.sort(glob.glob(directory)):
        data = fileread(filename, header, delim, col0, col1)[1]; yarray = data[:,1]
        filelist.append(extract_filename(filename)[1])
        #if needed, filter data x number of times; x = (stop - reps)
        reps = 0; stop = 3
        if stop > 0:
            while reps < stop:
                yarray = sig.savgol_filter(yarray, 3, 0, delta=16) #run noisy data through a filter
                reps += 1
            All_data[:,counter] = yarray
        else:
            All_data[:,counter] = yarray
        counter += 1
        
    return All_data, filelist#, init_array, yarray
    
def extract_filename(path):
    counter = -1
    if '/' in path:
        while -counter < len(path):
            if path[counter] == '/':
                dir = path[:counter+1]
                file = path[counter+1:]; filename = file[:-4]; extension = file[-4:]
                return [dir, file, filename, extension]
            else:
                counter -= 1
    else:
        return path
    
"""    
def fileread(filename, header, delim, col0, col1): #header should be '0' for SRI data.
    f = open(filename, 'rt')
    list = []#; colnames = []
    reader = csv.reader(f, delimiter = delim)
    for row in reader:
        list.append(row)
    f.close()
    array = list[header:] #remove file headers/excess data for array formatting
    MAT = []
    counter = 0
    while counter < len(array):
        a = float(array[counter][col0])
        b = float(array[counter][col1])
        MAT.append([a, b])
        counter += 1
    #print array, len(array)
    return None, np.array(MAT), array
"""

def fileread(filename, headers, delim, col0, col1):
    #generate organized dataframe with all products in a unique column; extract from col0 and col1, output as 2xn array.
    if '.log' in filename: #if the file is an SRI .log file:
        f = pd.read_csv(filename, sep = delim, header = headers)
        farray = np.array(f)
        fcleaned = np.zeros((len(f), 10)) #generate empty array of proper dimensions
        colcount  = 0; rowcount = 0
        while rowcount + 1 <= len(farray[:,0]):
            while colcount + 1 <= len(farray[0,:]): #if the column counter is <= the number of columns (length of a row) in farray:
                if farray[rowcount,colcount] == 'H2':
                    fcleaned[rowcount,0:2] = farray[rowcount,colcount+3:colcount+5]
                    colcount += 1
                if farray[rowcount,colcount] == 'CO':
                    fcleaned[rowcount,2:4] = farray[rowcount,colcount+3:colcount+5]
                    colcount += 1
                if farray[rowcount,colcount] == 'CH4':
                    fcleaned[rowcount,4:6] = farray[rowcount,colcount+3:colcount+5]
                    colcount += 1
                if farray[rowcount,colcount] == 'Ethylene':
                    fcleaned[rowcount,6:8] = farray[rowcount,colcount+3:colcount+5]
                    colcount += 1
                if farray[rowcount,colcount] == 'Ethane':
                    fcleaned[rowcount,8:10] = farray[rowcount,colcount+3:colcount+5]
                    colcount += 1
                else:
                    colcount += 1
            rowcount += 1; colcount = 0
        return fcleaned, fcleaned[:,col0:col1+1], pd.DataFrame(fcleaned)
    else: #assume data is .DTA file for potentionstatic data
        f = open(filename, 'rt')
        list = []#; colnames = []
        reader = csv.reader(f, delimiter = delim)
        for row in reader:
            list.append(row)
        f.close()
        array = list[headers:] #remove file headers/excess data for array formatting
        MAT = []
        counter = 0
        while counter < len(array):
            a = float(array[counter][col0])
            b = float(array[counter][col1])
            MAT.append([a, b])
            counter += 1
        #print array, len(array)
        return None, np.array(MAT), array

#return the column headers in a file.
def read_header(filename, header, delim):
    f = open(filename, 'rt')
    colnames = []
    reader = csv.reader(f, delimiter = delim)
    colnames.append(reader[header-1]) #return the last row of header (assumes these are the column names)
    return colnames

#given any value, returns the location (index) of which number in a list of numbers is closest to the value specified.
def approx_index(value, list):
    residual = []
    if value == None or value in list:
        return indexer(value, list)
    else:
        #make vector of the form [X1-value, X2-value,..., Xn-value...]
        #find index of min value of vector (location of the number closest to value)
        for every in list:
            residual.append(abs(value - every))
        return indexer( min(residual), residual)

#gives the location of an item in a list of items.
def indexer(value, list):
    index = 0
    indices = []
    listvalues = []
    if value in list:
        while (index + 1) <= len(list):
            if value == list[index]:
                indices.append(index)
                listvalues.append(list[index])
            index += 1
        return indices, listvalues
    return 'no such value in list'

def baseline(array, low, high): #use to average baseline given some array of values.
    counter = 0
    sum = 0
    adjusted_array = []
    subarray = array[low : high]
    for every in subarray:
        sum += every
        counter += 1
    average = sum / float(counter)
    for every in array:
        value = every - average
        adjusted_array.append(value)
    return average , np.array(adjusted_array)
    
#use as input to the variable 'directory' if current directory contains the files
def currentdir():
    return os.getcwd() + '/input_data/*.*'

def listfiles():
    print 'Files to be analyzed must be in the current working directory.'
    files = np.sort(glob.glob(currentdir()))
    counter = 0
    print 'File number\t', 'Filename'
    for each in files:
        print counter, '\t', each
        counter += 1
    return files