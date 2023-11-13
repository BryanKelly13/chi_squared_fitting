# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 21:56:15 2023

@author: grpat
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.optimize import curve_fit
import fnmatch


#lab_angles = [15,20,25,30,35,40,45,50]

fresco_output_path = '/home/bkelly/Programs/chi_squared_fitting/fresco_output_files'

def custom_sort_key(file_name):
    order = {}
    for char in 'spdf':
        order[char] = file_name.find(char) if char in file_name else float('inf')
    return [order[char] for char in 'spdf']

def sort_lists(list1, list2):
    return [x for _, x in sorted(zip(list2, list1))]

def ADWA_Fit(x,SF):
    return SF*x

def ADWA_Mixed_Fit(x, SF1, SF2):
    return SF1*x[0] + SF2*x[1]


def plot_data(MIXED_STRENGTH, SF, covariance, lab_angles, x_sec_vals, error_vals, fresco_vals, energy):
    os.chdir(fresco_output_path)
    totalFiles = os.listdir(fresco_output_path)
    fileFound = False
    if MIXED_STRENGTH == False:
        for filename in totalFiles:
            if filename.endswith('.sorted'):
                if energy in filename:
                    fileFound = True
                    file = filename
                    break
        if fileFound:
            data = np.loadtxt(file)

            # Step 2: Separate x and y values from the data
            x = data[:, 0]
            y = data[:, 1]
            scaled_y = y * SF
            std_dev = np.sqrt(covariance[0])
            print(f"Spectroscopic factor: {SF[0]:.5f} with a std deviation: {std_dev[0]:.5f}")

            # plt.plot(x,y, marker='', linestyle='--', color='k', label="DWBA UnScaled")
            plt.plot(x, scaled_y, marker='', linestyle='--', color='darkblue', label='Scaled DWBA',alpha=0.5)
            plt.fill_between(x, scaled_y - std_dev, scaled_y + std_dev, color='dodgerblue', alpha=0.3, label='Std. Dev.') # change to use the error from cov matrix from chisquared
            plt.errorbar(lab_angles, x_sec_vals, error_vals, color='black', marker='x', label='experimental data', ecolor='red', capsize=2, linestyle='None')
            plt.xlabel('CM Angle')
            plt.ylabel('Cross-Section [mb/sr]')
            plt.title(f'DWBA Analysis for {energy} keV state')
            plt.yscale('log')

            plt.minorticks_on()

            plt.legend()
            plt.show()
            
        else:
            print(f"No matching files for {energy} kev state found in fresco_output directory!")
            return
    else:
        # there are two separate files for the mixed states here then, unpack both separately and plot with corresponding SF's
        os.chdir(fresco_output_path)
        totalFiles = sorted(os.listdir(fresco_output_path))
        matching_files = []
        for file in totalFiles:
            if file.endswith('.sorted'):
                if fnmatch.fnmatch(file, f'*{energy}*'):
                    matching_files.append(file)
        if len(matching_files) == 2:
            SF0 =SF[0]
            SF1 = SF[1]
            std_err = np.sqrt(np.diag(covariance))
            std_dev0 = std_err[0]
            std_dev1 = std_err[1]
            print(f"Low l-transfer Spectroscopic factor: {SF0:.5f} with a std deviation: {std_dev0:.5f}")
            print(f"High l-transfer Spectroscopic factor: {SF1:.5f} with a std deviation: {std_dev1:.5f}")
            diag = np.diag(covariance)
            total_std = np.sqrt(diag[0] + diag[1])
            print("Total std: ", total_std)
            print(matching_files)
            configs = []
            for file in matching_files:
                tempconf = file.split('_')[2]
                configs.append(tempconf)
            sorted_files = sorted(configs, key=custom_sort_key) # make sure that the lowest transfer config is listed 1st, for creating the .inp file that is created
            sorted_list1 = sort_lists(matching_files, sorted_files)

            data0 = np.loadtxt(sorted_list1[0])
            x0 = data0[:,0]
            y0 = data0[:,1]

            y0_scaled = SF0*y0

            data1 = np.loadtxt(sorted_list1[1])
            x1 = data1[:,0]
            y1 = data1[:,1]

            y1_scaled = SF1*y1

            adwa_curve = y0_scaled + y1_scaled

            plt.plot(x0, adwa_curve, marker='', linestyle='--', color='darkblue', label='Scaled DWBA',alpha=0.5)
            plt.fill_between(x0, adwa_curve - total_std, adwa_curve + total_std, color='dodgerblue', alpha=0.3, label='Std. Dev.') # change to use the error from cov matrix from chisquared
            plt.errorbar(lab_angles, x_sec_vals, error_vals, color='black', marker='x', label='experimental data', ecolor='red', capsize=2, linestyle='None')
            plt.xlabel('CM Angle')
            plt.ylabel('Cross-Section [mb/sr]')
            plt.title(f'DWBA Analysis for {energy} keV state')
            plt.yscale('log')

            plt.minorticks_on()

            plt.legend()
            plt.show()

        else:
            print(f'There were not two fresco outputs for {energy} keV state, which is a mixed l-transfer')
        


def get_chi_squared(MIXED_STRENGTH, x_sec_vals, error_vals, fresco1, fresco2):
    if MIXED_STRENGTH == False: ##if peak is single L value
        SF = 0.001 ##sets intital SF
        best_fit = 1000.0
        best_SF = 1.0
        while SF < 0.9: ##sets the upper limit of SF
            iter = 0
            chi_sum = 0.0
            for data in fresco1: ##loops through the array of points
                SF_fit = data*SF
                chi_val = ((SF_fit - x_sec_vals[iter])/error_vals[iter])**2  ##individual chi val for each point
                chi_sum += chi_val   ##sums up the chi values
                iter += 1
            if chi_sum < best_fit: ##checks to see if the current SF is the lowest chi^2 val and sets it as the best fit.
                best_fit = chi_sum
                best_SF = SF
            #print ("SF is " + str(SF) + " with a chi^2 value of "+str(chi_sum)) ##can select to print to see each individual SF & chi^2 value
            SF += 0.001
        print ("Best SF is " + str(best_SF) + " with a chi^2 value of "+str(best_fit))
        return best_SF, best_fit
    else:
        SF1 = 0.01
        SF2 = 0.01 #sets intital SF
        best_fit = 1000000.
        best_SF1 = 1.0
        best_SF2 = 1.0 
        while SF2 < 0.5: ##setes upper limit 
            SF1 = 0.01
            while SF1 < 0.5: ##sets upper limit
                iter = 0
                chi_sum = 0.0
                for data in fresco1:
                    SF_fit1 = data*SF1
                    SF_fit2 = fresco2[iter]*SF2
                    SF_fit_combined = SF_fit1+SF_fit2
                    chi_val = ((SF_fit_combined - x_sec_vals[iter])/error_vals[iter])**2  ##individual chi val for each point
                    #print(chi_val)
                    chi_sum += chi_val                    ##sums up the chi values
                    iter += 1    
                if chi_sum < best_fit:
                    best_fit = chi_sum
                    best_SF1 = SF1
                    best_SF2 = SF2
                ##print ("SF is L=0: " + str(SF1) + " and L=2: " + str(SF2)+" with a chi^2 values of "+str(chi_sum)) ##can select to print the individual SF & chi^2 values
                SF1 += 0.01
            SF2 += 0.01
        #print ("Best SF is L=0: " + str(best_SF1) + " and L=2: " + str(best_SF2)+" with a chi^2 value of "+str(best_fit))
        return best_SF1, best_SF2


def main():
    dir = os.getcwd()
    minimizing_folder = f'{dir}/minimizing_folder'
    sorted_files = sorted(os.listdir(minimizing_folder))
    MIXED_STRENGTH = None
    for file in sorted_files:
        energy = file.split('_')[1]
        print(f"Starting Chi-Squared for {energy} keV state!")
        data = np.loadtxt(f'{minimizing_folder}/{file}')
        if data.shape[1] == 4:
            MIXED_STRENGTH = False
            angles, x_sec_vals, error_vals, fresco_vals = np.transpose(data)
            fresco = np.asarray(fresco_vals)
            SF, covariance = curve_fit(ADWA_Fit, fresco, x_sec_vals, sigma=error_vals)
            plot_data(MIXED_STRENGTH, SF, covariance, angles, x_sec_vals, error_vals, fresco_vals, energy)
            print('\n')
        else:
            MIXED_STRENGTH = True
            angles, x_sec_vals, error_vals, fresco_vals1, fresco_vals2 = np.transpose(data)
            combined_fresco_arr = np.array([fresco_vals1, fresco_vals2])
            SF, covariance = curve_fit(ADWA_Mixed_Fit, combined_fresco_arr, x_sec_vals, sigma=error_vals)
            print(SF)
            plot_data(MIXED_STRENGTH, SF, covariance, angles, x_sec_vals, error_vals, combined_fresco_arr, energy)
            print('\n')



if __name__ == '__main__':
    main()
