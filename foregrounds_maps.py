'''
Created on September 29, 2024
author: Gabriel S. Costa

Description: This script generates intensity maps (I, Q, U) at a user-defined
frequency range using PySM3 (Python Sky Model https://pysm3.readthedocs.io/en/latest/ )
for selected emission models. Maps are saved as .png images and data is stored 
in a .fits file, separated into directories for each frequency.


'''

import pysm3
import pysm3.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os

####################################################################
#### SET SIMULATION PARAMETERS AND OUTPUT DIRECTORY
####################################################################

# Define the simulation resolution (HEALPix nside)
nside = 128
preset_strings=["d1","s1"] #Default models ["d1","s1"] (dust and synchrotron). Summary of Models available in https://pysm3.readthedocs.io/en/latest/models.html

# Set the observation frequency range (in MHz)
Frequency_range = [980,1260]  # [initial,final] MHz
Frequency_steps = 10 # 10 MHz

# Set the output directory for saving the results
output_dir = "/path/to/output/directory/" #ex.: "/home/gcosta/input_maps/teste_spectra/"

####################################################################
#### DEFINE FUNCTIONS FOR SAVING PLOTS AND FITS DATA
####################################################################

# Function to save the map plot as .png, with overwrite
def save_figure(filename, figure, output_dir):
    '''
    Saves the map plot as a .png file. If the file already exists,
    it will be overwritten.

    filename: Name of the .png file.
    figure: Matplotlib figure object to save.
    '''
    file_path = os.path.join(output_dir, filename)
    if os.path.exists(file_path):
        print(f"File {filename} already exists, it will be overwritten.")
    figure.savefig(file_path)
    plt.close()

# Function to save data as a .fits file, with overwrite
def save_fits_data(data, filename, output_dir):
    '''
    Saves the emission data (HEALPix map) as a .fits file. If the file
    already exists, it will be overwritten.

    data: HEALPix map data to save.
    filename: Name of the .fits file (without extension).
    '''
    fits_path = os.path.join(output_dir, filename + ".fits")
    if os.path.exists(fits_path):
        print(f"File {fits_path} already exists, it will be overwritten.")
    hp.write_map(fits_path, data, overwrite=True)

#######################################################################
#### RUN SKY MODEL SIMULATION AND GENERATE MAPS FOR EACH FREQUENCY
#######################################################################

    
for Frequency in np.arange(Frequency_range[0],Frequency_range[1]+1,Frequency_steps):
    # Create a sky model with emission components in preset_strings
    sky_model = pysm3.Sky(nside=nside, preset_strings=preset_strings)

    # Obtain the emission map and convert units to uK_CMB
    frequency = Frequency* u.MHz
    emission_map = sky_model.get_emission(frequency)
    # emission_map = emission_map.to(u.uK_CMB,  equivalencies=u.cmb_equivalencies(frequency)) #uncomment this line to convert to uK_CMB units
    units = emission_map.unit
    output_path = output_dir + f"{Frequency}"
    os.makedirs(output_path, exist_ok=True)

    ####################################################################
    #### PLOT AND SAVE INTENSITY (I, Q, U) MAPS
    ####################################################################

    # Plot and save the intensity maps in equatorial coordinates

    rotator = hp.Rotator(coord=['G', 'C'])  # Rotates from Galactic ('G') to Equatorial ('C')

    # Plot and save the intensity (I) map
    hp.mollview(rotator.rotate_map_pixel(emission_map[0]), norm='hist', title=f"Intensity Map I at {Frequency}MHz {preset_strings}", unit=units, coord='C')
    save_figure("intensity_I.png", plt, output_path)

    # Plot and save the intensity (Q) map
    hp.mollview(rotator.rotate_map_pixel(emission_map[1]), norm='hist', title=f"Intensity Map Q at {Frequency}MHz {preset_strings}", unit=units, coord='C')
    save_figure("intensity_Q.png", plt, output_path)

    # Plot and save the intensity (U) map
    hp.mollview(rotator.rotate_map_pixel(emission_map[2]), norm='hist', title=f"Intensity Map U at {Frequency}MHz {preset_strings}", unit=units, coord='C')
    save_figure("intensity_U.png", plt, output_path)

    ####################################################################
    #### SAVE THE EMISSION DATA TO A .FITS FILE
    ####################################################################

    # Save the HEALPix emission data as a .fits file
    filename = f"intensity_data_at_{Frequency}MHz_" + "_".join(preset_strings)
    save_fits_data(emission_map, filename, output_path)
