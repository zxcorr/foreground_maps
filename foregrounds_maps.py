'''
Created on September 29, 2024
author: Gabriel S. Costa

Description: This script generates intensity maps (I, Q, U) at a user-defined
frequency range using PySM3 (Python Sky Model https://pysm3.readthedocs.io/en/latest/ )
for selected emission models. The maps for all frequencies are combined into a
single FITS file with a matrix data structure, and plots are saved as .png images.
The output FITS file follows the convention of the example file, storing 30
data channels and a separate HDU with 31 frequency bin edges.

'''

import pysm3
import pysm3.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits

####################################################################
#### SET SIMULATION PARAMETERS AND OUTPUT DIRECTORY
####################################################################

# Define the simulation resolution (HEALPix nside)
nside = 128
preset_strings = ["d1", "s1"] # Default models ["d1","s1"] (dust and synchrotron)

# Set the observation frequency range (in MHz)
Frequency_range = [980, 1260]  # [initial, final] MHz
Frequency_nbins = 30           # Number of frequency channels

# Set the output directory for saving the results
output_dir = "/path/to/output/directory/"

####################################################################
#### DEFINE FUNCTIONS FOR SAVING PLOTS
####################################################################

# Function to save the map plot as .png, with overwrite
def save_figure(filename, figure, output_dir):
    '''
    Saves the map plot as a .png file.
    '''
    file_path = os.path.join(output_dir, filename)
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    if os.path.exists(file_path):
        print(f"File {filename} already exists, it will be overwritten.")
    figure.savefig(file_path)
    plt.close()

#######################################################################
#### RUN SKY MODEL SIMULATION AND GENERATE MAPS FOR ALL FREQUENCIES
#######################################################################

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Generate frequency bin edges (31 values)
frequency_edges = np.linspace(Frequency_range[0], Frequency_range[1], Frequency_nbins + 1)

# Generate frequency bin centers (30 values) for map generation
frequency_centers = (frequency_edges[:-1] + frequency_edges[1:]) / 2

# List to store the intensity maps for each frequency channel
map_I_data = []

print(f"Generating maps for {len(frequency_centers)} frequencies from {Frequency_range[0]}MHz to {Frequency_range[1]}MHz...")

for i, freq_center in enumerate(frequency_centers):
    # Create a sky model with emission components in preset_strings
    sky_model = pysm3.Sky(nside=nside, preset_strings=preset_strings)

    # Obtain the emission map and get the intensity component (Stokes I)
    # The frequency unit is required by PySM3
    emission_map = sky_model.get_emission(freq_center * u.MHz)
    map_I = emission_map[0].value # We only store the Stokes I component

    # Append the map data to our list
    map_I_data.append(map_I)

    # Plot and save the intensity (I) map for the first frequency as an example
    if i == 0:
        output_path = os.path.join(output_dir, f"plots/")
        rotator = hp.Rotator(coord=['G', 'C'])
        hp.mollview(rotator.rotate_map_pixel(map_I), norm='hist', title=f"Intensity Map I at {freq_center}MHz", unit=emission_map.unit, coord='C')
        save_figure("intensity_I.png", plt, output_path)

print("Map generation complete. Creating single FITS file...")

#######################################################################
#### ASSEMBLE DATA INTO A SINGLE FITS FILE
#######################################################################

# Convert the list of maps into a NumPy array
# The shape will be (nbins, npix) -> (30, 196608)
map_data_matrix = np.array(map_I_data)

# Transpose the matrix to match the example format (npix, nbins) -> (196608, 30)
map_data_matrix = map_data_matrix.T

# Create the primary HDU (Header and Data Unit) for the main data matrix
primary_hdu = fits.PrimaryHDU(data=map_data_matrix)

# Add custom header keywords to the primary HDU
header = primary_hdu.header
header['FIELD'] = ('FG', 'Field name')
header['STOKES'] = ('I', 'Stokes parameter')
header['NSIDE'] = (nside, 'HEALPix nside resolution')
header['FREQ_MIN'] = (Frequency_range[0], 'Minimum frequency in MHz')
header['FREQ_MAX'] = (Frequency_range[1], 'Maximum frequency in MHz')
header['NBINS'] = (Frequency_nbins, 'Number of frequency channels')
header['COMMENT'] = 'This file contains the full sky map for a range of frequencies.'
header['COMMENT'] = 'The data matrix has shape (pixels, frequencies).'
header['COMMENT'] = 'Frequencies used for data generation are the centers of the bins.'

# Create a second HDU to store the frequency vector (edges)
freq_vector_hdu = fits.ImageHDU(data=frequency_edges)
freq_vector_hdu.header['EXTNAME'] = 'FREQUENCIES'
freq_vector_hdu.header['COMMENT'] = 'This HDU contains the frequency bin edges in MHz for the main data matrix.'

# Create an HDUList to combine both HDUs
hdul = fits.HDUList([primary_hdu, freq_vector_hdu])

# Define the output FITS filename
output_fits_filename = os.path.join(output_dir, f"FG_I_N{nside}_{Frequency_range[0]}mhz{Frequency_range[1]}mhz_{Frequency_nbins}bins.fits")

# Write the HDUList to a single FITS file
hdul.writeto(output_fits_filename, overwrite=True)

print(f"Data saved to a single FITS file: {output_fits_filename}")