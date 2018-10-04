# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 13:36:26 2017

This script opens and processes raw HYSCORE data from a Bruker .DTA file. It
baseline corrects, phases, apodizes, and plots the data. The apodization
function used is a customized, tapered diagonal Blackman window
(Carson Mize and Joe Butler) that prioritizes signal on the diagonal.

Once you choose a file, you will be asked "optimize?" on the command line. If
you input y for yes, then you will be able to add a custom contour min/max for
the output plot. If you input n for no, then the default contour levels will
be used.

@author: mmlockart


opens a searchable file window
imports all necessary packages
"""
from math import sqrt, sin, cos, pi
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
#import scipy.interpolate as si
from scipy.optimize import minimize


def hyscore():
    """ SAMPLE DOCSTRING but really, its just hyscore..."""
    def choosefile():
        """this pops up a tkinter window to choose a file. The file must be
        a .DTA or .DSC; otherwise you get an error and the code stops here"""
        root = tk.Tk()
        filename = filedialog.askopenfilename(parent=root)
        if not filename.endswith('.DTA'):
            if filename.endswith('.DSC'):
                filename = filename.replace('.DSC', '.DTA')
            else:
                messagebox.showinfo("Error", "must be a .DTA file")
        return filename


    def create_dict(filename):
        """ SAMPLE DOCSTRING this makes a searchable dictionary out of the
        description file(.DSC file) """
        if filename.endswith(".DSC"):
            filename = filename
        else:
            try:
                filename = filename.replace(".DTA", ".DSC")
                # just in case the .DTA file is chosen instead
            except:
                print('this is not a valid .DSC file')

        dictionary = {}

        file_open = open(filename, 'r')
        for line in file_open:
            line = line.strip()
            lookup = 'PlsSPELPrg'
            # the last line that needs to be in the dictionary
            # before the pulse sequence starts

            if lookup in line:
                break
            else:
                if ("#" not in line and ".DVC" not in line and "begin"
                        not in line and "end" not in line
                        and not line.startswith("\\n\\")
                        and not line.startswith("\\n")
                        and not line.startswith("*")
                        and not line.startswith(";")
                        and line != ''):
                    line = line.split()
                    if "=" in line:
                        dictionary[line[0]] = line[2]
                    else:
                        dictionary[line[0]] = line[1:]
        return dictionary


    def get_from_dict(key):
        """ SAMPLE DOCSTRING this grabs information from the dictionary """
        value = file_dictionary[key]
        if ((key != 'XPTS') and (key != 'XMIN') and (key != 'XWID') \
            and (key != 'YPTS') and (key != 'YMIN') and (key != 'YWID')):
        #if (key not in('XTPS', 'XMIN', 'XWID', 'YPTS', 'YMIN', 'YWID')
            value = " ".join(value)
            value = value.strip("'")
        return value


    def read_dta_file(filename):
        """ SAMPLE DOCSTRING this reads the .DTA file and
        creates the data matrix"""
        rawbin = np.fromfile(filename, dtype='>f8')
        # big endian, f8 is 64-bit floating point number
        rawbin = np.reshape(rawbin, (len(rawbin), 1, ))  # transopse
        value = get_from_dict('IKKF')  # if complex data

        if value == 'CPLX':
            size = rawbin.size/2
            size = int(size)
            raw_data = np.zeros((size, 1), dtype=complex)
        for i in range(0, len(raw_data)):
            raw_data[i, 0] = np.complex(rawbin[((2*i)-1)+1], rawbin[(2*i)+1])

        if value != 'CPLX':
            raw_data = raw_data
        dim_matrix = int(sqrt(len(raw_data)))
        # gets number of points to reshape data into a square matrix
        raw_data = np.reshape(raw_data, (dim_matrix, dim_matrix))

        return raw_data, dim_matrix
        # this returns both in one tuple; must index the
        # tuple to pull out either data or dim_matrix


    def phase(data):
        """ SAMPLE DOCSTRING this phases the data"""
        def f_phase(x_a):
            return np.sum(abs((sin(x_a)*data.real-cos(x_a)*data.imag)))
        start_pos = np.array(0)
        res = minimize(f_phase, start_pos, method='Powell')
        phi = res.x
        hyscorecomp = np.exp(-phi*1j)*data
        phased_data = hyscorecomp.real
        phased = phased_data
        # these variables will all be cleaned up once the script is complete
        return phased


    def baseline(phased):
        """SAMPLE DOCSTRING this baseline corrects both the x and y axis"""
        x_a = np.arange(0, len(phased))
        # vector the same length as the data in one direction

        base = np.zeros((dim_matrix, dim_matrix))
        for i in range(0, len(phased)):
            phased_baseline = np.polyfit(x_a, phased[i, :], 3)
            # baseline correct in y direction
            baseline = np.polyval(phased_baseline, x_a)
            phased[i, :] = phased[i, :]-baseline
            base[i, :] = baseline

        for i in range(0, len(phased)):
            phased_baseline = np.polyfit(x_a, phased[:, i], 3)
            # baseline correct in x direction
            baseline = np.polyval(phased_baseline, x_a)
            phased[:, i] = phased[:, i]-baseline
            base[:, i] = base[:, i]+baseline
            baselined = phased  # this will all be cleaned up later
        return baselined


    def butlermizewindow(baselined):  # Default
        """SAMPLE DOCSTRING this is windowing the data from an improved
        diagonal blackman window that can change the width of the window"""
        wind = np.ones([1, dim_matrix*2])
        a_0 = (1-alpha_0)/2
        a_1 = 1/2
        a_2 = a_1-a_0
        for i in range(0, dim_matrix*2):
            wind[0, i] = (a_0-a_1*cos(2*pi*(i)/(2*dim_matrix-1)) +
                          a_2*cos(4*pi*(i)/(2*dim_matrix-1)))
            # normal blackman function
        for i in range(1+x_m+2*x_0+dim_matrix, x_m+2*x_0+2*dim_matrix):
            # creates the right half
            wind[0, i] = wind[0, i]
        for i in range(0, x_m):  # define left side
            wind[0, i] = x_m+1
        wind[0, x_m+dim_matrix] = 1
        wind[0, x_m+dim_matrix-1] = 1
        # these two lines define the center;
        # they make positions 127, 128 both 1

        wind[0, x_m] = 0  # makes left side zero
        wind[0, 2*dim_matrix-1] = 0  # makes right side zero
        dwind = np.ones([dim_matrix, dim_matrix])
        # create the array for the next step
        for i in range(0, dim_matrix):
            dwind[:, i] = wind[0, (dim_matrix-i):2*dim_matrix-i]
        wind2 = np.ones([1, dim_matrix])
        for i_2 in range(dim_matrix-round(dim_matrix/4), dim_matrix):
            wind2[0, i_2] = abs(sin(pi/2*((i_2)/(round(dim_matrix/4)))))
            # Taper
        wind2[0, dim_matrix-1] = 0
        wind3 = (wind2*(np.ones([dim_matrix, dim_matrix])))
        windowed = baselined*wind3*np.transpose(wind3)*dwind
        return windowed


    def fourier_transform2d(windowed):
        """ SAMPLE DOCSTRING 2D Fourier transform. Data is zerofilled
        (a matrix of all zeros is created and data is dropped in it) first.
        The transform is shifted so that zero frequency is in the middle."""
        zerofill = np.zeros([2048, 2048])
        zerofill[:len(windowed), :len(windowed)] = windowed
        transform = np.fft.fft2(zerofill)
        transform = np.fft.fftshift(transform)
        transformed = np.absolute(transform)
        tmax = transformed.max()
        znormalized = (transformed - 0)/(tmax - 0)
        # if you want a minimum at zero for your normalized plot,
        # then replace zero with the min.
        return znormalized


    def buildxy(file_dictionary):
        """ SAMPLE DOCSTRING this constructs x and y axis based
        off of description file"""
        x_dim = list(map(float, get_from_dict('XPTS')))
        x_dim = float(x_dim[0])
        xmin = list(map(float, get_from_dict('XMIN')))
        xmin = float(xmin[0])
        xrange = list(map(float, get_from_dict('XWID')))
        xrange = float(xrange[0])
        d_x = xrange/(x_dim-1)
        x_axis = (np.arange(xmin, xmin+x_dim*d_x, d_x))

        y_dim = list(map(float, get_from_dict('YPTS')))
        y_dim = float(y_dim[0])
        ymin = list(map(float, get_from_dict('YMIN')))
        ymin = float(ymin[0])
        yrange = list(map(float, get_from_dict('YWID')))
        yrange = float(yrange[0])
        d_y = yrange/(y_dim-1)
        y_axis = np.arange(ymin, ymin+y_dim*d_y, d_y)

        # converts time to frequency scale

        frwidth = 1000/(x_axis[0])
        frinc = frwidth/(len(znormalized))
        freq = np.arange(-frwidth, frwidth, frinc*2)
        x_axis = freq
        y_axis = freq
        return x_axis, y_axis


    def plot1(x_axis, y_axis, z_coord):
        """ SAMPLE DOCSTRING"""
        plt.contourf(newx, newy, znormalized, 100,
                     cmap=plt.cm.rainbow, vmax=0.12, vmin=0.07)
        # more color contours = lower vmax
        # more noise = lower vmin
        # Changing the vmin here, create a better color range since
        # the highest z values cannot be ignored.
        plt.cm.rainbow.set_under(color='white')
        # Changes all values under vmin to white.
        #Makes most of the noise go away.
        title = get_from_dict('TITL')
        #plt.title(title)
        plt.gca().set_aspect('equal')
        plt.ylabel('Frequency (MHz)', fontsize=12)
        plt.xlabel('Frequency (MHz)', fontsize=12)
       
        plt.xlim([8,18])
        plt.ylim([8,18])
        plt.grid(False)
        save = [title, ".tiff"]
        plt.savefig("".join(save), dpi=900, bbox_inches = 'tight', format = "tiff")

#
#        def format_coord(x_axis, y_axis):
#            """SAMPLE DOCSTRING"""
#            z_coord = (np.take(si.RectBivariateSpline(x_1, y_1, znormalized)
#                               (x_axis, y_axis), 0))
#            return ('x = {x_axis:.5f} y = {y_axis:.5f} z = {z_coord:.3f}'
#                    .format(x_axis=x_axis, y_axis=y_axis, z_coord=z_coord))
#        plt.gca().format_coord = format_coord
#        plt.gca().set_aspect('equal')
#        plt.show()
#        return plot1


    def plot2(x_axis, y_axis, z_coord):
        """ SAMPLE DOCSTRING"""
        plt.contourf(newx, newy, znormalized, 100,
                     cmap=plt.cm.rainbow, vmax=nmax, vmin=nmin)
        # Changing the vmin here, create a better color range since the
        # highest z values cannot be ignored. -CM
        plt.cm.rainbow.set_under(color='white')
        # Changes all values under vmin to white. Makes most of the
        # noise go away. -CM
        #title = get_from_dict('TITL')
        #plt.title(title)
        plt.ylabel('Frequency (MHz)', fontsize=12)
        plt.xlabel('Frequency (MHz)', fontsize=12)
        plt.grid(False)
        # plt.ylim((8, 18)) #can be used to automatically zoom the x-axis -CM
        # plt.xlim((8, 18)) #can be used to automatically zoom the y-axis -CM

#
#        def format_coord(x_axis, y_axis):
#            """SAMPLE DOCSTRING"""
#            z_coord = (np.take(si.RectBivariateSpline(x_1, y_1, znormalized)
#                               (x_axis, y_axis), 0))
#            return ('x = {x_axis:.5f} y = {y_axis:.5f} z = {z_coord:.3f}'
#                    .format(x_axis=x_axis, y_axis=y_axis, z_coord=z_coord))
#        plt.gca().format_coord = format_coord
#        plt.gca().set_aspect('equal')
#        plt.show()
#        return plot2


    filename = choosefile()
    file_dictionary = create_dict(filename)
    data = (read_dta_file(filename)[0])
    dim_matrix = (read_dta_file(filename))[1]
    phased = phase(data)
    baselined = baseline(phased)
    dim_matrix = len(baselined)
    x_m = 0
    x_0 = 0
    alpha_0 = 0.16
    windowed = butlermizewindow(baselined)
    # can change to regular hamming and blackamn window if need be. -CM
    znormalized = fourier_transform2d(windowed)
    z_n = znormalized
    x_1 = buildxy(file_dictionary)[0]
    y_1 = buildxy(file_dictionary)[1]
    # x_build= np.transpose(buildxy(file_dictionary)[0])
    # y_build= buildxy(file_dictionary)[1]
    # unused? Can't tell what it does. Commented out
    ones = np.ones((2048, 2048))
    newy = y_1*ones
    x_2 = np.reshape(x_1, (-1, 2048))*ones
    newx = np.swapaxes(x_2, 0, 1)
    

    plotting = input('Optimize?\n')
    # automatic is vmax = 0.25 and vmin = 0.06 -CM
    if plotting in ['yes', 'Yes', 'y', 'Y']:
        nmax = float(input("enter a number for vmax:"))  # 0 to 1
        nmin = float(input("enter a number for vmin:"))  # 0 to 1
        plot2(newx, newy, z_n)
    else:
        plot1(newx, newy, z_n)
    field_pos = get_from_dict('FieldPosition')
    frequency_val = get_from_dict('MWFQ')
    power = get_from_dict('Power')
    tau_val = get_from_dict('d1')
    print('Field:', field_pos, 'Frequency:', frequency_val, 'Tau:', tau_val, 'Power', power)

    return (newx, newy, z_n)


[X, Y, Z] = hyscore()

