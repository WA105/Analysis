#!/usr/bin/env python

import os
from glob import glob
import numpy as np
import ROOT
import root_numpy
from datetime import datetime
from dateutil import tz
import scipy
import warnings

import sys
import argparse

def file_exists(filename): #returns True if file exists, False if it does not
    return bool(glob(filename))

def get_number_of_subruns(run): #returns number of subruns for a given run
    number_of_subruns = 0
    subrun = 0
    while file_exists('/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'):
        number_of_subruns += 1
        subrun += 1
    return int(number_of_subruns)

def get_number_of_events_subrun(run, subrun): #returns number of events for a given run and subrun
    filename = '/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'
    file_size = os.path.getsize(filename) #get file size in bytes
    total_event_number = file_size - 5.0 - 4.0 #First 5 bytes contain run header, last 4 bytes contain run footer.
    total_event_number /= 35.0 + 1280.0 * 1667.0 * 1.5 #For each event, the first 35 bytes are the event header, then come the 1667 ADC counts for the 1280 channels stored in 12 bit (1.5 bytes) format.
    if not total_event_number.is_integer():
        print('Number of events is not an integer.')
    else:
        return int(total_event_number)

def get_number_of_events_run(run): #returns number of events for a given run for all its subruns
    total_event_number_run = 0
    subrun = 0
    while file_exists('/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'):
        total_event_number_run += get_number_of_events_subrun(run=run, subrun=subrun)
        subrun += 1
    return int(total_event_number_run)

def read_one_event(run, subrun, event): # Reads 3x1x1 binary files and returns ADC counts for each channel in a 1280(channels) x 1667(ticks) array
    filename = '/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofull/' + str(run) + '/' + str(run) + '-' + str(subrun) + '-RecoFull-Parser.root'
    root_file_entries_list = root_numpy.list_branches(filename, 'analysistree/anatree')
    reco_file_values = root_numpy.root2array(filename, 'analysistree/anatree', start=event, stop=event+1, step=1)

    dictionary_reco_file_values = {}

    for root_file_index in range(len(root_file_entries_list)):
        dictionary_reco_file_values[root_file_entries_list[root_file_index]] = reco_file_values[0][root_file_index]

    all_channel_waveform_adc = []
    count = 0
    for i in range(1280):
        if i in dictionary_reco_file_values['RecoWaveform_Channel']:
            all_channel_waveform_adc.append(dictionary_reco_file_values['RecoWaveform_ADC'][count*1667:(count+1)*1667])
            count += 1
        else:
            all_channel_waveform_adc.append(np.zeros(1667))

    return np.array(all_channel_waveform_adc).reshape((1280, 1667))

def subtract_pedestal(ADC_values, method='median'): # subtracts the base level of the readout, setting the true 0:
    ADC_values_minped = np.zeros((1280, 1667))
    if method == 'median':
        for i in range(1280):
            ADC_values_minped[i] = ADC_values[i] - np.median(ADC_values[i])
    elif method == 'mean':
        for i in range(1280):
            ADC_values_minped[i] = ADC_values[i] - np.mean(ADC_values[i])
    else:
        print 'Method not recognized.'
    return ADC_values_minped

def get_reconstruction_variables(run, subrun, event):
    root_file_entries_list = root_numpy.list_branches('/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofast/'
                                                  + str(run) + '/' + str(run) + '-' + str(subrun) + '-RecoFast-Parser.root', 'analysistree/anatree')
    reco_file_values = root_numpy.root2array('/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofast/'
                                                 + str(run) + '/' + str(run) + '-' + str(subrun) + '-RecoFast-Parser.root', 'analysistree/anatree', start=int(event), stop=int(event)+1, step=1)
    dictionary_reco_file_values = {}

    for root_file_index in range(len(root_file_entries_list)):
        dictionary_reco_file_values[root_file_entries_list[root_file_index]] = reco_file_values[0][root_file_index]

    if dictionary_reco_file_values['NumberOfTracks'] >= 1:
        track_number_of_hits_index_position = [dictionary_reco_file_values['Track_NumberOfHits'][0]]
        for i in range(1, dictionary_reco_file_values['NumberOfTracks']):
            track_number_of_hits_index_position.append(track_number_of_hits_index_position[i-1] + dictionary_reco_file_values['Track_NumberOfHits'][i])

        for key in dictionary_reco_file_values.keys():
            if key[:10] == 'Track_Hit_':
                dictionary_reco_file_values[key] = np.split(dictionary_reco_file_values[key], track_number_of_hits_index_position)

    return dictionary_reco_file_values

def point_to_line_dist(m, b, x0, y0):
    #Formulae derived analytically.
    # minimum distance satisfies (slope <dot> (x1-x0, y1-y0) = 0)
    #     where (x0,y0) is the point, and (x1,y1) is the point on the line that is closest the point in question
    linePointX = (x0 + y0*m -b*m)/(m*m + 1)
    linePointY = (m*m*y0 + m*x0 +b)/(m*m + 1)
    distance = ((linePointX-x0)**2 + (linePointY-y0)**2)**0.5
    return distance

#uses the numpy polyfit to make a fit for a list of [x] points and a list of [y] points
# returns a dictionary of the fit parameters
def get_line_of_best_fit(x, y):
    fit_parameters = np.polyfit(x, y, 1)  # polyfit returns highest-degree parameters first. y=mx+b -> [m,b]
    angle_radians = np.arctan(fit_parameters[0])
    angle_degrees = angle_radians * 180.0 / np.pi
    # note: this angle is measured from the 'horizontal' to the slope.
    #       /|
    #      / |
    #     /  | m    a = arctan(m/1)
    #    /   |
    #   /a___|
    #     1
    return {'slope': fit_parameters[0], 'y_intercept': fit_parameters[1], 'angle_radians': angle_radians, 'angle_degrees': angle_degrees}

#takes a list of [x] points and [y] points, plus the width of a rectangle to put around the line
#         change- feed the fit into this function since it's calculated right before calculated again here
def get_points_rectangle_around_line(x, y, width, line_fit):
#    start_point_x = np.min(x)
#    end_point_x = np.max(x)

    max_x = np.max(x)
    min_x = np.min(x)
    max_y = np.max(y)
    min_y = np.min(y)


    slope = line_fit['slope']
    y_intercept = line_fit['y_intercept']
    angle_radians = line_fit['angle_radians']

    #find the y-intercepts of top and bottom lines of the box, given the width of the box
    #       to visualize, expand a narrow rectangle around an angled line, draw a line straight up from the middle of a narrow end and extend the top of the box to the drawn line
    #       the extra distance up is what we're adding to the y_intercept (the height of the box-side-midpoint
    parallel_line_1_y_intercept = y_intercept + width / 2.0 / np.cos(angle_radians)
    parallel_line_2_y_intercept = y_intercept - width / 2.0 / np.cos(angle_radians)

    #finding the slope of the other sides of the rectangle. SINCE these are perpendicular, thie is the slope
    perpendicular_lines_slope = - 1.0 / slope

    # these y-intercepts are calculated by knowing the slope of the perpendicular lines and that they intersect the OG fit at the smallest x-value of the fit.
    #       basicaly::   mx_{0} + b = (-1/m)x_{0} + b_{perp}, and solve for b_{perp}
    #perpendicular_line_1_y_intercept = (slope * start_point_x + y_intercept) - perpendicular_lines_slope * start_point_x
    #perpendicular_line_2_y_intercept = (slope * end_point_x + y_intercept) - perpendicular_lines_slope * end_point_x

    if slope > 0.0:
        perpendicular_line_1_y_intercept = min_y - perpendicular_lines_slope * min_x
        perpendicular_line_2_y_intercept = max_y - perpendicular_lines_slope * max_x
    else:
        perpendicular_line_1_y_intercept = max_y - perpendicular_lines_slope * min_x
        perpendicular_line_2_y_intercept = min_y - perpendicular_lines_slope * max_x

    # these are similar. A parallel line and a perpendicular line meet at these points, and so you solve for (x) when the parametrizations are equal, then use that (x) to find (y).
    point_1 = np.array([(perpendicular_line_1_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_1_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_1_y_intercept])
    point_2 = np.array([(perpendicular_line_2_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_2_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_1_y_intercept])
    point_3 = np.array([(perpendicular_line_2_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_2_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_2_y_intercept])
    point_4 = np.array([(perpendicular_line_1_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_1_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_2_y_intercept])
    area_rectangle = np.sqrt((point_1[0] - point_2[0])**2.0 + (point_1[1] - point_2[1])**2.0) * np.sqrt((point_2[0] - point_3[0])**2.0 + (point_2[1] - point_3[1])**2.0)
    return {'rectangle_points': np.array([point_1, point_2, point_3, point_4]), 'fit_slope': slope, 'fit_y_intercept': y_intercept, 'parallel_y_intercepts': np.array([parallel_line_1_y_intercept, parallel_line_2_y_intercept]), 'perpendicular_slope': perpendicular_lines_slope, 'perpendicular_y_intercepts': np.array([perpendicular_line_1_y_intercept, perpendicular_line_2_y_intercept]), 'area_rectangle': area_rectangle }

#calculates the area of a triangle given these three points. **deprecated**
#no longer used since now using a different method for checking if a point is in a rectangle
def area_triangle(point_1, point_2, point_3):
    area = 0.5 * np.abs(((point_2[0] * point_3[1] - point_3[0] * point_2[1]) - (point_1[0] * point_3[1] - point_3[0] * point_1[1]) + (point_1[0] * point_2[1] - point_2[0] * point_1[1])))
    return area

# vector product, used in is_point_in_rectangle
def dotProd(p0, p1):
    return(p0[0]*p1[0] + p0[1]*p1[1])

#vector difference
def vSub(p0, p1):
    p = [p1[0]-p0[0], p1[1]-p0[1]]
    return p

def sort_rectangle_points(points_rectangle):
    # sort the rectangle points so they go in order around.
    # this is important, becase the analytical way requires three rect-points ABC, such that AB is perp to BC
    xmax = 0
    xmin = 3
    for i in range(3):
        if points_rectangle[i+1][0]>points_rectangle[xmax][0]:
            xmax=i+1
        if points_rectangle[2-i][0]<points_rectangle[xmin][0]:
            xmin=2-i
    # this give points on extremes of the rectangle, which will be on opposite sides of the rectangle.
    ordRectangle = [points_rectangle[xmin], [0,0],points_rectangle[xmax], [0,0]]

    #found the minimum and maximum ones, so now we just need to choose one of the remainders to go in the middle
    #can go either clockwise around (ABCD) or (ADBC). These don't make a difference!

    # I'm sure there's a better way...
    added = False
    for i in range(4):
        if (i!=xmax and i!=xmin):
            if (added):
                ordRectangle[3]=points_rectangle[i]
            else:
                ordRectangle[1]=points_rectangle[i]
                added = True

    return(ordRectangle)

def get_slope_between(p1, p2):
    return( (p2[1]-p1[1])/(p2[0]-p1[0]))

#calculates the cross product of vectors v1 and v2
def xProd(v1, v2):
    if len(v1) != len(v2):
        print("error. Vectors need be of the same dimensionality")
        return([0,0,0])
    if len(v1)==2:
        #assuming 3rd component is zero.
        return([0.0, 0.0, 	v1[0]*v2[1] - v2[0]*v1[1]])
    if len(v1)==3:
        print("3d vectors not yet supported")
        return([0,0,0])
    else:
        return([0,0,0])

# given a an array of three points defining a triangle, determines whether a provided point lies within the triangle. Assumes 2D euclidean geometry.
def is_point_in_triangle( triangle, point):
    # tri points ABC, point P
    AB = vSub(triangle[0], triangle[1])
    BC = vSub(triangle[1], triangle[2])
    CA = vSub(triangle[2], triangle[0])
    AP = vSub(triangle[0], point)
    BP = vSub(triangle[1], point)
    CP = vSub(triangle[2], point)

    x1 = xProd(AB, AP)[2]
    x2 = xProd(BC, BP)[2]
    x3 = xProd(CA, CP)[2]

    x1x2 = x1*x2
    x2x3 = x2*x3

    if (x1x2>0 and x2x3>0):
        return True
    else:
        return False

#formula a fun little byproduct of analytical geometry!
def is_point_in_rectangle(points_rectangle, point, arg):

    #Starting from the min-x one, going sequentially, label the points ABC and D. Our point of curiosity is M
	# not sure if I can explain the logic behind this in the comments, sorry.
	#
	# I verified this by generating random rectangles, ordering the points with the sort_rectangle_points function
	# 	   plotting the points, and using this function to fill in the space contained by the points.

    AB = vSub(points_rectangle[0], points_rectangle[1])
    AM = vSub(points_rectangle[0], point)
    BC = vSub(points_rectangle[1], points_rectangle[2])
    BM = vSub(points_rectangle[1], point)
    dotABAM = dotProd(AB, AM)
    dotABAB = dotProd(AB, AB)
    dotBCBM = dotProd(BC, BM)
    dotBCBC = dotProd(BC, BC)

    if( (0 < dotABAM) and (dotABAM < dotABAB) and (0 < dotBCBM) and (dotBCBM < dotBCBC) ):
        return True
    else:
        return False

# ---- begin unit conversions ----

def convert_tick_to_x_position(tick, drift_velocity=0.00160562889065):
#    drift_velocity_cm_per_us = ((-0.0464*(87.0 - 105.749) + 1.0) * (1.88125 * 0.5 * np.log(1.0 + 0.99408 / 0.5) + 0.01172 * 0.5**4.20214) + 0.01712 * (87.0 - 105.749)) / 10.0 #in cm/us
#    drift_velocity_m_per_us = drift_velocity_cm_per_us * 0.01
    x_position = - float(tick) * 0.4 * drift_velocity + 0.5
    return x_position

def convert_x_position_to_tick(x_position, drift_velocity=0.00160562889065):
    tick = int(np.floor((0.5-x_position)/(0.4*drift_velocity)))
    return tick

def convert_x_position_to_us(x_position, drift_velocity=0.00160562889065):
#    drift_velocity_cm_per_us = ((-0.0464*(87.0 - 105.749) + 1.0) * (1.88125 * 0.5 * np.log(1.0 + 0.99408 / 0.5) + 0.01172 * 0.5**4.20214) + 0.01712 * (87.0 - 105.749)) / 10.0 #in cm/us
#    drift_velocity_m_per_us = drift_velocity_cm_per_us * 0.01
    time_us = (0.5 - x_position / 100.0) / drift_velocity
    return time_us

def convert_y_position_to_channel_number_view0(y_position):
    channel_number = int(np.floor((y_position+0.48)/0.003))
    return channel_number

def convert_z_position_to_channel_number_view1(z_position):
    channel_number = int(np.floor(z_position/ 0.003))
    return channel_number

def convert_channel_number_view0_to_y_position(channel_number):
    y_position = float(channel_number) * 0.003 - 0.48
    return y_position

def convert_channel_number_view1_to_z_position(channel_number):
    z_position = float(channel_number) * 0.003
    return z_position

# ---- end unit conversions ----

def get_track_2D_positions_in_each_view(run, subrun, event, track):
    track_variables = get_reconstruction_variables(run, subrun, event)
    reco_hits_x_view0 = []
    reco_hits_y_view0 = []
    reco_hits_x_view1 = []
    reco_hits_z_view1 = []
    for i in range(len(track_variables['Track_Hit_X'][track])):
        if track_variables['Track_Hit_View'][track][i] == 0:
            reco_hits_x_view0.append(track_variables['Track_Hit_X'][track][i])
            reco_hits_y_view0.append(track_variables['Track_Hit_Y'][track][i])
        elif track_variables['Track_Hit_View'][track][i] == 1:
            reco_hits_x_view1.append(track_variables['Track_Hit_X'][track][i])
            reco_hits_z_view1.append(track_variables['Track_Hit_Z'][track][i])
    reco_hits_x_view0 = np.array(reco_hits_x_view0) / 100.0
    reco_hits_y_view0 = np.array(reco_hits_y_view0) / 100.0
    reco_hits_x_view1 = np.array(reco_hits_x_view1) / 100.0
    reco_hits_z_view1 = np.array(reco_hits_z_view1) / 100.0
    return {'reco_hits_x_view0': reco_hits_x_view0, 'reco_hits_y_view0': reco_hits_y_view0, 'reco_hits_x_view1': reco_hits_x_view1, 'reco_hits_z_view1': reco_hits_z_view1}

def get_rectangle_parameters_for_track(run, subrun, event, track, rectangle_width):
    reco_hits_track = get_track_2D_positions_in_each_view(run, subrun, event, track)
    lin_fit_view0 = get_line_of_best_fit(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'])
    lin_fit_view1 = get_line_of_best_fit(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'])
    rectangle_points_view0 = get_points_rectangle_around_line(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'], rectangle_width, lin_fit_view0)
    rectangle_points_view1 = get_points_rectangle_around_line(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'], rectangle_width, lin_fit_view1)
    return {'lin_fit_view0': lin_fit_view0, 'lin_fit_view1': lin_fit_view1, 'rectangle_points_view0': rectangle_points_view0, 'rectangle_points_view1': rectangle_points_view1, 'reco_hits_x_view0': reco_hits_track['reco_hits_x_view0'], 'reco_hits_y_view0': reco_hits_track['reco_hits_y_view0'], 'reco_hits_x_view1': reco_hits_track['reco_hits_x_view1'], 'reco_hits_z_view1': reco_hits_track['reco_hits_z_view1']}

def segment_track_in_y(x, y, number_of_segments):
    max_y = np.max(y)
    min_y = np.min(y)
    segments_y = []
    for i in range(number_of_segments):
        segments_y.append([min_y + i * (max_y - min_y) / number_of_segments, min_y + (i + 1) * (max_y - min_y) / number_of_segments])
    x_segmented = []
    y_segmented = []
    for i in segments_y:
        x_1segment = []
        y_1segment = []
        for j in range(np.size(x)):
            if y[j] + 1.0e-6 >= i[0] and y[j] < i[1] + 1.0e-6:
                x_1segment.append(x[j])
                y_1segment.append(y[j])
        x_segmented.append(np.array(x_1segment))
        y_segmented.append(np.array(y_1segment))
    return np.array(x_segmented), np.array(y_segmented)

def get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, rectangle_width, reco_hits_track):
    lin_fit_view0 = get_line_of_best_fit(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'])
    lin_fit_view1 = get_line_of_best_fit(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'])
    rectangle_points_view0 = get_points_rectangle_around_line(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'], rectangle_width, lin_fit_view0)
    rectangle_points_view1 = get_points_rectangle_around_line(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'], rectangle_width, lin_fit_view1)
    return {'lin_fit_view0': lin_fit_view0, 'lin_fit_view1': lin_fit_view1, 'rectangle_points_view0': rectangle_points_view0, 'rectangle_points_view1': rectangle_points_view1, 'reco_hits_x_view0': reco_hits_track['reco_hits_x_view0'], 'reco_hits_y_view0': reco_hits_track['reco_hits_y_view0'], 'reco_hits_x_view1': reco_hits_track['reco_hits_x_view1'], 'reco_hits_z_view1': reco_hits_track['reco_hits_z_view1']}

def sum_up_charge_in_rectangle(run, subrun, event, track, rectangle_width, include_all_points_in_rectangle=False):
    rectangle_parameters_track = get_rectangle_parameters_for_track(run, subrun, event, track, rectangle_width)
    raw_info_minped = read_one_event(run, subrun, event)
#    raw_info_minped = subtract_pedestal(raw_info)
    max_chan_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    min_chan_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    max_tick_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    min_tick_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    max_chan_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    min_chan_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    max_tick_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    min_tick_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    charge_rectangle_view0 = 0.0
    charge_rectangle_view1 = 0.0
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = []
        points_in_rectangle_view0_x = []
        points_in_rectangle_view1_z = []
        points_in_rectangle_view1_x = []
    for channel_view0 in range(320):
        y_position_channel = convert_channel_number_view0_to_y_position(channel_view0)
        if min_chan_view0_rectangle < y_position_channel and y_position_channel < max_chan_view0_rectangle:
            for tick_view0 in range(1667):
                x_position_tick = convert_tick_to_x_position(tick_view0)
                if min_tick_view0_rectangle < x_position_tick and x_position_tick < max_tick_view0_rectangle:
                    if is_point_in_rectangle(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'], np.array([y_position_channel, x_position_tick]), rectangle_parameters_track['rectangle_points_view0']['area_rectangle']):
                        charge_rectangle_view0 += raw_info_minped[channel_view0][tick_view0]
                        if include_all_points_in_rectangle:
                            points_in_rectangle_view0_y.append(y_position_channel)
                            points_in_rectangle_view0_x.append(x_position_tick)
    for channel_view1 in range(960):
        z_position_channel = convert_channel_number_view1_to_z_position(channel_view1)
        if min_chan_view1_rectangle < z_position_channel and z_position_channel < max_chan_view1_rectangle:
            for tick_view1 in range(1667):
                x_position_tick = convert_tick_to_x_position(tick_view1)
                if min_tick_view1_rectangle < x_position_tick and x_position_tick < max_tick_view1_rectangle:
                    if is_point_in_rectangle(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'], np.array([z_position_channel, x_position_tick]), rectangle_parameters_track['rectangle_points_view1']['area_rectangle']):
                        charge_rectangle_view1 += raw_info_minped[channel_view1+320][tick_view1]
                        if include_all_points_in_rectangle:
                            points_in_rectangle_view1_z.append(z_position_channel)
                            points_in_rectangle_view1_x.append(x_position_tick)
    charge_fC_rectangle_view0 = charge_rectangle_view0 / 54.77942
    charge_fC_rectangle_view1 = charge_rectangle_view1 / 67.08538
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = np.array(points_in_rectangle_view0_y)
        points_in_rectangle_view0_x = np.array(points_in_rectangle_view0_x)
        points_in_rectangle_view1_z = np.array(points_in_rectangle_view1_z)
        points_in_rectangle_view1_x = np.array(points_in_rectangle_view1_x)
    if include_all_points_in_rectangle:
        dictionary = {'ADC_counts_minped': raw_info_minped, 'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'points_in_rectangle_view0_y': points_in_rectangle_view0_y, 'points_in_rectangle_view0_x': points_in_rectangle_view0_x, 'points_in_rectangle_view1_z': points_in_rectangle_view1_z, 'points_in_rectangle_view1_x': points_in_rectangle_view1_x}
    else:
        dictionary = {'ADC_counts_minped': raw_info_minped, 'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1}
    dictionary_tot = dict(dictionary.items() + rectangle_parameters_track.items())
    return dictionary_tot

#using our defined rectangles around the tracks, we find the charge enclosed.
def sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, rectangle_width, rectangle_parameters_track, include_all_points_in_rectangle=False):
    raw_info_minped = read_one_event(run, subrun, event)
#    raw_info_minped = subtract_pedestal(raw_info)
    max_chan_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    min_chan_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    max_tick_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    min_tick_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    max_chan_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    min_chan_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    max_tick_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    min_tick_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    charge_rectangle_view0 = 0.0
    charge_rectangle_view1 = 0.0
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = []
        points_in_rectangle_view0_x = []
        points_in_rectangle_view1_z = []
        points_in_rectangle_view1_x = []

    # get the channel range around the rectangle
    # I add the +/- 1 as a buffer.
    # start:
    channel_view0 = convert_y_position_to_channel_number_view0(min_chan_view0_rectangle) - 1
    #channel_view0 = 0
    # end:
    last_channel =  convert_y_position_to_channel_number_view0(max_chan_view0_rectangle) + 1

    # here I call a function to sort the rectangle points into a proper sequential format needed by my
    # is_point_in_rectangle function
    sorted_rectangle_points_0= sort_rectangle_points( rectangle_parameters_track['rectangle_points_view0']['rectangle_points'] )
    sorted_rectangle_points_1= sort_rectangle_points( rectangle_parameters_track['rectangle_points_view1']['rectangle_points'] )

    while(channel_view0<last_channel and channel_view0 < 320 ):
        y_position_channel = convert_channel_number_view0_to_y_position(channel_view0)
        #now, for these, we're looking for the range in the box. need fit
        # I'm using the fit, which is used to define the box as well, to find the bottom and top ticks of the box at this y-address
        # please note the weird notation. The "y_intercept" of this line is actually an x-intercept
        # ALSO, the conversion between tick and x-position has a sign change, meaning small tick correspond to high position and vice versa
        #      what thsi means is that the max position is the smallest tick and the min position is the biggest tick!
        # similar to before, I add the  +/-1 as a buffer
        tick_view0=convert_x_position_to_tick(y_position_channel*rectangle_parameters_track['lin_fit_view0']['slope'] + rectangle_parameters_track['lin_fit_view0']['y_intercept'] + 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view0']['angle_radians']) ) -1
        last_tick= convert_x_position_to_tick(y_position_channel*rectangle_parameters_track['lin_fit_view0']['slope'] + rectangle_parameters_track['lin_fit_view0']['y_intercept'] - 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view0']['angle_radians']) ) +1

        # the buffer and the fit can give nonsense negative ticks (since I parametrizing the above conversion for every contingency is not worth it),
        # so I instead just make sure it doesn't go too low or too high
        if tick_view0<0:
            tick_view0 = 0

        while(tick_view0<1667 and tick_view0<last_tick):
            # scanning across all ticks in the tick range
            # convert tick to position
            x_position_tick = convert_tick_to_x_position(tick_view0)

            if is_point_in_rectangle( sorted_rectangle_points_0 , [y_position_channel, x_position_tick], rectangle_parameters_track['lin_fit_view0']  ):
                #print("returned true")
                charge_rectangle_view0 += raw_info_minped[int(channel_view0)][int(tick_view0)]
                if include_all_points_in_rectangle:
                    points_in_rectangle_view0_y.append(y_position_channel)
                    points_in_rectangle_view0_x.append(x_position_tick)
            tick_view0 += 1
        channel_view0+=1

    # doing the same for channel 1! Note there is a different dimensionality here!
    # start:
    channel_view1 = convert_z_position_to_channel_number_view1(min_chan_view1_rectangle) - 1
    #channel_view1 =  0
    # end:
    last_channel =  convert_z_position_to_channel_number_view1(max_chan_view1_rectangle) + 1
    while(channel_view1<last_channel and channel_view1 < 960 ):
        z_position_channel = convert_channel_number_view1_to_z_position(channel_view1)
        #again, we use the fits (now in view1) to get the top and bottom tick at this channel.
        # I know, the y-intercept part is weird. It's more like, the fit generator is general purpose
        # and in most cases your vertical axis is the y-axis. So the fit generator returns a "y_intercept" in its dictionary
        # and our geometry is weird, so that "y-axis" actually points along the x-axis: the channel axis.


        tick_view1=convert_x_position_to_tick(z_position_channel*rectangle_parameters_track['lin_fit_view1']['slope'] + rectangle_parameters_track['lin_fit_view1']['y_intercept'] + 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view1']['angle_radians']) ) -1
        last_tick= convert_x_position_to_tick(z_position_channel*rectangle_parameters_track['lin_fit_view1']['slope'] + rectangle_parameters_track['lin_fit_view1']['y_intercept'] - 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view1']['angle_radians']) ) -1

        if tick_view1<0:
            tick_view1 = 0


        while(tick_view1<1667 and tick_view1<last_tick):
            x_position_tick = convert_tick_to_x_position(tick_view1)
            #print("outer if state")
            if is_point_in_rectangle(sorted_rectangle_points_1, [z_position_channel, x_position_tick], rectangle_parameters_track['lin_fit_view0'] ):
                #print("returned true")
                charge_rectangle_view1 += raw_info_minped[int(channel_view1+320)][int(tick_view1)]
                if include_all_points_in_rectangle:
                    points_in_rectangle_view1_z.append(z_position_channel)
                    points_in_rectangle_view1_x.append(x_position_tick)
            tick_view1 += 1
        channel_view1 += 1

    #scale the charge properly using the calibration constants:
    charge_fC_rectangle_view0 = charge_rectangle_view0 / 54.77942
    charge_fC_rectangle_view1 = charge_rectangle_view1 / 67.08538
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = np.array(points_in_rectangle_view0_y)
        points_in_rectangle_view0_x = np.array(points_in_rectangle_view0_x)
        points_in_rectangle_view1_z = np.array(points_in_rectangle_view1_z)
        points_in_rectangle_view1_x = np.array(points_in_rectangle_view1_x)
    if include_all_points_in_rectangle:
        dictionary = {'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'points_in_rectangle_view0_y': points_in_rectangle_view0_y, 'points_in_rectangle_view0_x': points_in_rectangle_view0_x, 'points_in_rectangle_view1_z': points_in_rectangle_view1_z, 'points_in_rectangle_view1_x': points_in_rectangle_view1_x, 'run': run, 'subrun': subrun, 'event': event, 'track': track, 'rectangle_width': rectangle_width}
    else:
        dictionary = {'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'run': run, 'subrun': subrun, 'event': event, 'track': track, 'rectangle_width': rectangle_width}
    dictionary_tot = dict(dictionary.items() + rectangle_parameters_track.items())
    return dictionary_tot

def get_charge_in_multiple_rectangles(run, subrun, event, track, number_of_segments, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle=False):
    dictionaries_charge_sum_small_rectangle = []
    dictionaries_charge_sum_large_rectangle = []
    reco_hits_tot_track = get_track_2D_positions_in_each_view(run=run, subrun=subrun, event=event, track=track)
    reco_hits_segmented_track_view0_y, reco_hits_segmented_track_view0_x = segment_track_in_y(x=reco_hits_tot_track['reco_hits_y_view0'], y=reco_hits_tot_track['reco_hits_x_view0'], number_of_segments=number_of_segments)
    reco_hits_segmented_track_view1_z, reco_hits_segmented_track_view1_x = segment_track_in_y(x=reco_hits_tot_track['reco_hits_z_view1'], y=reco_hits_tot_track['reco_hits_x_view1'], number_of_segments=number_of_segments)
    for i in range(number_of_segments):
        reco_hits_track_segment = {'reco_hits_x_view0': reco_hits_segmented_track_view0_x[i], 'reco_hits_y_view0': reco_hits_segmented_track_view0_y[i], 'reco_hits_x_view1': reco_hits_segmented_track_view1_x[i], 'reco_hits_z_view1': reco_hits_segmented_track_view1_z[i]}
        small_rectangle_parameters_track_segment = get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, small_rectangle_width, reco_hits_track_segment)
        large_rectangle_parameters_track_segment = get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, large_rectangle_width, reco_hits_track_segment)
        charge_in_small_rectangle = sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, small_rectangle_width, small_rectangle_parameters_track_segment, include_all_points_in_rectangle)
        charge_in_large_rectangle = sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, large_rectangle_width, large_rectangle_parameters_track_segment, include_all_points_in_rectangle)
        dictionaries_charge_sum_small_rectangle.append(charge_in_small_rectangle)
        dictionaries_charge_sum_large_rectangle.append(charge_in_large_rectangle)
    return dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle

def get_charge_in_multiple_rectangles_handle_warnings(run, subrun, event, track, number_of_segments, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle=False):
    for i in np.arange(number_of_segments, 0, -1):
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle = get_charge_in_multiple_rectangles(run, subrun, event, track, i, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle)
                break
            except (TypeError, np.RankWarning):
                print 'Not enough hits per box for ' + str(i) + ' segments (run: ' + str(run) + ', subrun: ' + str(subrun) + ', event: ' + str(event) + ', track: ' + str(track) + '). Trying again with fewer boxes.'
    return dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle

def select_tracks_based_on_length_position( run, subrun , max_number_of_tracks=10, min_track_length=90.0 ):
    center_lem_z_position = [50.0, 238.0]
    chosen_tracks = []
    number_of_events_in_subrun = get_number_of_events_subrun(run=run, subrun=subrun)
    for event in range(number_of_events_in_subrun):
        reco_info = get_reconstruction_variables(run=run, subrun=subrun, event=event)
        if reco_info['NumberOfTracks'] <= max_number_of_tracks:
            for track in range(reco_info['NumberOfTracks']):
                if np.abs(reco_info['Track_StartPoint_X'][track] - reco_info['Track_EndPoint_X'][track]) >= min_track_length:
#                   if reco_info['Track_StartPoint_Z'][track] > center_lem_z_position[0] and reco_info['Track_StartPoint_Z'][track] < center_lem_z_position[1]:
#                   if reco_info['Track_EndPoint_Z'][track] > center_lem_z_position[0] and reco_info['Track_EndPoint_Z'][track] < center_lem_z_position[1]:
                    chosen_tracks.append(np.array([run, subrun, event, track]))
                    print 'Possible track for Run: ' + str(run) + ', SubRun: ' + str(subrun) + ', Event: ' + str(event) + ', Track: ' + str(track)
    return np.array(chosen_tracks)

def get_box_charge_ratios_multiple_boxes(run, initial_track_selection, small_rectangle_width, large_rectangle_width, number_of_boxes=1, min_box_length=-1.0):
    charge_ratios_tracks = []
    for i in initial_track_selection:
        if min_box_length > 0.0:
            number_of_boxes = int(get_reconstruction_variables(i[0], i[1], i[2])['Track_Length_StraightLine'][int(i[3])] / min_box_length)
            print number_of_boxes
        if number_of_boxes == 0:
            number_of_boxes = 1
        try:
            dictionary_small_box, dictionary_large_box = get_charge_in_multiple_rectangles_handle_warnings(i[0], i[1], i[2], i[3], number_of_boxes, small_rectangle_width, large_rectangle_width)
            charge_ratios_view0 = []
            charge_ratios_view1 = []
            for j in range(len(dictionary_small_box)):
                charge_small_box_view0 = dictionary_small_box[j]['charge_fC_rectangle_view0']
                charge_small_box_view1 = dictionary_small_box[j]['charge_fC_rectangle_view1']
                charge_large_box_view0 = dictionary_large_box[j]['charge_fC_rectangle_view0']
                charge_large_box_view1 = dictionary_large_box[j]['charge_fC_rectangle_view1']
                charge_ratio_view0 = (charge_large_box_view0 - charge_small_box_view0) / charge_small_box_view0
                charge_ratio_view1 = (charge_large_box_view1 - charge_small_box_view1) / charge_small_box_view1
                charge_ratios_view0.append(charge_ratio_view0)
                charge_ratios_view1.append(charge_ratio_view1)
            charge_ratios_tracks.append({'run': i[0], 'subrun': i[1], 'event': i[2], 'track': i[3], 'charge_ratios_view0': charge_ratios_view0, 'charge_ratios_view1': charge_ratios_view1})
            print 'Charge ratios for Run ' + str(i[0]) + ', SubRun ' + str(i[1]) + ', Event ' + str(i[2]) + ', Track ' + str(i[3]) + ':   View 0 ' + str(charge_ratios_view0) + ', View 1 ' + str(charge_ratios_view1)
        except UnboundLocalError:
            print 'Could not fit track properly with one segment: Run ' + str(i[0]) + ', SubRun ' + str(i[1]) + ', Event ' + str(i[2]) + ', Track ' + str(i[3])
    return charge_ratios_tracks

def splitPath(filename):

    path = filename.split( "/" )
    name = path[-1]

    return path, name

def isroot(filename):

    isroot = False
    if filename.endswith('.root'):
        isroot = True
    return isroot

def isreco( path, recoversion ):
    
    found = False
    for entry in path:
        if entry == recoversion:
            found = True
    return found

#TODO: remove hardcoded filename in functions.
#Temporary fix: check on the reco version chosen
def getRunAndSubrun(filename, recoversion = '2018_June_24' ):
    """
    Assume filename given in the form /path/to/file/run-subrun-RecoFull-Parser.root
    """

    path, name = splitPath( filename )

    if isreco( path, recoversion ):
        run = name.split( "-" )[0]
        subrun = name.split( "-" )[1]
    else:
        print 'Invalid reconstruction version! Valid is: ' + recoversion

    return run, subrun

#here starts the main code
def main():

    #/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofull/' + str(run) + '/' + str(run) + '-' + str(subrun) + '-RecoFull-Parser.root'
    parser = argparse.ArgumentParser(description='3x1x1 track selection.')
    parser.add_argument('-i', '--input', help="Input file", default='')
    parser.add_argument('-o', '--outdir', help="Output directory", default='./')
    args = parser.parse_args()

    filename = args.input
    opath = args.outdir

    if isroot( filename ):
        run, subrun = getRunAndSubrun( filename )
    else:
        print "File is not a .root file"
        sys.exit()

    #select tracks using only length and position
    initial_track_selection = select_tracks_based_on_length_position(run, subrun , max_number_of_tracks=10000000, min_track_length=20.0 )
    np.save(opath+'initial_track_selection_run' + str(run) + '_subrun' + str(subrun) + '.npy', initial_track_selection)

    #select tracks using highway algorithm
    initial_track_selection_file = opath+'initial_track_selection_run' + str(run) + '_subrun' + str(subrun) + '.npy'
    box_charge_ratios = get_box_charge_ratios_multiple_boxes(run, initial_track_selection, small_rectangle_width=0.035, large_rectangle_width=0.1, number_of_boxes=1, min_box_length=50.0)
    np.save(opath+'box_charge_ratios_run' + str(run) + '_subrun' + str(subrun) + '.npy', box_charge_ratios)

    print 'all done'

if __name__ == "__main__":
    main()
