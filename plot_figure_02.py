#!/usr/bin/python


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import sys
import pickle
import time


def calcIntegral(delaysVoltageNumpy):
  start = 100 / 0.025
  end = 300 / 0.025
  res = []
  minVolt = -64
  delays = sorted(delaysVoltageNumpy.keys())[1:]
  controlDelay = sorted(delaysVoltageNumpy.keys())[0]
  minVolt = np.min(delaysVoltageNumpy[controlDelay][start:end])
  maxRes = -np.inf
  minRes = np.inf
  controlRes = sum(delaysVoltageNumpy[controlDelay][start:end] + abs(minVolt)) * 0.025
  for delay in delays:
    delayTime = start + (delay / 0.025)
    res.append(sum(delaysVoltageNumpy[delay][start:end] + abs(minVolt)) * 0.025)
    if (res[-1] > maxRes): maxRes = res[-1]
    if (res[-1] < minRes): minRes = res[-1]
  for ind in range(len(res)):
    res[ind] = res[ind] / (controlRes)
  return res

def plot_figure_02_c():
  delayDiff = 1
  delays = np.arange(-20, 60, delayDiff)
  
  res = (pickle.load(open('figure_02_data.pickle','rb')))
  resClose = calcIntegral(res['close'])
  resMiddle = calcIntegral(res['middle'])
  resFar = calcIntegral(res['far'])
  
  delayDiff = 1
  delays = np.arange(-20, 60, delayDiff)
  font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
  matplotlib.rc('font', **font)
  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[20,1]}, sharex=True)
  for ax in [ax1, ax2]:
    ax.plot(delays, resFar[:len(delays)], linewidth = 1, color='red')
    ax.plot(delays, resMiddle[:len(delays)], linewidth = 1, color='black')
    ax.plot(delays, resClose[:len(delays)], linewidth = 1, color='blue')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.tick_params(direction='out', width=1, size=5)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')   
    ax.xaxis.set_ticks([-15,0, 15, 30, 45, 60])
    ax.xaxis.set_ticklabels([-15,0, 15, 30, 45, 60], fontsize = 7)
  ax1.text(delays[np.where(resFar == np.min(resFar))[0]] - 2, np.min(resFar) + 0.2, r'$\mathregular{x\/ =\/ 35\mu m}$', fontsize=6, color="red")
  ax1.text(delays[np.where(resMiddle == np.min(resMiddle))[0]], np.min(resMiddle) - 0.07, r'$\mathregular{x\/ =\/ 0\mu m}$', fontsize=6, color="black")
  ax1.text(delays[np.where(resClose == np.min(resClose))[0]] - 5, np.min(resClose) + 0.12, r'$\mathregular{x\/ =\/ -35\mu m}$', fontsize=6, color="blue")
  ax2.set_xlabel(r'$\mathregular{\Delta t}$ Inhibition vs. excitation (msec)', fontsize = 8, fontweight = 'regular')
  ax1.set_ylabel('Time Integral of \nNMDA Spike (Norm.)', fontsize = 8, fontweight = 'regular')
  ax1.axis([delays[0], 60, 0.3625, 1.01])
  ax2.axis([delays[0], 60, 0, 0.1])
  kwargs = dict(linewidth = 1, transform=ax1.transAxes, color='k', clip_on=False)
  ax1.plot((-0.015, 0.015), (0 + 0.005, 0 - 0.005), **kwargs)  
  kwargs = dict(linewidth = 1, transform=ax2.transAxes, color='k', clip_on=False)
  ax2.plot((-0.015, 0.015), (1 + 0.1, 1 - 0.1), **kwargs)  
  ax1.spines['top'].set_visible(False)
  ax2.spines['top'].set_visible(False)
  ax1.spines['bottom'].set_visible(False)
  ax1.xaxis.set_ticks_position('none')
  plt.subplots_adjust(hspace = 0.1) 
  ax1.yaxis.set_ticks([round(0.4,2), round(0.55,3), round(0.7,2), round(0.85,3), round(1.0,2)])
  ax1.yaxis.set_ticklabels([round(0.4,2), round(0.55,3), round(0.7,2), round(0.85,3), round(1.0,2)], fontsize=7)
  ax2.yaxis.set_ticks([round(0.0,2)])
  ax2.yaxis.set_ticklabels([round(0.0,2)], fontsize=7)
  ax1.yaxis.set_label_coords(- 0.255, 0.4)
  ax.tick_params(direction='out', width=1, size=5)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_linewidth(1)
  ax.spines['left'].set_linewidth(1)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')   
  fig.set_size_inches(1.7, 1.7)
  fig.savefig('figure_02_c.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300)
  plt.close()


def plot_figure_02_b():
  delayDiff = 5
  delays = [-100000, 20] 
  gabaAmps = [0.001]
  nmdaAmps = [0.0004]
  
  font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
  matplotlib.rc('font', **font)
  colors = {-100000: 'black', -10: 'brown', 0 : 'green', 10 : 'blue', 20 : 'red'}
  plots = {}
  res = (pickle.load(open('figure_02_data.pickle','rb')))
  resClose = res['close']
  resMiddle = res['middle']
  resFar = res['far']
  
  timeVectorNumpy = np.arange(0, 250.025, 0.025)
  for d in delays[1:]:
    location_texts = [r'$\mathregular{x\/ =\/ 35\mu m}$' + '\n   %d msec' % d,
    r'$\mathregular{x\/ =\/ 0\mu m}$' + '\n   %d msec' % d,
    r'$\mathregular{x\/ =\/ -35\mu m}$' + '\n   %d msec' % d]
    location_colors = ['red','black','blue']
    location_files = ['figure_02_b_top_far_%d.pdf' % d,
    'figure_02_b_top_middle_%d.pdf' % d,
    'figure_02_b_top_close_%d.pdf' % d]
    location_res = [resFar,resMiddle,resClose]
    for location in range(3):
      fig = plt.figure(figsize=(8, 8))
      ax = fig.add_subplot(111)
      ax.plot([100 + d, 100 + d],[location_res[location][-100000][(100 + d) / 0.025], -4],linewidth=0.75,color='black',linestyle=':')
      ax.plot(timeVectorNumpy, location_res[location][-100000], linewidth = 0.75, color="black", linestyle='--', dashes=(2,2));
      ax.plot(timeVectorNumpy, location_res[location][d], linewidth = 0.75, color='black');
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
      ax.spines['left'].set_visible(False)
      ax.spines['bottom'].set_visible(False)  
      ax.yaxis.set_ticks_position('left')
      ax.xaxis.set_ticks_position('bottom')       
      ax.axis([80, 200, -80, 10], width=4)
      ax.set_xticks([])
      ax.set_yticks([])
      plt.text(100 + d - 22, 10, location_texts[location], color = location_colors[location], fontsize=6)
      ax.arrow(100 + d, 5, 0, -7.5, head_width=2, width=0.5, head_length=1, fc=location_colors[location], ec=location_colors[location])
      ax.fill_between(timeVectorNumpy, location_res[location][d], -70,color=location_colors[location], alpha=0.3)
      fig.set_size_inches(1, 1)
      fig.savefig(location_files[location], transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
      plt.close()
  
  location_delays = [20, 16, 40]
  location_texts = [r'$\mathregular{x\/ =\/ 35\mu m}$' + '\n   20 msec',
  r'$\mathregular{x\/ =\/ 0\mu m}$' + '\n   16 msec',
  r'$\mathregular{x\/ =\/ -35\mu m}$' + '\n   40 msec']
  location_colors = ['red','black','blue']
  location_files = ['figure_02_b_bottom_far_20.pdf',
  'figure_02_b_bottom_middle_16.pdf',
  'figure_02_b_bottom_close_40.pdf']
  location_res = [resFar,resMiddle,resClose]
  for location in range(3):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot([100 + location_delays[location], 100 + location_delays[location]],[location_res[location][-100000][(100 + location_delays[location]) / 0.025], -4],linewidth=0.75,color='black',linestyle=':')
    ax.plot(timeVectorNumpy, location_res[location][-100000], linewidth = 0.75, color="black", linestyle='--', dashes=(2,2));
    ax.plot(timeVectorNumpy, location_res[location][location_delays[location]], linewidth = 0.75, color='black');
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)  
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')       
    ax.axis([80, 200, -80, 10], width=4)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.text(100 + location_delays[location] - 22, 10, location_texts[location], color = location_colors[location], fontsize=6)
    ax.arrow(100 + location_delays[location], 5, 0, -7.5, head_width=2, width=0.5, head_length=1, fc=location_colors[location], ec=location_colors[location])
    ax.fill_between(timeVectorNumpy, location_res[location][location_delays[location]], -70,color=location_colors[location], alpha=0.3)
    fig.set_size_inches(1, 1)
    fig.savefig(location_files[location], transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
    plt.close()      

plot_figure_02_b()
plot_figure_02_c()