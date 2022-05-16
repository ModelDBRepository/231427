#!/usr/bin/python 

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle
import time

timeVectorNumpy = np.arange(0,250.025, 0.025)
font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
matplotlib.rc('font', **font)

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

def plot_figure_01_b():
  res = pickle.load(open('figure_1_data.pickle','rb'))[0.001]
  font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
  matplotlib.rc('font', **font)
  colors = {-100000: 'black', -10: 'brown', 0 : 'green', 10 : 'blue', 20 : 'red'}
  
  for d in [-10, 0, 10, 20]:
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(np.array([100 + d, 100 + d]), np.array([res[-100000][(100 + d) / 0.025], -4]),linewidth=0.75,color='black',linestyle=':')
    ax.plot(np.arange(0, len(res[0]) * 0.025, 0.025), res[-100000], linewidth = 0.75, color="black", linestyle='--', dashes=(3, 3));
    ax.plot(np.arange(0, len(res[0]) * 0.025, 0.025), res[d], linewidth = 0.75,label=r'$\mathregular{\Delta t\/ =\/ %d\/ (msec)}$' % d, color='black');
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')   
    ax.axis([80, 200, -80, 10], width=1)
    ax.set_xticks([])
    ax.set_yticks([])
    if (d == 20):
      ax.plot([150, 160],[-5, -5],linewidth=0.75,color='black')
      ax.plot([160, 160],[-5, 5],linewidth=0.75,color='black')
      ax.text(148, -15, "10 msec", fontsize = 10)
      ax.text(163, -1, "10 mV", fontsize = 10)
    # plots[d][0].set_color(colors[d])
    plt.text(100 + d - 10, 10, '%d msec' % d, color = colors[d], fontsize=10)
    ax.arrow(100 + d, 5, 0, -7.5, head_width=2, width=0.75, head_length=1, fc=colors[d], ec=colors[d])
    ax.fill_between(np.arange(0, len(res[0]) * 0.025, 0.025), res[d], -70, where = res[d] > -70, color=colors[d], alpha=0.3)
    fig.set_size_inches(1.7, 1.7)
    fig.savefig('figure_01_b_%d.pdf' % d, transparent=True, bbox_inches='tight', format='pdf', dpi=300)
    plt.show(block = 0)  
    

def plot_figure_01_c():
  res = calcIntegral(pickle.load(open('figure_1_data.pickle','rb'))[0.001])

  delayDiff = 1
  delays = np.arange(-20, 60, delayDiff)
  font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
  matplotlib.rc('font', **font)
  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[20,1]}, sharex=True)
  for ax in [ax1, ax2]:
    ax.plot(delays, res, linewidth = 1, color="black"); 
    ax.plot([-10],[res[list(delays).index(i)] for i in [-10]], "o", markersize = 7, color="brown"); 
    ax.plot([0],[res[list(delays).index(i)] for i in [0]], "o", markersize = 7, color="green"); 
    ax.plot([10],[res[list(delays).index(i)] for i in [10]], "o", markersize = 7, color="blue"); 
    ax.plot([20],[res[list(delays).index(i)] for i in [20]],"o", markersize = 7, color="red"); 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.tick_params(direction='out', width=1, size=5)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')   
    ax.xaxis.set_ticks([-15,0, 15, 30, 45, 60])
    ax.xaxis.set_ticklabels([-15,0, 15, 30, 45, 60], fontsize=10)
  ax1.fill_between(delays[:(res.index(min(res)))] , 0 , res[:(res.index(min(res)))], color="lightblue")
  ax1.fill_between(delays[(res.index(min(res)) - 1):] , 0 , res[(res.index(min(res)) - 1):], color="lightpink")
  ax1.text(delays[0] + 3, 0.415, 'Regeneration', fontsize=9, fontweight='regular')
  ax1.text(30, 0.415, 'Termination', fontsize=9, fontweight='regular')  
  ax2.set_xlabel(r'$\mathregular{\Delta t}$ Inhibition vs excitation (msec)', fontsize = 10, fontweight = 'regular')
  ax1.set_ylabel('Time Integral of \nNMDA Spike (Norm.)', fontsize = 10, fontweight = 'regular')
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
  ax1.yaxis.set_ticklabels([round(0.4,2), round(0.55,3), round(0.7,2), round(0.85,3), round(1.0,2)], fontsize=10)
  ax2.yaxis.set_ticks([round(0.0,2)])
  ax2.yaxis.set_ticklabels([round(0.0,2)], fontsize=10)
  ax1.yaxis.set_label_coords(- 0.205, 0.4)
  ax.tick_params(direction='out', width=1, size=5)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_linewidth(1)
  ax.spines['left'].set_linewidth(1)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')   
  fig.set_size_inches(3, 3)
  fig.savefig('figure_01_c.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300)
  plt.close('all')

def plot_figure_01_d():
  res = (pickle.load(open('figure_1_data.pickle','rb')))
  resWeak = calcIntegral(res[0.0005])
  resNormal = calcIntegral(res[0.001])
  resStrong = calcIntegral(res[0.0015])
  delayDiff = 1
  delays = np.arange(-20, 60, delayDiff)
  font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
  matplotlib.rc('font', **font)
  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[20,1]}, sharex=True)
  for ax in [ax1, ax2]:
    ax.plot(delays, resWeak[:len(delays)], linewidth = 1, color='red')
    ax.plot(delays, resNormal[:len(delays)], linewidth = 1, color='black')
    ax.plot(delays, resStrong[:len(delays)], linewidth = 1, color='blue')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.tick_params(direction='out', width=1, size=5)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')   
    ax.xaxis.set_ticks([-15,0, 15, 30, 45, 60])
    ax.xaxis.set_ticklabels([-15,0, 15, 30, 45, 60], fontsize = 10)
  ax1.text(delays[np.where(resWeak == np.min(resWeak))[0]] - 8, np.min(resWeak) + 0.2, r'$\mathregular{g_{GABA_{A}}}$ = 0.5 (nS)', fontsize=9, color="red")
  ax1.text(delays[np.where(resNormal == np.min(resNormal))[0]] - 2, np.min(resNormal) + 0.15, '1.0 (nS)', fontsize=9, color="black")
  ax1.text(delays[np.where(resStrong == np.min(resStrong))[0]] + 10, np.min(resStrong) + 0.025, '1.5 (nS)', fontsize=9, color="blue")
  ax2.set_xlabel(r'$\mathregular{\Delta t}$ Inhibition vs excitation (msec)', fontsize = 10, fontweight = 'regular')
  ax1.set_ylabel('Time Integral of \nNMDA Spike (Norm.)', fontsize = 10, fontweight = 'regular')
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
  ax1.yaxis.set_ticklabels([round(0.4,2), round(0.55,3), round(0.7,2), round(0.85,3), round(1.0,2)], fontsize=10)
  ax2.yaxis.set_ticks([round(0.0,2)])
  ax2.yaxis.set_ticklabels([round(0.0,2)], fontsize=10)
  ax1.yaxis.set_label_coords(- 0.205, 0.4)
  ax.tick_params(direction='out', width=1, size=5)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_linewidth(1)
  ax.spines['left'].set_linewidth(1)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')   
  fig.set_size_inches(3, 3)
  fig.savefig('figure_01_d.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300)
  plt.close()

plot_figure_01_b()
plot_figure_01_c()
plot_figure_01_d()