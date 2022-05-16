#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import sys
import pickle
import time

font = {'weight':'regular','size':24.0}
matplotlib.rc('font', **font)

def calculateFixedPoints(valuesOverVoltage):
  fixedPoints = []
  keys = valuesOverVoltage.keys()
  timeVector = sorted(valuesOverVoltage[keys[0]].keys())
  for tInd in timeVector:
    IV_PAS = (-1) * np.array(valuesOverVoltage['PAS'][tInd])
    IV_NMDA = np.array(valuesOverVoltage['NMDA'][tInd])
    IV_NMDA_rolled = np.append(IV_NMDA, IV_NMDA[-1])
    IV_NMDA_rolled = np.roll(IV_NMDA_rolled, 1)[:-1]
    IV_PAS_rolled = np.append(IV_PAS, IV_PAS[-1])
    IV_PAS_rolled = np.roll(IV_PAS_rolled, 1)[:-1]
    fixedPoints.append(np.where((((IV_NMDA - IV_PAS) > 0) == ((IV_NMDA_rolled - IV_PAS_rolled) > 0)) == False)[0])
    if fixedPoints[-1][0] == 0:
        fixedPoints[-1] = fixedPoints[-1][1:]
  return fixedPoints


def add_arrow(line, position=None, direction='right', size=30, color=None):
    if color is None:
        color = line.get_color()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 4
    else:
        end_ind = start_ind - 4
    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="-|>", color=color),
        size=size
    )

def plot_figure_03_a(valuesOverVoltage, valuesOverTime):
  fixedPoints = calculateFixedPoints(valuesOverVoltage)
  keys = valuesOverVoltage.keys()
  timeVector = sorted(valuesOverVoltage[keys[0]].keys())
  fig = plt.figure(figsize=(6, 6))
  ax = fig.add_subplot(111)
  ax.axis([-90, 0, -0.05, 0.001], fontsize=24.0)
  plt.xticks([-80, -45, -10], fontsize=24.0)    
  plt.yticks([-0.01, -0.03, -0.05], [-10, -30, -50], fontsize=24.0)
  ax.set_xlabel('mV', weight="regular", fontsize=24.0)
  ax.set_ylabel('pA', weight="regular", fontsize=24.0)
  for recordingTiming in [130 / dt]:
    timeInMillis = recordingTiming * dt
    NMDA_current = np.array(valuesOverVoltage['NMDA'][recordingTiming])
    AMPA_current = np.array(valuesOverVoltage['AMPA'][recordingTiming])
    GABA_current = - np.array(valuesOverVoltage['GABA'][recordingTiming])
    PAS_current =  - np.array(valuesOverVoltage['PAS'][recordingTiming])
    CAP_current =  - np.array(valuesOverVoltage['CAP'][recordingTiming])
    ABS_CAP_current =  - np.abs(np.array(valuesOverVoltage['CAP'][recordingTiming]))
    inhibitory_current = GABA_current + PAS_current + CAP_current
    excitatory_current = NMDA_current
    voltagePoint = valuesOverTime['VOLTAGE'][recordingTiming]
    vclampValues = np.array(valuesOverVoltage['vclamp'][recordingTiming])
    excitatoryLine = ax.plot(voltages, excitatory_current, color ='red', linewidth = 4.0)[0]
    inhibitoryLine = ax.plot(voltages, PAS_current, color ='blue', linewidth = 4.0, alpha = 1.0, linestyle='-')[0]
    add_arrow(excitatoryLine, position=-85, direction='right', size=25, color='red')
    add_arrow(excitatoryLine, position=-80, direction='right', size=25, color='red')
    add_arrow(inhibitoryLine, position=-70, direction='left', size=25, color='blue')
    add_arrow(inhibitoryLine, position=-60, direction='left', size=25, color='blue')
    add_arrow(inhibitoryLine, position=-50, direction='left', size=25, color='blue')
    add_arrow(excitatoryLine, position=-40, direction='right', size=25, color='red')
    add_arrow(excitatoryLine, position=-35, direction='right', size=25, color='red')
    add_arrow(excitatoryLine, position=-30, direction='right', size=25, color='red')
    add_arrow(excitatoryLine, position=-25, direction='right', size=25, color='red')
    add_arrow(excitatoryLine, position=-17, direction='right', size=25, color='red')
    add_arrow(inhibitoryLine, position=-10, direction='left', size=25, color='blue')
    add_arrow(inhibitoryLine, position=-5, direction='left', size=25, color='blue')
    # ax.arrow(-20, excitatory_current[np.where(voltages == -20)], 5, excitatory_current[np.where(voltages == -15)] - excitatory_current[np.where(voltages == -20)], shape='full', lw=0, length_includes_head=True, head_width=.5, color='blue')
    fixedPoint = fixedPoints[timeVector.index(recordingTiming)]
    for point in fixedPoint:
      ax.plot(voltages[point], excitatory_current[point], marker = 'o', markeredgecolor=fixedPointColors[np.where(fixedPoint == point)[0][0]],markeredgewidth=8, markerfacecolor='none', markersize = 26.0)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['left'].set_linewidth(2)
  ax.spines['bottom'].set_linewidth(2)
  ax.tick_params(direction='out', width=2, size=10)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')  
  plt.savefig('figure_03_a.pdf' , transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
  plt.close()


def plot_figure_03_b(controlValuesOverTime, controlValuesOverVoltage): 
  fixedPoints = calculateFixedPoints(controlValuesOverVoltage)
  keys = controlValuesOverVoltage.keys()
  T = sorted(controlValuesOverVoltage[keys[0]].keys())
  fixedPointVoltages = []
  fixedPointTimes = []
  numOfPaths = [len(fixedPoint) for fixedPoint in fixedPoints]
  indices = [i for i, x in enumerate(np.diff(numOfPaths)) if x <> 0]
  indices.append(len(fixedPoints) - 1)
  voltageGraph = []
  for tInd in range(len(T)):
    voltageGraph.append(controlValuesOverTime['VOLTAGE'][T[tInd]])
  startPath = 0
  fixedPointTrajectories = {}
  for p in indices:
    fixedPointVoltages.append([])
    fixedPointTimes.append([])
    for path in range(len(fixedPoints[p])):
      if (path not in fixedPointTrajectories.keys()):
        fixedPointTrajectories[path] = []
      fixedPointVoltages[-1].append([])
      fixedPointTimes[-1].append([])
      for t in range(startPath, p + 1):
        fixedPointTrajectories[path].append((voltages[fixedPoints[t][path]]))
        fixedPointVoltages[-1][-1].append((voltages[fixedPoints[t][path]]))
        fixedPointTimes[-1][-1].append(T[t] * dt)
    startPath = p + 1
  fig = plt.figure(figsize=(6, 6))
  ax = fig.add_subplot(111)
  ax.axis([90, 210, -90, 10], fontsize=24.0)
  plt.xticks([100, 150, 200], [0, 50, 100], fontsize=24.0)
  plt.yticks([-80, -45, -10], fontsize=24.0)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['left'].set_linewidth(2)
  ax.spines['bottom'].set_linewidth(2)
  ax.tick_params(direction='out', width=2, size=10)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')   
  plt.xlabel('time (msec)', weight='regular', fontsize=24.0)
  plt.ylabel('voltage (mV)', weight='regular', fontsize=24.0)
  ax.plot(controlValuesOverTime['TIME'], controlValuesOverTime['VOLTAGE'], color ='black', linewidth = 2, linestyle='-', alpha=1)
  voltageGraph = np.array(voltageGraph)
  ax.fill_between(np.array(T) * dt, voltageGraph, fixedPointTrajectories[0], alpha = 0.25, color = 'blue')
  for difference in range(len(fixedPointTimes)):
    for path in range(len(fixedPointTimes[difference])):
      ax.plot(fixedPointTimes[difference][path], fixedPointVoltages[difference][path], color=fixedPointColors[path], linewidth = 5, linestyle='-', alpha=1)
    if (len(fixedPointTimes[difference]) == 3):
      beginTimeFixedPoints = np.where(np.round(T,3) == (round(fixedPointTimes[difference][0][0] / dt, 3)))[0][0]
      endTimeFixedPoints = np.where(np.round(T,3) == (round(fixedPointTimes[difference][0][-1] / dt, 3)))[0][0]
      startTimeTermination = np.where(np.round(controlValuesOverTime['TIME'],3) == (round(fixedPointTimes[difference][0][-1], 3)))[0][0]
      ax.fill_between(fixedPointTimes[difference][0], fixedPointVoltages[difference][2], fixedPointVoltages[difference][1], alpha = 1, color = 'white')
      ax.fill_between(fixedPointTimes[difference][0], fixedPointVoltages[difference][2], fixedPointVoltages[difference][1], alpha = 0.25, color = 'red')

  plt.savefig('figure_03_b.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
  plt.close()

def plot_figure_03_c_d(controlValuesOverTime, controlValuesOverVoltage, inhibitionValuesOverTime, inhibitionValuesOverVoltage, Dt, figure_letter):
  timeInds = np.array([99 + Dt, 105 + Dt, 125 + Dt]) / dt
  fixedPoints = calculateFixedPoints(controlValuesOverVoltage)
  keys = controlValuesOverVoltage.keys()
  timeVector = sorted(controlValuesOverVoltage[keys[0]].keys())
  for recordingTiming in timeInds:    
    timeInMillis = recordingTiming * dt
    NMDA_current = np.array(inhibitionValuesOverVoltage['NMDA'][recordingTiming])
    AMPA_current = np.array(inhibitionValuesOverVoltage['AMPA'][recordingTiming])
    GABA_current = - np.array(inhibitionValuesOverVoltage['GABA'][recordingTiming])
    PAS_current =  - np.array(inhibitionValuesOverVoltage['PAS'][recordingTiming])
    CAP_current =  - np.array(inhibitionValuesOverVoltage['CAP'][recordingTiming])
    ABS_CAP_current =  - np.abs(np.array(inhibitionValuesOverVoltage['CAP'][recordingTiming]))
    inhibitory_current = GABA_current + PAS_current 
    excitatory_current = NMDA_current
    voltagePoint = inhibitionValuesOverTime['VOLTAGE'][recordingTiming]
    vclampValues = np.array(inhibitionValuesOverVoltage['vclamp'][recordingTiming])
    
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(voltages, excitatory_current, color ='red', linewidth = 4)
    ax.plot(voltages, PAS_current, color ='blue', linewidth = 4, alpha = 0.5, linestyle='--')
    ax.plot(voltages, inhibitory_current, color ='blue', linewidth = 4, alpha = 1, linestyle='-')
    
    fixedPoint = fixedPoints[timeVector.index(recordingTiming)]
    for point in fixedPoint:
      ax.plot(voltages[point], excitatory_current[point], marker = 'o', markeredgecolor=fixedPointColors[np.where(fixedPoint == point)[0][0]],markeredgewidth=6, markerfacecolor='none', markersize = 18.0)
    ax.plot(voltages[np.where(voltages == round(voltagePoint))[0][0]], excitatory_current[np.where(voltages == round(voltagePoint))[0][0]], marker = 'o', markerfacecolor='white',markeredgewidth=6, markeredgecolor='black', markersize = 18)
        
    ax.axis([-90, 0, -0.06, 0.001], fontsize=24.0)
    plt.xticks([-80, -45, -10], [-80, -45, -10], fontsize=24.0)
    plt.yticks([-0.01, -0.03, -0.05], [-10, -30, -50], fontsize=24.0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(direction='out', width=2, size=10)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')      
    ax.set_xlabel('mV', weight="regular", fontsize=24.0)
    ax.set_ylabel('pA', weight="regular", fontsize=24.0)
    ax.text(-80, 0.01,'t = %d msec' % (timeInMillis - 100), bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'), fontsize=24)
    plt.savefig('figure_03_%s_top_%f_Dt_%f.pdf' % (figure_letter, timeInMillis, Dt) , transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
    plt.close()
    
  fixedPoints = calculateFixedPoints(controlValuesOverVoltage)
  keys = controlValuesOverVoltage.keys()
  T = sorted(controlValuesOverVoltage[keys[0]].keys())
  fixedPointVoltages = []
  fixedPointTimes = []
  numOfPaths = [len(fixedPoint) for fixedPoint in fixedPoints]
  indices = [i for i, x in enumerate(np.diff(numOfPaths)) if x <> 0]
  indices.append(len(fixedPoints) - 1)
  voltageGraph = []
  for tInd in range(len(T)):
    voltageGraph.append(max(voltages[fixedPoints[tInd]]))
  startPath = 0
  for p in indices:
    fixedPointVoltages.append([])
    fixedPointTimes.append([])
    for path in range(len(fixedPoints[p])):
      fixedPointVoltages[-1].append([])
      fixedPointTimes[-1].append([])
      for t in range(startPath, p + 1):
        fixedPointVoltages[-1][-1].append((voltages[fixedPoints[t][path]]))
        fixedPointTimes[-1][-1].append(T[t] * dt)
    startPath = p + 1
  for recordingTiming in timeInds:
    timeInMillis = recordingTiming * dt
    timeInInds = np.where(np.round(inhibitionValuesOverTime['TIME'],3) == (round(recordingTiming * dt, 3)))[0][0]
    fig, ax = plt.subplots(figsize=(4,4))
    ax.axis([90, 210, -90, 10], fontsize=24.0)
    plt.xticks([100, 150, 200], [0, 50, 100], fontsize=24.0)
    plt.yticks([-80, -45, -10], [-80, -45, -10], fontsize=24.0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(direction='out', width=2, size=10)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')  
    plt.xlabel('time (msec)', weight='regular', fontsize=24.0)
    plt.ylabel('voltage (mV)', weight='regular', fontsize=24.0)
    voltageGraph = np.array(voltageGraph)
    for difference in range(len(fixedPointTimes)):
      for path in range(len(fixedPointTimes[difference])):
        ax.plot(fixedPointTimes[difference][path], fixedPointVoltages[difference][path], color=fixedPointColors[path], linewidth = 5, linestyle='-', alpha=1)
    ax.plot(controlValuesOverTime['TIME'], controlValuesOverTime['VOLTAGE'], color ='black', linewidth = 2, linestyle='--', alpha=0.75)
    ax.plot(inhibitionValuesOverTime['TIME'][0:timeInInds], inhibitionValuesOverTime['VOLTAGE'][0:timeInInds], color ='black', linewidth = 2, linestyle='-', alpha=1)
    if (recordingTiming == timeInds[-1]):
      ax.plot(inhibitionValuesOverTime['TIME'], inhibitionValuesOverTime['VOLTAGE'], color ='black', linewidth = 2, linestyle='-', alpha=1)
    ax.plot(inhibitionValuesOverTime['TIME'][timeInInds], inhibitionValuesOverTime['VOLTAGE'][timeInInds], marker = 'o', markerfacecolor='white',markeredgewidth=6, markeredgecolor='black', markersize = 18)
    ax.arrow(100 + Dt, 5, 0,  -8, head_width=4, width=1.5, head_length=2, fc="black", ec="black")
    plt.savefig('figure_03_%s_bottom_%f_Dt_%f.pdf' % (figure_letter, timeInMillis, Dt), transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
    plt.close()

    

dt = 0.025
voltages = (np.array(np.linspace(-90,0,901)) / 1.0);
recordingTimings = range(int(80 / dt), int(250 / dt), int(0.025 / dt))
fixedPointColors = ['blue','green','red']
finalStop = 300
timeAfterClamp = 0.1 / dt

controlFilename = 'results_for_depolarStrength_0.007000_hyperStrength_0.000000_shift_-1000.000000.pickle'
(controlValuesOverTime, controlValuesOverVoltage) = pickle.load(open(controlFilename, 'rb'))

plot_figure_03_a(controlValuesOverVoltage, controlValuesOverTime)
plot_figure_03_b(controlValuesOverTime, controlValuesOverVoltage)

inhibitionFilename_for_Dt_10 = 'results_for_depolarStrength_0.007000_hyperStrength_0.001500_shift_10.000000.pickle'
(inhibitionValuesOverTime_for_Dt_10, inhibitionValuesOverVoltage_for_Dt_10) = pickle.load(open(inhibitionFilename_for_Dt_10, 'rb'))

plot_figure_03_c_d(controlValuesOverTime, controlValuesOverVoltage, inhibitionValuesOverTime_for_Dt_10, inhibitionValuesOverVoltage_for_Dt_10, 10, 'c')

inhibitionFilename_for_Dt_40 = 'results_for_depolarStrength_0.007000_hyperStrength_0.001500_shift_40.000000.pickle'
(inhibitionValuesOverTime_for_Dt_40, inhibitionValuesOverVoltage_for_Dt_40) = pickle.load(open(inhibitionFilename_for_Dt_40, 'rb'))

plot_figure_03_c_d(controlValuesOverTime, controlValuesOverVoltage, inhibitionValuesOverTime_for_Dt_40, inhibitionValuesOverVoltage_for_Dt_40, 40, 'd')



