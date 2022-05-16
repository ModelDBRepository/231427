#!/usr/bin/python

## A point neuron model for testing inhibition on NMDA spikes

import matplotlib
matplotlib.use('Agg')
from neuron import h
from neuron import gui
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle
import time
import os

dt = 0.025
voltages = (np.array(np.linspace(-90,0,901)) / 1.0);
depolarStrength = 0.007
depolarLength = 1
hyperStrength = 0
hyperLength = 1
shift = -1000
recordingTimings = range(int(80 / dt), int(250 / dt), int(0.025 / dt))
fixedPointColors = ['blue','green','red']
finalStop = 300
timeAfterClamp = 0.1 / dt

startTime = time.time()
h.load_file("nrngui.hoc")

h.tstop = finalStop
h.dt = dt;
h.v_init = -80 
timeDiff = 40;

h("create soma")
h("access soma")
h("nseg = 1")
h("L = 20")
h("diam = 20")
h("insert pas")
h("cm = 1")
h("Ra = 100")
h("cm = 1")
h("forall nseg = 1")
h("forall e_pas = -80")
h("objref resistanceVector")
h("objref eSynlist, ePreconlist, iPreconlist, eStimlist, iStimlist, iSynlist, iclamp, voltageClamp, resistanceVector")
h("g_pas = 0.00005")


h.voltageClamp = h.SEClamp(0.5)
h.iclamp = h.IClamp(0.5)
eSynlist = []
ePreconlist = []
iPreconlist = []
eStimlist = []
iStimlist = []
iSynlist = []
h.resistanceVector = h.Vector(int(h.tstop/h.dt))


def placeGABA(location):
  iStimlist.append(h.NetStim())
  iStimlist[-1].interval = 1
  iStimlist[-1].number = 1
  iStimlist[-1].start = 100
  iStimlist[-1].noise = 0                   
  
  iSynlist.append(h.ProbUDFsyn2_lark(location))
  iSynlist[-1].tau_r = 0.18
  iSynlist[-1].tau_d = 5
  iSynlist[-1].e = - 80
  iSynlist[-1].Dep = 0
  iSynlist[-1].Fac = 0
  iSynlist[-1].Use = 0.6
  iSynlist[-1].u0 = 0
  iSynlist[-1].gmax = 0
  
  iPreconlist.append(h.NetCon(iStimlist[-1], iSynlist[-1]))
  iPreconlist[-1].weight[0] = 1
  iPreconlist[-1].delay = 0


def placeNMDA(location):
  eStimlist.append(h.NetStim())
  eStimlist[-1].interval = 1
  eStimlist[-1].number = 1
  eStimlist[-1].start = 100
  eStimlist[-1].noise = 0
  eSynlist.append(h.ProbAMPANMDA2_RATIO(location))
  eSynlist[-1].gmax = 0
  eSynlist[-1].mgVoltageCoeff = 0.08
  
  ePreconlist.append(h.NetCon(eStimlist[-1], eSynlist[-1]))
  ePreconlist[-1].weight[0] = 1
  ePreconlist[-1].delay = 0
  

placeNMDA(0.5)
placeGABA(0.5)


def IVCurveAtTime(depolarStrength, depolarLength, hyperStrength, hyperLength, shift, recordingTiming):
  recordingTimingInMilliseconds = recordingTiming * h.dt
  recordingTimingInInds = recordingTiming
  h.voltageClamp.dur1 = 10000000
  h.voltageClamp.rs = 0.0000001
  res = {'vclamp' : {},'NMDA' : {},'AMPA' : {},'PAS' : {},'GABA' : {},'CAP' : {}}
  tVector = np.zeros(h.tstop / h.dt)
  timeStart = time.time()
  res['vclamp'][recordingTimingInInds] = {}
  res['CAP'][recordingTimingInInds] = {}
  res['NMDA'][recordingTimingInInds] = {}
  res['AMPA'][recordingTimingInInds] = {}
  res['PAS'][recordingTimingInInds] = {}
  res['GABA'][recordingTimingInInds] = {}
  for j in xrange(int(h.tstop / h.dt)-1):   
    h.resistanceVector.x[j] = 1e9
  
  for j in np.arange(recordingTimingInInds, len(tVector), 1):   
    h.resistanceVector.x[int(j)] = 0.0000001
  
  eSynlist[-1].gmax = depolarStrength
  iSynlist[-1].gmax = hyperStrength
  iStimlist[-1].start = eStimlist[-1].start + shift
  
  h("resistanceVector.play_remove()")
  h("resistanceVector.play(&voltageClamp.rs, dt)")      
  for voltage in voltages:
    neuronNMDAVectors = {'vclamp' : h.Vector(), 'NMDA' : h.Vector(), 'AMPA' : h.Vector(), 'PAS' : h.Vector(), 'GABA' : h.Vector(), 'CAP' : h.Vector()}
    vecNMDA = {'vclamp' : [], 'NMDA' : [], 'AMPA' : [], 'PAS' : [], 'GABA' : [], 'CAP' : []}
    neuronNMDAVectors['vclamp'].record(h.voltageClamp._ref_i)
    neuronNMDAVectors['GABA'].record(iSynlist[-1]._ref_i)
    neuronNMDAVectors['AMPA'].record(eSynlist[-1]._ref_i_AMPA)
    neuronNMDAVectors['NMDA'].record(eSynlist[-1]._ref_i_NMDA)
    neuronNMDAVectors['PAS'].record(h.soma(0.5)._ref_i_pas)
    neuronNMDAVectors['CAP'].record(h.soma(0.5)._ref_i_cap)
    h.voltageClamp.amp1 = voltage
    h.tstop = ((recordingTimingInInds + timeAfterClamp + 1) * h.dt)
    h.run()
    res['vclamp'][recordingTimingInInds][voltage] = {}
    res['GABA'][recordingTimingInInds][voltage] = {}
    res['NMDA'][recordingTimingInInds][voltage] = {}
    res['AMPA'][recordingTimingInInds][voltage] = {}
    res['PAS'][recordingTimingInInds][voltage] = {}
    res['CAP'][recordingTimingInInds][voltage] = {}
    for key in vecNMDA.keys():
      if key == 'vclamp': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])
        res['vclamp'][recordingTimingInInds][voltage] = np.array(vecNMDA[key][-1])[recordingTimingInInds + timeAfterClamp]
      elif key == 'NMDA': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])   
        res['NMDA'][recordingTimingInInds][voltage] =  (np.array(vecNMDA[key][-1]))[recordingTimingInInds + timeAfterClamp]
      elif key == 'AMPA': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])   
        res['AMPA'][recordingTimingInInds][voltage] =  (np.array(vecNMDA[key][-1]))[recordingTimingInInds + timeAfterClamp]
      elif key == 'PAS': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])   
        res['PAS'][recordingTimingInInds][voltage] =  (np.array(vecNMDA[key][-1]) * h.area(0.5) / 100.0)[recordingTimingInInds + timeAfterClamp]
      elif key == 'GABA': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])   
        res['GABA'][recordingTimingInInds][voltage] =  (np.array(vecNMDA[key][-1]))[recordingTimingInInds + timeAfterClamp]
      elif key == 'CAP': 
        vecNMDA[key].append(h.Vector())
        vecNMDA[key][-1].copy(neuronNMDAVectors[key])   
        res['CAP'][recordingTimingInInds][voltage] =  (np.array(vecNMDA[key][-1]) * h.area(0.5) / 100.0)[recordingTimingInInds - 1]
             
  h.tstop = finalStop
  for j in xrange(int(h.tstop / h.dt)):   
    h.resistanceVector.x[j] = 1e9
  
  h("resistanceVector.play_remove()")
  h("resistanceVector.play(&voltageClamp.rs, dt)")
  
  neuronRealVectors = {'voltage' : h.Vector(), 'NMDA_total' : h.Vector(), 'AMPA_total' : h.Vector(), 'pas_total' : h.Vector(), 'cap_total' : h.Vector(), 'GABA_total' : h.Vector(), 'time' : h.Vector()}
  vecReal = {'voltage' : [], 'NMDA_total' : [], 'AMPA_total' : [], 'GABA_total' : [], 'pas_total' : [], 'cap_total' : [], 'time' : []}
  neuronRealVectors['voltage'].record(h.soma(0.5)._ref_v)
  neuronRealVectors['GABA_total'].record(iSynlist[-1]._ref_i)  
  neuronRealVectors['NMDA_total'].record(eSynlist[-1]._ref_i_NMDA)  
  neuronRealVectors['AMPA_total'].record(eSynlist[-1]._ref_i_AMPA)  
  neuronRealVectors['pas_total'].record(h.soma(0.5)._ref_i_pas)
  neuronRealVectors['cap_total'].record(h.soma(0.5)._ref_i_cap)
  neuronRealVectors['time'].record(h._ref_t)  
  h.run()
  vecReal['voltage'].append(h.Vector())
  vecReal['voltage'][-1].copy(neuronRealVectors['voltage'])
  res['voltage'] = np.array(vecReal['voltage'][-1])
  
  vecReal['NMDA_total'].append(h.Vector())
  vecReal['NMDA_total'][-1].copy(neuronRealVectors['NMDA_total'])
  res['NMDA_total'] = np.array(vecReal['NMDA_total'][-1])
  
  vecReal['AMPA_total'].append(h.Vector())
  vecReal['AMPA_total'][-1].copy(neuronRealVectors['AMPA_total'])
  res['AMPA_total'] = np.array(vecReal['AMPA_total'][-1])
  
  vecReal['GABA_total'].append(h.Vector())
  vecReal['GABA_total'][-1].copy(neuronRealVectors['GABA_total'])
  res['GABA_total'] = np.array(vecReal['GABA_total'][-1])
    
  vecReal['pas_total'].append(h.Vector())
  vecReal['pas_total'][-1].copy(neuronRealVectors['pas_total'])
  res['pas_total'] = np.array(vecReal['pas_total'][-1]) * h.area(0.5) / 100.0
  
  vecReal['cap_total'].append(h.Vector())
  vecReal['cap_total'][-1].copy(neuronRealVectors['cap_total'])
  res['cap_total'] = np.array(vecReal['cap_total'][-1]) * h.area(0.5) / 100.0
    
  vecReal['time'].append(h.Vector())
  vecReal['time'][-1].copy(neuronRealVectors['time'])
  res['time'] = np.array(vecReal['time'][-1])
  
  valuesOverTime = {'VOLTAGE' : res['voltage'], 'NMDA' : res['NMDA_total'], 'AMPA' : res['AMPA_total'], 'GABA' : res['GABA_total'], 'PAS' : res['pas_total'], 'CAP' : res['cap_total'], 'TIME' : res['time']}
  valuesOverVoltage = {'NMDA' : res['NMDA'], 'AMPA' : res['AMPA'], 'GABA' : res['GABA'], 'PAS' : res['PAS'], 'CAP' : res['CAP'], 'vclamp' : res['vclamp']}
  timeEnd = time.time()
  print 'recorded values for time %f in %f seconds' % (recordingTimingInMilliseconds, timeEnd - timeStart)
  return (valuesOverTime, valuesOverVoltage)



valuesOverTime = {}
valuesOverVoltage = {}

for recordingTiming in recordingTimings:
  (overTime, overVoltage) = IVCurveAtTime(depolarStrength, depolarLength, hyperStrength, hyperLength, shift, recordingTiming)
  for key in overTime:
    valuesOverTime[key] = overTime[key]
  for key in overVoltage:
    if (key not in valuesOverVoltage.keys()):
      valuesOverVoltage[key] = {}
    valuesOverVoltage[key][recordingTiming] = []
    for voltage in sorted(voltages):
      valuesOverVoltage[key][recordingTiming].append(overVoltage[key][recordingTiming][voltage])

filename = 'results_for_depolarStrength_%f_hyperStrength_%f_shift_%f.pickle' % (depolarStrength, hyperStrength, shift)
pickle.dump((valuesOverTime, valuesOverVoltage), open(filename, 'wb'),protocol=2)
