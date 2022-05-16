#!/usr/bin/python

from neuron import h
from neuron import gui
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import sys
import pickle
import time

startTime = time.time()
h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")

h.tstop = 250
h.v_init = -64

sectionNum = 44;
h.load_file("L5PCbiophys3_noActive.hoc")
h.load_file("TTC.hoc")
h("objref L5PC")
h.L5PC = h.TTC("cell1.asc")
h("forall nseg = 1")
h("forall e_pas = -70")
h("access L5PC.apic[" + str(sectionNum) + "]")
eSynlist = []
ePreconlist = []
iPreconlist = []
eStimlist = []
iStimlist = []
iSynlista = []

def placeNMDA(location):
  eStimlist.append(h.NetStim())
  eStimlist[-1].interval = 1
  eStimlist[-1].number = 1
  eStimlist[-1].start = 100
  eStimlist[-1].noise = 0
  eSynlist.append(h.ProbAMPANMDA2_RATIO(location))
  eSynlist[-1].gmax = 0.0004
  eSynlist[-1].mgVoltageCoeff = 0.08
  ePreconlist.append(h.NetCon(eStimlist[-1], eSynlist[-1]))
  ePreconlist[-1].weight[0] = 1
  ePreconlist[-1].delay = 0

def placeGABA(location):
  iStimlist.append(h.NetStim())
  iStimlist[-1].interval = 1
  iStimlist[-1].number = 1
  iStimlist[-1].start = 100
  iStimlist[-1].noise = 0       
  
  iSynlista.append(h.ProbUDFsyn2_lark(location))
  iSynlista[-1].tau_r = 0.18
  iSynlista[-1].tau_d = 5
  iSynlista[-1].e = - 80
  iSynlista[-1].Dep = 0
  iSynlista[-1].Fac = 0
  iSynlista[-1].Use = 0.6
  iSynlista[-1].u0 = 0
  iSynlista[-1].gmax = 0.0007
  
  iPreconlist.append(h.NetCon(iStimlist[-1], iSynlista[-1]))
  iPreconlist[-1].weight[0] = 1
  iPreconlist[-1].delay = 0


voltageVector = h.Vector()
timeVector = h.Vector()
tempVector = h.Vector()
delaysVoltageVector = {}
delaysVoltageNumpy = {}
voltageVector.record(h.L5PC.apic[sectionNum](0.5)._ref_v)
timeVector.record(h._ref_t)
h("access L5PC.apic[" + str(sectionNum) + "]")
h("nseg = 50")
for i in range(1):
  for location in np.linspace(0.2,0.7,20):
    placeNMDA(location)
  for location in np.linspace(0.50,0.50,1):      
    placeGABA(location)

h.finitialize()
h.run()
vec = h.Vector()
vec.copy(voltageVector)

delayDiff = 1
delays = np.array([-100000] + list(np.arange(-20, 60, delayDiff)))
gabaAmps = [0.0005, 0.001, 0.0015]

delaysVoltageNumpy = {}
for gabaAmp in gabaAmps:
  for syn in iSynlista:
    syn.gmax = gabaAmp
  
  delaysVoltageNumpy[gabaAmp] = {}
  for delay in delays:
    for stim in iStimlist:
      stim.start = 100 + delay
    h.run()
    delaysVoltageVector[delay] = h.Vector()
    delaysVoltageVector[delay].copy(voltageVector)
    delaysVoltageNumpy[gabaAmp][delay] = np.array(delaysVoltageVector[delay])

pickle.dump(delaysVoltageNumpy, open('figure_01_data.pickle', 'wb'),protocol=2)

