#!/usr/bin/python

from neuron import h
from neuron import gui
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
mpl.rc('font', **font)
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

import numpy as np
import sys
import pickle
import time
import itertools
import numpy as np

def indexes(y, thres=0.3, min_dist=1):
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")
    thres = thres * (np.max(y) - np.min(y)) + np.min(y)
    min_dist = int(min_dist)
    dy = np.diff(y)
    zeros,=np.where(dy == 0)
    
    while len(zeros):
        zerosr = np.hstack([dy[1:], 0.])
        zerosl = np.hstack([0., dy[:-1]])

        dy[zeros]=zerosr[zeros]
        zeros,=np.where(dy == 0)

        dy[zeros]=zerosl[zeros]
        zeros,=np.where(dy == 0)

    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]
    return peaks

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
  iSynlista[-1].gmax = 0.001
  
  iPreconlist.append(h.NetCon(iStimlist[-1], iSynlista[-1]))
  iPreconlist[-1].weight[0] = 1
  iPreconlist[-1].delay = 0

def plot_figure_06_b_c():
  startTime = time.time()
  h.tstop = 300
  h.v_init = -70
  h.dt = 0.025

  delays = [0, 90, 100, 110, 120, 130]

  ampaCond = 0.008
  gabaCond = 0.0015
  seg = 0.9

  secNames = []
  secRins = []
  for sec in h.L5PC.all:
    a = h.SectionRef(sec=sec)
    if (a.nchild() == 0) and ('dend' in sec.name()):
      im = h.Impedance()
      h("access " + sec.name())
      im.loc(seg)
      im.compute(0)
      Ri = im.input(seg) 
      if (Ri > 400):
        secNames.append(sec.name())
        secRins.append(Ri)

  secRins, excitable_secs = zip(*sorted(zip(secRins, secNames)))
  num_of_branches = 16

  sectionNum = "L5PC.soma[0]"

  no_zero = True

  excitable_secs = ['TTC[0].dend[82]', 'TTC[0].dend[69]', 'TTC[0].dend[80]',
         'TTC[0].dend[66]', 'TTC[0].dend[10]', 'TTC[0].dend[66]',
         'TTC[0].dend[82]', 'TTC[0].dend[66]', 'TTC[0].dend[37]',
         'TTC[0].dend[54]', 'TTC[0].dend[31]', 'TTC[0].dend[72]',
         'TTC[0].dend[54]', 'TTC[0].dend[49]', 'TTC[0].dend[76]',
         'TTC[0].dend[13]']

  for syn in eSynlist:
    syn.gmax = 0
    del syn
  for syn in iSynlista:
    syn.gmax = 0
    del syn
  for stim in iStimlist:
    del stim 
  for sec in h.L5PC.all:
    h("nseg = 1")
  for sec in excitable_secs:
    h("access " + str(sec))
    h("nseg = 10")
    placeNMDA(seg)
    placeGABA(seg)
    h.pop_section()
  voltageNpVec = {}
  soma_voltageNpVec = {}
  for delay in delays:
    inhibitionTiming = delay
    voltageVec = h.Vector()
    soma_voltageVec = h.Vector()
    timeVec = h.Vector()
    for syn in eSynlist:
      syn.gmax = ampaCond
    for syn in iSynlista:
      syn.gmax = gabaCond
    for stim in iStimlist:
      stim.start = inhibitionTiming
    voltageVec.record(eval("h." + excitable_secs[-14] + "(" + str(seg) + ")._ref_v"))
    soma_voltageVec.record(eval("h." + sectionNum + "(" + str(seg) + ")._ref_v"))
    timeVec.record(h._ref_t)
    h.init()
    h.run()
    soma_voltageNpVec[delay] = np.array(soma_voltageVec)
    voltageNpVec[delay] = np.array(voltageVec)
    print len(indexes(soma_voltageNpVec[delay], thres=np.abs(10 - np.min(soma_voltageNpVec[delay])) / (np.abs(np.max(soma_voltageNpVec[delay]) - np.min(soma_voltageNpVec[delay]))), min_dist=10))
  len_orinial_spike_indices = len(indexes(soma_voltageNpVec[0], thres=np.abs(10 - np.min(soma_voltageNpVec[0])) / (np.abs(np.max(soma_voltageNpVec[0]) - np.min(soma_voltageNpVec[0]))), min_dist=10))
  len_spike_indices = len(indexes(soma_voltageNpVec[120], thres=np.abs(10 - np.min(soma_voltageNpVec[120])) / (np.abs(np.max(soma_voltageNpVec[120]) - np.min(soma_voltageNpVec[120]))), min_dist=10))
  if (len_spike_indices == 0 and len_orinial_spike_indices == 2): 
    no_zero = False

  colors = {90: 'black', 100: 'green', 110 : 'brown', 120 : 'blue', 130 : 'red'}

  for delay in delays[1:]:
    fig = plt.figure(figsize=(4, 6), frameon = False)
    ax = plt.Axes(fig, [0., 0., 1., 1.], )
    fig.add_axes(ax)

    ax.plot(timeVec, soma_voltageNpVec[delay], linewidth = 0.75, color='black');
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)  
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')       
    ax.axis([80, 200, -90, 60])
    ax.set_xticks([])
    ax.set_yticks([])
    if (delay == 130):
      ax.plot([150, 170],[-18, -18],linewidth=0.75,color='black')
      ax.plot([170, 170],[-18, 2],linewidth=0.75,color='black')
      ax.text(148, -44, "20 msec", fontsize = 8)
      ax.text(173, -16, "20 mV", fontsize = 8)
    fig.set_size_inches(0.6, 0.6)
    fig.savefig('figure_06_b_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=3000, pad_inches=0)
    plt.show(block = 0)

  for delay in delays[1:]:
    fig = plt.figure(figsize=(4, 6), frameon = False)
    ax = plt.Axes(fig, [0., 0., 1., 1.], )
    fig.add_axes(ax)
    ax.plot([delay, delay],[voltageNpVec[0][delay / 0.025], 0],linewidth=0.5,color='black',linestyle=':')
    ax.plot(timeVec, voltageNpVec[delay], linewidth = 0.75, color='green');
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)  
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')       
    ax.axis([80, 200, -90, 30])
    ax.set_xticks([])
    ax.set_yticks([])
    if (delay == 130):
      ax.plot([150, 170],[-18, -18],linewidth=0.75,color='black')
      ax.plot([170, 170],[-18, 2],linewidth=0.75,color='black')
      ax.text(148, -44, "20 msec", fontsize = 8)
      ax.text(173, -16, "20 mV", fontsize = 8)
    plt.text(delay - 14, 20, '%d msec' % (delay - 100), color='green', fontsize = 8)
    ax.arrow(delay, 10, 0, -10, head_width=1.25, width=0.75, head_length=0.75, fc='green', ec='green')
    fig.set_size_inches(0.6, 0.6)
    fig.savefig('figure_06_c_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=3000, pad_inches=0)
  # plt.close('all')
    


def plot_figure_06_d():
  eps = np.finfo(float).eps
  (avg, std) = pickle.load(open('spikes_num.pickle', 'rb'))
  delays = [0] + range(90, 161, 1)

  fig = plt.figure(figsize=(1.7,1.7))
  ax = plt.Axes(fig, [0., 0., 1., 1.], )
  fig.add_axes(ax)

  ax.plot([delays[1], delays[-1]], [avg[0],avg[0]], linewidth = 1, color="black", linestyle='--', dashes=(3,3)); 
  ax.plot(delays[1:], avg[1:], linewidth = 1, color="black"); 
  ax.plot(delays[1:], avg[1:] + std[1:], linewidth = 1, color="green"); 
  ax.plot(delays[1:], avg[1:] - std[1:], linewidth = 1, color="green"); 
  ax.fill_between(delays[1:], avg[1:] + std[1:], avg[1:] - std[1:],facecolor='green', alpha=0.2)

  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['left'].set_linewidth(1)
  ax.spines['bottom'].set_linewidth(1)
  ax.tick_params(direction='out', width=1, size=5)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')   
  ax.set_xlabel(r'${\Delta t}$ Inhibition vs excitation (msec)', fontsize = 10, fontweight = 'regular')
  ax.set_ylabel('Number of spikes', fontsize = 10, fontweight = 'regular')
  ax.axis([90, 160, -0.1, 3.25], fontsize = 24)
  kwargs = dict(linewidth = 1, transform=ax.transAxes, color='k', clip_on=False)
  ax.xaxis.set_ticks([100, 120, 140, 160])
  ax.xaxis.set_ticklabels([0, 20, 40, 60], fontsize=8)
  ax.yaxis.set_ticks([0,1,2,3])
  ax.yaxis.set_ticklabels([0,1,2,3], fontsize=8)
  plt.subplots_adjust(hspace = 0.1) 
  fig.set_size_inches(1.5, 1.5)
  fig.savefig('figure_06_d.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)


h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")
h.load_file("L5PCbiophys3.hoc")
h.load_file("TTC.hoc")

h("objref L5PC")
h("celsius = 34.5")
h.L5PC = h.TTC("cell1.asc")
h("forall nseg = 1")
h("forall vo = 1")
h("objref eSynlist, ePreconlist, iPreconlist, eStimlist, iStimlist, iSynlista, voltageClamp, resistanceVector")

plot_figure_06_b_c()
plot_figure_06_d()