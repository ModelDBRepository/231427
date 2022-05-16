#!/usr/bin/python

import matplotlib
matplotlib.use('Agg') 
from neuron import h
from neuron import gui
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle
import time

h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")
h.load_file("cellmorphology.hoc")
h.load_file("iniparameter.hoc")

h('''objref all_spines
all_spines = new SectionList()''')

h("forall insert cadiffus")


def place_NMDA(nmdaCond):
  nmda_eStimlist.append(h.NetStim())
  nmda_eStimlist[-1].interval = 1
  nmda_eStimlist[-1].number = 1
  nmda_eStimlist[-1].start = 100
  nmda_eStimlist[-1].noise = 0
  nmda_eSynlist.append(h.nmda(0.5))
  nmda_ePreconlist.append(h.NetCon(nmda_eStimlist[-1], nmda_eSynlist[-1]))
  nmda_ePreconlist[-1].weight[0] = nmdaCond
  nmda_ePreconlist[-1].delay = 0

def place_AMPA(ampaCond):
  ampa_eStimlist.append(h.NetStim())
  ampa_eStimlist[-1].interval = 1
  ampa_eStimlist[-1].number = 1
  ampa_eStimlist[-1].start = 100
  ampa_eStimlist[-1].noise = 0
  ampa_eSynlist.append(h.ampa(0.5))
  ampa_ePreconlist.append(h.NetCon(ampa_eStimlist[-1], ampa_eSynlist[-1]))
  ampa_ePreconlist[-1].weight[0] = ampaCond
  ampa_ePreconlist[-1].delay = 0

def place_GABA(gabaCond):
  iStimlist.append(h.NetStim())
  iStimlist[-1].interval = 1
  iStimlist[-1].number = 1
  iStimlist[-1].start = 100000
  iStimlist[-1].noise = 0                   
  iSynlist.append(h.ProbUDFsyn2_lark( (branch_length - (number_of_spines / 2)) / branch_length))
  iSynlist[-1].tau_r = 0.18
  iSynlist[-1].tau_d = 5
  iSynlist[-1].e = - 80
  iSynlist[-1].Dep = 0
  iSynlist[-1].Fac = 0
  iSynlist[-1].Use = 0.6
  iSynlist[-1].u0 = 0
  iSynlist[-1].gmax = gabaCond
  iPreconlist.append(h.NetCon(iStimlist[-1], iSynlist[-1]))
  iPreconlist[-1].weight[0] = 1
  iPreconlist[-1].delay = 0


locs = ['basal','apical']
branches = ['dendA2_01100','dendA5_01111111111111111010']

number_of_spines = 20

for loc_ind in range(len(locs)):
  branch = branches[loc_ind]

  h('access ' + branch)
  branch_length = h.L
  h('nseg = %d ' % branch_length)
  center_of_distribution =  (branch_length - (number_of_spines / 2))

  h.pop_section()
  spine_segs = np.arange((center_of_distribution) - (number_of_spines / 2), (center_of_distribution) + (number_of_spines / 2), 1) / branch_length

  for seg in spine_segs:
    h('''{create %s_spine_neck_%s}'''  % (locs[loc_ind], str(seg).replace('.','')))
    h('''{%s connect %s_spine_neck_%s(0),%f}''' % (branch, locs[loc_ind], str(seg).replace('.',''),seg))
    h('''{access %s_spine_neck_%s}''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''necklength = 1  /*spine neck length in um*/''')
    h('''neckdiam = 0.1 /*spine neck diameter*/''')
    h('''%s_spine_neck_%s {nseg = 2
                pt3dclear()
                for j = 0, nseg-1 {
                    ty = (j*necklength)/(nseg-1)
                    pt3dadd(97.58,ty,20.87,neckdiam)
                    }
                }''' % (locs[loc_ind], str(seg).replace('.','')))

    h('''{create %s_spine_head_%s}''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''{%s_spine_neck_%s connect %s_spine_head_%s(0), 1}''' % (locs[loc_ind], str(seg).replace('.',''), locs[loc_ind], str(seg).replace('.','')))
    h('''{access %s_spine_head_%s}''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''spineradius = 0.297''')
    h('''%s_spine_head_%s {nseg = 7
                pt3dclear()
                for i = 0, nseg-1 {
                    ty = -(i*2*spineradius)/(nseg-1)
                    td = 2*sqrt(spineradius^2-(ty+spineradius)^2)
                    if (td<neckdiam){
                        td = neckdiam
                        } 
                        pt3dadd(97.58,ty+1,20.87,td)
                    }
                }''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''%s_spine_head_%s all_spines.append''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''%s_spine_neck_%s all_spines.append''' % (locs[loc_ind], str(seg).replace('.','')))
    h('''
    %s_spine_head_%s {
        insert pas  e_pas=Vleak  g_pas=spinefactor/Rm  Ra=global_ra  cm=spinefactor*Cm 
        insert canmda
  //      insert car
        insert cadiffus
    }
    %s_spine_neck_%s{
        insert pas  e_pas=Vleak  g_pas=0.00005  Ra=global_ra  cm=spinefactor*Cm 
        insert cadiffus
        
    }''' % (locs[loc_ind], str(seg).replace('.',''),locs[loc_ind], str(seg).replace('.','')))
    h('access %s_spine_neck_%s' % (locs[loc_ind], str(seg).replace('.','')))
    h('%s_spine_neck_%s.diam = 0.0394' % (locs[loc_ind], str(seg).replace('.','')))
    h.pop_section()


  h('''initchannels()

  proc init() {

       /* add initchannels() to init(), so parameter changes show up */
       initchannels()       
       finitialize(v_init)
       fcurrent()
  }
  ''')

  h.tstop = 300
  h.dt = 0.025;
  h.v_init = -65
  h.celsius = 34

  nmda_eSynlist = []
  nmda_ePreconlist = []
  nmda_eStimlist = []

  ampa_eSynlist = []
  ampa_ePreconlist = []
  ampa_eStimlist = []

  iPreconlist = []
  iStimlist = []
  iSynlist = []

  single_spine_AMPA_nS = 0.8
  single_spine_NMDA_nS = 0.8
  single_spine_GABA_nS = 0.5

  nmdaCond =  (single_spine_NMDA_nS / (1000 * 0.000045))
  ampaCond =  (single_spine_AMPA_nS / (1000 * 0.00001))
  gabaCond = (single_spine_GABA_nS * 0.001)

  middle =  (spine_segs)[(number_of_spines / 2)]

  for seg in spine_segs:
    h('access %s_spine_head_%s' % (locs[loc_ind], str(seg).replace('.','')))
    place_NMDA(nmdaCond)
    place_AMPA(ampaCond)
    im = h.Impedance()
    im.loc(0.5)
    im.compute(0)
    Rin = im.input(0.5)
    h.pop_section()

  inhibitory_locations = ['branch','spine']

  h('access ' + branch)
  h('nseg = %d ' % branch_length)
  place_GABA(gabaCond)
  h.pop_section()

  h('access %s_spine_head_%s' % (locs[loc_ind], str(middle).replace('.','')))
  place_GABA(gabaCond)
  h.pop_section()

  time_traces = {}

  delays = [-100000] + range(80, 160, 1)
  # spine_neck_diams = np.linspace(0.0394, 0.1, 3)
  spine_neck_diams = [0.0394]

  spine_head_voltage_traces = {}
  spine_head_current_traces = {}
  spine_head_cai_traces = {}

  colors = {-100000 : 'black', 110 : 'green', 120 : 'blue', 130 : 'red'}


  y_labels = {}
  y_labels['spine_head_voltage_traces'] = 'voltage (mV)'
  y_labels['spine_head_current_traces'] = 'pA'
  y_labels['spine_head_cai_traces'] = r'$\mathregular{[Ca^{2+}}$] (mM)'

  font = {'weight':'regular', 'size':8}
  matplotlib.rc('font', **font)
  for spine_neck_diam in spine_neck_diams:
    for spine_seg in spine_segs:
      h('%s_spine_neck_%s.diam=%f' % (locs[loc_ind], str(spine_seg).replace('.',''), spine_neck_diam))
    
    spine_head_voltage_traces[spine_neck_diam] = {}
    spine_head_current_traces[spine_neck_diam] = {}
    spine_head_cai_traces[spine_neck_diam] = {}
    
    mean_spine_head_voltage_traces_integral_full = []
    mean_spine_head_cai_traces_integral_full = []
    local_spine_head_voltage_traces_integral_full = []
    local_spine_head_cai_traces_integral_full = []

    for inhibitory_location in inhibitory_locations:
      spine_head_voltage_traces[spine_neck_diam][inhibitory_location] = {}
      spine_head_current_traces[spine_neck_diam][inhibitory_location] = {}
      spine_head_cai_traces[spine_neck_diam][inhibitory_location] = {}

      for delay in delays:
        spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay] = []
        spine_head_current_traces[spine_neck_diam][inhibitory_location][delay] = []
        spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay] = []
        
        time_trace = h.Vector()
        for spine_seg in spine_segs:
          spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay].append(h.Vector())
          spine_head_current_traces[spine_neck_diam][inhibitory_location][delay].append(h.Vector())
          spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay].append(h.Vector())
        
        time_trace.record(h._ref_t)
        for ind in range(len(spine_segs)):
          spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay][ind].record(eval('h.%s_spine_head_%s' % (locs[loc_ind], str((spine_segs)[ind]).replace('.','')))(0.5)._ref_v)
          spine_head_current_traces[spine_neck_diam][inhibitory_location][delay][ind].record(nmda_eSynlist[ind]._ref_i)
          spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay][ind].record(eval('h.%s_spine_head_%s' % (locs[loc_ind], str((spine_segs)[ind]).replace('.','')))(0.5)._ref_cai)
        
        for stim in iStimlist:
          stim.start = delay
        for syn in iSynlist:
          syn.gmax = 0
        
        iSynlist[inhibitory_locations.index(inhibitory_location)].gmax = gabaCond
        
        h.init()
        h.run()
        
        spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay] = np.array(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay])
        spine_head_current_traces[spine_neck_diam][inhibitory_location][delay] = np.array(spine_head_current_traces[spine_neck_diam][inhibitory_location][delay])
        spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay]     = np.array(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay])
      
      mean_spine_head_voltage_traces_integral = []
      mean_spine_head_cai_traces_integral = []
      local_spine_head_voltage_traces_integral = []
      local_spine_head_cai_traces_integral = []
      

      for delay in delays[1:]:
        mean_spine_head_voltage_traces_integral.append(np.mean((np.sum(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay] + 65, axis=1)) / (np.sum(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delays[0]] + 65, axis=1)), axis=0))
        mean_spine_head_cai_traces_integral.append(np.mean((np.sum(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay], axis=1)) / (np.sum(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delays[0]], axis=1)), axis=0))
        local_spine_head_voltage_traces_integral.append((np.sum(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)] + 65)) / (np.sum(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delays[0]][(number_of_spines / 2)] + 65)))
        local_spine_head_cai_traces_integral.append((np.sum(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)])) / (np.sum(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delays[0]][(number_of_spines / 2)])))

      mean_spine_head_voltage_traces_integral_full.append(mean_spine_head_voltage_traces_integral)
      mean_spine_head_cai_traces_integral_full.append(mean_spine_head_cai_traces_integral)
      local_spine_head_voltage_traces_integral_full.append(local_spine_head_voltage_traces_integral)
      local_spine_head_cai_traces_integral_full.append(local_spine_head_cai_traces_integral)

      if (locs[loc_ind] == 'apical') and (spine_neck_diam == spine_neck_diams[0]):
        for delay in [-100000, 120, 130]:
          fig = plt.figure(figsize=(1.7, 1.7))
          ax = fig.add_subplot(111)
          
          plt.plot(time_trace, spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delays[0]][(number_of_spines / 2)], linewidth=1, color='black',linestyle='--', dashes=(3,3))
          plt.plot(time_trace, spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)], linewidth=1, color='black')
          
          ax.fill_between(time_trace, spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)], -65, where = (np.mean(spine_head_voltage_traces[spine_neck_diam][inhibitory_location][delay], axis=0) > -65), color=colors[delay], alpha=0.3)
          ax.yaxis.set_ticks([-70, -50, -30, -10])
          ax.yaxis.set_ticklabels([-70, -50, -30, -10])
          plt.text(delay + 10, 2.5, '%d msec' % (delay - 100), color = colors[delay], fontsize=8)
          ax.arrow(delay, 10, 0, -7.5, head_width=6, width=2, head_length=2, fc=colors[delay], ec=colors[delay])      
          plt.axis([time_trace[0], time_trace[-1], -70, 10])
          
          ax.spines['top'].set_visible(False)
          ax.spines['right'].set_visible(False)
          ax.spines['left'].set_visible(1)
          ax.spines['bottom'].set_visible(1)  
          ax.spines['left'].set_linewidth(1)
          ax.spines['bottom'].set_linewidth(1)
          ax.tick_params(direction='out', width=1, size=5)
          ax.yaxis.set_ticks_position('left')
          ax.xaxis.set_ticks_position('bottom')       
          ax.xaxis.set_ticks([0,100,200,300])
          ax.xaxis.set_ticklabels([-100, 0, 100, 200])
          ax.set_xlabel('time (msec)', fontsize=8)
          ax.set_ylabel('voltage (mV)', fontsize=8)
          fig.set_size_inches(1.2, 1.2)
          if inhibitory_location == 'branch':
            plt.savefig('figure_05_b_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
          if inhibitory_location == 'spine':
            plt.savefig('figure_05_e_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
          plt.close('all')
        
          fig = plt.figure(figsize=(1.7, 1.7))
          ax = fig.add_subplot(111)
          
          plt.plot(time_trace, spine_head_cai_traces[spine_neck_diam][inhibitory_location][delays[0]][(number_of_spines / 2)], linewidth=1, color='black',linestyle='--', dashes=(3,3))
          plt.plot(time_trace, spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)], linewidth=1, color='black')
          
          ax.fill_between(time_trace, spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay][(number_of_spines / 2)], 0, where = (np.mean(spine_head_cai_traces[spine_neck_diam][inhibitory_location][delay], axis=0) > -65), color=colors[delay], alpha=0.3)
          
          plt.axis([time_trace[0], time_trace[-1], 0, 1])
          
          ax.spines['top'].set_visible(False)
          ax.spines['right'].set_visible(False)
          ax.spines['left'].set_visible(1)
          ax.spines['bottom'].set_visible(1)  
          ax.spines['left'].set_linewidth(1)
          ax.spines['bottom'].set_linewidth(1)
          ax.tick_params(direction='out', width=1, size=5)
          ax.yaxis.set_ticks_position('left')
          ax.xaxis.set_ticks_position('bottom')       
          ax.xaxis.set_ticks([0,100,200,300])
          ax.xaxis.set_ticklabels([-100, 0, 100, 200])
          ax.set_xlabel('time (msec)', fontsize=8)
          ax.set_ylabel(r'$\mathregular{[Ca^{2+}}$] (mM)', fontsize=8)
          fig.set_size_inches(1.2, 1.2)
          if inhibitory_location == 'branch':
            plt.savefig('figure_05_c_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
          if inhibitory_location == 'spine':
            plt.savefig('figure_05_f_%d.pdf' % (delay), transparent=True, bbox_inches='tight', format='pdf', dpi=300, pad_inches=0)
          plt.close('all')    

  if (locs[loc_ind] == 'apical'):
    font = {'weight':'regular', 'size':8, 'family':'sans-serif', 'sans-serif':'Arial'}
    matplotlib.rc('font', **font)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={'height_ratios':[20,1]}, sharex=True)
    for ax in [ax1, ax2]:
      ax.plot(delays[1:], local_spine_head_voltage_traces_integral_full[0], linewidth = 1, color="red", label='Inhibition on dendrite'); 
      ax.plot(delays[1:], local_spine_head_voltage_traces_integral_full[1], linewidth = 1, color="green", label='Inhibition on spine'); 
      
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
      ax.spines['left'].set_linewidth(1)
      ax.spines['bottom'].set_linewidth(1)
      ax.tick_params(direction='out', width=1, size=5)
      ax.yaxis.set_ticks_position('left')
      ax.xaxis.set_ticks_position('bottom')   
      ax.xaxis.set_ticks(np.array([-15,0, 15, 30, 45, 60]) + 100)
      ax.xaxis.set_ticklabels([-15,0, 15, 30, 45, 60], fontsize=10)
    ax2.set_xlabel(r'$\mathregular{\Delta t}$ Inhibition vs excitation (msec)', fontsize = 10, fontweight = 'regular')
    ax1.set_ylabel('Time Integral of \nNMDA Spike (Norm.)', fontsize = 10, fontweight = 'regular')
    ax1.axis([delays[1], 160, 0.3625, 1.01])
    ax2.axis([delays[1], 160, 0, 0.1])
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
    ax1.yaxis.set_label_coords(- 0.2, 0.4)
    fig.set_size_inches(3, 3)
    fig.savefig('figure_05_g.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300)
    plt.close('all')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), gridspec_kw={'height_ratios':[20,1]}, sharex=True)
    for ax in [ax1, ax2]:
      ax.plot(delays[1:], local_spine_head_cai_traces_integral_full[0], linewidth = 1, color="red", label='Inhibition on dendrite'); 
      ax.plot(delays[1:], local_spine_head_cai_traces_integral_full[1], linewidth = 1, color="green", label='Inhibition on spine'); 
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
      ax.spines['left'].set_linewidth(1)
      ax.spines['bottom'].set_linewidth(1)
      ax.tick_params(direction='out', width=1, size=5)
      ax.yaxis.set_ticks_position('left')
      ax.xaxis.set_ticks_position('bottom')   
      ax.xaxis.set_ticks(np.array([-15,0, 15, 30, 45, 60]) + 100)
      ax.xaxis.set_ticklabels([-15,0, 15, 30, 45, 60], fontsize=10)
    ax2.set_xlabel(r'$\mathregular{\Delta t}$ Inhibition vs excitation (msec)', fontsize = 10, fontweight = 'regular')
    ax1.set_ylabel('Time Integral of \n' + r'$\mathregular{[Ca^{2+}]}$ (Norm.)', fontsize = 10, fontweight = 'regular')
    ax1.axis([delays[1], 160, 0.3625, 1.01])
    ax2.axis([delays[1], 160, 0, 0.1])
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
    ax1.yaxis.set_label_coords(- 0.2, 0.4)
    fig.set_size_inches(3, 3)
    fig.savefig('figure_05_h.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=300)
    plt.close('all')
