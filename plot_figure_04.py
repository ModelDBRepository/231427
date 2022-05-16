
#!/usr/bin/python

## A point neuron model for testing inhibition on NMDA spikes

from neuron import h
from neuron import gui
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import copy
import numpy as np
import sys
import pickle
import time
from scipy.stats.stats import pearsonr
from scipy.stats.stats import linregress
from scipy.optimize import curve_fit

def get_threshold(results):
  nmdaConds = sorted(results.keys())
  gabaConds = sorted(results[np.float(np.round(nmdaConds[0], 6))].keys())
  keySecNames = sorted(results[np.float(np.round(nmdaConds[0], 6))][np.float(np.round(gabaConds[0],6))].keys())
  delays = sorted(results[np.float(np.round(nmdaConds[0], 6))][np.float(np.round(gabaConds[0],6))][keySecNames[0]].keys())
  thresholds = {}
  for sec in keySecNames: 
    integrals = []
    for nmdaAmp in nmdaConds:
      integrals.append(results[nmdaAmp][gabaConds[0]][sec][delays[0]])
    diffs = np.diff(integrals)
    threshold = (nmdaConds[np.where(diffs ==  np.max(diffs))[0][0]]) if len(np.where(diffs == np.max(diffs))[0]) > 0 else -1
    threshold = threshold if (threshold > nmdaConds[0] and threshold < nmdaConds[-2]) else -1
    thresholds[sec] = np.round(threshold * 20, 5)
  return thresholds


def get_termination_threshold_for_nmda_cond(results, sec, nmda_cond, delay):
  if nmda_cond == -1: return -1
  inhibition_integrals = []
  gabaConds = sorted(results[nmda_cond].keys())
  for gabaCond in gabaConds:
    inhibition_integrals.append(results[nmda_cond][gabaCond][sec][delay])
  inhibition_diffs = np.diff(inhibition_integrals)    
  inhibition_threshold = (np.where((inhibition_diffs) == np.min(inhibition_diffs))[0][-1]) if len(np.where((inhibition_diffs) == np.min(inhibition_diffs))[0]) > 0 else -1
  inhibition_threshold = inhibition_threshold if (inhibition_threshold < (len(inhibition_diffs) - 1) and inhibition_threshold > 0) else -1  
  return inhibition_threshold
  



def plot_figure_04_a():
    startTime = time.time()
    h.load_file("nrngui.hoc")
    h.load_file("import3d.hoc")

    h.tstop = 250
    h.v_init = -64

    h("forall nseg = 1")
    h("forall e_pas = -70")
    h("objref eSynlist, ePreconlist, iPreconlist, eStimlist, iStimlist, iSynlist")
    eSynlist = {}
    ePreconlist = {}
    iPreconlist = {}
    eStimlist = {}
    iStimlist = {}
    iSynlist = {}

    def placeNMDA(location, sec):
      eStimlist[sec].append(h.NetStim())
      eStimlist[sec][-1].interval = 1
      eStimlist[sec][-1].number = 1
      eStimlist[sec][-1].start = 100
      eStimlist[sec][-1].noise = 0
      eSynlist[sec].append(h.ProbAMPANMDA2_RATIO(location))
      eSynlist[sec][-1].gmax = 0
      eSynlist[sec][-1].mgVoltageCoeff = 0.08
      ePreconlist[sec].append(h.NetCon(eStimlist[sec][-1], eSynlist[sec][-1]))
      ePreconlist[sec][-1].weight[0] = 1
      ePreconlist[sec][-1].delay = 0

    def placeGABA(location, sec):
      iStimlist[sec].append(h.NetStim())
      iStimlist[sec][-1].interval = 1
      iStimlist[sec][-1].number = 1
      iStimlist[sec][-1].start = 100
      iStimlist[sec][-1].noise = 0       
      iSynlist[sec].append(h.ProbUDFsyn2_lark(location))
      iSynlist[sec][-1].tau_r = 0.18
      iSynlist[sec][-1].tau_d = 5
      iSynlist[sec][-1].e = - 80
      iSynlist[sec][-1].Dep = 0
      iSynlist[sec][-1].Fac = 0
      iSynlist[sec][-1].Use = 0.6
      iSynlist[sec][-1].u0 = 0
      iSynlist[sec][-1].gmax = 0
      iPreconlist[sec].append(h.NetCon(iStimlist[sec][-1], iSynlist[sec][-1]))
      iPreconlist[sec][-1].weight[0] = 1
      iPreconlist[sec][-1].delay = 0


    delaysVoltageVector = {}
    delaysVoltageNumpy = {}

    resultsDict = {}
    resultsList = []

    secs = {'low' : 'L5PC.dend[2]', 'high' : 'L5PC.apic[52]'}
    states = ['low', 'high']
    nmdaAmp = 0.0004

    for state in states:
        h("access " + secs[state])
        eSynlist[secs[state]] = []
        ePreconlist[secs[state]] = []
        iPreconlist[secs[state]] = []
        eStimlist[secs[state]] = []
        iStimlist[secs[state]] = []
        iSynlist[secs[state]] = []
        excitatory_locations = np.arange(h.L - 20, h.L, 1) / h.L
        inhibitory_locations = [excitatory_locations[10]]
        h("nseg = " + str(h.L))

        for i in range(1):
          for location in excitatory_locations:
            placeNMDA(location, secs[state])
          for location in inhibitory_locations:
            placeGABA(location, secs[state])


    for state in states:
        for other_state in states:
            for sec in iSynlist[secs[other_state]]:
                sec.gmax = 0.001
            for sec in eSynlist[secs[other_state]]:
                sec.gmax = 0.0004
        
        for sec in iSynlist[secs[state]]:
            sec.gmax = 0.001
        for sec in eSynlist[secs[state]]:
            sec.gmax = 0.0004
        for stim in iStimlist[secs[state]]:
            stim.start = 110

        h.finitialize()

        delaysVoltageNumpy[secs[state]] = {}

        im = h.Impedance()
        h("access " + secs[state])
        im.loc(0.93)
        im.compute(0)

        voltageTraces = []
        for gabaAmp in [0, 0.001]:
          for syn in iSynlist[secs[state]]:
            syn.gmax = gabaAmp
          voltageVector = h.Vector()
          voltageVector.record(eval("h." + secs[state] + "(0.93)._ref_v"))
          timeVector = h.Vector()
          timeVector.record(h._ref_t)
          h.init()
          h.run()
          voltageTraces.append(np.array(voltageVector))
          del voltageVector

        fig = plt.figure(figsize=(8,8)) 
        ax = fig.add_subplot(111)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)  
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')       
        ax.axis([95, 200, -80, 10], width=4)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.plot(timeVector, voltageTraces[0], linewidth = 6, color='white', linestyle='--', dashes = (3,3))
        plt.plot(timeVector, voltageTraces[1], linewidth = 6, color='white')
        plt.text(100, voltageTraces[0][110 / 0.025] + 25, '10 msec', color='white', fontsize = 50)
        plt.arrow(110, voltageTraces[0][110 / 0.025] + 20, 0, -10, head_width=6, width=2, head_length=3, fc='white', ec='white')
        plt.fill_between(timeVector, voltageTraces[1], -70, facecolor='white', where=voltageTraces[1] > -70, alpha=1)
        plt.fill_between(timeVector, voltageTraces[1], -70, facecolor='red', where=voltageTraces[1] > -70, alpha=0.8)
        ax.plot([170, 180],[-7, -7],linewidth=6,color='white')
        ax.plot([180, 180],[-7, 3],linewidth=6,color='white')
        ax.text(168, -17, "10 msec", fontsize = 50, color='white')
        ax.text(183, -6, "10 mV", fontsize = 50, color='white')
        fig.savefig('figure_04_a_%s.pdf' % state, transparent=True, bbox_inches='tight', format='pdf', dpi=1000, facecolor='black')
        plt.close()


def plot_figure_04_b():
    nmdaConds = sorted(results.keys())
    gabaConds = sorted(results[nmdaConds[0]].keys())
    colors=['black','red','blue']
    colorInd = 0
    secNames = sorted(results[nmdaConds[0]][gabaConds[0]].keys())
    inputResistance_dict = {}

    fig = plt.figure(figsize=(8,8)) 
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(direction='out', width=2, size=10)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    colors = {10:'blue', 25:'red', 40:'green'}
    for delay in [10, 25, 40]: 
        inh_strong_inputResistance = []
        inhibitory_strong_thresolds = []
        inh_diams = []
        inh_surfaces = []
        inh_lengths = []
        inh_integrals = []
        inh_distanceFromSoma = []
        im = h.Impedance()
        h("access L5PC.soma[0]")
        h.distance(0)
        h.pop_section()  
        for sec in secNames: 
            secModified = sec.replace('TTC[0]','L5PC')
            h("access " + secModified)
            h("nseg = " + str(h.L))
            excitatory_locations = np.arange(h.L - 20, h.L, 1) / h.L
            seg = excitatory_locations[10]
            inhibitory_thresold = get_termination_threshold_for_nmda_cond(results, secModified, nmda_cond = 0.008 / 20, delay = delay)
            if(inhibitory_thresold < 0): continue
            Rins = []
            for seg in excitatory_locations:
                im.loc(seg)
                im.compute(0)
                Rins.append(im.input(seg))
            im.loc(excitatory_locations[10])
            im.compute(0)
            Rins = [im.input(excitatory_locations[10])]
            inh_strong_inputResistance.append(1000.0 / np.mean(Rins))
            inhibitory_strong_thresolds.append(gabaConds[inhibitory_thresold] * 1000)
            h.pop_section()
        
        y_vector = np.array(copy.deepcopy(inhibitory_strong_thresolds))
        x_vector = np.array(copy.deepcopy(inh_strong_inputResistance))
        x_vector, y_vector = zip(*sorted(zip(x_vector, y_vector)))
        plt.plot(x_vector, y_vector,  linewidth=0, marker='o', color=colors[delay], markeredgewidth=0, markersize = 10, label=r'$\mathregular{\Delta}$t = ' + str(delay) + ' ms')
        regress = linregress(x_vector, y_vector)
        slope_line = regress.slope * np.array(x_vector) + regress.intercept
        plt.plot(x_vector, slope_line, '-', color=colors[delay], label='_nolegend_')
    
    ax.set_xlabel('Input conductance\nat branch (nS)', fontsize = 34, fontweight = 'regular')
    ax.set_ylabel(r'GABA$\mathregular{_A}$ necessary' + '\nfor NMDA termination (nS)', fontsize = 34, fontweight = 'regular')
    ax.set_xticks([0.6, 1.0, 1.4])
    ax.set_xticklabels([0.6, 1.0, 1.4], fontsize=28)
    ax.set_yticks([0.0, 1.0, 2.0])
    ax.set_yticklabels([0.0, 1.0, 2.0], fontsize=28)
    plt.show(block = 0)
    fig.savefig('figure_04_b.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
    plt.close()


def plot_figure_04_c():
    h("access L5PC.soma[0]")
    im = h.Impedance()
    h.distance(0)
    h.pop_section()
    delays = np.arange(1,100,1)

    nmdaConds = sorted(results.keys())
    gabaConds = sorted(results[nmdaConds[0]].keys())

    ratio_inputResistance = []
    ratio_lengths = []
    ratio_diams = []
    ratio_distanceFromSoma = []
    ratio_surfaces = []
    ratio = []
    pearson = []
    for nmda_cond in nmdaConds:
        for delay in delays:
            ratio_inputResistance.append([])
            ratio_lengths.append([])
            ratio_diams.append([])
            ratio_distanceFromSoma.append([])
            ratio_surfaces.append([])
            ratio.append([])
            pearson.append([])
            for sec in secNames:
                inhibitory_thresold = get_termination_threshold_for_nmda_cond(results, sec, nmda_cond = nmda_cond, delay = delay)
                # inhibitory_thresold = get_termination_threshold_for_nmda_cond(results, sec, nmda_cond = 0.008 / 20, delay = delay)
                if (inhibitory_thresold < 0): continue
                secModified = sec.replace('TTC[0]','L5PC')
                h("access " + secModified)
                h("nseg = " + str(h.L))
                excitatory_locations = np.arange(h.L - 20, h.L, 1.0) / h.L
                inhibitory_locations = [excitatory_locations[10]]
                seg = excitatory_locations[10]
                
                im.loc(seg)
                im.compute(0)
                ratio_inputResistance[-1].append(1000.0 / im.input(seg))
                ratio_diams[-1].append(eval("h." + secModified + ".diam"))
                ratio_lengths[-1].append(eval("h." + secModified + ".L"))
                ratio_distanceFromSoma[-1].append(h.distance(seg))
                ratio_surfaces[-1].append((2 * np.pi * eval("h." + secModified + ".diam") * eval("h." + secModified + ".L")))
                ratio[-1].append(gabaConds[inhibitory_thresold] * 1000)
                h.pop_section()
            pearson[-1].append(pearsonr(ratio_inputResistance[-1], ratio[-1])[0])
            pearson[-1].append(pearsonr(ratio_diams[-1], ratio[-1])[0])
            pearson[-1].append(pearsonr(ratio_lengths[-1], ratio[-1])[0])
            pearson[-1].append(pearsonr(ratio_surfaces[-1], ratio[-1])[0])
            pearson[-1].append(pearsonr(ratio_distanceFromSoma[-1], ratio[-1])[0])

    nonan_pearson = []
    for p in pearson:
        if (np.isnan(p) == False).all():
            nonan_pearson.append(p)

    pearson = np.array(nonan_pearson)
    x_axis_labels = ['\nInput\nconductance', 'Diameter', 'Length', '\nSurface\narea', '\nDistance\nfrom soma']
    mean_pearson = np.mean(np.array(pearson), axis=0)
    std_pearson = np.std(np.array(pearson), axis=0)
    mean_pearson, x_axis_labels = zip(*sorted(zip(mean_pearson, x_axis_labels)))  

    fig = plt.figure(figsize=(8,8)) 
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(direction='out', width=2, size=10)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    # y_vector = np.array(copy.deepcopy(y_vector_arg))
    # x_vector = np.array(copy.deepcopy(x_vector_arg))
    # x_vector, y_vector = zip(*sorted(zip(x_vector, y_vector)))
    a = plt.bar(left = range(0, len(mean_pearson) * 2, 2), height = mean_pearson, width=1, linewidth = 3, edgecolor= 'black', color='white', yerr=std_pearson, error_kw={'capsize' : 5, 'ecolor' : 'black', 'elinewidth' : 3, 'capthick' : 3})
    # a = plt.errorbar(x = range(0, len(mean_pearson) * 2, 2), y = mean_pearson, ecolor = 'black', capsize = 5, color='white', yerr=std_pearson)
    ax.xaxis.set_ticks(np.arange(0.5, len(mean_pearson) * 2 + 0.5, 2))
    ax.xaxis.set_ticklabels(x_axis_labels, fontsize = 24)
    ax.yaxis.set_ticks([-1, -0.5, 0, 0.5, 1])
    ax.yaxis.set_ticklabels([-1, -0.5, 0, 0.5, 1], fontsize = 24)
    ax.axis([-1, 10, -1.5, 1.5])
    ax.set_ylabel('Pearson Coefficient', fontsize = 32, fontweight = 'regular')
    # ax.set_title(r'Correlation between NMDA termination threshold' + '\nand the passive properties of the branch', fontsize = 32, fontweight = 'regular')
    plt.show(block = 0)
    fig.savefig('figure_04_c.pdf', transparent=True, bbox_inches='tight', format='pdf', dpi=1000)
    plt.close()  


results = pickle.load(open('data_same_excitation.pickle','rb'))

h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")
h.load_file("L5PCbiophys3_noActive.hoc")
h.load_file("TTC.hoc")
h("objref L5PC")
h.L5PC = h.TTC("cell1.asc")

secNames = []
for sec in h.L5PC.all:
    a = h.SectionRef(sec=sec)
    if (a.nchild() == 0.0) and ((int(eval("h." + sec.name() + ".L"))) > 20):
        secNames.append(sec.name().replace('TTC[0]','L5PC'))

secNames = sorted(secNames)

plot_figure_04_a()
plot_figure_04_b()
plot_figure_04_c()