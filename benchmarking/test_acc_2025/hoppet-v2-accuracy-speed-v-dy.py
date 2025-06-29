#!/usr/bin/env python3
"""Matplotlib template generated automatically with 
/Users/gsalam/scripts/mptemplate.py accuracy-speed-v-dy.py
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages
from   matplotlib.ticker import ScalarFormatter
from   matplotlib.ticker import FuncFormatter
from   matplotlib import cycler
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import re
import copy
import hfile as h
from glob import glob

styles = [
    {'color':'#f00000'},
    {'color':'#0000e0'},
    {'color':'#00c000'},
    {'color':'#404040'},
    {'color':'#e0a040'},
    ]
colors = cycler('color', [style['color'] for style in styles])
# see options for things to set with 
#     python -c "import matplotlib as mpl; print(mpl.rcParams)"
plt.rc('axes', prop_cycle=colors)
mpl.rcParams.update({"axes.grid" : True, "grid.linestyle": ":"})
plt.rc('figure', figsize=(5,5))

dir="../../build/pp/"

def main(pdf):
    #res = h.get_array_plus_comments("filename",regexp='',columns={'x':1,'y':(3,4)}) 
    run_stats_pre   = RunStats('nloop3-preev-dy*.dat')
    run_stats_nopre = RunStats('nloop3-nopreev-dy*.dat')

    fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=True)
    fig.subplots_adjust(hspace=0.04)
    ax2.set_xlabel(r'dy')

    mask = run_stats_pre.dy > 0.035

    ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_allf_xlt09[mask])
    ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_guds_xlt07[mask])

    ax2.plot(run_stats_pre  .dy[mask], run_stats_pre  .t_ev_s[mask], label='M2Pro, cached evolution'   )
    ax2.plot(run_stats_nopre.dy[mask], run_stats_nopre.t_ev_s[mask], label='M2Pro, one-off evolution')
    ax1.set_ylabel("accuracy [all flavs, $x<0.9$]")
    ax2.set_ylabel("time [s]")

    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.tick_params(axis='both', which='both', left=True, right=True, direction='in')
    ax2.tick_params(axis='both', which='both', left=True, right=True, direction='in')


    xticks_major = ax2.get_xticks().tolist()
    extra_xticks = [0.05, 0.2]
    all_xtics = sorted(xticks_major + extra_xticks)
    ax2.set_xticks(all_xtics)

    ax2.set_xticklabels(f"{xt}" for xt in all_xtics)
    ax2.set_xlim(0.029,0.31)
    ax2.yaxis.set_major_formatter(FuncFormatter(h.log_formatter_fn))

    ax1.text(0.03,0.93, "Hoppet v2.0.0, NNLO\nymax = 12, dlnlnQ = dy/4", va='top', transform=ax1.transAxes)
    #ax1.text(0.03,0.86, "", va='top', transform=ax1.transAxes)

    # ax.set_xlim(0,2.5)
    # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    # ax.tick_params(top=True,right=True,direction='in',which='both')
    # ax.set_xscale('log')


    # ax.grid(True,ls=":")
    # ax.set_title("title")
    # ax.text(x,y,'hello',transform=ax.transAxes)
    #ax.plot(res.x, res.y, label='label', **styles[0])
    ax2.legend(loc='lower left')
    pdf.savefig(fig,bbox_inches='tight')
    #pdf.savefig(fig,bbox_inches=Bbox.from_extents(0.0,0.0,7.5,4.8))
    plt.close()

class RunStats(object):
    '''Class to extra run stats for all files matching a certain glob pattern (over dy values)'''
    def __init__(self,glob_str):
        self.files_timing = sorted(glob(f"{dir}/{glob_str}"))
        self.files_prec   = sorted(glob(f"{dir}/{glob_str}.acc.smry"))
        self.n = len(self.files_timing)
        self.dy = np.array([
          float(re.findall(r"dy(0\.[0-9]+)",f)[0]) for f in self.files_timing
        ])

        # get the timing numbers
        self.t_init_s = get_number(self.files_timing,"Initialisation time")
        self.t_preev_s = get_number(self.files_timing,"Pre-evolution time")
        self.t_ev_s = get_number(self.files_timing,"Evolution time")
        self.t_eval_all_flav_ns = get_number(self.files_timing,"Evaluation..all")
        self.t_eval_one_flav_ns = get_number(self.files_timing,"Evaluation..one")

        # get the accuracy numbers from the prec files
        self.accuracies = np.zeros([self.n,4])
        for i, f in enumerate(self.files_prec):
            self.accuracies[i,:] = h.get_array(f,"max.rel.error")[0,:]
        self.acc_guds_xlt07 = self.accuracies[:,0]
        self.acc_allf_xlt07 = self.accuracies[:,1]
        self.acc_guds_xlt09 = self.accuracies[:,2]
        self.acc_allf_xlt09 = self.accuracies[:,3]
    

        #print(self.accuracies)

def get_number(filename,regex):
    if isinstance(filename,list):
        return np.array([get_number(f,regex) for f in filename])
    
    with open(filename,'r') as f:
        for l in f:
            matches = re.findall(regex+r".*?([0-9.]+)",l)
            if matches: return float(matches[0])

def line_and_band(ax,x,val_and_err,**extra):
    extra_no_label = copy.copy(extra)
    if ('label' in extra_no_label): del extra_no_label['label']
    ax.fill_between(x,
                    val_and_err.value-val_and_err.error,
                    val_and_err.value+val_and_err.error,
                    alpha=0.2,
                    **extra_no_label
                    )
    ax.plot(x,val_and_err.value, **extra)

if __name__ == "__main__": 
    pdfname = __file__.replace('.py','.pdf')
    with PdfPages(pdfname) as pdf: main(pdf)


