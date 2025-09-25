#!/usr/bin/env python3
"""Matplotlib template generated automatically with 
/Users/gsalam/scripts/mptemplate.py accuracy-speed-v-dy.py
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages
from   matplotlib.ticker import LogLocator, ScalarFormatter
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
    {'color':'#404040'},
    {'color':'#0000e0'},
    {'color':'#00c000'},
    {'color':'#e0a040'},
    ]
colors = cycler('color', [style['color'] for style in styles])
# see options for things to set with 
#     python -c "import matplotlib as mpl; print(mpl.rcParams)"
plt.rc('axes', prop_cycle=colors)
mpl.rcParams.update({"axes.grid" : True, "grid.linestyle": ":"})
plt.rc('figure', figsize=(5,5))

dirM2pro      ="../../../2025-prec-and-timing/M2Pro-gfortran15.1-O3-2025-09-18"
dirM2proLHAPDF="../../../2025-prec-and-timing/M2Pro-gfortran15.1-O3-2025-09-18-LHAPDF"
nloop_names = {1: "LO", 2: "NLO", 3: "NNLO", 4: "N3LO"}

def main(pdf,nloop):
    
    #res = h.get_array_plus_comments("filename",regexp='',columns={'x':1,'y':(3,4)}) 
    run_stats_pre   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-dy*.dat')
    run_stats_nopre = RunStats(f'{dirM2pro}/nloop{nloop}-nopreev-dy*.dat')

    fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=True,height_ratios=[3,2])
    fig.subplots_adjust(hspace=0.04)
    ax2.set_xlabel(r'dy')

    mask = run_stats_pre.dy > 0.035

    ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_guds_xlt07[mask], label='guds, $x<0.7$'    , **styles[0], ls="-")
    ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_allf_xlt07[mask], label='all-flav, $x<0.7$', **styles[0], ls="--")
    #ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_guds_xlt09[mask], label='guds, $x<0.9$'    , **styles[0], ls=":")

    ax2.plot(run_stats_nopre.dy[mask], run_stats_nopre.t_init_s[mask], **styles[1], label='setup ($t_s$)')
    ax2.plot(run_stats_nopre.dy[mask], run_stats_nopre.t_ev_s  [mask], **styles[2], label='one-off evol. ($t_i$)')
    ax2.plot(run_stats_pre  .dy[mask], run_stats_pre  .t_ev_s  [mask], **styles[3], label='cached evol. ($t_c$)'   )
    ax1.set_ylabel("rel. accuracy")
    ax2.set_ylabel("time [s]")

    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.tick_params(axis='both', which='both', left=True, right=True, direction='in')
    ax2.tick_params(axis='both', which='both', left=True, right=True, direction='in')

    ax2.text(0.95,0.95,"M2Pro, gfortran 15.1 (-O3)", va='top',ha='right', transform=ax2.transAxes)



    ax2.set_xlim(0.050,0.31)
    ax2.yaxis.set_major_formatter(FuncFormatter(h.log_formatter_fn))


    ax1.text(0.03,0.93, f"Hoppet v2.0.0, {nloop_names[nloop]} evolution\nymax = 12, dlnlnQ = dy/4", va='top', transform=ax1.transAxes)
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
    ax1.legend(loc='lower right',reverse=True, markerfirst=False)
    ax2.legend(loc='lower left')

    standard_ticks(ax1,ax2)
    
    pdf.savefig(fig,bbox_inches='tight')
    #pdf.savefig(fig,bbox_inches=Bbox.from_extents(0.0,0.0,7.5,4.8))

    # add on the all-falv < 0.90 line
    ax1.plot(run_stats_pre.dy[mask], run_stats_pre.acc_allf_xlt09[mask], label='all-flav, $x<0.9$', **styles[0], ls=":")
    ax1.legend(loc='lower right',reverse=True, markerfirst=False)



    pdf.savefig(fig,bbox_inches='tight')
    plt.close()

    if nloop != 3: return

    #--------------------------------------------------------------------------------
    # Plot over different interpolation orders
    run_stats_oQ2_oY2   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-oQ2-oY2-dy*.dat')
    run_stats_oQ3_oY3   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-oQ3-oY3-dy*.dat')
    run_stats_oQ4_oY4   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-oQ4-oY4-dy*.dat')
    run_stats_oQ4_oY5   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-dy*.dat')
    run_stats_oQ4_oY6   = RunStats(f'{dirM2pro}/nloop{nloop}-preev-oQ4-oY6-dy*.dat')

    run_stats_LHAPDF    = RunStats(f'{dirM2proLHAPDF}/nloop{nloop}-preev-dy*.lhapdf.dat')


    fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=True,height_ratios=[3,2])
    fig.subplots_adjust(hspace=0.04)
    ax2.set_xlabel(r'dy')

    mask = run_stats_oQ2_oY2.dy > 0.035

    ax1.plot(run_stats_LHAPDF .dy[mask], run_stats_LHAPDF .acc_guds_xlt07[mask], label='LHAPDF', **styles[4], ls="-")
    ax1.plot(run_stats_oQ2_oY2.dy[mask], run_stats_oQ2_oY2.acc_guds_xlt07[mask], label='oY=2, oQ=2', **styles[3], ls="-")
    ax1.plot(run_stats_oQ3_oY3.dy[mask], run_stats_oQ3_oY3.acc_guds_xlt07[mask], label='oY=3, oQ=3', **styles[1], ls="-")
    ax1.plot(run_stats_oQ4_oY4.dy[mask], run_stats_oQ4_oY4.acc_guds_xlt07[mask], label='oY=4, oQ=4', **styles[2], ls="-")
    ax1.plot(run_stats_oQ4_oY5.dy[mask], run_stats_oQ4_oY5.acc_guds_xlt07[mask], label='default: oY=5, oQ=4', **styles[0], ls="-")

    #ax1.plot(run_stats_LHAPDF .dy[mask], run_stats_LHAPDF .acc_allf_xlt07[mask], label='LHAPDF', **styles[4], ls="-")
    #ax1.plot(run_stats_oQ2_oY2.dy[mask], run_stats_oQ2_oY2.acc_allf_xlt07[mask], label='oY=2, oQ=2', **styles[3], ls="-")
    #ax1.plot(run_stats_oQ3_oY3.dy[mask], run_stats_oQ3_oY3.acc_allf_xlt07[mask], label='oY=3, oQ=3', **styles[1], ls="-")
    #ax1.plot(run_stats_oQ4_oY4.dy[mask], run_stats_oQ4_oY4.acc_allf_xlt07[mask], label='oY=4, oQ=4', **styles[2], ls="-")
    #ax1.plot(run_stats_oQ4_oY5.dy[mask], run_stats_oQ4_oY5.acc_allf_xlt07[mask], label='default: oY=5, oQ=4', **styles[0], ls="-")

    #ax1.plot(run_stats_oQ4_oY6.dy[mask], run_stats_oQ4_oY6.acc_allf_xlt07[mask], label='oQ=4, oY=6', **styles[4], ls="-")

    #
    print("LHAPDF accuracy stats (allf_xlt07):")
    print(h.reformat(run_stats_LHAPDF .dy[mask],run_stats_LHAPDF .acc_guds_xlt07[mask] ))

    ax2.plot(run_stats_oQ2_oY2.dy[mask], run_stats_oQ2_oY2.t_interp_ns[mask], **styles[3], label='')
    ax2.plot(run_stats_oQ3_oY3.dy[mask], run_stats_oQ3_oY3.t_interp_ns[mask], **styles[1], label='')
    ax2.plot(run_stats_oQ4_oY4.dy[mask], run_stats_oQ4_oY4.t_interp_ns[mask], **styles[2], label='')
    ax2.plot(run_stats_oQ4_oY5.dy[mask], run_stats_oQ4_oY5.t_interp_ns[mask], **styles[0], label='')
    ax2.plot(run_stats_LHAPDF .dy[mask], run_stats_LHAPDF .t_interp_ns[mask], **styles[4], label='LHAPDF')
    #ax2.plot(run_stats_oQ4_oY6.dy[mask], run_stats_oQ4_oY6.t_interp_ns[mask], **styles[4], label='')

    ax1.set_ylabel("rel. accuracy")
    ax2.set_ylabel("interp. time ($t_{xQ}$) [ns]")

    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_ylim(47.0, 120.0)
    #ax2.set_yscale('log')
    ax1.tick_params(axis='both', which='both', left=True, right=True, direction='in')
    ax2.tick_params(axis='both', which='both', left=True, right=True, direction='in')

    ax2.text(0.95,0.95,"M2Pro, gfortran 15.1 (-O3)", va='top',ha='right', transform=ax2.transAxes)

    #ax2.yaxis.set_major_formatter(FuncFormatter(h.log_formatter_fn))
    ax1.text(0.03,0.95, f"Hoppet v2.0.0, {nloop_names[nloop]} evolution\n" +
             "ymax = 12, dlnlnQ = dy/4\n[guds, $x<0.7$]", va='top', transform=ax1.transAxes)
    ax1.legend(loc='lower right',markerfirst=False)
    standard_ticks(ax1,ax2)

    pdf.savefig(fig,bbox_inches='tight')
    plt.close()

def standard_ticks(ax1,ax2):

    ax1.set_ylim(2e-8,1.5e-1)
    #ax1.yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0]))

    xticks_major = ax2.get_xticks().tolist()
    xticks_minor = ax2.get_xticks(minor=True).tolist()
    extra_xticks = [0.05, 0.2]
    all_xtics = sorted(xticks_major + extra_xticks)
    all_xtics_minor = sorted(xticks_minor + [0.15, 0.25])
    ax2.set_xticks(all_xtics, [f"{xt}" for xt in all_xtics])
    ax2.set_xticks(all_xtics_minor, ["" for xt in all_xtics_minor], minor=True)
    ax2.set_xlim(0.050,0.31)

class RunStats(object):
    '''Class to extra run stats for all files matching a certain glob pattern (over dy values)'''
    def __init__(self,glob_str):
        self.files_timing = sorted(glob(f"{glob_str}"))
        self.files_prec   = sorted(glob(f"{glob_str}.acc.smry"))
        self.n = len(self.files_timing)
        if self.n == 0: raise RuntimeError(rf"Found no files in glob {glob_str}")
        self.dy = np.array([
          float(re.findall(r"dy(0\.[0-9]+)",f)[0]) for f in self.files_timing
        ])

        # get the timing numbers
        self.t_init_s = get_number(self.files_timing,"Initialisation time")
        self.t_preev_s = get_number(self.files_timing,"Pre-evolution time")
        self.t_ev_s = get_number(self.files_timing,"Evolution time")
        self.t_interp_ns = get_number(self.files_timing,"Interpolation time")
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
    with PdfPages(pdfname) as pdf: 
        #main(pdf,nloop=2) -> results messed up by weird properties of (small) NLO charm close to threshold that the mask doesn't handle well
        main(pdf,nloop=3)
        main(pdf,nloop=4)


