#!/usr/bin/env python3
"""Matplotlib template generated automatically with 
/Users/gsalam/scripts/mptemplate.py n3lo-scale-tests.py
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages
from   matplotlib.ticker import ScalarFormatter
from   matplotlib.ticker import FuncFormatter
from   matplotlib import cycler
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import copy
import hfile as h

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
plt.rc('figure', figsize=(9,3))

def main(pdf):
    
    nloop = 4
    file00125 = f"nloop{nloop}-asQ0-0.0125-range100000000.dat"
    file00250 = f"nloop{nloop}-asQ0-0.025-range10000.dat"
    file00500 = f"nloop{nloop}-asQ0-0.05-range100.dat"
    file01000 = f"nloop{nloop}-asQ0-0.1-range10.dat"
    files = [file01000, file00500, file00250, file00125, ]
    alphas_values = [0.1, 0.05, 0.025, 0.0125, ]
    re_muR10 = f"xmuR *= *1\.0"
    re_muR05 = f"xmuR *= *0\.5"
    re_muR20 = f"xmuR *= *2\.0"

    fig,(ax1,ax2) = plt.subplots(ncols=2)
    for ax in ax1, ax2:
        ax.set_xscale('log')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(rf'$[f(x_{{\mu_R}}=2)/f(x_{{\mu_R}}=1) - 1] / \alpha_s^{nloop}$')

    ax1.set_title("gluon")
    ax2.set_title("d quark")
    #ax1.text(0.1, 0.1, r'Q/Q_0 = 10^(0.1/\alpha_s(Q_0))$', transform=ax.transAxes)


    for ffile, alphas in zip(files, alphas_values):

        columns = {'x':0, 'g':7, 'd':8}
        res_muR10 = h.get_array_plus_comments(ffile,regexp=re_muR10,columns=columns)
        res_muR05 = h.get_array_plus_comments(ffile,regexp=re_muR05,columns=columns)
        res_muR20 = h.get_array_plus_comments(ffile,regexp=re_muR20,columns=columns)

        #res = h.get_array_plus_comments("filename",regexp='',columns={'x':1,'y':(3,4)}) 
        # ax.set_xlim(0,2.5)
        # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        # ax.yaxis.set_minor_locator(MultipleLocator(0.02))
        # ax.tick_params(top=True,right=True,direction='in',which='both')
        # ax.xaxis.set_major_formatter(FuncFormatter(h.log_formatter_fn))
        # ax.grid(True,ls=":")
        # ax.set_title("title")
        ax1.plot(res_muR10.x, (res_muR20.g / res_muR10.g - 1)/alphas**(nloop), label=fr'$\alpha_s(Q_0) = {alphas}$')
        ax2.plot(res_muR10.x, (res_muR20.d / res_muR10.d - 1)/alphas**(nloop), label=fr'$\alpha_s(Q_0) = {alphas}$')

    ax2.text(0.1, 0.1, r'$Q/Q_0 = 10^{0.1/\alpha_s(Q_0)}$', transform= ax.transAxes)
    ax1.legend()
    pdf.savefig(fig,bbox_inches='tight')
    #pdf.savefig(fig,bbox_inches=Bbox.from_extents(0.0,0.0,7.5,4.8))
    plt.close()

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
    with PdfPages('n3lo-scale-tests.pdf') as pdf: main(pdf)


