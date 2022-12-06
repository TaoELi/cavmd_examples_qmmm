# QM/MM CavMD simulations of a dilute Fe(CO)<sub>5</sub> in n-dodecane

This folder contains all necessary files for reproducing the following paper:

- Li, T. E., Hammes-Schiffer, S.  QM/MM Modeling of Vibrational Polariton Induced Energy Transfer and Chemical Dynamics. [J. Am. Chem. Soc. 2022, accepted](https://arxiv.org/abs/2212.02322).

## A quick reproducing of all figures in the paper

Please go to folder **plotting/**, and then run, e.g., *python plot_fig2_polariton_spectrum.py* in your terminal to reproduce Fig. 2 in the paper.

Please use a recent version of python 3 (e.g., python 3.8). Necessary python libraries for plotting: numpy, scipy, matplotlib, pandas, [MDAnalysis](https://www.mdanalysis.org/pages/installation_quick_start/).

## Want to rerun the simulations in your Linux environment?

- Because QM/MM CavMD simulations need a parallel scheme to model a lot of molecules, **it can only be efficiently run on a Linux cluster**. The files in folders **qmmm/** and **qm_solvent_free/** include everything you need to start a job.

- A tutorial for QM/MM simulations is available at https://github.com/TaoELi/cavity-md-ipi/tree/master/tutorials. In brief, this QM/MM CavMD simulation is implemented by calling many instances of [Q-Chem](https://www.q-chem.com/) at each time step in parallel. In practice, one can also use free electronic structure packages such as [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) instead of Q-Chem (but in this case you need to write your own glue code from i-pi to ORCA). Please contact me if you want to use the CavMD code to start a large-scale QM/MM calculation.  

- Before trying this large-scale QM/MM CaMD simulations, why not first try to run a first-principles CavMD simulations for a single molecule inside the cavity **in your laptop**? The CavMD code can do it, and this feature is fully open source. Please check https://github.com/TaoELi/cavity-md-ipi/tree/master/tutorials for an example! 


##### If you find any error on reproducing the results, please feel free to send me an email.

My contact info is available at my personal website: https://taoeli.github.io/.
