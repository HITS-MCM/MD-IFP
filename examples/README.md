#  This is a collection of example recipes showing how to carry out particular tasks using MD-IFP scripts


### 1. Analysis of tauRAMD interaction fingerprints:
####    For running on a google account
From the [tauRAMD repository](https://github.com/HITS-MCM/tauRAMD) use the first part to get the residence times [1sttutorial_tauRAMD_residencetime.ipynb](https://github.com/HITS-MCM/tauRAMD/blob/master/1sttutorial_tauRAMD_residencetime.ipynb)
Then use the MD-IFP analysis to get the dissociation mechanism: [2nd_tutorial_tauRAMD_dissociationmechanism.ipynb](./2nd_tutorial_tauRAMD_dissociationmechanism.ipynb)
Please upload this jupyter notebook to google colab.
####    For running on an [ebrains](https://wiki.ebrains.eu/bin/view/Main/) account
Please start a [lab session](https://lab.ebrains.eu/) and upload the above jupyter notebook in your space.
From the [tauRAMD repository](https://github.com/HITS-MCM/tauRAMD) use the first part to get the residence times [1st_tutorial_tauRAMD-residencetime-ebrains.ipynb](https://github.com/HITS-MCM/tauRAMD/blob/master/1st_tutorial_tauRAMD-residencetime-ebrains.ipynb)
Then use the MD-IFP analysis to get the dissociation mechanism: [2nd_tutorial_tauRAMD_egress-routes-ebrains.ipynb](./2nd_tutorial_tauRAMD_egress-routes-ebrains.ipynb)
Please upload this jupyter notebook to you [ebrains lab account](https://lab.ebrains.eu).

### 2. Generation of the Interaction Fingerprint table for a set of MD trajectories obtained either from standard MD simulations or RAMD simulations and saving IFPs in a pkl file
    IFP.py
### 3. Visualization  of computed IFPs (from pkl file with IFPs)
    IFP_contacts_quickView.py

