# FASERCAL - electronic calorimeter for FASER Run4

FASERcal code to simulate and analyse events in the FASERCAL detector

# BatchReco

Code to read FASERCAL GEANT4 output and batch reconstruct events, filling histograms, ...

# EvDisplay

Interactive Event Display of FASERCAL GEANT4 output

# FASERTuple

Convert official FASER MC files into FASERCAL PO files (generator level) 

# TauSearch
A generator level tau search analysis code

- t.C : code to convert FASER ntuple into event summary tuples
- s.C : analyse event summary tuples for each tau decay channel and create sig/background tuples
- a.C : read sig/bkg tuples for each decay channel and perform BDT analysis

# Installation

- Get the source code:

$ git clone https://github.com/rubbiaa/FASER.git

- Set up ROOT and GEANT4 environment in the setup.sh file:

setup.sh:
    source <ROOTINSTAL>/bin/thisroot.sh
    source <GEANT4INSTALL>/bin/geant4.sh

$ source setup.sh

- On lxplus use the following command instead:

$ source lxplus_setup.csh




