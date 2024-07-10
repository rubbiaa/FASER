# FASERCAL - electronic calorimeter for FASER Run4

FASERCAL code to simulate and analyse events in the FASERCAL detector

# Event data flow

This project follows a structured workflow for simulating and processing neutrino interactions within the FASERCAL experiment using GEANT4. The process begins with the FASERMC official Monte Carlo (MC) simulation. These interactions are converted into the TPOEvent class format by **ConvertFASERMC**. Next, the **FASERG4** module uses GEANT4 to simulate events from the TPOEvent class, generating TcalEvent objects. These simulated events can be visualized using the **EvDisplay** module or processed through the **BatchReco** module for batch reconstruction, providing a comprehensive analysis pipeline for neutrino interaction events.

![Diagram of the project](images/eventchainflow.png)

# Batch (BatchReco)

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

# Installation preliminaries

- Get the source code:

$ git clone https://github.com/rubbiaa/FASER.git

- Set up ROOT and GEANT4 environment in the setup.sh file:

setup.sh:
    source <ROOTINSTAL>/bin/thisroot.sh
    source <GEANT4INSTALL>/bin/geant4.sh

$ source setup.sh

- On lxplus use the following command instead:

$ source lxplus_setup.csh

# Install event display

 - move to the evDisplay directory and compile with "make"

$ cd evDiplay

$ make

 - if compilation and linking was successful, the executable is "evDisplay.eve"

 - make sure G4 FASERCAL simulated files are linked at the "input" subdirectory
 
$ ln -fs </path_to_g4_simulated_data> input

 - on lxplus.cern.ch, there is data available

 $ ln -fs /eos/home-r/rubbiaa/FASER/G4DATA input

 - run the event display

 $ evDisplay.exe

 