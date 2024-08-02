# Instructions for Reading ROOT Files using PyROOT

This folder provides instructions and necessary files to read ROOT files using PyROOT. The software versions used are:

- ROOT version: `6.30/06`
- Python version: `3.10.9`
- NumPy version: `2.0.0`

## Directory Structure

The repository contains the following directory:

- `lib`: Contains the necessary header files (`Classes.h`, `LinkDef.h`) to create a shared object (.so) library for PyROOT.

## Classes Information

The classes defined in `Classes.h` are `DigitizedTrack` and `TcalEvent:GEOM_DETECTOR`. These classes are based on the implementations found in `FASER/CoreUtils/TcalEvent.hh`. If there are any changes in the source, these classes should be updated accordingly.

## Steps to Create the .so Library

Follow the steps below to create the shared object (.so) library:

1. Navigate to the `lib` directory:

    ```bash
    cd lib
    ```

2. Generate the dictionary source files using `rootcint`:

    ```bash
    rootcint -f ClassesDict.cxx -c Classes.h LinkDef.h
    # This will create the files: ClassesDict.cxx and ClassesDict_rdict.pcm
    ```

3. Compile the dictionary source file into a shared object (.so) library:

    ```bash
    g++ -o ClassesDict.so -shared -fPIC `root-config --cflags` ClassesDict.cxx `root-config --libs`
    # This will create the file: ClassesDict.so
    ```

4. Navigate back to the root directory:

    ```bash
    cd ..
    ```

## Running the Python Script

To run the `generate_dataset.py` script, use the following command:

```bash
python generate_dataset.py --lib lib/ClassesDict.so --input ROOT_DIRECTORY --output OUTPUT_DIRECTORY
```

Where:

- `ROOT_DIRECTORY`: The directory where the ROOT files are stored.
- `OUTPUT_DIRECTORY`: The directory where you want to save the NumPy converted files.

Example of usage:

```bash
python generate_dataset.py --lib lib/ClassesDict.so --input /scratch/faser/FASERCALDATA_v2.0 --output /scratch/faser/np_events
```

