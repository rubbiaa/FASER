# Instructions for Reading ROOT Files using PyROOT

This folder provides instructions and necessary files to read ROOT files using PyROOT. The software versions used are:

- ROOT version: `6.30/06`
- Python version: `3.10.9`
- NumPy version: `2.0.0`

Make sure that the root version you use to compile this library with, is the same as you will use in PyROOT.

## Directory Structure

The repository contains the following directory:

- `lib`: Contains the necessary header file (`LinkDef.h`) to create a shared object (.so) library for PyROOT.

## Classes Information

The actual class implementations in `FASER/CoreUtils/TcalEvent.hh` have to be used, as we need to compile against the source, to have a properly compiled shared object library.

## Automatic Script to create the .so file

Execute the `generateLib.sh` script, using the absolute path to the `FASER` directory as argument.
This compiles the `ClassesDict.so` file and stores it in `FASER/Python_io/ClassesDict.so`.

## Manual Steps to Create the .so Library

Follow the steps below to create the shared object (.so) library:

1. Navigate to the `CoreUtils` directory:

    ```bash
    cd CoreUtils
    ```
2. Copy the `lib/LinkDef.h` file to the `CoreUtils` directory:

   ```bash
   cp ../Python_io/lib/LinkDef.h .
   ```

3. Generate the dictionary source files using `rootcint`:

    ```bash
    rootcint -f ClassesDict.cxx -c *.hh LinkDef.h
    # This will create the files: ClassesDict.cxx and ClassesDict_rdict.pcm
    ```

4. Compile the dictionary source file into a shared object (.so) library:

    ```bash
    g++ -o ClassesDict.so -shared -fPIC `root-config --cflags` *.cc ClassesDict.cxx `root-config --libs`
    # This will create the file: ClassesDict.so
    ```

5. Navigate back to the root directory:

    ```bash
    cd ..
    ```

## Usage in python

To use the python IO-file, you need to import it using

	```py
	ROOT.gSystem.Load("/Absolute/Path/FASER/Python_io/lib/ClassesDict.so")
	```

Then you can use the individual functions and classes stored in the CoreUtils source files.

	```py
	import ROOT
	ROOT.gSystem.Load("/Absolute/Path/FASER/Python_io/lib/ClassesDict.so")
	example = ROOT.TPORec()
	example.FUNCTIONOFCLASS()
	result = ROOT.STANDALONEFUNCTION(a, b, c)
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
