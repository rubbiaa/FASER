To run HL-LHC neutron studies in FASERCalProtoG4, use the existing single-particle mode with the new spectrum option.

What changed:
- You can now choose the single-particle species with /generator/singleParticleName.
- You can turn on spectrum sampling with /generator/useEnergySpectrum true.
- You can point to a two-column file with /generator/energySpectrumFile.

Expected spectrum file format:
- Column 1: neutron kinetic energy in GeV.
- Column 2: dN/dE at that energy.
- Lines starting with # are ignored.
- Energies must be positive and sorted order is not required.

Important physics note:
- The black line in your plot is a differential flux, not a normalized probability distribution.
- The code samples energies with probability proportional to dN/dE times the inferred bin width.
- This preserves the spectrum shape for transport studies.
- If you need an absolute rate for one HL-LHC year at 250 fb^-1, normalize after the run using the integral of your digitized spectrum over the illuminated area and solid angle assumed by the source paper.

How to prepare the black curve:
- Digitize the black curve into a text file. WebPlotDigitizer works well for this.
- Save the points as kineticEnergy_GeV and dNdE.
- Replace the placeholder values in input/neutron_spectrum_template.dat.

How to run:
- Build from the build directory with cmake --build . -j4.
- Run batch mode with ./faserps ../RunFASER_neutron_spectrum.mac from the build directory.

What the current macro does:
- Uses prototypeMode true and the tilted geometry settings already used in this directory.
- Fires neutrons from the existing single-particle start point at x=240 mm, y=240 mm, z=-800 mm.
- Sends them along +z with kinetic energy sampled from your spectrum file.

Current limitation:
- The source position and angular distribution are still mono-directional, following the current single-particle implementation.
- If the HL-LHC source model also requires angular spread or a spatial source plane, extend PrimaryGeneratorAction in the same way by adding sampled x, y, theta, and phi distributions.