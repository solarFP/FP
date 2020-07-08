-------------
SUMMARY
-------------
During solar flares numerous particles (electrons, protons, and ions) are accelerated to high energies at the tops of magnetic flux loops in the solar corona. They travel down the loops interacting with and heating the ambient plasma over the course of their transport. This heating drives upflows of dense material into the magnetic loops, resulting in the intense bursts of radiation that characterize flares. The interaction of these high energy particles with the ambient solar atmosphere is well-represented by the Fokker-Planck kinetic theory, which includes Coulomb force interactions coupled with magnetic and external electric field forces. Mathematically, the Fokker-Planck theory is represented by a nonlinear partial differential equation that must be solved numerically. The particular solution of the Fokker-Planck equation depends on the atmospheric conditions within the magnetic flux loop (e.g., temperature, density, ionization fractions, magnetic flux density) as well as the nonthermal distribution function of the flare-accelerated particles. Thus, particular solutions must be obtained for each flare that is to be analyzed. We have developed the FP package to efficiently perform this computational task for arbitrary loop conditions and nonthermal particle distributions.

FP is a Fortran package developed with MPI for parallel processing and includes a Python and IDL wrapper to facilitate ease of use. It is designed for use by solar and stellar physics researchers. It has been used as input to the radiation hydrodynamic model, RADYN, to model the solar atmosphere's response to the heating and nonthermal ionization produced by impacting flare-accelerated particles. It has also been used as a forward-modeling tool to invert X-rays observed during flares by NASA's Ramaty High Energy Solar Spectroscopic Imager (RHESSI) to determine the nonthermal particle distribution function leaving the acceleration site.

-------------
Documentation
-------------
Documentation describing how to compile and run FP as well as examples is located in the file, "FP.pdf", in the "doc" subdirectory of this package. 

