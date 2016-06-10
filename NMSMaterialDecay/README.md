### NMSMaterialDecay

The NMSMaterialDecay library provides a particle source that can be placed in any Geant4 application. It can creates starting α particles or electrons according to either alpha- and beta-decays, as well as neutrons or gammas from spontaneous fission of given isotopes. The isotopes are taken from a source material.

The library forms a part of the application NMS (Neutron Multiplicity Simulation)

#### Requirements

  - CMake
  - Geant4, version >= 9.6
  - Fission, from the Physics Simulation Package (http://nuclear.llnl.gov/simulation/main.html)

#### Build
     1. Extract / clone repository in a directory ("A"). 
     2. Create a new directory ("B") and go to B. 
     3. Edit CMakeLists.txt to reflect your location of fission
     4. Execute `cmake`
     5. Execute `make` / `make install`