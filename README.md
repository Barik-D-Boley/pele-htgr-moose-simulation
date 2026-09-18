# Project Pele HTGR Core Block Thermal-Mechanical Simulation

A coupled thermo-mechanical finite element simulation of a transportable High-Temperature Gas-Cooled Reactor (HTGR) microreactor core block based on Project Pele, implemented using the MOOSE (Multiphysics Object-Oriented Simulation Environment) Framework.

---

## Technical Overview

This repository models the thermal distribution and von Mises stress state in a gas-cooled TRISO fuel compact / graphite matrix core block. The model accounts for volumetric heat generation, thermal expansion, temperature-dependent material properties, and convective cooling channels under steady-state and transient conditions.

* **Physics Modules:** `HeatConduction`, `TensorMechanics`
* **Mesh Engine:** Coreform Cubit (Automated `.jou` Journal Scripting)
* **Target Materials:** TRISO compacts, H-451 Nuclear Graphite, Inconel 625 structural containment
* **Post-Processing:** ParaView (`.pvsm` state pipelines)

---

## 4-Phase Development Workflow

### Phase 1: Concept Selection & Material Definition
* Ground core geometry and volumetric power density in published DoD / BWXT Project Pele baseline specifications.
* Formulate temperature-dependent material properties ($k(T)$, $C_p(T)$, $E(T)$, $\alpha(T)$) for nuclear graphite matrix and TRISO fuel compacts.
* Define operating boundary conditions: nominal coolant temperature, convective heat transfer coefficient ($h$), and thermal output.

### Phase 2: CAD Geometry & Coreform Cubit Meshing
* Model a representative 3D hex-block core unit cell in SolidWorks containing fuel compact channels, coolant channels, and structural tie rods.
* Export geometry to `.STEP` format and develop a Coreform Cubit journal script (`mesh_generator.jou`) to automate meshing.
* Generate a high-fidelity HEX8 element mesh with refined node density along channel interfaces.
* Export `.e` (Exodus II) mesh file with defined Block IDs and Boundary Side Sets.

### Phase 3: MOOSE Physics Setup & Solver Execution
* Construct the primary MOOSE input file (`pele_core_sim.i`).
* Implement non-linear heat conduction coupled with isotropic linear elasticity (`TensorMechanics`).
* Apply a 3-2-1 boundary constraint scheme to prevent rigid body motion without introducing artificial thermal stress concentration.
* Execute steady-state solve followed by transient thermal shock simulations (Loss of Coolant / Loss of Heat Sink scenarios).

### Phase 4: Post-Processing & GitHub Documentation
* Process Exodus II results in ParaView to visualize temperature contours, hydrostatic pressure, and von Mises stress distributions.
* Extract line-plots across critical fuel-graphite interface nodes to verify peak thermal gradient locations.
* Document structural safety margins against graphite failure limits under transient temperature spikes.

---

## Project Structure

```text
.
├── cad/                  # SolidWorks .STEP files
├── images/               # Images
├── meshes/               # Coreform Cubit .jou scripts and Exodus .e files
├── moose/                # MOOSE input files (.i) and material parameters
├── post-processing/      # ParaView state files (.pvsm) and high-res renders
└── README.md
```

## How to Run

1. **Generate Mesh:**
   ```bash
   cubit -nographics mesh/generate_mesh.jou
   ```

2. **Run MOOSE Simulation:**
   ```bash
   combined-opt -i moose/pele_core_sim.i
   ```

3. **Visualize Results:**
   Open `post_processing/pele_results.e` in ParaView and load `post_processing/pele_view.pvsm`.