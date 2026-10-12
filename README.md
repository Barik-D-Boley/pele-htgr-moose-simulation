# Coupled OpenMC–MOOSE Multiphysics Simulation of HTGR Microreactor RVE
## Status: Active Development — Phase 3 (Data Coupling & MOOSE Multiphysics Solve)

A 2-step coupled neutronics and thermo-mechanical finite element simulation of a transportable High-Temperature Gas-Cooled Reactor (HTGR) microreactor core unit cell based on Project Pele. The workflow couples Monte Carlo neutron transport (**OpenMC**) with continuum thermal-stress analysis (**MOOSE Framework**).

## Technical Overview

This repository models the coupled reactor physics, spatial heat generation, temperature distribution, and von Mises stress state in a TRISO fuel compact / H-451 graphite matrix Representative Volume Element (RVE). 

Rather than assuming uniform heat generation, this project calculates the spatially dependent volumetric fission power density $Q(x,y,z)$ from first principles using OpenMC neutron transport, maps the resulting tally grid into MOOSE, and evaluates local stress peaking factors under steady-state and transient conditions.

* **Physics Modules:** Monte Carlo Neutron Transport, `HeatConduction`, `TensorMechanics`
* **Software Stack:** OpenMC, MOOSE Framework, Coreform Cubit, Python (`h5py`, `scipy`), ParaView
* **Target Materials:** HALEU TRISO compacts (19.75% U-235), H-451 Nuclear Graphite, Helium Coolant
* **Coupling Paradigm:** Off-line 2-step data transfer (Neutronics Mesh Tally $\rightarrow$ Spatial Heat Source Kernel)

## 4-Phase Development Workflow

### Phase 1: OpenMC Neutronics & Power Tallying [*COMPLETE*]
* Define isotopic compositions for HALEU TRISO compacts, H-451 graphite matrix, and helium channels.
* Build a Constructive Solid Geometry (CSG) model of the $18.80\text{ mm} \times 9.40\text{ mm}$ RVE in OpenMC.
* Execute Monte Carlo eigenvalue solve ($k_{\text{eff}}$) with a 3D RegularMesh tally to extract volumetric fission heating density $Q(x,y,z)$ in $\text{W/m}^3$.

### Phase 2: CAD Geometry & Coreform Cubit Meshing [*COMPLETE*]
* Model the 3D unit cell in SolidWorks and export to `.STEP` format.
* Execute automated Coreform Cubit journal script (`create_mesh.jou`) to generate a 10,152-element HEX8 mesh.
* Define explicit surface sidesets (`symm_x_min`, `symm_y_min`, `bottom_face`, `fuel_wall`, `coolant_wall`) and scaling to SI units.

### Phase 3: Data Coupling & MOOSE Multiphysics Solve [*IN PROGRESS - VERIFICATION*]
* Run Python bridge script (`coupling/map_power.py`) to parse OpenMC `statepoint.*.h5` tally outputs and format a 3D spatial grid function for MOOSE.
* Construct primary MOOSE input file (`moose/pele_core_sim.i`) using `PiecewiseMultilinear` functions driven by OpenMC power tallies.
* Enforce planar symmetry roller constraints and solve coupled non-linear heat conduction and isotropic linear elasticity.

### Phase 4: Comparative Stress Analysis & Documentation [*PLANNED*]
* Compare thermal gradients and von Mises stress fields between realistic neutronics power profiles vs. uniform volumetric heat assumptions.
* Process Exodus II results in ParaView to identify peak stress concentrations along the fuel-graphite interface.
* Document structural safety margins and local peaking factors for microreactor core licensing considerations.

## Verification & Quality Assurance (QA)

To ensure numerical stability, physical conservation, and accurate multivariable coupling between OpenMC and MOOSE, this workflow executes specific verification checks across four modeling domains:

### 1. Spatial Geometry & Mesh Quality Checks (SolidWorks & Coreform Cubit)
* **Domain Bounding:** Verified unit cell dimensions ($18.80\text{ mm} \times 9.40\text{ mm}$) against Project Pele RVE design specs to prevent volume overlaps.
* **Element Metrics:** Inspected the 10,152-element HEX8 mesh in Coreform Cubit to enforce minimum scaled Jacobian values $> 0.6$ and aspect ratios $< 3.0$ across high-gradient fuel-matrix interfaces.
* **Unit Consistency:** Automated scale transformations in `create_mesh.jou` to convert CAD dimensions ($\text{mm}$) to SI units ($\text{m}$) required by MOOSE.
* **Sideset Integrity:** Validated explicit surface sideset IDs (`symm_x_min`, `symm_y_min`, `bottom_face`, `fuel_wall`, `coolant_wall`) against Exodus block exports to prevent unconstrained degrees of freedom.

### 2. Neutronics & Cross-Section Verification (OpenMC)
* **Isotopic Composition:** Verified HALEU enrichment ($19.75\text{\% } ^{235}\text{U}$) and H-451 graphite composition against nominal HTGR material specifications.
* **Eigenvalue Convergence:** Monitored $k_{\text{eff}}$ convergence across Monte Carlo batches to ensure active source distribution equilibrium before tally generation.
* **Spatial Alignment:** Confirmed 3D RegularMesh tally bounds ($x, y, z$) match the spatial domain and coordinate origin of the Cubit FEA mesh.

### 3. Data Transfer & Energy Conservation (Python Bridge)
* **Spatial Mapping:** Verified `map_power.py` grid interpolation routines using `h5py` and `SciPy` to eliminate energy leakage or spatial distortion during HDF5-to-CSV data transfer.
* **Integrated Heat Check:** Verified total integrated volumetric fission heat ($\iiint Q(x,y,z) \, dV$) calculated from OpenMC matches the integrated heat flux input assigned across MOOSE material blocks.

### 4. Continuum Mechanics & Boundary Constraints (MOOSE)
* **Rigid-Body Mode Elimination:** Applied planar symmetry roller constraints (`symm_x_min`, `symm_y_min`, `bottom_face`) to eliminate 6 rigid-body modes without inducing artificial thermal stresses.
* **Material Property Mapping:** Cross-referenced temperature-dependent thermal conductivity, coefficient of thermal expansion (CTE), and Young's modulus across graphite and fuel blocks against nuclear materials handbook baselines.

## Project Structure

```text
.
├── cad/                  # SolidWorks .STEP geometry files
├── coupling/             # Python scripts to map OpenMC HDF5 tallies to MOOSE functions
├── images/               # High-res renders and comparison plots
├── meshes/               # Coreform Cubit .jou scripts and exported Exodus (.e) meshes
├── moose/                # MOOSE input files (.i) and material definitions
├── neutronics/           # OpenMC Python model, cross-section links, and tally definitions
├── post-processing/      # ParaView state files (.pvsm) and CSV line plots
└── README.md
```

## How to Run

1. **Run OpenMC Neutronics Solve:**
   ```bash
   python neutronics/openmc_rve_model.py
   ```
   *Outputs `statepoint.100.h5` containing $k_{\text{eff}}$ and 3D fission power tallies.*

2. **Map OpenMC Power Density to MOOSE Grid:**
   ```bash
   python coupling/map_power.py --input neutronics/statepoint.100.h5 --output moose/power_density.csv
   ```

3. **Generate Exodus FEA Mesh (Optional if `.e` exists):**
   ```bash
   cubit -nographics meshes/create_mesh.jou
   ```

4. **Execute MOOSE Multiphysics Solver:**
   ```bash
   mpiexec -n 4 combined-opt -i moose/pele_core_sim.i
   ```

5. **Visualize Results:**
   Open `moose/pele_core_sim_out.e` in ParaView and load `post-processing/pele_multiphysics_view.pvsm`.