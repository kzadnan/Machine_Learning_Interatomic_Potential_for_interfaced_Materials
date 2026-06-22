# Machine Learning Interatomic Potentials-semiconductor-interfaces
Welcome to the official repository for the machine learning interatomic potentials developed and utilized in our computational materials research. 

This repository hosts highly accurate potentials designed explicitly to simulate complex heterogeneous interfaces. Thermal management in advanced semiconductor and electronic devices often relies on integrating high-thermal-conductivity substrates, such as diamond. However, interfacial thermal resistance frequently serves as a major bottleneck. 

These potentials have been specifically trained to accurately capture complex atomic interactions, phonon dispersion, and thermal boundary conductance (TBC) across these critical interfaces using Non-Equilibrium Molecular Dynamics (NEMD) simulations in LAMMPS.

---

## 📂 Included Potentials & Interface Descriptions

### 1. Silicon/Diamond (Si/Diamond) Interface
Direct integration of silicon with diamond is highly desirable for advanced thermal management, but lattice mismatches and weak bonding often degrade thermal transport. The deep learning potential and moment tensor potentials are developed to accurately model the Si/diamond interface, including the intricate effects of intermediate bonding layers. 
* **Applications:** Specifically optimized to evaluate and predict the enhancement of Thermal Boundary Conductance (TBC) through the introduction of ultrathin interlayers, such as $\mathrm{SiN}_{x}$ and amorphous carbon. 
* **Compatibility:** Designed for use with the DeePMD-kit and MTP backend in LAMMPS.

### 2. Gold/Diamond (Au/Diamond) Interface
Metal-dielectric interfaces present unique challenges in computational modeling due to the stark differences in their vibrational densities of states. This Moment Tensor Potential (MTP) bridges that gap, providing near-DFT accuracy at molecular dynamics scales.
* **Applications:** Developed to predict the thermal boundary conductance of metal-diamond interfaces. It accurately captures the phonon dynamics and atomic-level structural relaxations required to evaluate heat dissipation from metallic contacts into diamond substrates.
* **Compatibility:** Designed for use with the MTP package in LAMMPS.

---

## 📝 Citations & Usage

If you utilize these potentials, training datasets, or associated simulation methodologies in your research, please cite the corresponding publications below:

### Si/Diamond Potential (DeePMD & MTP)
**Citation:**
Adnan, K. Z., Abedien, T., & Feng, T. (2026). Si/diamond thermal boundary conductance enhanced by $\mathrm{SiN}_{x}$ and amorphous carbon interlayers. *Physical Review Applied*, 25(4), 044013. https://doi.org/10.1103/84jf-svtg

**BibTeX:**
```bibtex
@article{adnan2026sidiamond,
  title = {$\mathrm{Si}$/diamond thermal boundary conductance enhanced by ${\mathrm{Si}\mathrm{N}}_{x}$ and amorphous carbon interlayers},
  author = {Adnan, Khalid Zobaid and Abedien, Tanvirul and Feng, Tianli},
  journal = {Phys. Rev. Appl.},
  volume = {25},
  issue = {4},
  pages = {044013},
  numpages = {14},
  year = {2026},
  month = {Apr},
  publisher = {American Physical Society},
  doi = {10.1103/84jf-svtg},
  url = {[https://link.aps.org/doi/10.1103/84jf-svtg](https://link.aps.org/doi/10.1103/84jf-svtg)}
}
```


### Metal(Al,Mo,Zr,Au)/Diamond Potential (MTP)
Please cite the following paper when using the Au/diamond Moment Tensor Potential (MTP):

**Citation:**
Adnan, K. Z., Neupane, M. R., & Feng, T. (2024). Thermal boundary conductance of metal–diamond interfaces predicted by machine learning interatomic potentials. International Journal of Heat and Mass Transfer. https://doi.org/10.1016/j.ijheatmasstransfer.2024.126227

**BibTeX:**
```bibtex
@article{adnan2024thermal,
  title = {Thermal boundary conductance of metal–diamond interfaces predicted by machine learning interatomic potentials},
  author = {Adnan, Khalid Zobaid and Neupane, Mahesh R. and Feng, Tianli},
  journal = {International Journal of Heat and Mass Transfer},
  year = {2024},
  publisher = {Elsevier},
  doi = {10.1016/j.ijheatmasstransfer.2024.126227},
  url = {[https://doi.org/10.1016/j.ijheatmasstransfer.2024.126227](https://doi.org/10.1016/j.ijheatmasstransfer.2024.126227)}
}