# Third-party tools

This project orchestrates **open-source** binaries and libraries. Install and comply with each upstream license when you enable the corresponding backend:

| Component | Typical license | Notes |
|-----------|-----------------|--------|
| **GNINA** | GPL-2.0 / project license | GPU docking; not bundled here |
| **AutoDock Vina** | Apache-2.0 | CPU docking |
| **fpocket** | MIT | Pocket prediction |
| **OpenMM** | MIT / LGPL parts | See OpenMM docs |
| **PDBFixer** | MIT | Receptor repair |
| **RDKit** | BSD-like | Cheminformatics |
| **Open Babel** | GPL-2.0 | Formats / protonation |
| **PubChem** | Public data; see NCBI terms of use | PUG-REST API |

The **mock docking backend** ships only in Python and does not embed third-party docking code.
