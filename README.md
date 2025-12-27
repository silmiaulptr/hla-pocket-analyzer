# hla-pocket-analyzer
Automated Hydrogen Bond detection and HLA Pocket mapping for Vaccine Design using Python and PyMOL.
# HLA Pocket Analyzer

## Overview
This Python script was developed to automate the analysis of molecular interactions between receptors (HLA/MHC) and peptides in *in silico* vaccine design studies. It is specifically designed to overcome the efficiency limitations of manual analysis in PyMOL by automatically mapping hydrogen bonds to specific binding pockets.

## Key Features
* **Automated H-Bond Detection:** Automatically identifies hydrogen bonds based on standard physical criteria (distance of $2.5 - 3.5 \text{ \AA}$) and specific donor/acceptor atom types (N/O).
* **Functional Pocket Mapping:** Maps residue positions to the six major **HLA Pockets (A, B, C, D, E, and F)**.
* **High Throughput Data Export:** Compiles interaction details (atom paths, distances, and pocket identities) directly into a **CSV** format for further statistical analysis.
* **Visual Automation:** Automatically configures PyMOL displays, including pocket color schemes and stick/cartoon representations for qualitative review.

## Technical & Personal Context
This script is an integral part of my undergraduate thesis at **UIN Sunan Gunung Djati Bandung**. As the first student to initiate *in silico* vaccine design research in my department, I developed this tool to prove that high quality molecular research can be achieved through computational automation, even when physical laboratory facilities are limited. 

This project reflects my journey as a self taught bioinformatician, demonstrating my ability to bridge biological concepts with practical programming solutions.

## Prerequisites
* **Software:** PyMOL 3.1.5.1.
* **Language:** Python 3.13.
* **Libraries:** `pymol`, `csv`, `os`.

## How to Use
1. Load your receptor and peptide structures into PyMOL.
2. Ensure your selections are named `receptor` and `peptide`.
3. Run the script using the command: `hla-pocket-analyzer.py`.
4. Find your results in the `Documents/PyMOL_Results` folder.
