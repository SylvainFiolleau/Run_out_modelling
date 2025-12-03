# Run_out_modelling

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/SylvainFiolleau/Run_out_modelling)](https://github.com/SylvainFiolleau/Run_out_modelling/issues)

**Run_out_modelling** is a Python GUI application for running **SLBL**, **AvaFrame**, and **wave propagation** simulations.

---

## Table of Contents

- [Prerequisites](#prerequisites)  
- [Installation](#installation)  
- [Usage](#usage)  
- [Building an Executable](#building-an-executable)  
- [Simpler Alternative](#simpler-alternative-in-bash-or-windows-batch)  
- [Troubleshooting](#troubleshooting)  
- [License](#license)  

---

## Prerequisites

- **Anaconda** or **Miniconda** (Python 3.8+ recommended)  
- **Git**

---

## Installation

1. **Clone the repository:**

    ```bash
    git clone https://github.com/SylvainFiolleau/Run_out_modelling.git
    cd Run_out_modelling
    ```

2. **Create the Conda environment:**

    ```bash
    conda env create -f environment.yml
    conda activate Avaframe_env
    ```

3. **Install AvaFrame**  

   Follow the official guide without installing numpy again: [AvaFrame Installation](https://docs.avaframe.org/en/1.1_rc1/developinstallwin.html)  

   a. Install `pip` and `cython` in the environment:

    ```bash
    conda install pip cython
    ```

   b. Clone the AvaFrame repository (choose your directory `[YOURDIR]`):

    ```bash
    cd [YOURDIR]
    git clone https://github.com/avaframe/AvaFrame.git
    cd AvaFrame
    ```

   c. Compile Cython parts and install:

    ```bash
    python setup.py build_ext --inplace
    pip install -e .
    ```

4. **(Optional) Install PyInstaller**  

    ```bash
    pip install pyinstaller
    ```

---

## Usage

1. **Activate the environment:**

    ```bash
    conda activate Avaframe_env
    ```

2. **Launch the GUI:**

    ```bash
    python Main_gui.py
    ```

---

## Building an Executable

Create a standalone executable:

```bash
pyinstaller --onefile --windowed Main_gui.py
