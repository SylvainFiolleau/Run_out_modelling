# Run_out_modelling

**Run_out_modelling** is a Python GUI application ruuning SLBL, Avaframe and Wave propagation.

## Prerequisites

- Anaconda or Miniconda (Python 3.8+ recommended)
- Git

## Installation

1. **Clone the repository:**
    ```
    git clone https://github.com/SylvainFiolleau/Run_out_modelling.git
    cd Run_out_modelling
    ```

2. **Create the Conda environment:**
        ```
        conda env create -f environment.yml
        conda activate runout
        ```
3. **(Optional) Install PyInstaller**  
   If you want to build a standalone executable:
    ```
    pip install pyinstaller
    ```

## Usage

1. **Activate the Conda environment:**
    ```
    conda activate runout
    ```

2. **Launch the GUI:**
    ```
    python main.py
    ```
    *(Replace `main.py` with the actual entry point if different.)*

## Building an Executable

To create a standalone executable:
