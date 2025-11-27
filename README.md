# Run_out_modelling

**Run_out_modelling** is a Python GUI application running SLBL, Avaframe and Wave propagation.

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
        conda activate Avaframe_env
        ```
2b. **install Avaframe in the same environment (https://docs.avaframe.org/en/1.1_rc1/developinstallwin.html):**
    - install pip, numpy and cython in this environment:
    ```
    conda install pip numpy cython
    ```
    
   - Clone the AvaFrame repository (in a directory of your choice: [YOURDIR]) and change into it:

    ```
    cd [YOURDIR]
    git clone https://github.com/avaframe/AvaFrame.git
    cd AvaFrame
   ```
   - Compile the cython com1DFA part and install Avaframe:
    ```
   python setup.py build_ext --inplace

    pip install -e .
    ```

    
4. **(Optional) Install PyInstaller**  
   If you want to build a standalone executable:
    ```
    pip install pyinstaller
    ```

## Usage

1. **Activate the Conda environment:**
    ```
    conda activate Avaframe_env
    ```

2. **Launch the GUI:**
    ```
    python Main_gui.py
    ```

## Building an Executable

To create a standalone executable:

```
pyinstaller --onefile --windowed Main_gui.py
```

The executable will be in the `dist` directory.

## Troubleshooting

- Make sure your Conda environment is activated before running the GUI.
- For issues with raster or DEM files, check file paths and formats[1].

## License

Distributed under the MIT License. See `LICENSE` for details.
