Prerequisites

Anaconda or Miniconda (Python 3.8+ recommended)

Git

Installation

Clone the repository:

git clone https://github.com/SylvainFiolleau/Run_out_modelling.git
cd Run_out_modelling


Create the Conda environment:

conda env create -f environment.yml
conda activate Avaframe_env


Install AvaFrame in the same environment
Follow the official instructions: AvaFrame Installation Guide

a. Install pip and cython:

conda install pip cython


b. Clone the AvaFrame repository (choose a directory [YOURDIR]):

cd [YOURDIR]
git clone https://github.com/avaframe/AvaFrame.git
cd AvaFrame


c. Compile the Cython com1DFA part and install AvaFrame in the Conda environment:

python setup.py build_ext --inplace
pip install -e .


(Optional) Install PyInstaller
If you want to build a standalone executable:

pip install pyinstaller

Usage

Activate the Conda environment:

conda activate Avaframe_env


Launch the GUI:

python Main_gui.py

Building an Executable

To create a standalone executable:

pyinstaller --onefile --windowed Main_gui.py


The executable will be created in the dist directory.

Simpler Alternative in Bash / Windows Batch
@echo off
REM === Activate conda and run script ===
CALL "C:\ProgramData\anaconda3\Scripts\activate.bat" "C:\ProgramData\Conda_envs\Avaframe_env"
python "C:\ProgramData\Run_out_modelling\Main_gui.py"
pause

Troubleshooting

Make sure your Conda environment is activated before running the GUI.

For issues with raster or DEM files, check file paths and formats.

License

Distributed under the MIT License. See LICENSE for details.
