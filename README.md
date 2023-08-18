# BIO-482 Miniproject  

Repository for the data analysis miniproject of the course _Neuroscience: cellular and circuit mechanisms_ (BIO-482), EPFL. 

This repository contains a MATLAB and Python version of this project, as well as associated functions used for computations.

## Setting up
### Download
- Clone `git clone ...` or download as zip file this repository (green button).
- Download `Data_Bio482.mat` dataset at the provided link.

The repository should contain the following folders:
- `doc`: documentation related to the project and to the data.
- `matlab`: code and helper functions to do the project with MATLAB.
- `python`: code and helper functions to do the project using Python.

**Note**: add the `.mat` file in a `bio482_miniproject/data` folder.

### Installation and usage

1. **If you're using MATLAB**:
  - Add the `/matlab` folder to your MATLAB path to run the code.
2. **If you're using Python**:
  - You need to have Anaconda installed for your system: [install anaconda here](https://docs.anaconda.com/anaconda/install/index.html). 
  - Once Anaconda is installed, open a terminal and go to the `python` directory of the project.
  - Install the provided "bio482" conda environment: `conda env create -f env.yml`.
  - Close terminal to make the conda environment effective.
  
    Then, to work on the project:
  - Go to `python` and open a terminal
  - Activate the environment: `conda activate bio482`.
  - Then run: `jupyter lab`.
  - Open notebooks to start working on the project.


    **Note about paths**:
    - **MATLAB**: You must add folders to MATLAB's path: Right-click on folder -> Add to path...
    - **Python**: In jupyter notebooks, you must change paths based on where you cloned this repository.
  
  _fin_ 
 

