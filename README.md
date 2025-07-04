# fNIRS Image Reconstruction and Motion Correction

This academic project explores motion correction strategies and image reconstruction techniques for functional Near-Infrared Spectroscopy (fNIRS) data collected during a motor task involving texting with the right or left hand. The analysis focuses on identifying reliable methods for artifact removal and reconstructing hemodynamic responses from DOT (Diffuse Optical Tomography) measurements using MATLAB.

## Project Description
This project was completed as part of the Imaging for Neuroscience course. The goal was to process and analyze DOT data from a single subject and compare cortical activation patterns during two motor conditions: texting with the right hand and texting with the left hand.

## Aims of the Study
- Identify and remove "bad" fNIRS channels based on intensity and SNR criteria.
- Apply and evaluate motion correction techniques on optical density data.
- Preprocess the fNIRS signal and compute block-averaged responses.
- Reconstruct spatial maps of HbO and HbR concentration changes on the GM surface.
- Compare activation patterns between right-hand and left-hand texting tasks.

## Technologies & Tools
- Language: MATLAB
- Toolboxes:
  - Homer2 – for fNIRS data preprocessing and motion correction
  - iso2mesh – for 3D mesh processing and visualization
  - Custom MATLAB scripts for image reconstruction and visualization  

The .jac file used in this study is not included in this repository due to size restrictions but it is available [here](https://www.dropbox.com/s/7guwkrj52781b08/S07_texting.jac?dl=0).

## Folder Structure
```bash
dot-motion-correction/
│
├── MNI/                             # MNI anatomical structures and mesh files
│   ├── 10-5_Points.txt              # 10-5 positions of the asymmetric MNI152 atlas
│   ├── GMSurfaceMesh.mat            # Grey matter surface mesh
│   ├── HeadVolumeMesh.mat           # Head volume surface mesh
│   ├── LandmarkPoints.txt           # Cranial landmark coordinates
│   ├── ScalpSurfaceMesh.mat         # Scalp surface mesh
│   ├── TissueMask.mat
│   ├── TissueTypes.txt
│
├── main_script.m                    # Main MATLAB script for channel QC and preprocessing
├── S07_texting.nirs                 # NIRS data file of one subject
├── vol2gm.mat                       # Mapping matrix from volumetric to GM surface mesh
├── removeNoisyChannels.m            # Custom function to remove bad channels
├── dodWavelet.mat                   # Output of wavelet motion correction from HOMER
├── results_HW2_gruppo2.mat          # Final results as required by the assignment
├── solution_report.pdf              # Final project report with plots and results
├── assignment.pdf                   # Original assignment instructions
└── README.md
```

## Main Steps of the Analysis

### Channel Quality Check
- 3D plot of source-detector array configuration
- Histogram of source-detector distances
- Channel exclusion criteria:
  - Average intensity < 500 or > 1e10
  - SNR < 0  
- Bad channels flagged in SD.MeasListAct

### Preprocessing Pipeline
- Conversion to optical density (OD)
- Motion correction
  - Multiple methods tested: PCA, spline, wavelet
  - Wavelet filtering selected as best trade-off between aggressiveness and artifact removal
- Band-pass filtering (0.01 – 0.5 Hz)
- Block averaging from -2 to 40 s around stimulus onset

### Sensitivity and Image Reconstruction
- Computation of array sensitivity using the Jacobian matrix
- Mapping from volumetric nodes to GM surface using vol2gm
- HbO and HbR image reconstruction using the modified Beer-Lambert Law
- Regularization parameter: λ = 0.1
- Reconstructed maps visualized at 0 s, 10 s, and 18 s for both texting conditions

## Key Findings
- Wavelet filtering effectively reduced motion artifacts without suppressing neural signals
- Clear hemodynamic responses observed after block averaging
- Reconstructed activation maps showed distinct patterns for right vs. left hand texting
- As expected, HbO and HbR signals showed inverse trends over time and space
- Some channels (e.g., 39 and 57) appeared inactive across both conditions
