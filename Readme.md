# Distribution-Free Testing for Sample Comparison
## Directory
```
.
├── code
│   └── analysis
│       └── plots
├── data
└── readings
```

Simulations can be found in the analysis folder. The data folder contains the liver cancer data which will be used in the `CuMiDa analysis.R` file.  

Code for the painted turtles can be found in the `baseline.R` file in the *code* folder.

## Important files
- `test_fn.R` - File which contains a switch between MCM, MMCM and FR tests. This is used extensively in the analysis and simulations.
- `analysis/graphical power simulations.R` - File which contains the graphs of power for the different tests with increasing dimensions.

## Key R Packages used
- [MCM from multicross](https://rdrr.io/cran/multicross/man/mcm.html)  
- [MMCM from multicross](https://rdrr.io/cran/multicross/man/mmcm.html)  
- [FR-Match](https://github.com/JCVenterInstitute/FRmatch?tab=readme-ov-file)

## Additional notes
For ease of use: When pulling from the github repository, use your own directory when calling source functions in some scripts for better navigation.

## Acknowledgements
Special thanks to Somabha Mukherjee for guidance on graph-based feature selection and high-dimensional data analysis. I also thank the SBCB Lab for access to the CuMiDa dataset and the Moismann 1958 Department of Biology, University of Montreal, for the painted turtles data.  
