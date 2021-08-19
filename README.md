

#### PRUC:
This is the source code and dataset for the paper PRUC: P-Regions with User-Defined Constraint.

#### Environment:
The code is written in Java 14. The original developing platform is IntelliJ IDEA. The maven version is 3.6.3. pom.xml stores all the denpencies for this project.

#### Quick Start:
To run the code, locate the TestGeneral Class under GSLO package and run the main method. There is a one-to-one relationship between the static methods in this Class and the experiments in the paper, which is documented in the annotation of the methods.

#### Explanation:
The util package stores the common Class that are used among all the methods.
The GSLO package stores all the Class related to Global Search and Local Optimization.
The baseline package stores our baseline competitors, it is further divided into greedy, SKATER and SKATERCON.

#### Dataset:
The datasets used in experiment come from two sources,  [TIGER/Line Shapefile](https://catalog.data.gov/dataset/tiger-line-shapefile-2016-series-information-for-the-current-census-tract-state-based-shapefile "TIGER/Line Shapefile") dataset and [2000 Health, Income and Diversity](https://geodacenter.github.io/data-and-lab/co_income_diversity_variables/ "2000 Health, Income and Diversity") dataset. The original datasets are slightly modified by removing the island areas. The [diversity](https://github.com/Yongyi-Liu/PRUC/tree/master/DataFile/30K "diversity") folder under DataFile folder refers to the 2000 Health, Income and Diversity dataset. All other datasets are from TIGER/Line Shapefile dataset. Different states in TIGER/Line Shapefile dataset are merged to form larger datasets.

#### Contact:
If you have any question with this project, please feel free to reach me by yliu786 [at] ucr.edu
