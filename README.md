# Imagine regression analysis using Excel, Geogebra, Graph, etc. except it gets done much faster (sciplot-edition)
+ [Getting started](https://github.com/Saltworker/RegressionAnalysis_CAS-plotting#to-run-regressionanalysis_sciplotcpp)
![Screenshot](Screenshots/Screenshot_1.0_test_sciplot-edition.jpg) 

**NOTE: The program Graph is simply used for testing the results of the regression model & coefficient of determination R^2**

## RegressionAnalysis_sciplot
Regression analysis program utilizing sciplot to plot the results. Using .txt files as input, the following can be made:
  - Regression model(s) and coefficient(s) of determination R^2
     1. Linear
     2. Exponential
     3. Logarithmic
     4. Power
     5. Second-degree polynomium [degree of polynomium adjustablue in code]
 - Plotting
    1. Regression models
    2. Residual plot hereof
 - Other statistical descriptors (x- and y-axis respectively
    1. Domain & range; quarters of dataset
    2. Average/mean
    3. Variance, standard deviation, sample standard deviation
 - Correlation coefficient r of ONLY linear regression model

**RegressionAnalysis_sciplot.cpp imports the following files:**
- ...sciplot-master/sciplot/sciplot.h
- [FILENAME].txt

![Screenshot](Screenshots/Screenshot_1.0-repository_sciplot-edition.jpg) 

The format of this [FILENAME].txt can be seen with the available example datasets. Compatible data (in the current state of the program) can be boiled down to the following datapoints:

| **x**  | **y** |
| --- | --- |
| 1  | 1  |
| 1,  | 1  |
| 1.0  | 1.0  |
| 1.0,  | 1.0  |
| 1,0  | 1,0  |
| 1,0,  | 1,0  |

*NOTE: Do not forget to enter SPACEBAR or TAB in between the x- and y-axis for each line in [FILENAME].txt* 
## PRO TIP regarding datasets from excel, geogebra, graph, etc
1. DATASET -> CTRL + C
2. REPOSITORY -> Make/insert .txt file
3. .txt file -> CTRL + V

You're welcome
## To run regressionAnalysis_sciplot.cpp:
### Install gnuplot
Once this is installed, you'll be able to run the executable attached to the release "1.0 Sciplot-edition"

https://sourceforge.net/projects/gnuplot/

**NOTE: Make sure to set {...\gnuplot\bin} as an environment variable**
### Install sciplot
https://github.com/sciplot/sciplot

Include the downloaded "sciplot-master" following under "include Directories": {...\sciplot-master}
