# Improved testing inferences for beta regressions with parametric mean link function

The computer code can be used to perform data analysis and Monte Carlo simulations.

## Authors

Cristine Rauber, Francisco Cribari-Neto, Fábio M. Bayer

## Programming language 

[Ox](https://www.doornik.com/) The console version of Ox is freely distributed for academic use. It is currently available for Windows, Linux, macOS, and several Unix platforms.

## Description

### Data analysis 

The program is used to estimate the beta regression model with parametric mean link function presented in Section 5.1. 

### Monte Carlo simulation 1 

The program is used to compute the null rejection rates of the tests of H0: \gamma_2 = \gamma_3 = 0 (fixed precision) in the varying precision beta regression model with parametric mean link function; \gamma_1 = \log(30).

### Monte Carlo simulation 2 

The program is used to compute the null rejection rates of the tests of H0: \lambda = 1 in the fixed precision beta regression model with parametric mean link function; \phi = 30.

### Monte Carlo simulation 3 

The program is used to compute the null rejection rates of the tests of H0: \beta_4 = 0 in the fixed precision beta regression model with parametric mean link function; \phi = 30.

## Files

### Data analysis

* gasoline.mat: data
* appPrater.ox: Ox code for reproducing some of the results in Section 5.1
* appPrater.out: Ox output

### Monte Carlo simulation 1

* log30fixedprec.ox: Ox code
* log30fixedprec.out: Ox output

### Monte Carlo simulation 2

* phi30lambda.ox: Ox code
* phi30lambda.out: Ox output

### Monte Carlo simulation 3

* phi30beta4.ox: Ox code
* phi30beta4.out: Ox output

## Execution 

### Data analysis 

1. Command line: ```oxl appPrater.ox > appPrater.out``` (On some systems, one must use oxl64 instead of oxl).  

2. IDE, via OxEdit: Open the file appPrater.ox in OxEdit (Ox IDE) and then click on: Run -> Ox

### Monte Carlo simulation 1 

1. Command line: ```oxl log30fixedprec.ox > log30fixedprec.out``` (On some systems, one must use oxl64 instead of oxl).  

2. IDE, via OxEdit: Open the file log30fixedprec.ox in OxEdit (Ox IDE) and then click on: Run -> Ox

### Monte Carlo simulation 2

1. Command line: ```oxl phi30lambda.ox > phi30lambda.out``` (On some systems, one must use oxl64 instead of oxl).  

2. IDE, via OxEdit: Open the file log30fixedprec.ox in OxEdit (Ox IDE) and then click on: Run -> Ox

### Monte Carlo simulation 3

1. Command line: ```oxl phi30beta4.ox > phi30beta4.out``` (On some systems, one must use oxl64 instead of oxl).  

2. IDE, via OxEdit: Open the file log30fixedprec.ox in OxEdit (Ox IDE) and then click on: Run -> Ox

## Contact 

Cristine Rauber (rauberoliveira at gmail.com)

## License
[MIT](https://choosealicense.com/licenses/mit/)
