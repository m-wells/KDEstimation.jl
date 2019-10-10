# KDEstimation (Kernel Density Estimation)
[![Build Status](https://travis-ci.com/m-wells/KDEstimation.jl.svg?branch=master)](https://travis-ci.com/m-wells/KDEstimation.jl)
[![codecov](https://codecov.io/gh/m-wells/KDEstimation.jl/branch/master/graph/badge.svg?branch=master)](https://codecov.io/gh/m-wells/KDEstimation.jl)
[![Coverage Status](https://coveralls.io/repos/github/m-wells/KDEstimation.jl/badge.svg?branch=master)](https://coveralls.io/github/m-wells/KDEstimation.jl?branch=master)

The purpose of this package is to provide a general framework for implementing Kernel Density Estimation methods.

## Univariate KDE
The density estimator

![equation](https://latex.codecogs.com/svg.latex?&space;\hat{f}(x)=\frac{1}{n}\sum_{i=1}^nK\left(\frac{x-x_i}{h}\right))

where

* ![equation](https://latex.codecogs.com/svg.latex?\inline&space;\hat{f}(x)) is the estimator
* ![equation](https://latex.codecogs.com/svg.latex?\inline&space;K(u)) is the kernel function
* ![equation](https://latex.codecogs.com/svg.latex?\inline&space;h) is the bandwidth

can be evaluated using one of three implemented methods.

* `Direct() # O(N^2) where N is the sample size`
* `Binned() # O(M^2) where M is the number of evaluation points (M=4096 by default)`
* `Direct() # O(MlogM) where M is the number of evaluation points (M=4096 by default)`

## Multivariate KDE (work in progress)

## Kernels implemented
Here is a link to the [relevant wikipedia article](https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use)
<table>
    <tr>
        <th>Kernel</th>
        <th><img src="https://latex.codecogs.com/svg.latex?&space;K(u)" /></th>
        <th>Support</th>
    </tr>
    <tr>
        <td>Biweight</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{15}{16}(1-u^2)^2" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /></td>
    </tr>
    <tr>
        <td>Cosine</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{\pi}{4}\cos(\frac{\pi}{2}u)" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /></td>
    </tr>
    <tr>
        <td>Epanechnikov</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{3}{4}(1-u^2)" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /></td>
    </tr>
    <tr>
        <td>Logistic</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{1}{e^u+2+e^{-u}}" /></td>
        <td></td>
    </tr>
    <tr>
        <td>Normal</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{1}{\sqrt{2\pi}}\exp\left(-\frac{1}{2}u^2\right)" /></td>
        <td></td>
    </tr>
    <tr>
        <td>SymTriangularDist</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=1-|u|" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /></td>
    </tr>
    <tr>
        <td>Triweight</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{35}{32}(1-u^2)^3" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /> </td>
    <tr>
        <td>Uniform</td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;K(u)=\frac{1}{2}" /></td>
        <td><img src="https://latex.codecogs.com/svg.latex?\inline&space;|u|\le1" /></td>
    </tr>
</table>

This package uses [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) to suppy kernels such that ![equation](https://latex.codecogs.com/svg.latex?\inline&space;K_h\left(x-x_i\right)=\textnormal{pdf}(D(x_i,h),x)) where ![equation](https://latex.codecogs.com/svg.latex?\inline&space;K_h(u)=\tfrac{1}{h}K(\tfrac{u}{h})) and ![equation](https://latex.codecogs.com/svg.latex?\inline&space;D) is one of the kernels listed in the table above.

__Note:__ for the Uniform distribution, [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) defines `(loc,scale) = (a, b-a))` where `a` and `b` are the bounds lower and upper bounds, respectively.
This package accounts for this inconsistancy by evaluating the Uniform kernel as ![equation](https://latex.codecogs.com/svg.latex?\inline&space;\textnormal{pdf}(\textnormal{Uniform}(x_i-\frac{h}{2},x_i+\frac{h}{2}),x))


## Bandwidth selection via Least Squares Cross Validation
The objective function to minimize is given by

![equation](https://latex.codecogs.com/svg.latex?&space;LSCV(h)=\int\hat{f}^2_h(x)dx-2\frac{1}{n}\sum_i\hat{f}_{h,-i}(X_i))

where

![equation](https://latex.codecogs.com/svg.latex?&space;\hat{f}_{h,-i}(X_i)=\frac{1}{(n-1)h}\sum_{j\ne&space;i}K\left(\frac{X_i-X_j}{h}&space;\right&space;))

This has also been implemented using `Direct`, `Binned`, and `FFT` methods.

# Example usage

    julia> using KDEstimation, Distributions

    # generate random data
    julia> using Random: seed!

    julia> seed!(1234);

    julia> x = randn(1000);

    julia> lscv_res = lscv(Normal,x,FFT())
    LSCV{Normal,FFT(4096),1}
    Results of Optimization Algorithm
     * Algorithm: Golden Section Search
     * Search Interval: [0.128239, 0.195882]
     * Minimizer: 1.688927e-01
     * Minimum: -2.788450e-01
     * Iterations: 34
     * Convergence: max(|x - x_upper|, |x - x_lower|) <= 2*(1.5e-08*|x|+2.2e-16): true
     * Objective Function Calls: 35

    julia> using Plots; pyplot()

    julia> plot(lscv_res)
    
![plot(lscv_res)](docs/lscv_normal_fft4096.png)
    
    julia> h = minimizer(lscv_res)
    0.16761124952746273

    julia> fkde = kde(Normal, h, x, FFT())
    KDE via FFT Evaluation Method (4096 bins)
      Kernel: Normal
      Bandwidth: 0.16761124952746273

    julia> fkde(0.3)
    0.3694698807401497

    julia> plot(fkde)

![plot(fkde)](docs/kde_normal_fft4096.png)

# Visualization using `Plots.jl`

# Further Reading
This work has been heavily influenced by Artur Gramacki's book "Nonparametric Kernel Density Estimation and Its Computational Aspects" 
