# ECON628 - Fall 2018

## Topics in Applied Econometrics I



- **Instructors:** jesse.perla@ubc.ca and schrimpf@mail.ubc.ca
- **Office Hours:** TBD
- **Teaching Assistant:** Jasmine Hao, haojasmine@gmail.com

**Course Description**
This is a graduate topics course in applied econometrics and computational economics.


A key purpose of this class is to teach specific techniques, algorithms, and tools to ensure that students write robust, correct, and tested code - and hopefully open the research opportunities for students to move to the cutting edge of quantitative economics.  Beyond the necessary algorithms and new programming languages, another goal is to ensure that economists are using modern software engineering tools to allow collaboration -  as most projects involve multiple coauthors and research assistants.  Finally, all of the practice in this class will be done with the goal of showing how code used in research can be shared as open-source with the economics research community (and the scientific computing community as a whole).

**Grading**
While students will all have significant experience with Matlab, the only way to learn how to apply new programming languages and methods to economic problems is practice.  To aid in this, a significant portion of the grade will regular problem sets.  The remainder of the grade will be a computational project.

- Nearly weekly problem sets: 50%
- Final Project: 50%

While the problem sets will be frequent, many will be short to force practice - and hence they will not be weighting uniformly in the grade.

## Course Parts and Programming Languages

The course will be taught in 3 parts, starting with some econometric theory and assignments, and then moving to 
1. Econometric topics (Paul Schrimpf)
2. Introduction to Julia and scientific computing (Jesse Perla)
3. Dynamic Programming applications (Jesse Perla)
4. Solving the econometric models (Paul Schrimpf)

In part (1), we will use R for the majority of examples (where we may give you the option to use Matlab and/or Stata for some examples).

For parts (2) to (4) of the class, we will be using Julia.  While it would be nice to start the class by learning Julia, the package ecosystem and development tools are in the process of stabilizing after its 1.0 release.  While we will be able to use it when the 2nd part of the class starts, we suggest that you **avoid using Julia** until we give you the green-light.


## Econometric Topics

1. Linear panel data
     - Fixed/Random effects, first differencing, strict exogeneity
     - Potentially: dynamic panel

2. Extremum estimators & optimization
     - Review of extremum estimators 
     - Application to limited dependent variables
     - Introduction to optimization algorithms
       
3. Inference for extremum estimators
     - Review of usual asymptotic distribution
     - Identification robust inference from inverting LR-tests
     - Bootstrap
       
4. Unobserved heterogeneity and simulation based inference
     - Numeric integration, especially monte carlo
     - BLP
       
5. Nonparametric methods
     - Kernel, sieve
     - Semiparametric models
       
6. Machine learning methods in econometrics
     - As an alternative to traditional nonparametric methods      
     - For high dimensional data 
     - Instrument selection in BLP
       
7. (Hopefully) Dynamic structural models 
     - Dynamic discrete choice, dynamic games
     
## Computational Topics

1. Introduction to Julia 
   - Learning the Julia programming language, with simple applications
  
2. Software engineering tools: source-code control, unit testing, and continuous integration
   - Git and Github version tracking, diffs, collaboration, Pull Requests, etc.
    - Reproducible environments: package managers, and virtual environments
    - Unit and regression testing frameworks, benchmarking, and continuous-integration
3. Translating "white-board" models to code
   - Generic and Functional programming
   - Muliple dispatch
   - As an application, will discuss the design of the https://github.com/JuliaStats/Distributions.jl package for generic specification of univariate distributions, and the https://github.com/econtoolkit/Expectations.jl package for computing expectations of functions over those distributions.
4. Function approximation and quadrature for solving dynamic economic models
    - Interpolation with https://github.com/JuliaMath/Interpolations.jl
    - Newton–Cotes, Gaussian, and adaptive Gauss-Kronod quadrature with https://github.com/ajt60gaibb/FastGaussQuadrature.jl, https://github.com/JuliaMath/QuadGK.jl, and https://github.com/econtoolkit/Expectations.jl
    - Spectral solutions to functional equations with https://github.com/JuliaApproximation/ApproxFun.jl/
5. Sparsity and auto-differentiation for scaling to large problems
   - Dealing with large dimensions
   - Sparse matrices, banded matrices, etc. in Julia
    - Forward/reverse/etc. auto-differentiation with https://github.com/JuliaDiff/ForwardDiff.jl, https://github.com/JuliaDiff/ReverseDiff.jl, and maybe https://github.com/JuliaDiff/Capstan.jl
    - Iterative and matrix-free solutions to linear systems, eigensystems, and sparse linear least squares: https://github.com/JuliaSmoothOptimizers/Krylov.jl and https://github.com/JuliaMath/IterativeSolvers.jl
    - Anderson Fixed-point iteration in https://github.com/JuliaNLSolvers/NLsolve.jl/


6. Markov-chains and asset-pricing applications
   - We will introduce concepts of finite markov chains and use them for some asset-pricing examples.
   - https://lectures.quantecon.org/jl/finite_markov.html and https://lectures.quantecon.org/jl/markov_asset.html
7. Dynamic programming
    - Select topics in solving dynamic models in Julia, such as https://lectures.quantecon.org/jl/mccall_model.html and https://lectures.quantecon.org/jl/optgrowth.html
