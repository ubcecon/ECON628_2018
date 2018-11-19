# Final Project
## Dates
* Final project **proposal due** on **November 28th**
  * The writeup for this can be very short (e.g. 1/2 page or whatever), just so we can all agree on whether it is acceptable or not.  It will not be graded.  If we think that the proposal is inappropriate, or want some clarifications, we will tell you that week
* The final project is **due December 19th**

## Requirements
We will accept a wide-range of proposals on many different topics.  The key requirements are that:
* There is some element of practicing your coding skills (in Julia or R)
* The practice is on a tool with the intention of learning something which would be useful for doing your own research (e.g. replicating an estimation procedure, adding some new examples to a public package for some econometric procedure, doing a web scraping exercise to generate new datasource, etc.)
* It does not need to be narrowly focused on econometrics as long as it fulfills the other criteria.

## Size of Project and Grading Criteria
You should target roughly a similar length as the amount of code in one of the Jupyter notebooks used for lectures.  If you are contributing to an open-source project, then the contribution can be significantly smaller as it would require more work to iterate.

You will be able to get a decent grade on the project by just providing a single Jupyter notebook or R markdown file.  However, if you want to get the highest possible marks on the project, you will need to:
 * Create a thorough set of tests for your project
 * Carefully creating tests for each individual function rather than writing out the code all the way to the solution.
 * Ensure that it runs with continuous integration (i.e. Travis) and with code-coverage
 * Ensure that the project is reproducible (i.e. someone else can run the code without any modifications)
 * Demonstrate some use of generic programming (for Julia) or functional programming (for R)
 
 ## Some Ideas
 
 If you are able to make a small contribution to an existing open-source project, it is preferred.  Quite often, adding in a careful test and/or example to the document is a major contribution for packages.  Feel free to look around for ideas, and we can potentially contact the maintainers if you are interested.
 
 We have asked the maintainers of a few projects who would be excited for some contributions:
 * There are a few places you could help out on [FixedEffectModels.jl](https://github.com/matthieugomez/FixedEffectModels.jl) - where the maintainer is willing to help guide you a little
   * i.e. comparing the output with another software, and checking weird combinations of collinear variables / missing observations / fixed effects / clustered standard errors, and if the tests fail, try to understand where the problem comes from
   * Add another kind of standard errors (like heteroscedastic standard errors) following the template for the existing standard errors.
   * We would really love to see this built out and I think it could support a few of you interested in panel data.
 * [InteractiveFixedEffectModels.jl](https://github.com/matthieugomez/InteractiveFixedEffectModels.jl) could use more examples
 * [LeastSquaresOptim](https://github.com/matthieugomez/LeastSquaresOptim.jl) could use some econometric examples, if you can think of any
 * [MomentBasedEstimators.jl](https://github.com/gragusa/MomentBasedEstimators.jl) could use a few examples, and the maintainer would be willing to guide you 
