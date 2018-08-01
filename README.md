# MetaRowEchelon #

[![Build Status](https://travis-ci.org/jgoldfar/MetaRowEchelon.jl.svg?branch=master)](https://travis-ci.org/jgoldfar/MetaRowEchelon.jl)
[![Coverage Status](https://coveralls.io/repos/github/jgoldfar/MetaRowEchelon.jl/badge.svg?branch=master)](https://coveralls.io/github/jgoldfar/MetaRowEchelon.jl?branch=master)


This is an experiment in metaprogramming Julia for the automated solution of (well-posed) multiple assignments in the form of linear equations. Basically, the question is: can Gaussian elimination be computed ahead-of-time, and if so, should it be?

While I have possible applications in mind, right now the amount of copy/paste and general effectiveness suggest it should be considered a very rough WIP.

## TODO ##

* Deduplicate simplification/arithmetic for symbol logic

* Along with the above TODO, define unit tests for all transformations on Exprs and Coefficients to ensure good coverage and flexibility.

* Fix macro hygiene issue: output should pick up local variable bindings (?) and/or emit a function of the right-hand side (?) that allows extraction of the independent variables on the left-hand side. Right now, none of that works...

### Who do I talk to? ###

* Jonathan Goldfarb <jgoldfar@gmail.com>
