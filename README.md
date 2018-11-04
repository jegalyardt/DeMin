# DeMin
A Gradient-Free, Self-Adaptive Single Objective Optimization Utility

This is a Python 2.x module designed to perform single-objective global optimization using a gradient-free approach based on the Differential Evolution algorithm (Storn and Price, 1997). It adapts the optimization parameters during the course of a run via a superimposed multi-objective optimization seeking a specific tradeoff between _exploration_ of the parameter space of the user's problem and _eploitation_ of promising regions of parameter space already found. DeDriver.py is the user's entry point; documentation for the various command line parameters is available via 

```./DeDriver.py --help``` 

Storn, R.; Price, K. (1997). "Differential evolution - a simple and efficient heuristic for global optimization over continuous spaces". Journal of Global Optimization. 11: 341 - 359.
