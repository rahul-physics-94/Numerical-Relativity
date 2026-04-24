**Numerical Evolution of Reissner-Nordström Internal Instability**

This repository contains a Julia implementation of the numerical framework established by Simpson and Penrose in their seminal paper, "Internal instability in a Reissner-Nordström black hole".The project simulates the behaviour of a scalar field $H$ within the interior of a charged black hole, specifically investigating the stability of the Cauchy horizon using a double-null coordinate system $(u, v)$.


**Overview**

The code evolves a scalar field across a numerical grid to observe how gravitational perturbations behave as they approach the inner horizon. It utilises Newton's Method to solve the initial value problem for the tortoise-like coordinate $R$ along the initial $v$-null surface. Runge-Kutta (RK4): To evolve the metric function R(u, v) across the grid. Coupled ODE System: To solve for the scalar field $H$ and its derivatives (q, f) based on the Einstein-Maxwell-Scalar equations.

**Features High-Performance Numerical Loops:**

Leverages Julia’s @inbounds and @inline macros for fast grid computations. Dynamic Grid Setup: Uses a $1,000,000 \times 1,000,000$ grid for high-resolution tracking of the scalar field. Root Finding: A custom implementation of the Newton-Raphson method with central difference derivatives. Visualisation: Integration with PlotlyJS for generating slices and surface plots of the horizon evolution.

**Mathematical Context**

The simulation solves for the evolution of the function $R$ and the scalar field $H$ using the following logic: Initial Surface: Define $R$ along $u = u_0$ using the transcendental relation $f(R, v)$.Evolution: Use the derivative relation $\frac{dR}{du} = -0.5 \cdot g(u, R)$. Scalar Field: Evolve $H$ through intermediate variables $q$ and $f$, which represent the physical perturbations.

**Results**

The graphs are plotted against $H$ and $R$, however, the plotting should be done against $log(H)$ and $log(v)$ in order to make the oscillations on. If you plot against $R$, the physical perturbation will appear less oscillatory; however, it will blow up near the Cauchy horizon. 

