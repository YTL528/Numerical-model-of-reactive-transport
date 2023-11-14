# Numerical-model-of-reactive-transport
**Objective**

In environmental engineering, there is a need to characterize contaminant transport in groundwater over space and time. Crank-Nicolson is a second-order finite difference method of solving partial differential equations, which can be applied to advection-diffusion equations to model movement and change in concentration of a contaminant. This can be used to determine the rate of travel of contaminants, as well as the concentration at a specific location over time. The objective of this project is to create a numerical model of reactive transport using the Crank-Nicolson method applied to advection-diffusion equations and compare the results with analytical solutions and (potentially) other numerical methods.

**Methodology**

This project will be based on the following key equations:
1) Advection-diffusion equation
   [https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation]

$$
\theta \frac{\partial C}{\partial t}=\theta D_{x} \frac{\partial^2 C}{\partial x^2}-q_{x}\frac{\partial C}{\partial x} - kC
$$
where,
C - solute concentration
$\theta$ - effective porosity
$D_{x}$ - dispersive coefficient
$q_{x}$ - darcy velocity (q_x=v_x θ, where v_x is average groundwater velocity)
k - reaction rate coefficient

2) Crank-Nicolson method
   [https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method]

$$
\theta \frac{C^{n}_{i}-C^{n-1}_{i}}{\Delta t} = D_{x} \frac{C^{n}_{i+1}-2C^{n}_{i}+C^{n}_{i+1}+C^{n-1}_{i+1}-2C^{n-1}_{i}+C^{n-1}_{i+1}}{2 \Delta x^2} - q_{x}\frac{C^{n}_{i+1}-C^{n}_{i-1}+C^{n-1}_{i+1}-C^{n-1}_{i-1}}{4 \Delta x} - \frac{k}{2}(C^{n}_{i}+C^{n-1}_{i})
$$
where,
n - change in time
i - change in space

[both equations modified from Numerical modeling of contaminant transformation in a permeable reactive barrier, A.Rahman and Anurag 2021. https://link.springer.com/chapter/10.1007/978-981-16-5547-0_43]

**Implementation**

This project is implemented in Python with added tests and CI.
* `src/_.py`:
Python codes to implement the Crank-Nicolson method to solve the partial differential advection-diffusion equation.

* `tests/`:
Tests to verify that the written PDE solver works.

* `pyproject.toml`:
This file contains configuration settings for the build system, project metadata (name, version, dependencies), and tool specifications for linting and formatting using tools like Flake8, isort, Ruff, and pyupgrade.

* `.pre-commit-config.yaml`:
Defines a list of Git hooks using the pre-commit framework. It sets up hooks for checking and fixing trailing whitespace, using Ruff-specific checks, and applying Black code formatting.

* `noxfile.py`:
Contains Python code defining sessions for the Nox automation tool. In this case, it sets up a session named tests to install test dependencies and run the pytest test suite.

* `.github/workflow/ci.yml`:
Defines the continuous integration (CI) workflow using GitHub Actions. It has two jobs: formats and tests. The formats job checks for code formatting issues using pre-commit, and the tests job runs tests using different Python versions specified in a matrix.

