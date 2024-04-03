![](/images/banner.png)
# Model Predictive Control with a Periodic Observer ($\Pi$-MPC)

This repository contains example code for [Perfecting Periodic Trajectory Tracking: Model Predictive Control with a Periodic Observer](https://arxiv.org/abs/2404.01550) by _Luis Pabon, Johannes Köhler, John Irvin Alora, Patrick Benito Eberhard, Andrea Carron, Melanie N. Zeilinger,_ and _Marco Pavone._

The repository is designed to facilitate understanding and implementation of $\Pi$-MPC using a straightforward example. We showcase the application of $\Pi$-MPC by tracking a periodic trajectory around a pendulum's unstable equilibrium. The nominal model is a linearized model that additionally has incorrect parameters (i.e. mass, length).

## Getting Started
Ensure Python 3 is installed on your system. The project dependencies are listed in requirements.txt.

### Installation
1. **Clone the repository** to your local machine.
2. **Create a virtual environment** (recommended) for the project to manage dependencies.
3. **Install dependencies** by navigating to the repository's root directory and executing:
```
pip install -r requirements.txt
```
This command installs all the necessary Python packages as specified in requirements.txt.

## Running the Example
Open the provided Jupyter notebook and run the cells sequentially. The notebook will guide you through:

- **Offline Design:** Preparation steps and design considerations for $\Pi$-MPC.
- **Online Operation and Simulation:** Simulation with the pendulum tracking a periodic trajectory, demonstrating the ability of $\Pi$-MPC to achieve minimal tracking errors despite model mismatch.

## Support
For questions or support, please file an issue on the GitHub repository page. Contributions and feedback are welcome.
