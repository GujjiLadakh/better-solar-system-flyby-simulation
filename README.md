# Stellar Flyby + Solar System simulation

## Overview

This is a simulation of a stellar flyby on the Solar System using the IAS15 integrator from the REBOUND package.

## How it works

When you run the code, you are presented with 2 options: a single simulation or a sweep.

### Single Simulation

The parameters for the flyby (impact parameter, start distance, mass) can be defined at the start of the file. The program then produces a CSV and a matplotlib plot of the change in orbital properties (eccentricity, semi-major axis, semi-minor axis and angular momentum) of Earth over 20,000 years. The CSV contains changes for all planets but the plots are for a single planet (default Earth). This can be altered to display change in orbital elements for any planet.

It will probably be better for you to remove the lines that use the logger to make output more readable.

The plots and CSV files will go into their respective folders.

### Sweep

The sweep option will allow you to run multiple simulations in parallel. The number of simulations run will be equal to the length of the SWEEP_BS array mutiplied by the length of the SWEEP_MASSES array. The default is 400 simulations. They will be run in parallel on all of your computer's threads so the CPU usage of your computer will skyrocket to 100%.

I recommend having a good CPU and turning all open programs off before running a sweep. Also ensure your device is on a hard surface and there is proper ventilation. Refer to Google for more information. I am not responsible if your computer becomes damaged as a result of running a sweep simulation.

If successful, a contour map showing the mass and impact parameter will be shown with a cyan dashed line indicating the threshold mass and impact parameter that will result in the ejection of a planet from the Solar System.

You can alter the flyby parameters in the code as you so wish.

Licence: Apache Licence 2.0
