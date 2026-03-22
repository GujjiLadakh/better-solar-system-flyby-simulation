# Stellar Flyby + Solar System simulation

## Overview

This is a simulation of a stellar flyby on the Solar System using the IAS15 integrator from the REBOUND package.

## How it works

The parameters for the flyby (impact parameter, start distance, mass) can be defined at the start of the file. The program then produces a CSV and a matplotlib plot of the change in orbital properties (eccentricity, semi-major axis, semi-minor axis and angular momentum) of Earth over 20,000 years. The CSV contains changes for all planets but the plots are for a single planet. This can be altered to display change in orbital elements for any planet.

It will probably be better for you to remove the lines that use the logger to make output more readable.

The plots and CSV files will go into their respective folders.

You can alter the flyby parameters in the code as you so wish.

Licence: Apache Licence 2.0
