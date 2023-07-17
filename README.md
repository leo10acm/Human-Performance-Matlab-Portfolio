Biomechanics Data Analysis Project
This project focuses on the analysis of biomechanics data collected from force plates. The objective is to extract valuable insights from force, acceleration, velocity, and position data during various movements such as landings, drop jumps, and jumps. The analysis involves calculating key variables, including forces during the landing phase, drop height, jump height, breaking and propulsion phases, Reactive Strength Index (RSI), modified Reactive Strength Index (mRSI), and more.

Project Overview
Biomechanics data analysis plays a crucial role in understanding human movement and optimizing performance. By analyzing force plate data, we can gain insights into the forces exerted on the body, timing of different phases of movement, and quantify variables related to performance and injury risk.
 
File Structure
The project directory includes the following files:

I have included two files entitled DJ_N_*.c3d which contain the raw biomechanics data collected from force plates, including variables such as force, acceleration, velocity, and position.
The data was calculated from two different force plates, i.e. one foot per force plate
Calc_force2.m: This MATLAB script is used to analyze the force plate data and extract key variables of interest.
README.md: This file provides an overview of the project and instructions for usage.

Usage
To analyze the biomechanics data using the Calc_force2.m script:
Dowload the btk toolkit as you need to add it to your path to be able to extract the data from the c3d files http://biomechanical-toolkit.github.io/
Ensure that the data contains the relevant raw data from the force plates.
Open MATLAB and navigate to the project directory.
Load the Calc_force2.m script in MATLAB.

Execute the script to perform the data analysis.
The script utilizes various algorithms and calculations to process the force plate data. It segments the data into different phases of movement, calculates variables of interest (e.g., landing forces, jump height), and generates visualizations to aid in data interpretation.

Customization
The  script can be customized to meet specific analysis requirements. Users can modify the script to include additional calculations, apply different algorithms, or generate personalized visualizations based on their research goals or specific data characteristics.

Dependencies
The project relies on the following dependencies:

MATLAB 2021R

Contact
For any questions or inquiries regarding this project, please contact me.

We hope this biomechanics data analysis project proves to be a valuable resource for researchers, coaches, and scientists interested in understanding human movement and optimizing performance using force plate data.

Enjoy exploring the world of biomechanics data analysis!
For any questions or inquiries regarding this project, please reach out to leutrim.mehmeti@tum.de
We hope this project proves to be a valuable resource for sports scientists, researchers, and enthusiasts in the pursuit of enhancing sports performance through data-driven diagnostics.

Happy coding and sports science exploration!
