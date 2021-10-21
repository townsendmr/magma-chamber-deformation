# magma-chamber-deformation
This MATLAB code links magma chamber evolution to surface deformation. 

To cite this code, please reference: 

Townsend M., Linking surface deformation to thermal and mechanical magma chamber processes, Earth and Planetary Science Letters, 2021

Degruyter W. and Huber C., A model for eruption frequency of upper crustal silicic magma chambers, Earth and Planetary Science Letters, 2014

Run the code using the script "runCode.m"
Within runCode, the user chooses among the following setups: 
1. switch on/off effects of viscoelastic relaxation and heat conduction
2. choose between one of two viscoelastic models: shell model (Segall, 2010) or full viscoelastic crust (Bonafede and Ferrari, 2009)
3. choose one of three magma recharge conditions: constant recharge, pulse-like recharge, or recharge from a deeper connected chamber
Other input parameters specified in runCode include: chamber depth, volume, water content, critical overpressure for eruption, initial temperature, magma and crust properties
Model outputs include: time series of pressure, temperature, volume, mass, crystal fraction, gas fraction, and surface deformation
