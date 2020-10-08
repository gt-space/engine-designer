# engine-designer
Summary: Designs YJSP regen engine contour/regen circuit using NASA CEA, optimizes around thermal and structural requirements for peak thrust.
Point of Contact: Ben Woodman (benwoodman@gatech.edu)

# Getting Python

You will need python 3 on your machine. I’ve tested with both 3.7 and 3.8 and they work fine.  


I’ve been running this on Mac and this is what worked for me:

Homebrew + Python Installation Guide: https://docs.python-guide.org/starting/install3/osx/

Run this command in terminal to get 3.8 working (optional):

$ export PATH="/usr/local/opt/python@3.8/bin:$PATH"


For windows, I’d start here if you need Python installed: https://www.python.org/downloads/  

# Download Packages

With that taken care of, you should be good to install the packages listed in this guide:

https://rocketcea.readthedocs.io/en/latest/quickstart.html

Note for Mac:

Everything went fine except when trying to pip3 install genericf2py

I was met with: Preparing wheel metadata ... error

Good news is that this doesn’t matter. I don’t know what causes this, but the code will still run normally.

If you run into any other issues with installation let me know and we can work through them.

# Design an Engine

Now just navigate to the directory where you’ve saved the code and run:

$ python engine_designer.py 

This will generate an engine geometry and combustion gas property array based on your chosen parameters. To edit those input parameters, you can change the input parameters in engine_designer.py in the text editor of your choice.

# Run Regen Analysis
This part of the script generates a channel profile for the regenerative cooling system. To run it:

$ python regen.py

This returns an array defining the channel contour (the fin and channel widths at many points along the engine), and prints out the coolant outlet temperature and required inlet pressure.

How it works:
1. The user selects an initial gas side wall temp, min fin width, and desired coolant outlet temp. To edit these values, pass them in as parameters to the jacket object. You can see their ordering in the regenAnalysis.py init method.
2. get_coolant_temps() finds an estimate for coolant bulk temperature at throat, exit, and injector plane.
   a. Assumes a linear increase in temp (based on prior runs this is reasonable)
3. get_critical_dimensions() takes these coolant bulk temperatures and uses them to converge on the channel width for steady state 1D heat transfer through the wall (at each of those three stations).
   a. This width is contstrained by minimium fin and channel widths specified in the jacket object's initialization
   b. If the channel width can't get any smaller, the wall temp must increase to compensate.
4. get_all_dimensions() stitches together the geometry to generate a full contour
5. full_sim() simulates this full contour and determines the outlet temp and inlet pressure
   a. Solve for coolant convection coefficient
   b. Pick gas side wall temp and converge on equilibrium
   c. Solve for temperature increase
   d. Solve for pressure loss
6. The main loop checks if coolant outlet temp is too high, and if so it increases the gas side wall temp and repeats steps 2-5 until cooolant outlet temp is low enough or wall temps approach the limit defined in init.

The strain at the inner wall is also calculated, but it's not used as a metric for optimizing the design since it is well within acceptable values for our expected range of operating cycles.

I've done my best to comment these scripts, but some things might not be too clear. For further explanation feel free to message me.
