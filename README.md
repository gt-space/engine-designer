# engine-designer
Summary: This codebase designs an engine contour and regenerative cooling circuit using NASA's CEA solver.

Point of Contact: Ben Woodman (benwoodman@gatech.edu)

NOTE: The program is currently under revision. The latest working version can be found under archive > engineDesigner_v0.3.

# Python Setup

I have run this code successfully on Mac and Linux.

If you are on a Windows machine, running this program will be a difficult process. I strongly reccomend a Linux partition as the Windows installation steps for rocketCEA are terrible. If you do not have one, this video is a good start: https://youtu.be/aKKdiqVHNqw. Be aware that things may not work first try depending on how your laptop is configured, so you may need to Google some additional things.

You will need python 3 on your machine. I’ve tested with both 3.7 and 3.8 and they work fine.  

FOR MAC:

The rocketCEA setup requires that you have Python setup with homebrew. The following guide can help you set this up.

Homebrew + Python Installation Guide: https://docs.python-guide.org/starting/install3/osx/

Run this command in terminal to get 3.8 working (optional):
$ export PATH="/usr/local/opt/python@3.8/bin:$PATH"

FOR LINUX:

Your Linux installation should already come with Python 3 installed. 
To check version: $ python3 --version
I recommend 3.7 or newer, as earlier versions of Python have not been tested. If your machine is running an older version, these steps should help you update it: https://dev.to/serhatteker/how-to-upgrade-to-python-3-7-on-ubuntu-18-04-18-10-5hab

# Download Packages

With that taken care of, you should be good to install the packages listed in this guide:

https://rocketcea.readthedocs.io/en/latest/quickstart.html

Note for Mac:

Everything went fine except when trying to pip3 install genericf2py

I was met with: Preparing wheel metadata ... error

To fix, try running with this option: pip3 install genericf2py==0.1.17

If you attempt a Windows installation and you end up figuring it out, please let me know and we can document useful information here.

# Design an Engine

You will now need to either clone or download the code. Cloning will give you git control, letting you update the code and pull updates when they are available, whereas downloading just gives you the latest files with no git control. Unless you know what you're doing and would like to set up a branch, I recommend just downloading the latest version and not worrying about accidentally pushing any changes.

Click the green "Clone or Download" button and then the "Download ZIP" button.

Finally you can run the code! Navigate to the directory where you’ve saved the code, then move into 
/engineDesigner_v1.0/contourDesigner

Now run the main script:
$ python3 main.py

This will generate an engine geometry and combustion gas property array based on your chosen parameters. To edit those input parameters, you can change the input parameters in main.py in the text editor of your choice.

# Run Regen Analysis
This part of the script generates a channel profile for the regenerative cooling system. Currently this code is being re-written, so the latest working version can be found in the archive folder. Go there and run:

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
