# engine-designer
Summary: This codebase designs an engine contour and regenerative cooling circuit using NASA's CEA solver.

Point of Contact: Ben Woodman (benwoodman@gatech.edu), James Jutras (jjutras3@gatech.edu)

NOTE: Version 1.0 (and now 2.0!) is finally up to date and working!

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

# Design an Engine (no Regen)

You will now need to either clone or download the code. Cloning will give you git control, letting you update the code and pull updates when they are available, whereas downloading just gives you the latest files with no git control. Unless you know what you're doing and would like to set up a branch, I recommend just downloading the latest version and not worrying about accidentally pushing any changes.

Click the green "Clone or Download" button and then the "Download ZIP" button.

Finally you can run the code! Navigate to the directory where you’ve saved the code, then move into 
/engineDesigner_v1.0/contourDesigner

Now run the main script:
$ python3 main.py

This will generate an engine geometry and combustion gas property array based on your chosen parameters. To edit those input parameters, you can change the input parameters in main.py in the text editor of your choice.

# Run Analysis with Regen, v1.0
This part of the script generates a channel profile for the regenerative cooling system. Go to the engineDesigner_v1.0 directory and run:

$ python3 main.py

This returns general engine design outputs as well as regenerative cooling ciruit results. To fine tune the parameters of the regen circuit like channel dimenisions, change these values in main.py. If everything is working properly you should see two consecutive plots show up when you run the code, one of the enigne contour and another of the copper wall temperatures.

Be aware that the current setup solves for the regen circuit assuming an abltative diverging section (to prevent excess coolant temperature increase).

I've done my best to comment these scripts, but some things might not be too clear. For further explanation feel free to message me (Ben Woodman).

# Optimize Design with Regen, v2.0
Version 2 restructures the regen design process. Instead of inputting channel geometry and solving for wall temperatures, the code inputs a maximum wall temperature and solves for the channel geometry required to satisfty that thermal load.

It can be run in the same way as v1.0, and all relevant parameters can be input into main.py 

For more detailed documentation, see the Nuclino: Subscale -> Flight Teams -> Propulsion -> Engines -> Kraken
