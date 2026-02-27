# engine-designer
Summary: This codebase designs an engine contour and regenerative cooling circuit using NASA's CEA solver.

Point of Contact:
Engine Designer V5+ Kieran Yarberry (kyarberry3@gatech.edu),  Saaketh Ramoju (sramoju@gatech.edu)
Engine Designer V1-4 Ben Woodman (benwoodman@gatech.edu), James Jutras (jjutras3@gatech.edu), Taylor Hampson (thampson3@gatech.edu)

NOTE: Version 1.0 (and now 2.0!) is finally up to date and working!

**UPDATE**: Modifications has been made in the documentation to enable code's functionality in Windows platform

# Setup

## FOR WINDOWS: 

If you ask programmers about their preferred OS, most would say Linux. But meh! Which operation system has the coolest looking GUI? Windows! Which OS has the most installable softwares? Windows! All hail Windows.

OKAY, so in a nutshell, modifications are now done to support the installations in Windows, thanks to Windows Subsystem for Linux (WSL) 2.0. The Windows Subsystem for Linux lets developers run a GNU/Linux environment -- including most command-line tools, utilities, and applications -- directly on Windows, unmodified, without the overhead of a traditional virtual machine or dualboot setup.

Here are the procedures to setup your Windows machine to run codes from this repo:

### First objective is to install WSL 2.0 in Windows "as administrator"
- To install WSL2 on Windows 10, open Command Prompt as admin and run `wsl –install` (or)  `wsl -install -d Ubuntu`
- However, if you have a specific distro (a Linux OS) in mind, use the command `wls -install -d <DISTRO_NAME>`, where <DISTRO_NAME> is to be replaced with the OS of your choice. (https://winaero.com/list-available-wsl-linux-distros-windows-10/)

### Next objective 
### Installing Python 3 in WSL

- Before installation of Python package, make sure that the build is up-to-date. Run the command 
`sudo apt -get update`
AND
`sudo apt -get upgrade`

- Now, install the Python packages using the command
`sudo apt install python3 python3-pip ipython3`

- To check it the installation was successful, you may run a quick command
`python3 --version`
The code should spit out `Python 3.10.8` as the output

Viewing Python scripts in the Linux interface is horrendous, trust me. Although a few advocate for using "vim" in Linux interface for easy operatibility, I prefer GUI (ofcourse! I use windows and not stone age OS!). So, for better GUI experience, neat indentations, easy debugging yata yata yata..., I would prefer installing Atom text editor, which is the "Hackable text editor for the 21st century"....or atleast that's how they market it.

### Installing Atom

Download Atom from the official website in Windows.

Add Atom to the windows environment variables.
- Run the following: win+r and type in SystemPropertiesAdvanced.exe
- Open: Environment Variables
- Add your Atom path:
`C:\Users\<user-name>\AppData\Local\atom\bin`

Go into the WSL and add an alias for Atom in your bashrc file:
- Open your bash configuration: `vim ~/.bashrc`
- Add to the end of the file and save/exit:
`alias atom="/mnt/c/Windows/System32/cmd.exe /c 'atom'"`
- Update your bash profile: `source ~/.bashrc`
Now you can use `atom . &` to open your python projects from WSL command line.

### Installing Jupyter Notebook [Optional - but advisable]

Install jupyter: `pip3 install jupyter`

Create alias to launch jupyter without browser from the WSL:
- Open your bash configuration: `vim ~/.bashrc`
- Add to the end of the file and save/exit:
`alias jupyter-notebook="~/.local/bin/jupyter-notebook --no-browser"`
- Update your bash profile: `source ~/.bashrc`
Now you can run a jupyter server jupyter-notebook and access the service with your browser from Windowslocalhost:8888.

You are all set now, spartans. Go to "[Download Packages](https://github-research.gatech.edu/YJSP/engine-designer/edit/master/README.md#download-packages)" now and follow the procedures there.

**NOTE** - The procedures for Linux OS now is applicable (everywhere) for Windows command prompt, provided that the above procedures are followed correctly. Henceforth, if you see a procedure stated for Ubuntu OS in the upcoming RocketCEA installation, it is also applicable for Windows OS.

## FOR MAC:

The rocketCEA setup requires that you have Python setup with homebrew. The following guide can help you set this up.

Homebrew + Python Installation Guide: https://docs.python-guide.org/starting/install3/osx/

Run this command in terminal to get 3.8 working (optional):
$ export PATH="/usr/local/opt/python@3.8/bin:$PATH"

## FOR LINUX:

Your Linux installation should already come with Python 3 installed. 
To check version: $ python3 --version
I recommend 3.7 or newer, as earlier versions of Python have not been tested. If your machine is running an older version, these steps should help you update it: https://dev.to/serhatteker/how-to-upgrade-to-python-3-7-on-ubuntu-18-04-18-10-5hab

# Download Packages

With that taken care of, you should be good to install the packages listed in this guide:

https://rocketcea.readthedocs.io/en/latest/quickstart.html

**NOTE** - When you reach `Try a quick test of the install by pasting the following into a command terminal:` statement under **Install RocketCEA** section in the RocketCEA installation documentation in the above specified website, use this command 
`python3 -c "from rocketcea.cea_obj import CEA_Obj; C=CEA_Obj(oxName='LOX', fuelName='LH2'); print(C.get_Isp())"` 

instead of 

`python -c "from rocketcea.cea_obj import CEA_Obj; C=CEA_Obj(oxName='LOX', fuelName='LH2'); print(C.get_Isp())"` 

Yes, just add "3" after `python` (i.e. `python3`), instead of just typing `python .......`

## Note for Mac:

Everything went fine except when trying to pip3 install genericf2py

I was met with: Preparing wheel metadata ... error

To fix, try running with this option: pip3 install genericf2py==0.1.17


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
