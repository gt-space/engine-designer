# engine-designer
Summary: Designs YJSP regen engine contour/regen circuit using NASA CEA, optimizes around thermal and structural requirements for peak thrust.
Point of Contact: Benjamin Woodman (benwoodman@gatech.edu)

# GETTING PYTHON

You will need python 3 on your machine. I’ve tested with both 3.7 and 3.8 and they work fine.  


I’ve been running this on Mac and this is what worked for me: 

Homebrew + Python Installation Guide: https://docs.python-guide.org/starting/install3/osx/ 

Run this command in terminal to get 3.8 working (optional): 

$ export PATH="/usr/local/opt/python@3.8/bin:$PATH" 


For windows, I’d start here if you need Python installed: https://www.python.org/downloads/  

# DOWNLOAD PACKAGES

With that taken care of, you should be good to install the packages listed in this guide: 

https://rocketcea.readthedocs.io/en/latest/quickstart.html 

Note for Mac: 

Everything went fine except when trying to pip3 install genericf2py 

I was met with: Preparing wheel metadata ... error 

Good news is that this doesn’t matter. I don’t know what causes this, but the code will still run normally. 

If you run into any other issues with installation let me know and we can work through them.

# USAGE

Now just navigate to the directory where you’ve saved the code and run: 

$ python dataCollectionScript.py 

This will generate an engine geometry and combustion gas property array based on your chosen parameters. To edit those input parameters, you will need to modify these scripts in the text editor of your choice.
