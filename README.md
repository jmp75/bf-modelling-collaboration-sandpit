
This repository is a way to exchange code snippets/samples

# Getting set up

## Windows

Outline to be fleshed out:

* Docker??
    * Entry point to deploying these.

if not docker:

* It is recommended you use an Anaconda distribution to speed up the installation process (packages);
* https://www.continuum.io/downloads. Select the installer for Windows. Download the "64-Bit Graphical Installer (537 MB)", Python 3.6. 
* Install Anaconda (TODO document options. Certainly pandas and numpy at the least)
* Add batch cmd files to have access to `pip`, `jupyter` and co from the prompt. 

install spotpy with some fixes:

```
cd ~/src/github_jm/spotpy
python3 setup.py develop --user
```

* Install/clone the "packages" for server and gr4j-sc.
* Attempt to get jupyter notebook server thing to work
* VS code

* Sample data to run through a calibration example (non-batch mode)
