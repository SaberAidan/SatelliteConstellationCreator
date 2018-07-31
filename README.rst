
===============================
Satellite-Constellation-Creator
===============================

Library for representing satellite constellations and ground stations


Installation
-------------
To install latest released version::
    pip install satellite-constellation
    
To install from github master branch::
    pip install https://github.com/SaberAidan/SatelliteConstellationCreator

For development::

    # fork https://github.com/SaberAidan/SatelliteConstellationCreator to YOUR_GITHUB
    # clone your repo locally
    git clone https://YOUR_GITHUB@github.com/YOUR_GITHUB/SatelliteConstellationCreator.git
    cd SatelliteConstellationCreator

    # add upstream remote
    git remote add upstream https://github.com/SaberAidan/SatelliteConstellationCreator.git

    # create a virtualenv
    python3 -m venv scc_venv
    source scc_venv/bin/activate

    # install for development
    pip install -r requirements.txt

Testing if installation is good::
    $ satellite-constellation --test

Features
--------

JSON serialisable orbital elements for satellite constellation creation. Individual Satellite element creation class.


Credits
-------

This package is a conversion of existing academic Matlab code created by Saber Astronautics Australia Pty. Ltd. for use in creating scenes in the Predictive Interactive Ground Interface software.



    
