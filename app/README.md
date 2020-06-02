# pychemcurv web application

This directory contains a dash application that aims to use the pychemcurv 
package and visualize the geometrical or chemical atomic quantities mapped on 
the chemical structure of your system.

The application is available at this address: https://pychemapps.univ-pau.fr/mosaica/

Demo video:

[![youtube demo video](https://img.youtube.com/vi/q7UO5Gou-lw/0.jpg)](https://www.youtube.com/watch?v=q7UO5Gou-lw)

## Run the application locally

You can run the application locally once you have installed pychemcurv and
the dash dependencies of the application. The easiest way to this is to use
the `requirements.txt` or the `environment.yml` files to set up a python 
environment.

Using pip

    pip install -r requirements.txt

Using conda (recommended)

    conda env create -f environment.yml

This will create an environment named `curv` and install all dependecies. 
Then, to run the application, change to `pychemcurv/app` folder and run the
`app.py` file.

    conda activate curv
    cd pychemcurv/app
    python app.py

This will output something like this. Open the url to use the application.

    [user@computer] (curv) > $ python app.py
    Running on http://127.0.0.1:8050/mosaica/
    Debugger PIN: 065-022-191
    * Serving Flask app "app" (lazy loading)
    * Environment: production
    WARNING: This is a development server. Do not use it in a production deployment.
    Use a production WSGI server instead.
    * Debug mode: on

You can swith off the debug mode by setting `debug=False` on the last line of 
the `app.py` file.

## TODO

* Manage file format
* Manage periodic structure
* Show/Hide atom names = species + index
* Ball and stick representation
* Zoom fit the box at the beginning
