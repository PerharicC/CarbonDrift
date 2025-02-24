from setuptools import setup, find_packages

setup(name="carbondrift",
      version="1.0.0",
      author="Crtomir E. Perharic Bailey",
      author_email="crtomir11@hotmail.com",
      url="https://github.com/PerharicC/CarbonDrift",
      packages=find_packages(),
      install_requires= [
          "haversine>=2.8.1",
          "numba>0.60.0",
          "numba-progress>1.1.0",
      ],
      entry_points = {
          "console_scripts": [
              "run_simulation=carbondrift.simulation.run:main",
              "plot_simulation=carbondrift.plotting.plot_run:main",
              "CD_gui=carbondrift.gui:main"
          ]
      }
      )
