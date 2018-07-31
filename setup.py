from distutils.core import setup

setup(name="satellite-constellation",
      version="0.1dev",
      description="Create Satellite Constellations around the Earth, Luna and planets in the solar system",
      long_description=open('README.rst').read(),
      author="Aidan O'Brien - Saber Astronautics Australia",
      author_email="aidan.obrien@saberastro.com",
      url="https://github.com/SaberAidan/SatelliteConstellationCreator",
      packages=["satellite_constellation"],
      package_dir={'satellite_constellation': 'satellite_constellation'},
      keywords="satellite-constellation satellite satellites constellation orbital mechanics",
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics'
      ])
