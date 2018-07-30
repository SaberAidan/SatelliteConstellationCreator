from distutils.core import setup

setup(name="SatelliteConstellationCreator",
      version="0.1",
      description="Create Satellite Constellations around the Earth, Luna and planets in the solar system",
      author="Aidan O'Brien - Saber Astronautics Australia",
      author_email="aidan.obrien@saberastro.com",
      url="https://github.com/SaberAidan/pyScene",
      packages=["satellite-constellation"],
      package_dir={'scenecreator': 'scenecreator'},
      keywords="satellite satellites constellation orbital mechanics",
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Licence :: MIT',
          'Natural Language :: English',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
      ],
      test_suite='tests')
