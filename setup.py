from setuptools import find_packages, setup

setup(name="satellite-constellation",
    version="0.2",
    description="Create Satellite Constellations around the Earth, Luna and planets in the solar system",
    long_description=open('README.rst').read(),
    author="Aidan O'Brien - Saber Astronautics Australia",
    author_email="aidan.obrien@saberastro.com",
    url="https://github.com/SaberAidan/SatelliteConstellationCreator",
    license='MIT',
    entry_points={
        'console_scripts': ['satellite-constellation=satellite_constellation.main:main'],
    },
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    package_dir={'satellite_constellation': 'satellite_constellation'},
    keywords="satellite-constellation satellite satellites constellation orbital mechanics",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics'
    ])
