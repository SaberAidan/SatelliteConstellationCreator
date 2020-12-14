from satellite_constellation.Constellation import *
from plotter import visualiser
from orbital_analyzer.analysis import *
from satellite_constellation.Errors import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(8, 8, 2, 80, 35785, 0, 10)
    visualiser.draw_walker_plotly(myWalker, sensor_regions=True, links=True)

    # myFlower = FlowerConstellation(8, 1, 9, 1, 9, 0, 0, 2500, 10)
    # visualiser.draw_flower_plotly(myFlower, links=True)

    # myStreets = SOCConstellation(1, 20, 15000, 20, [0], 0, 10)
    # visualiser.draw_soc_plotly(myStreets, links=False, sensor_regions=True)

    myAnalyzer = orbital_analysis(myWalker)
    print(myAnalyzer.calculate_revisit(0,0))
    print(myAnalyzer.calculate_constellation_coverage(1))
    print(myAnalyzer.find_links())


    # file_name, sat_json = myWalker.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)

