from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(8, 8, 2, 80, 35785, 0, 10)
    draw_walker_plotly(myWalker, sensor_regions=True, links=False)
    #
    # myFlower = FlowerConstellation(8, 1, 9, 1, 9, 0, 0, 2500, 10)
    # draw_flower_plotly(myFlower, links=True)

    # myStreets = SOCConstellation(1, 20, 15000, 20, [0, 60, 120], 0, 10)
    # draw_soc_plotly(myStreets, links=False, sensor_regions=True)

    # file_name, sat_json = myWalker.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)

    print(myWalker.calculate_constellation_coverage())
    # print(myStreets.calculate_constellation_coverage())
