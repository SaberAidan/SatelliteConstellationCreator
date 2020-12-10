from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import time
import json

if __name__ == '__main__':

    # myWalker = WalkerConstellation(1, 1, 1, 90, 35786, 0, 10)
    # draw_walker_plotly(myWalker, sensor_regions=True, links=False)
    # print(myWalker.calculate_revisit(0, 0))

    myFlower = FlowerConstellation(31, 11, 30, 7, 10, 0, 80, 9000, 60)
    draw_flower_plotly(myFlower, links=False)
    print(myFlower.calculate_revisit(0,0))

    # myStreets = SOCConstellation(1, 20, 36000, 10, [0,  60,  120], 0, 10)
    # draw_soc_plotly(myStreets, links=False, sensor_regions=True)
    # print(myStreets.calculate_revisit(0,0))

    # file_name, sat_json = myWalker.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)
