from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import time
import json

if __name__ == '__main__':

    myWalker = WalkerConstellation(4, 1, 1, 90, 35786, 0, 10)
    draw_walker_plotly(myWalker, sensor_regions=True, links=True)

    myFlower = FlowerConstellation(8, 1, 9, 1, 9, 0, 0, 2500, 10)
    draw_flower_plotly(myFlower, links=True)

    myStreets = SOCConstellation(1, 30, 36000, 10, [0,  60,  120], 0, 10)
    draw_soc_plotly(myStreets, links=False, sensor_regions=True)

    # file_name, sat_json = myWalker.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)
