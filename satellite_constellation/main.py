from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(12, 6, 1, 60, 36000, 0, 45)
    draw_walker_plotly(myWalker, links=True, sensor_regions=True)
    # file_name, sat_json = myWalker.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)

    myFlower = FlowerConstellation(3, 1, 4, 1, 4, 0, 0, 600, 0)
    draw_flower_plotly(myFlower, links = False)
    # file_name, sat_json = myFlower.as_pigi_output()
    # with open(file_name, "w") as outfile:
    #     json.dump(sat_json, outfile, indent=1)

    myStreets = SOCConstellation(1, 20, 12000, 15, [0,  60,  120], 0, 10)
    draw_soc_plotly(myStreets, links=False, sensor_regions=True)



