from satellite_constellation.Constellation import *
from plotter import visualiser

#
if __name__ == '__main__':
    myWalker = WalkerConstellation(30, 15, 1, 60, 35786, 0, 20)
    visualiser.draw_walker_plotly(myWalker, sensor_regions=True, links=True)

    myStreets = SOCConstellation(1, 5, 15000, 60, [0], 0, 100)
    visualiser.draw_soc_plotly(myStreets, sensor_regions=True, links=True)

