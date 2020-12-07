from satellite_constellation.visualiser import *

if __name__ == '__main__':
    myWalker = WalkerConstellation(3,3,1, 45, 20000, 0, 45)
    draw_walker(myWalker)
    # for satellite in myWalker.satellites:
    #     print(satellite.right_ascension, satellite.ta, satellite.perigee)
    # draw_walker(WalkerConstellation(8, 4, 2, 60, 36000, 0, 50))
    # draw_walker(WalkerConstellation(12, 6, 2, 60, 36000, 0, 50))
    # draw_flower(FlowerConstellation(4, 1, 4, 1, 4, 0, 0, 600, 0))
