from satellite_constellation.visualiser import *

if __name__ == '__main__':
    draw_walker(WalkerConstellation(8, 4, 2, 60, 36000, 0, 50))
    draw_flower(FlowerConstellation(4, 1, 4, 1, 4, 0, 0, 600, 0))
