from satellite_constellation.visualiser import *

def main():
    draw_walker(WalkerConstellation(8, 2, 1, 45, 1000, 0, 20))
    draw_flower(FlowerConstellation(37, 18, 57, 6, 19, 0, 0, 19702, 30))

if __name__ == '__main__':
    main()
