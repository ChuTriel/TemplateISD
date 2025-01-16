# File that holds the definition of an Instance and the derived ReducedInstance.
# Additionally actual global instances are defined here and added to global dictionaries.

class Instance:

    def __init__(self, n, k, w):
        self.isReduced = False
        self.n = n
        self.k = k
        self.w = w
        self.ncols = n
        self.nrows = n-k
        self.addRows = 1 # 1 pce is always possible

    def print_info(self):
        print("n: {}, k: {}, w: {}, ncols: {}, nrows: {}".format(self.n, self.k, self.w, self.ncols, self.nrows))

class ReducedInstance(Instance):

    def __init__(self, n, k, w, newN, addRows, weightDistro = ""):
        Instance.__init__(self, n, k, w)
        self.isReduced = True
        self.addRows = addRows
        self.newN = newN
        self.ncols = newN
        self.nrows = n-k
        self.weightDistro = weightDistro

    def print_info(self):
        print("n: {}, k: {}, w: {}, ncols: {}, nrows: {}, newN: {}, addRows: {}".format(self.n, self.k, self.w, self.ncols, self.nrows, self.newN, self.addRows))


inst156 = Instance(156, 125, 4)
inst240 = Instance(240, 192, 6)
inst381 = Instance(381, 305, 9)
inst482 = Instance(482, 386, 11)
inst640 = Instance(640, 512, 13)
inst751 = Instance(751, 601, 16)
inst808 = Instance(808, 657, 17)
inst923 = Instance(923, 739, 19)
inst982 = Instance(982, 786, 20)
inst1041 = Instance(1041, 833, 19)
inst1101 = Instance(1101, 881, 21)
inst1223 = Instance(1223, 979, 23)
inst1473 = Instance(1473, 1179, 27)
inst1665 = Instance(1665, 1332, 31)
inst1995 = Instance(1995, 1596, 37)
inst2129 = Instance(2129, 1704, 36)
inst2197 = Instance(2197, 1758, 37)
inst2265 = Instance(2265, 1812, 38)
inst2472 = Instance(2472, 1978, 42)
inst3488 = Instance(3488, 2720, 64)

MapFromNToStandardInstance = {
    156: inst156,
    240: inst240,
    381: inst381,
    482: inst482,
    640: inst640,
    751: inst751,
    808: inst808,
    923: inst923,
    982: inst982,
    1041: inst1041,
    1101: inst1101,
    1223: inst1223,
    1473: inst1473,
    1665: inst1665,
    1995: inst1995,
    2129: inst2129,
    2197: inst2197,
    2265: inst2265,
    2472: inst2472,
    3488: inst3488
}

inst156R = ReducedInstance(156, 125, 4, 92, 3, "11002")
inst240R = ReducedInstance(240, 192, 6, 128, 4, "00011220")
inst381R = ReducedInstance(381, 305, 9, 221, 7, "100112010102")
inst482R = ReducedInstance(482, 386, 11, 256, 8, "1100110010102030")
inst640R = ReducedInstance(640, 512, 13, 352, 11, "20121100000101011011")
inst751R = ReducedInstance(751, 601, 16, 416, 13, "210010120010010210011110")
inst808R = ReducedInstance(808, 657, 17, 384, 12, "01011212000201011100000300")
inst923R = ReducedInstance(923, 739, 19, 448, 14, "01012001031100120010102001100")
inst982R = ReducedInstance(982, 786, 20, 502, 16, "1001101000201111310011201000001")
inst1041R = ReducedInstance(1041, 833, 19, 576, 18, "101001011001001001011201011011110")
inst1101R = ReducedInstance(1101, 881, 21, 544, 17, "00011001021100210101110010200120010")
inst1223R = ReducedInstance(1223, 979, 23, 608, 19, "011201001000001001001111100012001311010")
inst1473R = ReducedInstance(1473, 1179, 27, 736, 23, "01100101000121110011201101001000011120100200010")
inst1665R = ReducedInstance(1665, 1332, 31, 832, 26, "01101012201100010113000100000010010011101121011001100")
inst1995R = ReducedInstance(1995, 1596, 37, 991, 31, "102111111200100000201200011110011000001010110000011122101000010")
inst2129R = ReducedInstance(2129, 1704, 36, 1056, 33, "1010110001100200011000001110001010011101201200110110110100100110110")

inst2197R = ReducedInstance(2197, 1758, 37, 1024, 32, "101011110200000001000011112010000010201011010100111200110100001012010")
inst2197R2 = ReducedInstance(2197, 1758, 37, 928, 29, "010001001010100000102000010002001210111101000010000011120110202110020")

inst2265R = ReducedInstance(2265, 1812, 38, 1120, 35, "00011001001010101100111111103001000000011110110111100000010001211110100")
inst2472R = ReducedInstance(2472, 1978, 42, 1224, 39, "110111010101011110210011101221000100101000011010101111100100101000000000110001")
inst3488R = ReducedInstance(3488, 2720, 64, 1920, 60, "0110001011100200011100011111100110111121110000101011100111111111111010010110101100010012001012100100001001000")

MapFromNToTemplateInstance = {
    156: inst156R,
    240: inst240R,
    381: inst381R,
    482: inst482R,
    640: inst640R,
    751: inst751R,
    808: inst808R,
    923: inst923R,
    982: inst982R,
    1041: inst1041R,
    1101: inst1101R,
    1223: inst1223R,
    1473: inst1473R,
    1665: inst1665R,
    1995: inst1995R,
    2129: inst2129R,
    2197: inst2197R, #inst2197R,
    2265: inst2265R,
    2472: inst2472R,
    3488: inst3488R
}
