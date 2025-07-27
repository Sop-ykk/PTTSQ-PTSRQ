#encoding=utf-8
#this is used to convert (x,y) point to (x) point
#a simple 2D to 1D convertion

def convert(x: int, y: int, width: int)->int:
    return y*width+x