from ROOT import *
from array import array

def setPalette(name=None, ncontours=999):
  #Set a color palette from a given RGB list

  if name == "gray" or name == "grayscale":
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [1.00, 0.84, 0.61, 0.34, 0.00]
    green = [1.00, 0.84, 0.61, 0.34, 0.00]
    blue  = [1.00, 0.84, 0.61, 0.34, 0.00]

  elif name == "rainbow":
    stops = [ 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 ]
    red   = [ 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 ]
    green = [ 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 ]
    blue  = [ 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 ]

  elif name == "negative":
    stops = [ 0.00, 0.50, 1.00 ]
    red   = [ 0.00, 1.00, 1.00 ]
    green = [ 0.00, 1.00, 0.00 ]
    blue  = [ 0.00, 1.00, 0.00 ]

  elif name == "negblue":
    stops = [ 0.00, 0.50, 1.00 ]
    red   = [ 0.20, 1.00, 1.00 ]
    green = [ 0.20, 1.00, 0.20 ]
    blue  = [ 1.00, 1.00, 0.20 ]

  elif name == "rainbow2":
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [0.00, 0.00, 0.87, 1.00, 0.80]
    green = [0.00, 0.81, 1.00, 0.20, 0.00]
    blue  = [0.80, 1.00, 0.12, 0.00, 0.00]

  else:
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [0.00, 0.00, 0.87, 1.00, 0.51]
    green = [0.00, 0.81, 1.00, 0.20, 0.00]
    blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

  s = array('d', stops)
  r = array('d', red)
  g = array('d', green)
  b = array('d', blue)

  npoints = len(s)
  TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
  gStyle.SetNumberContours(ncontours)
