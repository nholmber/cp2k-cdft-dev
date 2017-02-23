#!/usr/bin/python
# -*- coding: utf-8 -*-

# Generates the CP2K CDFT Branch Logo
# Adapted from gen_cp2k_logo.py
# author: Nico Holmberg

import os
import sys
from PIL import Image

#-------------------------------------------------------------------------------
def main():
    gen_povray()
    crop_fraction = [0.28, 0.31]
    # run povray
    for res in (200, 800):
        cmd = "povray -D +UA +H%d +W%d +Q11 +A +Ocp2k_cdft_logo_%d.png logo_cdft.pov"%(res,res, res)
        print("Running: "+cmd)
        os.system(cmd)
        im = Image.open("cp2k_cdft_logo_%d.png" %(res)) 
        im.crop((0, int(crop_fraction[0]*res), res, int((1.0-crop_fraction[1])*res))).save("cp2k_cdft_logo_%d_cropped.png" %(res))


#-------------------------------------------------------------------------------
def gen_povray():
    txt = ["XXXX XXXX      X  X   OOO OOO  OOOO OOOOO",
           "X    X  X XXXX X X   O    O  O O      O  ",
           "X    X  X    X XX    O    O  O O      O  ",
           "X    XXXX    X XX    O    O  O OOOO   O  ",
           "X    X    XXXX X X   O    O  O O      O  ",
           "XXXX X    X    X  X   OOO OOO  O      O  ",
           "          X                              ",
           "          XXXX                           "]

    coords = []
    types = []
    for z in range(3):
        for y, line in enumerate(txt):
            for x, c in enumerate(line):
                if(c=="X" or c=="O"):
                    coords.append((x-9,-y+3, -z))
                    types.append(c)

    output  = '#include "colors.inc"\n'
    output += '#include "textures.inc"\n'
    output += 'global_settings { assumed_gamma 1.0 }\n'

    output += 'light_source { <-100, -70, -300> White }\n'
    output += 'light_source { <40, -20, 50> White }\n'

    output += 'camera {\n'
    output += '  location <11.52, -9.1, -60.6>\n'
    output += '  sky   <0.15, 1, 0>\n'
    output += '  angle 40\n'
    output += '  look_at <11,0,0>\n'
    output += '}\n'

    for i, c in enumerate(coords):
        output += "sphere { <%f, %f, %f>, 0.73 \n"%c
        output += "  texture {\n"
        if types[i] == 'X':
            output += "    finish { Shiny ambient 0.3 }\n"
            output += "    pigment { color rgb <1, 0.2, 0> }\n"
        else:
            output += "    finish { Shiny ambient 0.1 diffuse 0.3 }\n"
            output += "    pigment { color rgb <8./255., 180./255., 112./255.> }\n"
        output += "  }\n"
        output += "}\n"

    f = open("logo_cdft.pov", "w")
    f.write(output)
    f.close()


#-------------------------------------------------------------------------------
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
