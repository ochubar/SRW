# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example View Data File: View data stored in a file
# v 0.02
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from uti_plot import *
import optparse

if __name__=='__main__':
    p = optparse.OptionParser()
    p.add_option('-f', '--infile', dest='infile', metavar='FILE', default='', help='input file name')
    p.add_option('-e', '--e', dest='e', metavar='NUMBER', type='float', default=0, help='photon energy')
    p.add_option('-x', '--x', dest='x', metavar='NUMBER', type='float', default=0, help='horizontal position')
    p.add_option('-y', '--y', dest='y', metavar='NUMBER', type='float', default=0, help='vertical position')
    p.add_option('-l', '--readlab', dest='readlab', metavar='NUMBER', type='int', nargs=0, default=0, help='read labels from the file header (1) or not (0)')
    p.add_option('-j', '--joined', dest='joined', metavar='NUMBER', type='int', nargs=0, default=0, help='place different graphs jointly into one figure (1) or into separate figures (0)')

    opt, args = p.parse_args()

    if opt.readlab != 0: opt.readlab = 1
    if opt.joined != 0: opt.joined = 1

    if len(opt.infile) == 0:
        print('File name was not specified. Use -f option to specify the file name with path.')
        quit()

    #print(opt.joined)
    uti_plot_init('TkAgg')
    uti_data_file_plot(opt.infile, opt.readlab, opt.e, opt.x, opt.y, opt.joined)
    uti_plot_show()
