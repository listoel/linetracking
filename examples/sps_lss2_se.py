# TODO: Change ZS to effective length with doubleappdrift at start and
# end for correct anode length

import sys
sys.path.insert(0, '/afs/cern.ch/user/l/listoel/Desktop/pythonpacks')

import linetracking as lt
import numpy as np

##################################################
#                                                #
#        LSS2 parameters                         #
#                                                #
##################################################
# Fixed hardware parameters
QD219holeRadius = 0.0575  # TODO, change to actual number,
                          # this is a guess from a plot...
QD219holeFieldAxis = 0.30094  # From MADX/Sicopex
QD219holeApAxis = 0.247  # Guess from Sicopex documentation
QDArad = 0.1065  # Tentative
QFArad = 0.1065  # Tentative
r_tce = 0.0645
tpst1_thick = 0.0046
tpst2_thick = 0.0052
tpst_extrDiam = 0.04
tpst_circulDiam = 0

zs_len = 3.13
zs12_blade = 6E-5
zs345_blade = 1E-4
zs_extrDiam = 0.02
zs_circulDiam = 0

mst_len = 2.38
mst_blade = 0.0042
mst_extrDiam = 0.0979
mst_circulDiam = 0

mse_len = 2.38
mse_blade = 0.01725
mse_extrDiam = 0.061
mse_circulDiam = 0

# Practical hardware settings
zs12_blade *= 2
zs345_blade *= 2

kQDA21710 = -0.01196559835
kQFA21810 = 0.01466201431*9/11
kQDA21910quad = -0.01196559835  # From my own MADX match, =~ m.gyr
kQDA21910hole = kQDA21910quad*(-.16000500)  # From MADX/sicopex
MPNH21732angle = 0.339718E-3  # From trim editor

zs_an = 8.327E-5
mst_an = 1.69520713e-3 / 3.0
mse_an = 9.74519477e-3 / 5.0

zs1_up = 0.06817
zs1_down = 0.063912049
zs2_up = 0.062799488
zs2_down = 0.058541537
zs3_up = 0.057448975
zs3_down = 0.053191025
zs4_up = 0.052078463
zs4_down = 0.047820512
zs5_up = 0.046707951
zs5_down = 0.04245

tpstPosUp = 0.03969
tpstPosBreak = 0.04026137
tpstPosDown = 0.041086756

mst1_up = 0.041509243
mst1_down = 0.04307643
mst2_up = 0.04307643
mst2_down = 0.045188215
mst3_up = 0.045732813
mst3_down = 0.0473

mse1_up = 0.056875
mse1_down = 0.056109803
mse2_up = 0.056397256
mse2_down = 0.058539328
mse3_up = 0.060040694
mse3_down = 0.065090036
mse4_up = 0.067805315
mse4_down = 0.075761927
mse5_up = 0.079691118
mse5_down = 0.090555

##################################################
#                                                #
#        Specific LSS2 elements                  #
#                                                #
##################################################

# TODO: add apertures for non-septum elements??

# ZS

zs1 = lt.Septum('ZS1', zs_len, zs_an, zs1_up, zs1_down, zs12_blade,
                zs_circulDiam, zs_extrDiam)

zs2 = lt.Septum('ZS2', zs_len, zs_an, zs2_up, zs2_down, zs12_blade,
                zs_circulDiam, zs_extrDiam)

zs3 = lt.Septum('ZS3', zs_len, zs_an, zs3_up, zs3_down, zs345_blade,
                zs_circulDiam, zs_extrDiam)

zs4 = lt.Septum('ZS4', zs_len, zs_an, zs4_up, zs4_down, zs345_blade,
                zs_circulDiam, zs_extrDiam)

zs5 = lt.Septum('ZS5', zs_len, zs_an, zs5_up, zs5_down, zs345_blade,
                zs_circulDiam, zs_extrDiam)

drift_zs12 = lt.Drift('DRIFT_ZS12', 0.78, 0)

drift_zs23 = lt.Drift('DRIFT_ZS23', 0.78, 0)

drift_zs34 = lt.Drift('DRIFT_ZS34', 0.78, 0)

drift_zs45 = lt.Drift('DRIFT_ZS45', 0.78, 0)

# ZS - MST

tce = lt.Drift('TCE', 2.7, r_tce)

qda21710 = lt.Quadrupole('QDA21710', 3.791, kQDA21710, QDArad)

mpnh21732 = lt.Kicker('MPNH21732', 2.04, MPNH21732angle, 0)

tpst1 = lt.DoubleApDrift('TPST1', 0.875, tpstPosUp, tpstPosBreak,
                         tpst1_thick, tpst_circulDiam, tpst_extrDiam)

tpst2 = lt.DoubleApDrift('TPST2', 2.251-0.875, tpstPosBreak, tpstPosDown,
                         tpst2_thick, tpst_circulDiam, tpst_extrDiam)

drift_zs_tce = lt.Drift('DRIFT_ZS_TCE', 2.491, 0)

drift_tce_qd = lt.Drift('DRIFT_TCE_QD', 2.5867, 0)

drift_qd_bump = lt.Drift('DRIFT_QD_BUMP', 0.777, 0)

drift_bump_tpst = lt.Drift('DRIFT_BUMP_TPST', 10.087, 0)

drift_tpst_mst = lt.Drift('DRIFT_TPST_MST', 0.692, 0)

# MST

mst1 = lt.Septum('MST1', mst_len, mst_an, mst1_up, mst1_down, mst_blade,
                 mst_circulDiam, mst_extrDiam)

mst2 = lt.Septum('MST2', mst_len, mst_an, mst2_up, mst2_down, mst_blade,
                 mst_circulDiam, mst_extrDiam)

mst3 = lt.Septum('MST3', mst_len, mst_an, mst3_up, mst3_down, mst_blade,
                 mst_circulDiam, mst_extrDiam)

drift_mst12 = lt.Drift('DRIFT_MST12', 0.854, 0)

drift_mst23 = lt.Drift('DRIFT_MST23', 0.854, 0)

# MST - MSE

qfa21810 = lt.Quadrupole('QFA21810', 3.791, kQFA21810, QFArad)

drift_mst_qf = lt.Drift('DRIFT_MST_QF', 3.5117, 0)

drift_qf_mse = lt.Drift('DRIFT_QF_MSE', 1.302, 0)

# MSE

mse1 = lt.Septum('MSE1', mse_len, mse_an, mse1_up, mse2_down, mse_blade,
                 mse_circulDiam, mse_extrDiam)

mse2 = lt.Septum('MSE2', mse_len, mse_an, mse2_up, mse2_down, mse_blade,
                 mse_circulDiam, mse_extrDiam)

mse3 = lt.Septum('MSE3', mse_len, mse_an, mse3_up, mse3_down, mse_blade,
                 mse_circulDiam, mse_extrDiam)

mse4 = lt.Septum('MSE4', mse_len, mse_an, mse4_up, mse4_down, mse_blade,
                 mse_circulDiam, mse_extrDiam)

mse5 = lt.Septum('MSE5', mse_len, mse_an, mse5_up, mse5_down, mse_blade,
                 mse_circulDiam, mse_extrDiam)

drift_mse12 = lt.Drift('DRIFT_MSE12', 0.854, 0)

drift_mse23 = lt.Drift('DRIFT_MSE23', 0.854, 0)

drift_mse34 = lt.Drift('DRIFT_MSE34', 0.854, 0)

drift_mse45 = lt.Drift('DRIFT_MSE45', 0.854, 0)

# MSE - TT20

qda21910 = lt.QuadHole('QDA21910', 3.791, kQDA21910quad, kQDA21910hole,
                       QDArad, QD219holeRadius, QD219holeFieldAxis,
                       QD219holeApAxis, QD219holeApAxis)

drift_mse_qd = lt.Drift('DRIFT_MSE_QD', 11.5887, 0)

##################################################
#                                                #
#   Definition of the line                       #
#                                                #
##################################################

line =  [zs1, drift_zs12, zs2, drift_zs23, zs3, drift_zs34, zs4,
        drift_zs45, zs5,
        drift_zs_tce, tce, drift_tce_qd, qda21710, drift_qd_bump,
        mpnh21732, drift_bump_tpst, tpst1, tpst2, drift_tpst_mst,
        mst1, drift_mst12, mst2, drift_mst23, mst3,
        drift_mst_qf, qfa21810, drift_qf_mse,
        mse1, drift_mse12, mse2, drift_mse23, mse3, drift_mse34, mse4,
        drift_mse45, mse5,
        drift_mse_qd, qda21910]

##################################################
#                                                #
#   Plotting helpers                             #
#                                                #
##################################################

def _rgb(r,g,b):
    return (r/255,g/255,b/255)

def colormap():
    from matplotlib.colors import LinearSegmentedColormap
    cmaplist = [(.5,.5,.5,0.05),
                _rgb(170, 68, 85), _rgb(221, 119, 136),
                _rgb(119, 17, 85), _rgb(170, 68, 136), _rgb(204, 153, 187),
                _rgb(68, 119, 170),
                _rgb(17, 119, 119), _rgb(68, 170, 170), _rgb(119, 204, 204),
                _rgb(17, 119, 68), _rgb(68, 170, 119), _rgb(136, 204, 170),
                _rgb(119, 170, 221),
                _rgb(119, 119, 17), _rgb(170, 170, 68), _rgb(221, 221, 119),
                _rgb(221, 170, 119)]
    return LinearSegmentedColormap.from_list('Default LSS2 colormap',cmaplist,18)

def colorcodes():
    return [['free aperture',['CIRCULATING'],[],[]],
            ['ZS wire front',['ZS','start'],[],['blade','coll']],
            ['ZS wire downstream',['ZS','down'],[],['blade','coll']],
            ['ZS cathode front',['ZS','start_extr'],[],[]],
            ['ZS cathode downstream',['ZS','down_extr'],[],[]],
            ['ZS other',['ZS'],['blade','coll','extr','DRIFT'],[]],
            ['TCE',['TCE'],['DRIFT'],[]],
            ['TPST blade front',['TPST','start_coll'],[],[]],
            ['TPST blade downstream',['TPST','down_coll'],[],[]],
            ['TPST other',['TPST'],['coll','DRIFT'],[]],
            ['MST blade front',['MST','start'],[],['blade','coll']],
            ['MST blade downstream',['MST','down'],[],['blade','coll']],
            ['MST other',['MST'],['blade','coll','DRIFT'],[]],
            ['Q218',['QFA21810'],[],[]],
            ['MSE blade front',['MSE','start'],[],['blade','coll']],
            ['MSE blade downstream',['MSE','down'],[],['blade','coll']],
            ['MSE other',['MSE'],['blade','coll','DRIFT'],[]],
            ['Other',[],[],[]]]

def beamcolormap():
    from matplotlib.colors import LinearSegmentedColormap
    cmaplist = [(1,0,0),(0,1,0),(0,0,1),(0,0,0)]
    return LinearSegmentedColormap.from_list('Default LSS2 beam colormap',cmaplist,4)

def beamcolorcodes():
    return [['left up',['0'],[],[]],
            ['left down',['1'],[],[]],
            ['right up',['2'],[],[]],
            ['right down',['3'],[],[]]]

def aperture(infty=0.1, s=0):
    res = []
    s0 = s
    for element in line:
        res.extend(element.aperture(infty, s))
        s += element.len
    return res
        

beam = [[0.06817,-0.00143],[0.06817,-0.00147],[0.08, -0.00173], [0.082, -0.00173]]
