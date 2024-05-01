import pandas as pd
import sys

seq = sys.argv[1]
charge_num = int(sys.argv[2])
output = sys.argv[3]


monomass = {'A': 71.037113805, 'R': 156.10111105, 'N': 114.04292747, 'D': 115.026943065, 'C': 103.009184505,
            'E': 129.042593135, 'Q': 128.05857754, 'G': 57.021463735, 'H': 137.058911875, 'I': 113.084064015,
            'L': 113.084064015, 'K': 128.09496305, 'M': 131.040484645, 'F': 147.068413945, 'P': 97.052763875,
            'S': 87.032028435, 'T': 101.047678505, 'W': 186.07931298, 'Y': 163.063328575, 'V': 99.068413945}


def gen_bion(seq):
    bion = []
    for i in range(1, len(seq) + 1):
        bion.append(seq[0:i])
    return bion


def gen_yion(seq):
    yion = []
    for i in range(0, len(seq)):
        yion.append(seq[i:])
    return yion


# print(gen_bion(seq))
# print(gen_yion(seq))


def mono_bion(bion, monomass):  # bion = residual + H+
    im = []
    for i in bion:
        m = 0
        for j in i:
            m += monomass[j]
        im.append(m + 1.00727646688)
    return im


def mono_yion(yion, monomass):  # yion = mono + H+
    im = []
    for i in yion:
        m = 0
        for j in i:
            m += monomass[j]
        im.append(m + 3 * 1.007825035 + 15.99491463)
    return im


def chargecolumn(im, cha):  # b,y 自带1个电荷
    chargemass = []
    for i in im:
        chargemass.append(float((i + (cha - 1) * 1.00727646688) / cha))
    return chargemass


def modify(im, mdn, mdv, cha):
    modmass = []
    for i in im:
        modmass.append(float((i + mdn * mdv + (cha - 1) * 1.00727646688) / cha))
    return modmass


def modify_nocharge(im, mdn, mdv):
    modmass_nocharge = []
    for i in im:
        modmass_nocharge.append(float((i + mdn * mdv)))
    return modmass_nocharge


# seq = 'SCQCGSKNGTC'  # input
# charge_num = 10  # input

h2o_num = seq.count('S') + seq.count('T')

data_bion = {}
bion_seq = gen_bion(seq)
im_b = mono_bion(bion_seq, monomass)
data_bion['bion_seq'] = bion_seq
data_bion['isomass'] = im_b

data_yion = {}
yion_seq = gen_yion(seq)
im_y = mono_yion(yion_seq, monomass)
data_yion['yion_seq'] = yion_seq
data_yion['isomass'] = im_y

for i in range(1, charge_num + 1):
    data_bion['+%d' % i] = chargecolumn(im_b, i)
    data_yion['+%d' % i] = chargecolumn(im_y,i)

del_H2O = -(2 * 1.007825035 + 15.99491463)
del_car = -(12 + 15.99491463 * 2 + 1.007825035 * 2)
add_2H = 2 * 1.007825035
add_CH3 = 12 + 2 * 1.007825035
for num_h2o in range(0, h2o_num + 1):
    for num_2H in range(0, num_h2o + 1):
        for num_CH3 in range(0, 5):
            for num_de_car in range(0, 2):
                for charge in range(1, charge_num + 1):
                    data_bion['-%dH2O_-%dCO2H2+%dH_+%dCH2+%d' % (num_h2o,num_de_car, 2 * num_2H, num_CH3, charge)] = \
                        modify(modify_nocharge(modify_nocharge(modify_nocharge(im_b, num_h2o, del_H2O),
                                                               num_de_car, del_car),num_2H, add_2H), num_CH3, add_CH3, charge)  # -H2O, -CO2, +2H, +CH3
                    data_yion['-%dH2O_-%dCO2H2+%dH_+%dCH2+%d' % (num_h2o,num_de_car, 2 * num_2H, num_CH3, charge)] = \
                        modify(modify_nocharge(modify_nocharge(modify_nocharge(im_y, num_h2o, del_H2O), num_de_car, del_car),
                                               num_2H, add_2H), num_CH3, add_CH3, charge)  # -H2O, -CO2, +2H, +CH3

df_b = pd.DataFrame(data_bion)
df_y = pd.DataFrame(data_yion)
df_b.to_csv('%s/bion.csv'%output, sep=',', index=False)  # output input
df_y.to_csv('%s/yion.csv'%output, sep=',', index=False)  # output input
