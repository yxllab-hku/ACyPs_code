import pandas as pd
import sys

seq = sys.argv[1]
charge_num = int(sys.argv[2])
output = sys.argv[3]

monomass = {'A': 71.037113805, 'R': 156.10111105, 'N': 114.04292747, 'D': 115.026943065, 'C': 103.009184505,
            'E': 129.042593135, 'Q': 128.05857754, 'G': 57.021463735, 'H': 137.058911875, 'I': 113.084064015,
            'L': 113.084064015, 'K': 128.09496305, 'M': 131.040484645, 'F': 147.068413945, 'P': 97.052763875,
            'S': 87.032028435, 'T': 101.047678505, 'W': 186.07931298, 'Y': 163.063328575, 'V': 99.068413945}


def gen_multiseq(seq):
    seqlist = []
    for i in range(len(seq)):
        for j in range(i, len(seq)):
            seqlist.append(seq[i:j + 1])
    return seqlist


def isomass(multiseq, monomass):
    im = []
    for i in multiseq:
        m = 0
        for j in i:
            m += monomass[j]
        im.append(m + 2 * 1.007825035 + 15.99491463)
    return im


def chargecolumn(im, cha):
    chargemass = []
    for i in im:
        chargemass.append(float((i + cha * 1.00727646688) / cha))
    return chargemass


def modify(im, mdn, mdv, cha):
    modmass = []
    for i in im:
        modmass.append(float((i + mdn * mdv + cha * 1.00727646688) / cha))
    return modmass


def modify_nocharge(im, mdn, mdv):
    modmass_nocharge = []
    for i in im:
        modmass_nocharge.append(float((i + mdn * mdv)))
    return modmass_nocharge


# seq = 'SCQCGSKNGTC'  # input
# charge_num = 10  # input
h2o_num = seq.count('S') + seq.count('T')
data = {}
mulseq = gen_multiseq(seq)
im = isomass(mulseq, monomass)
data['seq'] = mulseq
data['isomass'] = im

for i in range(1, charge_num + 1):
    data['+%d' % i] = chargecolumn(im, i)

del_H2O = -(2 * 1.007825035 + 15.99491463)
del_car = -(12 + 15.99491463 * 2 + 1.007825035 * 2)
add_2H = 2 * 1.007825035
add_CH3 = 12 + 2 * 1.007825035
for num_h2o in range(1, h2o_num + 1):
    for charge in range(1, charge_num + 1):
        for num_2H in range(0, num_h2o + 1):
            for num_CH3 in range(0, 5):
                data['-%dH2O_-CO2H2_+%dH_+%dCH2_+%d' % (num_h2o, 2 * num_2H, num_CH3, charge)] = \
                    modify(modify_nocharge(modify_nocharge(modify_nocharge(im, num_h2o, del_H2O), 1, del_car), num_2H,
                                           add_2H), num_CH3, add_CH3, charge)  # -H2O, -CO2, +2H, +CH3
#
# for num_h2o in range(1, h2o_num + 1):
#     for charge in range(1, charge_num + 1):
#         for num_CH3 in range(0, 5):
#             data['-%dH2O_-CO2H2+%dCH2+%d' % (num_h2o, num_CH3, charge)] = \
#                 modify(modify_nocharge(modify_nocharge(im, num_h2o, del_H2O), 1, del_car), num_CH3, add_CH3,
#                        charge)  # -H2O, -CO2, +CH3

df = pd.DataFrame(data)
df.to_csv(output, sep='\t', index=False)  # output input
