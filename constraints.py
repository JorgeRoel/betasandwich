def hbond_constraints(strands, orientation):
    '''Generates ALL possible hbond constraints via harmonic potential'''

    constraints = []

    if orientation[0] == "+":
        ox_resi = strands[0][1%1::2]
    elif orientation[0] == "-":
        ox_resi = list(reversed(strands[0]))
        if (len(ox_resi) % 2) != 0:
            ox_resi = ox_resi[1%2::2]
        else:
            ox_resi = ox_resi[1%2::2]

    if orientation[1] == "+":
        ni_resi = strands[1][1%2::2]
    if orientation[1] == "-":
        if (len(strands[1]) % 2) != 0:
            ni_resi = list(reversed(strands[1][1%1::2]))
        else:
            ni_resi = list(reversed(strands[1][1%2::2]))

    for i,ox in enumerate(ox_resi):
        try:
            cnt = "AtomPair O {} N {} HARMONIC 2.5 0.1".format(ox, ni_resi[i])
            constraints.append(cnt)
        except:
            pass

    return constraints
