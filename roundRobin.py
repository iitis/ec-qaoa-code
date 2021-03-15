# soure: https://code.activestate.com/recipes/65200-round-robin-pairings-generator/
# adapted for python 3

def roundRobin(units, sets=None):
    """ Generates a schedule of "fair" pairings from a list of units """
    if len(units) % 2:
        units.append(None)
    count    = len(units)
    sets     = sets or (count - 1)
    half     = count / 2
    schedule = []
    for turn in range(sets):
        pairings = []
        for i in range(int(half)):
            pairings.append(units[i])
            pairings.append(units[count-i-1])
        units.insert(1, units.pop())
        schedule.append(pairings)
    return schedule

""" test code """
if __name__ == '__main__':
    a=list(range(1,14))
    b=[]
    for pairings in roundRobin(a):
        b.append(list(zip(pairings[0::2],pairings[1::2])))
    print(b[0])

#     pairings = roundRobin(a)
#     print(pairings)

#     for pairings in roundRobin(a):
#         print(pairings)

#     for pairings in roundRobin(a):
#         print(list(zip(pairings[0::2],pairings[1::2])))
