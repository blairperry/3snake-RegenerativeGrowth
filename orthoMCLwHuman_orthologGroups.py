
full_1to1 = []
unique_hum = []
notUnique_hum = []

with open('groups_allCore.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        if len(line) == 4:
            pCount = 0
            rCount = 0
            gCount = 0
            hCount = 0
            for item in line:
                species = item.split('|')[0]
                if species == 'python':
                    pCount += 1
                elif species == 'rattlesnake':
                    rCount += 1
                elif species == 'gartersnake':
                    gCount += 1
                elif species == 'human':
                    hCount += 1
            if pCount == 1 and rCount == 1 and gCount == 1 and hCount == 1:
                full_1to1.append(line)
        elif len(line) > 4:
            pCount = 0
            rCount = 0
            gCount = 0
            hCount = 0
            for item in line:
                species = item.split('|')[0]
                if species == 'python':
                    pCount += 1
                elif species == 'rattlesnake':
                    rCount += 1
                elif species == 'gartersnake':
                    gCount += 1
                elif species == 'human':
                    hCount += 1
            if hCount == 1:
                unique_hum.append(line)
            elif hCount > 1:
                notUnique_hum.append(line)

with open('nonOne-to-One_notUniqueHuman_02.04.18.txt','w') as o:
    for line in notUnique_hum:
        print >> o, '\t'.join(line)
