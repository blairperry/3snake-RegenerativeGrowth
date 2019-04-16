

with open('Full1to1andUniqueHuman_OrthologConversions.txt','w') as o:
    with open('groups_allCore.txt') as a:
        for line in a.readlines():
            line = line.rstrip().split('\t')
            if len(line) == 4:
                p = ''
                r = ''
                g = ''
                h = ''
                for item in line:
                    if 'python' in item:
                        p = item
                    elif 'rattlesnake' in item:
                        r = item
                    elif 'gartersnake' in item:
                        g = item
                    elif 'human' in item:
                        h = item
                print >> o, p + '\t' + h
                print >> o, r + '\t' + h
                print >> o, g + '\t' + h

            elif len(line) > 4:
                hCount = 0
                for item in line:
                    if 'human' in item:
                        hCount += 1
                if hCount == 1:
                    p = ''
                    r = ''
                    g = ''
                    h = ''
                    for item in line:
                        if 'python' in item:
                            p = item
                        elif 'rattlesnake' in item:
                            r = item
                        elif 'gartersnake' in item:
                            g = item
                        elif 'human' in item:
                            h = item
                    print >> o, p + '\t' + h
                    print >> o, r + '\t' + h
                    print >> o, g + '\t' + h
