

ortho_check = {}

with open('nonOne-to-One_notUniqueHuman_02.04.18.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        human_IDs = []
        for item in line:
            if 'human' in item:
                human_IDs.append(item)
        for item in line:
            if 'human' not in item:
                ortho_check[item] = human_IDs


best_hits = {}

with open('allSnake_human.BestBlast_results.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        best_hits[line[0]] = line[1]

with open('allSnake_human_BestBlast_IDconversions.txt','w') as o:
    for ID in ortho_check:
        if ID in best_hits:
            best = best_hits[ID]
            print >> o, '\t'.join([ID, best])
