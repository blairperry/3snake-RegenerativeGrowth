

best_hits = {}

with open('/Users/perryb/Downloads/rattlesnake_human.blast_results.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        human = ''
        snake = ''
        if (line[0].split('|'))[0] == 'human':
            human = line[0]
            snake = line[1]
        else:
            human = line[1]
            snake = line[0]
        if snake not in best_hits:
            best_hits[snake] = [human,float(line[-2]),float(line[-1])]
        else:
            prev_eval = best_hits[snake][1]
            new_eval = float(line[-2])
            if new_eval < prev_eval:
                best_hits[snake] = [human,float(line[-2]),float(line[-1])]
            elif new_eval == prev_eval:
                prev_bit = best_hits[snake][2]
                new_bit = float(line[-1])
                if new_bit > prev_bit:
                    best_hits[snake] = [human,float(line[-2]),float(line[-1])]

with open('rattlesnake_human.BestBlast_results.txt','w') as o:
    for key in best_hits:
        print >> o, '\t'.join([key,best_hits[key][0]])





best_hits = {}

with open('/Users/perryb/Downloads/python_human.blast_results.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        human = ''
        snake = ''
        if (line[0].split('|'))[0] == 'human':
            human = line[0]
            snake = line[1]
        else:
            human = line[1]
            snake = line[0]
        if snake not in best_hits:
            best_hits[snake] = [human,float(line[-2]),float(line[-1])]
        else:
            prev_eval = best_hits[snake][1]
            new_eval = float(line[-2])
            if new_eval < prev_eval:
                best_hits[snake] = [human,float(line[-2]),float(line[-1])]
            elif new_eval == prev_eval:
                prev_bit = best_hits[snake][2]
                new_bit = float(line[-1])
                if new_bit > prev_bit:
                    best_hits[snake] = [human,float(line[-2]),float(line[-1])]

with open('python_human.BestBlast_results.txt','w') as o:
    for key in best_hits:
        print >> o, '\t'.join([key,best_hits[key][0]])





best_hits = {}

with open('/Users/perryb/Downloads/gartersnake_human.blast_results.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        human = ''
        snake = ''
        if (line[0].split('|'))[0] == 'human':
            human = line[0]
            snake = line[1]
        else:
            human = line[1]
            snake = line[0]
        if snake not in best_hits:
            best_hits[snake] = [human,float(line[-2]),float(line[-1])]
        else:
            prev_eval = best_hits[snake][1]
            new_eval = float(line[-2])
            if new_eval < prev_eval:
                best_hits[snake] = [human,float(line[-2]),float(line[-1])]
            elif new_eval == prev_eval:
                prev_bit = best_hits[snake][2]
                new_bit = float(line[-1])
                if new_bit > prev_bit:
                    best_hits[snake] = [human,float(line[-2]),float(line[-1])]

with open('gartersnake_human.BestBlast_results.txt','w') as o:
    for key in best_hits:
        print >> o, '\t'.join([key,best_hits[key][0]])
