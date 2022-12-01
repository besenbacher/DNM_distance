import sys
import random
from collections import Counter
#from DNM_distance.readers import eprint, BedReader
from readers import eprint, BedReader

def pos_to_int(chrom, pos, cum_length):
    return cum_length[chrom] + int(pos)

# METHOD
#Find ud af hvor mange n_cpg n_gc og n_at mut der er i hvert individ.
#For hvert individ traek n_cpg,n_gc, og n_at tilfaeldige floats i intervallet 0..1
#sorter de tilfaeldige tal L_cpg : (r, individ, sample_i), L_gc og L_at
#Loeb callability fil igennem og udregn hvilken chrom og pos hvert r svarer til


def print_chrom_pos(chrom, pos, next_tuple, outfile):
    if not outfile is None:
        _, pn, i = next_tuple
        print(f'{chrom} {pos} {pn} {i}', file = outfile)


def do_analysis(mutfile, bedfile, tb, n_random, obs_outfile = None, random_outfile = None, random_chrompos_outfile = None):
    all_pns = set()
    n_cpg = Counter()
    n_cg = Counter()
    n_at = Counter()
    intposses = []
    used_chroms = set()
    f = open(mutfile)
    for line in f:
        chrom, pos, ref, alt, pn = line.split()
        pos = int(pos)
        used_chroms.add(chrom)
    f.close()

    used_chroms = list(used_chroms)
    used_chroms.sort()
    chrom_lengths = tb.chroms()
    cum_length = {}
    cum_sum = 0
    for chrom in used_chroms:
        cum_length[chrom] = cum_sum
        cum_sum += chrom_lengths[chrom]

    f = open(mutfile)
    for line in f:
        chrom, pos, ref, alt, pn = line.split()
        if len(ref) != 1 or len(alt) != 1:
            continue
        pos = int(pos)
        context = tb.sequence(chrom, pos-2, pos+1).upper()
        if len(context) < 3:
            continue
        intposses.append((pos_to_int(chrom,pos, cum_length), pn))
        all_pns.add(pn)
        if context[:2] == 'CG' or context[1:] == 'CG':
            n_cpg[pn] += 1
        elif context[1] == 'C' or context[1] == 'G':
            n_cg[pn] += 1
        else:
            n_at[pn] += 1
    f.close()
    intposses.sort()

    if not obs_outfile is None:
        for pos, pn in intposses:
            print(pos, pn, file = obs_outfile)

    obs_any, obs_same, obs_notsame = get_dist_lists(intposses, all_pns, cum_sum)
    eprint("observed", len(obs_any), len(obs_same), len(obs_notsame))

    csum_cpg = 0.0
    csum_cg = 0.0
    csum_at = 0.0

    dreader = BedReader(bedfile)
    for (chrom, pos, callability) in dreader:
        pos = pos + 1 # bed returnerer 0 baserede positioner
        if chrom not in used_chroms:
            continue
        #context = genome[chrom][pos-2:pos+1].upper()
        try:
            context = tb.sequence(chrom, pos-2, pos+1).upper()
        except:
            continue
        if len(context) < 3 or context[1]=="N":
            continue
        if context[:2] == 'CG' or context[1:] == 'CG':
            csum_cpg += callability
        elif context[1] == 'C' or context[1] == 'G':
            csum_cg += callability
        else:
            csum_at += callability
    eprint("csum_cpg:", csum_cpg)
    eprint("csum_cg:", csum_cg)
    eprint("csum_at:", csum_at)
    eprint("n_cpg:", sum(n_cpg[x] for x in n_cpg))
    eprint("n_cg:", sum(n_cg[x] for x in n_cg))
    eprint("n_at:", sum(n_at[x] for x in n_at))

    L_cpg = []
    L_cg = []
    L_at = []

    for pn in all_pns:
        for i in range(n_random):
            for x in range(n_cpg[pn]):
                L_cpg.append((random.random()*csum_cpg, pn, i))
            for x in range(n_cg[pn]):
                L_cg.append((random.random()*csum_cg, pn, i))
            for x in range(n_at[pn]):
                L_at.append((random.random()*csum_at, pn, i))

    L_cpg.sort(reverse=True)
    L_cg.sort(reverse=True)
    L_at.sort(reverse=True)
    next_cpg = L_cpg.pop()
    next_cg = L_cg.pop()
    next_at = L_at.pop()

    random_intposses = [[] for x in range(n_random)]

    dreader = BedReader(bedfile)

    cur_cpg = 0.0
    cur_cg = 0.0
    cur_at = 0.0
    last_chrom = "NA"
    for (chrom, pos, callability) in dreader:
        pos = pos + 1 # bed returnerer 0 baserede positioner
        if chrom not in used_chroms:
            continue
        if chrom != last_chrom:
            eprint("at chrom", chrom)
            last_chrom = chrom
        #context = genome[chrom][pos-2:pos+1].upper()
        try:
            context = tb.sequence(chrom, pos-2, pos+1).upper()
        except:
            continue
        if len(context) <3 or context[1] == "N":
            continue
        if context[:2] == 'CG' or context[1:] == 'CG':
            cur_cpg += callability
            while cur_cpg > next_cpg[0] and next_cpg[2] != -1:
                random_intposses[next_cpg[2]].append((pos_to_int(chrom,pos, cum_length), next_cpg[1]))
                print_chrom_pos(chrom, pos, next_cpg, random_chrompos_outfile)
                if len(L_cpg) > 0:
                    next_cpg = L_cpg.pop()
                else:
                    next_cpg = (csum_cpg + 10000, -1, -1)
        elif context[1] == 'C' or context[1] == 'G':
            cur_cg += callability
            while cur_cg > next_cg[0] and next_cg[2] != -1:
                random_intposses[next_cg[2]].append((pos_to_int(chrom,pos, cum_length), next_cg[1]))
                print_chrom_pos(chrom, pos, next_cg, random_chrompos_outfile)
                if len(L_cg) > 0:
                    next_cg = L_cg.pop()
                else:
                    next_cg = (csum_cg + 10000, -1, -1)
        else:
            cur_at += callability
            while cur_at > next_at[0] and next_at[2] != -1:
                random_intposses[next_at[2]].append((pos_to_int(chrom,pos, cum_length), next_at[1]))
                print_chrom_pos(chrom, pos, next_at, random_chrompos_outfile)
                if len(L_at) > 0:
                    next_at = L_at.pop()
                else:
                    next_at = (csum_at +10000, -1, -1)

    if not random_outfile is None:
        for i in range(n_random):
            for pos, pn in random_intposses[i]:
                print(pos, pn, i, file = random_outfile)

    random_any_list = []
    random_same_list = []
    random_notsame_list = []
    for i in range(n_random):
        r_any, r_same, r_notsame = get_dist_lists(random_intposses[i], all_pns, cum_sum)
        eprint("random", i, len(r_any), len(r_same), len(r_notsame))
        random_any_list.append(r_any)
        random_same_list.append(r_same)
        random_notsame_list.append(r_notsame)

    print("obs mean.random median.random type")
    for x1, x2, x3 in zip(obs_any,
                          list_mean(random_any_list),
                          list_median(random_any_list)):
        print(x1, x2, x3, "any")
    for x1, x2, x3 in zip(obs_same,
                          list_mean(random_same_list),
                          list_median(random_same_list)):
        print(x1, x2, x3, "same")
    for x1, x2, x3 in zip(obs_notsame,
                          list_mean(random_notsame_list),
                          list_median(random_notsame_list)):
        print(x1, x2, x3, "notsame")


def dist_to_closest(intposses, index, possible_pn):
    bestres = 9e10
    for direction in [1,-1]:
        i = index
        while (i + direction) >= 0 and (i + direction) < len(intposses):
            i = i + direction
            intpos, pn = intposses[i]
            if pn in possible_pn:
                bestres = min(bestres, abs(intposses[i] - intposses[index]))
                break
            else: continue
    return bestres

def dist_to_right(intposses, index, possible_pn, cum_sum):
    i = index
    start_intpos, start_pn = intposses[index]
    while (i + 1) <len(intposses):
        i += 1
        intpos, pn = intposses[i]
        if pn in possible_pn:
            return intpos - start_intpos
    i = 0
    while i < len(intposses):
        intpos, pn = intposses[i]
        if pn in possible_pn:
            return intpos + (cum_sum - start_intpos)
        i += 1
    return 'NA'

def get_dist_lists(intposses, all_pns, cum_sum):
    dist_any = []
    dist_same = []
    dist_notsame = []
    for index in range(len(intposses)):
        intpos, pn = intposses[index]
        dr_any = dist_to_right(intposses, index, all_pns, cum_sum)
        dr_same = dist_to_right(intposses, index, set([pn]), cum_sum)
        dr_notsame = dist_to_right(intposses, index, all_pns - set([pn]), cum_sum)
        assert dr_any != 'NA'
        dist_any.append(dr_any)
        if dr_same != 'NA':
            dist_same.append(dr_same)
        if dr_notsame != 'NA':
            dist_notsame.append(dr_notsame)
    dist_any.sort()
    dist_same.sort()
    dist_notsame.sort()
    return dist_any, dist_same, dist_notsame

def list_mean(L):
    res = []
    for i_index in range(len(L[0])):
        i_sum = 0
        for i_trial in range(len(L)):
            i_sum += L[i_trial][i_index]
        res.append(float(i_sum)/len(L))
    return res

def list_median(L):
    res = []
    for i_index in range(len(L[0])):
        i_sum = 0
        tmp_L = []
        for i_trial in range(len(L)):
            tmp_L.append(L[i_trial][i_index])
        tmp_L.sort()
        res.append(tmp_L[len(tmp_L)//2])
    return res
