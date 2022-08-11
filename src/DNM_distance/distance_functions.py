import sys
import random
from collections import Counter
import py2bit

autosomes = [f"chr{i}" for i in range(1,23)]

mutfile = sys.argv[1]
#coveredfile = sys.argv[2]
#n_random = int(sys.argv[3])
#hgver = sys.argv[4]
#out_dir = sys.argv[4]





tb = py2bit.open("/Users/au149884/Data/2bit/hg38.2bit")
#genome = TwoBitFile(open(datapaths.TwoBit(hgver)))
#chrom_lengths = get_chrom_lengths(hgver)
chrom_lengths = tb.chroms()
autosomes.sort()
cum_length = {}
cum_sum = 0
for chrom in autosomes:
    cum_length[chrom] = cum_sum
    cum_sum += chrom_lengths[chrom]

def pos_to_int(chrom,pos):
    return cum_length[chrom] + int(pos)


all_pns = set()
n_pn = Counter()
intposses = []
f = open(mutfile)
for line in f:
    chrom, pos, ref, alt, pn = line.split()
    if chrom not in autosomes:
        continue
    intposses.append((pos_to_int(chrom,pos), pn))
    all_pns.add(pn)
    n_pn[pn] += 1
f.close()
intposses.sort()


def dist_to_closest(intposses, index, possible_pn):
    bestres = 9e10
    for direction in [1,-1]:
        i = index
        while (i + direction)>=0 and (i + direction)<len(intposses):
            i = i + direction
            intpos, pn = intposses[i]
            if pn in possible_pn:
                bestres = min(bestres, abs(intposses[i] - intposses[index]))
                break
            else: continue
    return bestres

def dist_to_right(intposses, index, possible_pn):
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

def get_dist_lists(intposses):
    dist_any = []
    dist_same = []
    dist_notsame = []
    for index in range(len(intposses)):
        intpos, pn = intposses[index]
        dr_any = dist_to_right(intposses, index, all_pns)
        dr_same = dist_to_right(intposses, index, set([pn]))
        dr_notsame = dist_to_right(intposses, index, all_pns - set([pn]))
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
    for i_index in xrange(len(L[0])):
        i_sum = 0
        for i_trial in xrange(len(L)):
            i_sum += L[i_trial][i_index]
        res.append(float(i_sum)/len(L))
    return res

def list_median(L):
    res = []
    for i_index in xrange(len(L[0])):
        i_sum = 0
        tmp_L = []
        for i_trial in xrange(len(L)):
            tmp_L.append(L[i_trial][i_index])
        tmp_L.sort()
        res.append(tmp_L[len(tmp_L)/2])
    return res


obs_any, obs_same, obs_notsame = get_dist_lists(intposses)

print("observed", len(obs_any), len(obs_same), len(obs_notsame), file=sys.stderr)

#print >> sys.stderr, "observed", len(obs_any), len(obs_same), len(obs_notsame)
#OLD METHOD:
#Find ud af hvor mange sites er coveret: n_cov
#Find antal mut i hvert individ: n_pn
#For hvert individ traek n_pn tilfaeldige tal i intervallet 1..n_cov
#Udregn hvilke intpos de n_pn tal svarer til

#NEW METHOD
#Find ud af hvor mange n_cpg n_gc og n_at mut der er i hvert individ.
#For hvert individ traek n_cpg,n_gc, og n_at tilfaeldige floats i intervallet 0..1
#sorter de tilfaeldige tal L_cpg : (r, individ, sample_i), L_gc og L_at
#Loeb callability fil igennem og udregn hvilken chrom og pos hvert r svarer til

# all_pns = set()
# n_cpg = Counter()
# n_cg = Counter()
# n_at = Counter()
# intposses = []
# f = open(mutfile)
# for line in f:
#     chrom, pos, ref, alt, pn = line.split()
#     if chrom not in autosomes:
#         continue
#     context = genome[chrom][pos-2:pos+1].upper()
#     intposses.append((pos_to_int(chrom,pos), pn))
#     all_pns.add(pn)
#     if context[:2] == 'CG' or context[1:] == 'CG':
#         n_cpg[pn] += 1
        
#     elif context[1] == 'C' or context[1] == 'G':
#         n_cg[pn] += 1
#     else:
#         n_at[pn] += 1
# f.close()
# intposses.sort()


# csum_cpg = 0.0
# csum_cp = 0.0
# csum_at = 0.0

# wigr = WigReader(callabilityfile)
# for (chrom,pos,callability) in wigr:
#     if chrom not in autosomes:
#         continue
#     context = genome[chrom][pos-2:pos+1].upper()
#     if context[:2] == 'CG' or context[1:] == 'CG':
#         csum_cpg += callability
#     elif context[1] == 'C' or context[1] == 'G':
#         csum_cg += callability
#     else:
#         csum_at += callability
# wigr.close()

# print >> sys.stderr, "n_cov:", n_cov

# L_cpg = []
# L_gc = []
# L_at = []

# for pn in all_pns:
#     for i in n_random:
#         for x in n_cpg[pn]:
#             L_cpg.append((random.random()*csum_cpg, pn, i))
#         for x in n_cg[pn]:
#             L_cg.append((random.random()*csum_cg, pn, i))
#         for x in n_at[pn]:
#             L_at.append((random.random()*csum_at, pn, i))

# L_cpg.sort()
# L_cg.sort()
# L_at.sort()
# next_cpg = L_cpg.pop()
# next_cg = L_cg.pop()
# next_at = L_at.pop()

# intposes = [[] for x in range(n_random)]
# #pos_cg = [[] for x in range(n_random)]
# #pos_at = [[] for x in range(n_random)]

# wigr = WigReader(callabilityfile)
# for (chrom,pos,callability) in wigr:
#     if chrom not in autosomes:
#         continue
#     context = genome[chrom][pos-2:pos+1].upper()
#     if context[:2] == 'CG' or context[1:] == 'CG':
#         csum_cpg += callability
#         if csum_cpg > next_cpg[0]:
#             intposes[next_cpg[2]].append((pos_to_int(chrom,pos), next_cpg[1]))
#             if len(L_cpg) > 0:
#                 next_cpg = L_cpg.pop()
#             else:
#                 next_cg = (csum_cpg +1, -1, -1)                
#     elif context[1] == 'C' or context[1] == 'G':
#         csum_cg += callability
#         if csum_cg > next_cg[0]:
#             intposes[next_cg[2]].append((pos_to_int(chrom,pos), next_cg[1]))
#             if len(L_cg) > 0:
#                 next_cg = L_cg.pop()
#             else:
#                 next_cg = (csum_cg +1, -1, -1)
#     else:
#         csum_at += callability
#         if csum_at > next_at[0]:
#             intposes[next_at[2]].append((pos_to_int(chrom,pos), next_at[1]))
#             if len(L_at) > 0:
#                 next_at = L_at.pop()
#             else:
#                 next_at = (csum_at +1, -1, -1)
# wigr.close()

        
# random_any_list = []
# random_same_list = []
# random_notsame_list = []
# for i in range(n_random):
#     r_any, r_same, r_notsame = get_dist_lists(intposses[i])
#     print >> sys.stderr, "random", i, len(r_any), len(r_same), len(r_notsame)
#     random_any_list.append(r_any)
#     random_same_list.append(r_same)
#     random_notsame_list.append(r_notsame)

# #f = open(out_dir + '/dist_between_mutations.dat', 'w')
# f = sys.stdout
# print >> f, "obs mean.random median.random type"
# for x1, x2, x3 in zip(obs_any,
#                       list_mean(random_any_list),
#                       list_median(random_any_list)):
#     print >> f, x1, x2, x3, "any"
# for x1, x2, x3 in zip(obs_same,
#                       list_mean(random_same_list),
#                       list_median(random_same_list)):
#     print >> f, x1, x2, x3, "same"
# for x1, x2, x3 in zip(obs_notsame,
#                       list_mean(random_notsame_list),
#                       list_median(random_notsame_list)):
#     print >> f, x1, x2, x3, "notsame"
# #f.close()
