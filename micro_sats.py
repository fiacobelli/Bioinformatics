import sys
import itertools
from pprint import pprint
import re,math


def dna_combo(n):
    nt = ['A', 'C', 'G', 'T']
    temp = itertools.product(nt, repeat=n)
    return ["".join(x) for x in temp]
#def flank_region(length,start_m, end_m):


def find_reps(seq, dna, start_idx):
    start = dna.find(seq, start_idx)
    if start != -1:
        end = start + len(seq)
        reps = 1
        while seq[0] == dna[end]:
            end += len(seq)
            reps+=1
        return start, end,reps
    return -1, -1,-1


def find_all_ms(seq, dna, num_reps,lflank=40):
    s = 0
    ms_dict = {}
    while s > -1:
        #print "Finding reps from",s
        s,e,reps = find_reps(seq, dna, s)
        #print s,e,reps
        if s > -1 and num_reps <= reps:
            ms_dict[seq] = ms_dict.get(seq, [])
            sa=s-40 if s>lflank else 0
            ep=e+40 if e<len(dna)-lflank else len(dna)
            cga=cg_content(dna[sa:s])
            cgp=cg_content(dna[e:ep])
            ms_dict[seq].append((s, e,reps,cga,cgp))
        s = e
    return ms_dict

def find_all_ms_re(seq,dna,num_reps,lflank=40,target_cg=0.5,m_length=16,max_length=24):
    pat = re.compile(r'('+seq+'){'+str(num_reps)+',}')
    ms_dict={}
    for m in pat.finditer(dna):
        s,e,gr= (m.start(),m.end(),m.group())
        reps = (e-s)/len(seq)
        ms_dict[seq] = ms_dict.get(seq, [])
        sa=s-40 if s>lflank else 0
        ep=e+40 if e<len(dna)-lflank else len(dna)
        cgfo=cg_content(dna[sa:s])
        cgre=cg_content(dna[e:ep])
        opt_fcg,opt_f=find_best_primer(dna[sa:s],m_length,max_length,target_cg)
        opt_rcg,opt_r=find_best_primer(dna[e:ep],m_length,max_length,target_cg)
        ms_dict[seq].append((s, e,reps,cgre,cgfo,opt_f,opt_fcg,opt_r,opt_rcg))
    return ms_dict


def find_combinations(n, dna, num_reps):
    result = {}
    com = dna_combo(n)
    for n_nuc in com:
        print "Looking for ",n_nuc
        result.update(find_all_ms_re(n_nuc, dna, num_reps))
    return result


def cg_content(region_string):
    cg = region_string.count('C') + region_string.count('G')
    #print "CG Content of",region_string,"is",cg,float(cg)/len(region_string)
    return float(cg)/len(region_string)
#where to report this value?


def file_save(l_dict,ofile):
    f = open(ofile, 'w')
    f.write("microsat\tstart\tend\treps\tcg_forward\tcg_reverse\topt_forward_seq\tcg_content\topt_reverse_seq\tcg_content")
    for key,val in l_dict.items():
        for s, e,reps,cgre,cgfo,opt_f,opt_fcg,opt_r,opt_rcg in val:
            f.write("\n%s\t%d\t%d\t%d\t%f\t%f\t%s\t%f\t%s\t%f"%(key,s, e,reps,cgre,cgfo,opt_f,opt_fcg,opt_r,opt_rcg))
    f.close()


def ngrams(txt, n):
    output = []
    for i in range(len(txt)-n+1):
        output.append(txt[i:i+n])
    return output


def dna_windows(seq,min_length,max_length):
    " Returns all the sequences of length l (min_length < l <max_length) in seq"
    seqs = []
    for i in range(min_length,max_length+1):
        seqs.extend(ngrams(seq,i))
    return seqs

def optimum_cgc(seqs,target_cg):
    "Returns the sequence in seqs that has the closest cg content to target_cg"
    seql = [((cg_content(s)-target_cg)**2,s) for s in seqs]
    cgc,seq=min(seql)
    return (math.sqrt(cgc)+target_cg,seq)


def find_best_primer(seq,min_length,max_length,target_cg):
    return optimum_cgc(dna_windows(seq,min_length,max_length),target_cg)


def open_deNovo(fname):
    with open(fname) as f:
        lines=[l.strip() for l in f if ">" not in l]
    return "~".join(lines)


def main(dna, n, reps,ofile):
    d = find_combinations(n, dna, reps)
    #pprint(d)
    file_save(d,ofile)

if __name__== "__main__":
    trash,n,reps,dnafile,ofile = sys.argv
    n=int(n)
    reps = int(reps)
    if len(sys.argv) > 3:
        dna = open_deNovo(dnafile)
        print "Looking for repetions of "+str(n)+" base pairs"
        if n!=0:
           main(dna, n, reps,ofile)
        else:
           for i in [2,3,4]:
               main(dna,i,reps,ofile+".%d"%i)
    print "Done"
                

