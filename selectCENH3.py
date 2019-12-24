from Bio import SeqIO

db = 'D:\Programming\PycharmProjects\CRISPR\selectedORF_CENH3.faa'
blast_res = 'D:\Programming\PycharmProjects\CRISPR\CENH3_vs_selectedORF_CENH3'

def buildHitPerQuery(blast_res):
    d = {'LC075743.1_HaCENH3_Helianthus_annuus':{}, 'histone':{}}
    for line in blast_res:
        sp = line.split("\t")
        score = sp[-1]
        identity = sp[2]
        hit = sp[1]
        query = sp[0]
        if hit not in d[query]:
            d[query][hit] = score
        elif score > d[query][hit]:
            d[query][hit] = score
    return d

with open(blast_res) as inFile, open('D:\Programming\PycharmProjects\CRISPR\CENH3_selected.faa', 'w') as outCENH3, open('D:\Programming\PycharmProjects\CRISPR\H3_selected.faa', 'w') as outH3:
    best_score = buildHitPerQuery(inFile)
    ind = SeqIO.index(db, 'fasta')
    h3_cnt = 0
    cenh3_cnt = 0
    for hits in best_score['LC075743.1_HaCENH3_Helianthus_annuus']:
        if hits not in best_score['histone']:
            SeqIO.write(ind[hits], outCENH3, 'fasta')
            cenh3_cnt += 1
        elif best_score['histone'][hits] < best_score['LC075743.1_HaCENH3_Helianthus_annuus'][hits]:
            SeqIO.write(ind[hits], outCENH3, 'fasta')
            cenh3_cnt += 1
        else:
            s = ind[hits]
            s.id = "H3_{}".format(h3_cnt)
            SeqIO.write(s, outH3, 'fasta')
            h3_cnt += 1

    print("Number of CENH3 selected is {}".format(cenh3_cnt))
    print("Number of H3 selected is {}".format(h3_cnt))