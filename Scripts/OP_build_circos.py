import operator
import numpy as np
import pandas as pd


class Wedge:
    def __init__(self,name,ct,age):
        self.name = name
        self.ct = ct
        self.age = age
        self.genes = []
        self.muts = 0
        self.axis = []
        self.pos = []
        
    def add_genes(self,gene):
        self.genes.append(gene)

    def set_som_muts(self,mutcnt):
        self.muts = mutcnt

    def add_location(self,axis,pos):
        self.axis.append(axis)
        self.pos.append(pos)


def file2list_col1(cfname,h):
    out = []
    with open(cfname, "r") as f:
        if h :
            cfhead = f.readline()
        for line in f:
            out.append(line[:-1].split("\t")[0])
    return out

def file2dict(afname,h):
    out = {}
    with open(afname,"r") as f:
        if h:
            header = f.readline() 
        for line in f:
            k,v = line[:-1].split("\t")
            out[k] = v
    f.close()
    return out

def link(CIRC):
    out = []
    for s in CIRC:
        w = CIRC[s]
        if len(w.genes) > 1: 
            out.append(w.name)
    return out

def orderWedges(CIRC,GENE,gene2path,path2rgb):
    #Make an ordered Karyotype file 
    linked_samps = link(CIRC)
    #sorted_genes = sorted(GENE.items(), key=operator.itemgetter(1),reverse=True)
    sorted_genes = sorted(GENE.items(), key = lambda x: (x[1][1],x[1][0]),reverse=True)
    with open("germline.v4.karyotype.txt","w") as o:
        cnt = 0
        for i in sorted_genes:
            print(i)
            cnt += 1
            ax = "axis"+str(cnt)
            path = gene2path[i[0]]
            print(path)
            out = ["chr","-",ax,i[0],str(1),str(i[1][0]+1),path2rgb[path]]
            o.write(" ".join(out))
            o.write("\n")
    o.close()
    gene_ord = list(zip(*sorted_genes))[0]
    
    #Here I'm going to make 2 histograms sorted by the first (age)\
    #And second is the total number of somatic mutations
    ax_num = 0
    with open("germline.v4.age.histogram.txt","w") as oage:
        with open("germline.v4.somatic.txt","w") as osom:
            for axis in gene_ord:
                SLICE = {}
                ax_num += 1
                ax = "axis"+str(ax_num) #Just for reporting
                for samp in CIRC:
                    w = CIRC[samp]
                    if axis in w.genes:
                        SLICE[w.name] = w
                
                sorted_samps = sorted(SLICE.values(), key=operator.attrgetter('age'), reverse=True)
                pos = 0
                for s in sorted_samps:
                    pos += 1
                    #append to Age histogram
                    s.add_location(ax_num,pos)
                    aout = [ax,str(pos),str(pos+1),str(s.age)]
                    oage.write(" ".join(aout))
                    oage.write("\n")
                    
                    #append to somatic count histogram
                    sout = [ax,str(pos),str(pos+1),str(s.muts)] 
                    osom.write(" ".join(sout))
                    osom.write("\n")
        osom.close()
    oage.close()

    #Here I'm going to make the files with the links
    with open("germline.v4.links.txt","w") as olink:
        for s in CIRC: 
            w = CIRC[s]
            if w.name in linked_samps:
                for i in range(0,len(w.axis)):
                    for j in range(i+1,len(w.axis)):
                        a1 = "axis"+str(w.axis[i])
                        a2 = "axis"+str(w.axis[j])
                        p1s = str(w.pos[i]) #Pos 1 start
                        p1e = str(w.pos[i]+1) #Pos 1 end
                        p2s = str(w.pos[j])
                        p2e = str(w.pos[j]+1) 
                        color = "color=black_a4"
                        search = "FAN"
                        rescol = next((True for fan in w.genes if search in fan), False)
                        if rescol:
                            color="color=orange"
                        lout = [a1,p1s,p1e,a2,p2s,p2e,color]
                        olink.write(" ".join(lout))
                        olink.write("\n")
    olink.close()

def make_circos(afname, bfname, goodsamps, gene2path, path2rgb):
    CIRCOS = {}
    GENES = {}
    with open(afname,"r") as f:
        aheader = f.readline()
        for line in f: 
            line = line[:-1]
            l = line.split("\t")
            n = l[0]
            g = l[1]
            c = l[2]
            ##### Deal with types and missing data  
            a = l[136]
            if a == "NA":
                a = 0
            else:
                a = int(a)
            
            if n in goodsamps:
                if n in CIRCOS:
                    this = CIRCOS[n]
                    this.add_genes(g)
    
                else:
                    w = Wedge(n,c,a)
                    w.add_genes(g)
                    CIRCOS[n] = w
    
                ###### Count genes #######
                if g in GENES:
                    GENES[g][0] += 1
                else:
                    GENES[g] = [1,gene2path[g]]
        f.close()
             
    ###### Add somatic Mutation Counts ########
    with open(bfname, "r") as f:
        bheader = f.readline()
        for line in f:
            line = line[:-1]
            l = line.split("\t")
            nm = l[0]
            mutcnt = l[1] 
        
            if nm in CIRCOS:
                this = CIRCOS[nm]
                this.set_som_muts(mutcnt)

    orderWedges(CIRCOS,GENES,gene2path,path2rgb)
    


if __name__ == "__main__":
    import sys
    goodsamps = file2list_col1(sys.argv[1],1)
    gene2path = file2dict(sys.argv[3],1)
    path2rgb = file2dict(sys.argv[4],0)
    make_circos(sys.argv[2],sys.argv[1],goodsamps,gene2path,path2rgb) 
