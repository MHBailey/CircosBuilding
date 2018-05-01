import numpy as np
import pandas as pd


class Wedge:
    def __init__(self,name): #for each gene
        self.name = name
        self.panfreq = 0
        self.canfreqs = {} #k=Cancer v=cancerspec freq
        self.cancers = []
        self.pangenes = []
        self.can_score = {} #k=Cancer v=orginial cgatscore
        self.newcan_score = {} #k=Cancer v=corredted cgatscore score (0.05)
        self.specificity = {} #k=Cancer v=p-value
        self.cgat_cans = []
        self.mark_cans = []
        self.axis = 0
        self.position = 0
        self.organs = []
        self.level = ""
        self.topscore = 0
        self.pancan = False
        self.lda = 0
        self.history = 0
        self.onco = False
        self.li = False
        self.rescuable = False
        self.cgatg = False
        self.new = []
        self.where = {} #K=cancer, v="string where"
        self.status = {} #K=cancer, v="string"
        self.cgc = False #Unless made true
        self.methods = {} #K=cancer, v="comma sep list"
        self.marker = {} #K=cancer, v="yes or no" in marker paper.
        self.mmiss = {} #K=cancer, v="yes or no" Marginally Missed
        self.tmiss = {} #K=cancer, v="yes or no" Totally missed
        self.oncts = {} #K=cancer, v="20/20plus descriptor"
        self.decision = ""
        

    def __lt__(self, other):
        return len(self.organs) < len(other.organs)

    def add_name(self,name):
        self.name = name

    def add_panfreq(self,freq):
        self.panfreq = freq

    def set_lda(self,lda):
        self.lda = lda

    def add_canfreq(self,cancer,cf,pf):
        if cancer != "PANCAN":
            cfreq = float(cf)
            self.canfreqs[cancer] = cfreq
        #else:
        #    self.canfreqs[cancer] = pf
        #    self.panfreq = pf

    def add_cancer(self,cancer):
        if cancer == "PANCAN":
            self.pancan = True
            #self.cancers.append(cancer)
        else:
            self.cancers.append(cancer)

    def add_new(self,M_NEW):
        if self.name in M_NEW:
            self.new.append(M_NEW[self.name])

    def add_organ(self,organ):
        self.organs.append(organ)

    def add_wsig(self,c,cgatc,m):
        if cgatc == "Yes":
            self.cgat_cans.append(c)
        if m == "Yes":
            self.mark_cans.append(m)

    def add_specificity(self,c,SPECIFICITY):
        if self.name in SPECIFICITY:
            if c in SPECIFICITY[self.name]:
                self.specificity[c] = SPECIFICITY[self.name][c]

    def add_where(self,c,where):
        self.where[c]=where

    def add_status(self,c,stat):
        self.status[c] = stat

    def add_methods(self,c,methods):
        if methods == "":
            self.methods[c] = "none"
        else:
            self.methods[c] = methods

    def add_marker(self,c,m):
        self.marker[c] = m

    def add_mmiss(self,c,mmiss):
        self.mmiss[c] = mmiss

    def add_tmiss(self,c,tmiss):
        self.tmiss[c] = tmiss

    def add_oncts(self,c,oncts):
        ot = oncts
        if oncts == "":
            ot = "nope" 
        elif oncts == "possible tsg":
            ot = "ptsg"
        elif oncts == "possible oncogene":
            ot = "ponco"
        #print(c,ot) 
        self.oncts[c] = ot

    def set_axis(self,axis):
        self.axis = axis

    def set_position(self,position):
        self.position = position

    def set_canScore(self,can,score):
        self.can_score[can] = score

    def set_canNewScore(self,can,newscore):
        self.newcan_score[can] = newscore

    def set_level(self,level):
        if self.level == "":
            self.level = level
        else:
            if level == "gold":
                self.level = "gold"
            elif self.level == "gold":
                self.level = "gold"
            elif self.level == "silver" and level == "bronze":
                self.level = "silver"
            elif self.level == "bronze" and level == "silver":
                self.level = "silver"
            elif self.level == "silver" and level == "silver":
                self.level = "silver"
            elif level == "bronze":
                self.level = "bronze"
            elif level == "white":
                self.level = "white"
            else:
                sys.exit("set_level:::something went wrong: " + level +" __ " + self.level)

    def set_topscore(self,score):
        if score == "NA":
            score = 0
        if self.topscore < float(score):
            self.topscore = float(score)
    

    def set_decision(self,decision):
        self.decision = decision
    
    def set_history(self,HISTORY):
        if self.name in HISTORY:
            self.history = HISTORY[self.name]
        else:
            self.history = 0

    def set_li(self,LI):
        if self.name in LI:
            self.li = True

    def set_rescuable(self,M_RESCUE):
        if self.name in M_RESCUE:
            self.rescuable = True

    def set_cgata(self,cgata):
        if cgata == "Yes":
            self.cgatg = True

    def set_cgc(self,cgc):
        if cgc == "Yes":
            self.cgc = True

    def print_me(self):
        print(self.name,str(self.freq),self.cancers)

    def get_cancers(self):
        return list(self.cancers)



def file2list(cfname):
    out = []
    with open(cfname, "r") as f:
        for line in f:
            out.append(line.strip())
    return out

def get_level(score):
    lev = ""

    if score == "NA":
        lev = "white"
    else:
        s = float(score)
        if s > 1.5:
            lev = "gold"
        elif s == 1.5:
            lev = "silver"
        elif s > 0:
            lev = "bronze"
        elif s == 0:
            lev = "white"
    return lev



def make_circos(afname):
    CIRCOS = {}
    with open(afname,"r") as f:
        header = f.readline()
        for line in f:
            line = line[:-1]
            l = line.split("\t")
            g = l[0]
            c = l[1]
            cf = l[4]
            pf = float(l[5])
            s = l[7]
            oncts = l[2]
            d = l[3]

            lev = get_level(s)
 
            if g in CIRCOS:
                this = CIRCOS[g]
                this.add_canfreq(c,cf,pf)
                this.add_cancer(c)
                #this.add_organ(CAN[c])
                this.set_level(lev)
                this.set_topscore(s)
                this.set_canScore(c,s)
                #this.set_canNewScore(c,ns)
                #this.add_new(M_NEW)
                #this.add_wsig(c,cgatc,m)
                #this.set_li(LI)
                #this.add_specificity(c,SPECIFICITY)
                #this.add_where(c,where)
                #this.add_status(c,status)
                #this.add_methods(c,tool)
                #this.add_marker(c,m)
                #this.add_mmiss(c,mmiss)
                #this.add_tmiss(c,tmiss)
                this.add_oncts(c,oncts)

            else:
                w = Wedge(g)
                w.add_panfreq(pf)
                w.add_canfreq(c,cf,pf)
                w.add_cancer(c)
                #w.add_organ(CAN[c])
                w.set_level(lev)
                w.set_topscore(s)
                w.set_canScore(c,s)
                #w.set_canNewScore(c,ns)
                #w.set_lda(lda)
                #w.set_history(HISTORY)
                #w.add_new(M_NEW)
                #w.add_wsig(c,cgatc,m)
                #w.set_rescuable(M_RESCUE)
                #w.set_li(LI)
                #w.add_specificity(c,SPECIFICITY)
                #w.set_cgata(cgata)
                #w.add_where(c,where)
                #w.set_cgc(cgc)
                #w.add_status(c,status)
                #w.add_methods(c,tool)
                #w.add_marker(c,m)
                #w.add_mmiss(c,mmiss)
                #w.add_tmiss(c,tmiss)
                w.add_oncts(c,oncts)
                CIRCOS[g] = w
    f.close()
    return CIRCOS





def make_circos_files(circos,can_order,gene_order):
    ORG_K = {}
    SINGLE = set()
    PLUS = set()
    AXIS = {}
    PANAXIS = []
    axis_ord = ["Multi"]

    for g in circos:
        w = circos[g]
        if len(w.cancers) > 1:
            PLUS.add(w.name)
        else:
            SINGLE.add(w.name)


#MAKE THE KARYOTPE FILES
    with open("circos.v8.karyotype.txt","w") as o:
        out = ["chr - axis1 MULTI 1",str(len(PLUS)+1),"pblue"]
        o.write(" ".join(out))
        o.write("\n")

        for g in circos:
            w = circos[g]
            if w.name in SINGLE and len(w.cancers) > 0: #This commented statement is taken care of in the gerenation of the wedge
                if len(w.cancers) > 1:
                    sys.exit("Should have only seen this gene once")
                can = w.cancers[0]
                if can in ORG_K:
                    ORG_K[can] += 1
                else:
                    ORG_K[can] = 1

        axis_num = 1
        for org in ORG_K: #Note that 'org' here refers to a cancer type 
            axis_num += 1
            ORG_K[org] += 1
            axis_str = "axis"+str(axis_num)
            axis_ord.append(org)
            out = ["chr -",axis_str,org,"1",str(ORG_K[org]),"pblue"]
            o.write(" ".join(out))
            o.write("\n")
    o.close()


#MAKE THE LABELS FILES
    with open("circos.v8.geneLabels.txt","w") as o:
        AXIS["Multi"] = []
        for z in gene_order:
            for g in PLUS: #This is where I can sort them... 
                if g == z:
                    AXIS["Multi"].append(circos[g])

        for g in SINGLE:
            w = circos[g]
            if len(w.cancers) > 0: #This step was 'w.pancan == False' taken care of in making the wedge
                org = w.cancers[0] #Note that 'org' here refers to a cancer type 
                if org in AXIS:
                    AXIS[org].append(w)
                else:
                    AXIS[org] = []
                    AXIS[org].append(w)
            else: #this is going to be the Pancan loop
                PANAXIS.append(w)


        axis_num = 0
        for a in axis_ord:
            axis_num += 1
            axis_str = "axis"+str(axis_num)
            position = 0
            for w in AXIS[a]:
                position += 1
                w.set_axis(axis_num)
                w.set_position(position)
                out = [axis_str,str(position),str(position+1),w.name]
                o.write(" ".join(out))
                o.write("\n")
    o.close()


##############################MAKE PANCAN###################################
#Make the Pancan circle karyotype file
    with open("circos.v9.karyotype.txt","w") as o:
        kar = len(PANAXIS)+1
        out = ["chr - axis1 PanCancer 1",str(kar),"pblue"]
        o.write(" ".join(out))
    o.close()

#Make the Pancan circle genelist 
    with open("circos.v9.geneLabels.txt","w") as o:
        position = 0
        for p in PANAXIS:
            position += 1
            axis_str = "axis1"
            out = [axis_str,str(position),str(position+1),p.name]
            o.write(" ".join(out))
            o.write("\n")
    o.close()
#Make the Pancan circle frequency tables
    pan_textfile = []
    with open("Histograms/circos.v9.histogram.pancan.txt","w") as o:
        position = 0
        for p in PANAXIS:
            position += 1
            axis_str = "axis1"
            out = [axis_str,str(position),str(position+1),str(p.panfreq)]
            if float(p.panfreq) > .03:
                tex = [axis_str,str(position),str(position+1),str("%.1f" % (float(p.panfreq)*100))+"%\n"]
                pan_textfile.append(tex)
            o.write(" ".join(out))
            o.write("\n")
    o.close()
#Make the Pancan text of frequencies
    with open("Text/circos.v9.text.greater3.txt","w") as o:
        for t in pan_textfile:
            o.write(" ".join(t))
    o.close()

#Make the Pancan tiles for gold and silver 
    with open("Tiles/circos.v9.goldsilver.txt","w") as o:
        position = 0
        for p in PANAXIS:
            position += 1
            axis_str = "axis1"
            olevel = "score"+str("%.1f" % p.topscore)
            out = [axis_str,str(position),str(position+1),olevel]
            o.write(" ".join(out))
            o.write("\n")
    o.close()
############################################################################

################################MAKE CANCERS###################################
#ADD Frequency and Track text to add to samples with greater than 10% mutation frequency 
    textfile = []
    with open("Histograms/circos.v8.histogram.all.txt","w") as o:
        axis_num = 0
        multi_ord = ["Multi"]
        for a in multi_ord:
            axis_num += 1
            axis_str = "axis"+str(axis_num)
            position = 0
            for w in AXIS[a]:
                position += 1
                w.set_axis(axis_num)
                w.set_position(position)
                out = [axis_str,str(position),str(position+1),str(w.panfreq)]
                if float(w.panfreq) > .1:
                    tex = [axis_str,str(position),str(position+1),str("%.1f" % (float(w.panfreq)*100))+"%\n"]
                    textfile.append(tex)
                o.write(" ".join(out))
                o.write("\n")
    o.close()

    with open("Histograms/circos.v8.histogram.all.txt","a") as o:
        axis_num = 1
        for a in axis_ord[1:]:
            axis_num += 1
            axis_str = "axis"+str(axis_num)
            position = 0
            for w in AXIS[a]:
                position += 1
                w.set_axis(axis_num)
                w.set_position(position)
                out = [axis_str,str(position),str(position+1),str(w.canfreqs[a])]
                if float(w.canfreqs[a]) > .1:
                    tex = [axis_str,str(position),str(position+1),str("%.1f" % (w.canfreqs[a]*100))+"%\n"]
                    textfile.append(tex)
                o.write(" ".join(out))
                o.write("\n")
    o.close()

    with open("Text/circos.v8.text.greater10.txt","w") as o:
        for t in textfile:
            o.write(" ".join(t))
    o.close()

#ADD Gold Silver tiles
    with open("Tiles/circos.v8.goldsilver.txt","w") as o:
        axis_num = 0
        for a in axis_ord:
            axis_num += 1
            axis_str = "axis"+str(axis_num)
            position = 0
            for w in AXIS[a]:
                position += 1
                w.set_axis(axis_num)
                w.set_position(position)
                olevel = "score"+str("%.1f" % w.topscore)
                out = [axis_str,str(position),str(position+1),olevel]
                o.write(" ".join(out))
                o.write("\n")
    o.close()

#ADD Tumor suppressor Oncogene tiles 
    with open("Tiles/circos.v8.tsonc.txt","w") as o:
        axis_num = 1
        axis_ord_can = axis_ord[1:]
        for a in axis_ord_can:
            axis_num += 1
            axis_str = "axis"+str(axis_num)
            position = 0
            for w in AXIS[a]:
                position += 1
                w.set_axis(axis_num)
                w.set_position(position)
                out = [axis_str,str(position),str(position+1),w.oncts[a]]
                o.write(" ".join(out))
                o.write("\n")
    o.close()


#MAKE HEATMAP TRACKS
    R1 = .8
    with open("Plot.v8.conf","w") as p:
        for c in can_order:
            #CANCERTYPES
            ofname = "Heats/heat."+c+".txt"
            with open(ofname, "w") as o:
                for w in AXIS["Multi"]:
                    score = 0
                    if c in w.can_score:
                        score = w.can_score[c]
                    out = ["axis1",str(w.position), str(w.position+1),str(score),"color=score"+str("%.1f" % float(score))]
                    o.write(" ".join(out))
                    o.write("\n")
            o.close()
            #CONFS
            p.write("\n<plot>\n")
            p.write("type    = heatmap\n")
            p.write("file    = "+ofname+"\n")
            p.write("r1      = "+str(round(R1,3))+"r\n")
            R1 -= 0.015
            p.write("r0      = "+str(round(R1,3))+"r\n")
            p.write("</plot>\n")

    p.close()



if __name__ == "__main__":
    import sys
    CIRCOS = make_circos(sys.argv[1])
    can_order = file2list(sys.argv[2])
    gene_order = file2list(sys.argv[3])
    make_circos_files(CIRCOS,can_order,gene_order)

