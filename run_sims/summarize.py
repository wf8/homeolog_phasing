import os,sys
from ete3 import Tree
import pandas as pd
# Requires MESS: https://github.com/tgvaughan/mess
# Requires ETE3: etetoolkit.org



treetypes = ["unlinked","shared"]
nindv = [1,2,10,50,100]
numsims = 100
tetra_labels = ["tetra_A","tetra_B"]
hexa_labels = ["hexa_A","hexa_B","hexa_C"]

sys.stdout.write("Type\tnIndiv\tSimNum\tTetraDist\tHexaABDist\tHexaACDist\tHexaBCDist\tTetraAMinDip\tTetraBMinDip\tHexaAMinDip\tHexaBMinDip\tHexaCMinDip\tESS\tAvgTetraPP\tAvgHexaPP\n")
for t in treetypes:
    os.chdir("{}_gene_trees".format(t))
    for n in nindv:
        os.chdir("nindv_{}".format(n))
        #Calculate ESS with MESS
        for s in range(1,numsims+1):
            os.chdir("data{}".format(s))
            #distances between polyploid subgenomes on species tree
            simtree = Tree("species_tree.tree",format=1)
            dip_names = [x for x in simtree.traverse() if x.name.startswith("dip_")]
            for node in simtree.traverse():
                if node.is_leaf:
                    if node.name.startswith(tetra_labels[0]):
                        tetraA = node
                        #distance between a subgenome tip and the nearest diploid
                        tetraA_mindip = min([x.get_distance(tetraA) for x in dip_names])
                    elif node.name.startswith(tetra_labels[1]):
                        tetraB = node
                        tetraB_mindip = min([x.get_distance(tetraB) for x in dip_names])
                    elif node.name.startswith(hexa_labels[0]):
                        hexaA = node
                        hexaA_mindip = min([x.get_distance(hexaA) for x in dip_names])
                    elif node.name.startswith(hexa_labels[1]):
                        hexaB = node
                        hexaB_mindip = min([x.get_distance(hexaB) for x in dip_names])
                    elif node.name.startswith(hexa_labels[2]):
                        hexaC = node
                        hexaC_mindip = min([x.get_distance(hexaC) for x in dip_names])
            tetra_dist = tetraA.get_distance(tetraB)
            hexaAB_dist = hexaA.get_distance(hexaB)
            hexaAC_dist = hexaA.get_distance(hexaC)
            hexaBC_dist = hexaB.get_distance(hexaC)
            
            
            
            os.chdir("../output{}".format(s))
            # ess of rb run
            mess_call = "mess phasing.log -o branch --burnin 10 | tail -1"
            mess_stream = os.popen(mess_call)
            mess_output = mess_stream.read()
            ess = mess_output.split()[1]
            # posterior probability of phase calls
            tetra_output_fn = "joint_phase_probs_tetraploid_1.csv"
            tetra_probs = pd.read_csv(tetra_output_fn)
            avgTetraPP = tetra_probs.groupby("locus").max().mean()["joint_prob"]
            hexa_output_fn =  "joint_phase_probs_hexaploid_1.csv"
            hexa_probs = pd.read_csv(hexa_output_fn)
            avgHexaPP = hexa_probs.groupby("locus").max().mean()["joint_prob"]

            # write all to output
            output_string = "{}\t{}\t{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{}\t{:0.2f}\t{:0.2f}\n".format(t,n,s,
                tetra_dist,
                hexaAB_dist,hexaAC_dist,hexaBC_dist,
                tetraA_mindip,tetraB_mindip,hexaA_mindip,hexaB_mindip,hexaC_mindip,
                ess,
                avgTetraPP,avgHexaPP)
            sys.stdout.write(output_string)
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")
    