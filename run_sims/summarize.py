import os,sys
from ete3 import Tree
import pandas as pd
import subprocess
# Requires MESS: https://github.com/tgvaughan/mess
# Requires ETE3: etetoolkit.org

# Run this from within the sims directory, with the types as separate directories

nindv = [0.0001, 0.001,0.1,0.5,1,2]
tetra_labels = ["tetra_A","tetra_B"]
tetra_species = ["tetra_X_1","tetra_Y_1"]
hexa_labels = ["hexa_A","hexa_B","hexa_C"]
hexa_species = ["hexa_X_1","hexa_Y_1","hexa_Z_1"]

sys.stdout.write("Type\tnIndiv\tSimNum\tTetraDist\tHexaABDist\tHexaACDist\tHexaBCDist\tTetraAMinDip\tTetraBMinDip\tHexaAMinDip\tHexaBMinDip\tHexaCMinDip\tESS\tAvgTetraPP\tAvgHexaPP\tMinTetraPP\tMinHexaPP\tTetraCorrect\tHexaCorrect\tKF\tRF\tseq_length\n")
for n in nindv:
    os.chdir("nindv_{}".format(n))
    if n == 0.0001:
        numsims = 500
    elif n == 0.001:
        numsims = 400
    else:
        numsims = 100

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
        
        # get KF and RF tree distances
        stream = os.popen("Rscript ../../get_tree_distances.R")
        out = stream.read()
        out = out.split(',')
        KF = out[0]
        RF = out[1]

        # get sequence length
        stream = os.popen('sed -n 2p alignment1.fasta | wc -c')
        seq_length = stream.read().strip()
        
        # summarize posteriors 
        os.chdir("../output{}".format(s))
        subprocess.call(["Rscript", "../../../src/summarize_posteriors.R", ".", "."])

        # ess of rb run
        mess_call = "mess phasing.log -o branch --burnin 10 | tail -1"
        mess_stream = os.popen(mess_call)
        mess_output = mess_stream.read()
        ess = mess_output.split()[1]
        # posterior probability of phase calls
        tetra_output_fn = "joint_phase_probs_tetraploid_1.csv"
        tetra_probs = pd.read_csv(tetra_output_fn)
        tetraBestLoc = tetra_probs.groupby("locus")["joint_prob"].idxmax()
        tetraJMAPs = tetra_probs.loc[tetraBestLoc]
        avgTetraPP = tetraJMAPs.mean()["joint_prob"]
        minTetraPP = tetraJMAPs.min()["joint_prob"]
        if len(tetraJMAPs[tetra_species].apply(lambda x: ''.join(x),axis=1).unique()) == 1:
            tetraCorrect="True"
        else:
            tetraCorrect="False"
        
        hexa_output_fn =  "joint_phase_probs_hexaploid_1.csv"
        hexa_probs = pd.read_csv(hexa_output_fn)
        hexaBestLoc = hexa_probs.groupby("locus")["joint_prob"].idxmax()
        hexaJMAPs = hexa_probs.loc[hexaBestLoc]
        avgHexaPP = hexaJMAPs.mean()["joint_prob"]
        minHexaPP = hexaJMAPs.min()["joint_prob"]
        if len(hexaJMAPs[hexa_species].apply(lambda x: ''.join(x),axis=1).unique()) == 1:
            hexaCorrect="True"
        else:
            hexaCorrect="False"

        
        
        t = 'unlinked'
        # write all to output
        output_string = "{}\t{}\t{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{}\t{}\t{}\t{}\t{}\n".format(t,n,s,
            tetra_dist,
            hexaAB_dist,hexaAC_dist,hexaBC_dist,
            tetraA_mindip,tetraB_mindip,hexaA_mindip,hexaB_mindip,hexaC_mindip,
            ess,
            avgTetraPP,avgHexaPP,minTetraPP,minHexaPP,
            tetraCorrect,hexaCorrect,
            KF, RF, seq_length)
        sys.stdout.write(output_string)
        os.chdir("..")
    os.chdir("..")
    
