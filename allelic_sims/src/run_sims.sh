
parallel -j 50 echo "'sim_num={1}; source(\"src/sims.Rev\");' | rb" ::: {0..999} 
time parallel --eta -j 50 echo "'sim_num={}; dummy_tip=TRUE; source(\"src/inference.Rev\");' | rb" ::: {0..999}
time parallel --eta -j 50 echo "'sim_num={}; dummy_tip=FALSE; source(\"src/inference.Rev\");' | rb" ::: {0..999}

