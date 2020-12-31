
#parallel -j 50 echo "'sim_num={1}; source(\"src/sims.Rev\");' | rb" ::: {0..999} 
#time parallel --eta -j 50 echo "'sim_num={}; bayes_factors=FALSE; dummy_tip=TRUE; source(\"src/inference.Rev\");' | rb" ::: {0..999}
#time parallel --eta -j 50 echo "'sim_num={}; bayes_factors=FALSE; dummy_tip=FALSE; source(\"src/inference.Rev\");' | rb" ::: {0..999}
#time parallel --eta -j 50 echo "'sim_num={}; bayes_factors=TRUE; dummy_tip=TRUE; source(\"src/inference.Rev\");' | rb" ::: {0..999}
#time parallel --eta -j 50 echo "'sim_num={}; bayes_factors=TRUE; dummy_tip=FALSE; source(\"src/inference.Rev\");' | rb" ::: {0..999}

time parallel --eta -j 40 echo "'sim_num={}; bayes_factors=TRUE; dummy_tip=FALSE; source(\"src/inference.Rev\");' | rb" ::: {0..49} {500..549}
time parallel --eta -j 40 echo "'sim_num={}; bayes_factors=TRUE; dummy_tip=TRUE; source(\"src/inference.Rev\");' | rb" ::: {0..49} {500..549}

