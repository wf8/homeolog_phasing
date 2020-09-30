
parallel -j 50 echo "'sim_num={1}; pop_size={2}; source(\"msc_sims.Rev\");' | rb" ::: {1..100} ::: 0.001 0.1 0.5 1 2
parallel -j 50 echo "'sim_num={1}; pop_size={2}; source(\"msc_sims.Rev\");' | rb" ::: {101..400} ::: 0.001
parallel -j 50 echo "'sim_num={1}; pop_size={2}; source(\"msc_sims_seq_size.Rev\");' | rb" ::: {1..500} ::: 0.0001

for i in 0.001 0.1 0.5 1 2
do
    cd nindv_"$i"
    time parallel --eta -j 50 echo "'sim_num={}; source(\"../msc_inference.Rev\");' | rb" ::: {1..100}
    cd ..
done
i=0.001
cd nindv_"$i"
time parallel --eta -j 50 echo "'sim_num={}; source(\"../msc_inference.Rev\");' | rb" ::: {101..400}
cd ..
i=0.0001
cd nindv_"$i"
time parallel --eta -j 50 echo "'sim_num={}; source(\"../msc_inference.Rev\");' | rb" ::: {1..500}
cd ..

