
parallel -j 50 echo "'sim_num={1}; pop_size={2}; seed({1}); source(\"msc_sims.Rev\");' | rb" ::: {1..100} ::: 0.001 0.1 0.5 1 2

for i in 0.001 0.1 0.5 1 2
do
    cd nindv_"$i"
    time parallel --eta -j 50 echo "'sim_num={}; seed({}); source(\"../msc_inference.Rev\");' | rb" ::: {1..100}
    cd ..
done

