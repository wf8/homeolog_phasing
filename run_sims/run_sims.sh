
parallel -j 12 echo "'sim_num={1}; pop_size={2}; seed({1}); source(\"../msc_sims.Rev\");' | rb" ::: {1..100} ::: 1 2 10 50 #100

for i in 1 2 10 50 #100
do
cd nindv_"$i"
time parallel --eta -j 25 echo "'sim_num={}; seed({}); source(\"../../msc_inference.Rev\");' | rb" ::: {1..100}
cd ..
done

