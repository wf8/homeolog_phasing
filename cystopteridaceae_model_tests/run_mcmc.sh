rm -rf logs
mkdir logs
for r in {1..10}
do
    echo "run=${r}; source(\"mcmc/src/cystopteridaceae.Rev\")" | \
        rb cysto_1 > logs/${r}.log &
    sleep 0.01
done
