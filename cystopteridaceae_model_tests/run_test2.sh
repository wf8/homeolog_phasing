rm -rf logs2
mkdir logs2
for r in {1..4}
do
    for m in {1..2}
    do
        echo "model=${m}; run=${r}; source(\"test2/src/cystopteridaceae_model_test.Rev\")" | \
            rb test_2 > logs2/${m}_${r}.log &
        sleep 0.01
    done
done
