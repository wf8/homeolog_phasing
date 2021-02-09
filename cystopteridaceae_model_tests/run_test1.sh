rm -rf logs1
mkdir logs1
for r in {1..4}
do
    for m in {1..2}
    do
        echo "model=${m}; run=${r}; source(\"test1/src/cystopteridaceae_model_test.Rev\")" | \
            rb test_1 > logs1/${m}_${r}.log &
        sleep 0.01
    done
done
