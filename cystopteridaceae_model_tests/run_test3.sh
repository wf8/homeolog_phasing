rm -rf logs3
mkdir logs3
for r in {1..4}
do
    for m in {1..3}
    do
        echo "model=${m}; run=${r}; source(\"test3/src/cystopteridaceae_model_test.Rev\")" | \
            rb test_3 > logs3/${m}_${r}.log &
        sleep 0.01
    done
done
