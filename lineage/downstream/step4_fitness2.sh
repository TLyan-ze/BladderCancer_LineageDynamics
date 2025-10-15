cd fitnesses
ls *tree.fitness.txt|awk -F "_tree" '{print $1}' > temp1
cat temp1|while read l1
do
    cat ${l1}_tree.fitness.txt|grep -v False|awk -F "\t" '{print $8"\t"$9}' > mean_fitness.${l1}.txt
done


ls *fitness.txt|awk -F "_tree.fitness.txt" '{print $1}'|while read l1
do
    cat ${l1}_tree.fitness.txt |grep -v False|awk -F'\t' 'BEGIN {print "\tmean_fitness"} NR==1 {for(i=1;i<=NF;i++) if($i=="mean_fitness") col=i} NR>1 {print $(col-1)"\t"$col}' > mean_fitness.${l1}.txt
done
