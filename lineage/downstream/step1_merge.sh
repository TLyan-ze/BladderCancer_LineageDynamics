cat sample_id |while read l1 l2 l3
do
    cat raw/${l1}-AP.re.csv|sed "s/^/${l2}\\t${l2}./g"|sed "s/${l2}\\t${l2}.cellBC/sampleID\\tcellBC/g"  > clean/${l2}.alleleTable.csv
done
echo -e "Tumor\tsampleID\tcellBC\tintBC\tallele\tr1\tr2\tr3\tlineageGrp\tUMI\treadCount" > AP.alleleTable.unfiltered.txt
cat sample_id |while read l1 l2 l3; do cat clean/${l2}.alleleTable.csv|grep -v sampleID |awk -F "\t" '{print $1"_"$8"\t"$0}' ; done >> AP.alleleTable.unfiltered.txt
echo -e "Tumor" > tumor_list.txt
awk -F "\t" '{print $1}' AP.alleleTable.unfiltered.txt|grep -wv Tumor|sort|uniq >> tumor_list.txt

#awk -F"\t" '{ if ((index($6, "None") && index($7, "None") && index($8, "None")) ) { } else { print $0 } }' AP.alleleTable.unfiltered.txt > temp_info

