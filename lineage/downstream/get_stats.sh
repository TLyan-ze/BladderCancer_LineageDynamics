
input=$1


name=`basename $1`
cellBC=`cat $1|grep -v cellBC|awk '{print $1}'|sort|uniq|wc -l`
intBC=`cat $1|grep -v cellBC|awk '{print $2}'|sort|uniq|wc -l`
lineageGrp=`cat $1|grep -v cellBC|awk '{print $7}'|sort|uniq|wc -l`
UMI=`cat $1|grep -v cellBC|awk '{print $8}'|sort|uniq|wc -l`
#echo "name\tcellBC\tintBC\tlineageGrp\tUMI"
echo "$name\t$cellBC\t$intBC\t$lineageGrp\t$UMI"