#!/bin/zsh

L=(0 5000 20000 100000 1000000)

# for i in {1..$((${#L}-1))}; do
#     #echo -n $(($L[$i]+1))--$L[$i+1]" "
#     echo -n "<"$L[$i+1]" "
#     cat $1 | fgrep -v Indel| awk 'length($3)==1 && length($4)==1' | gsort -k1,1 -k2,2n |python src/add_closest_in_same.py | awk '$NF<='$L[$i+1]'{ n+=1 } END {print n}'
# done

for i in {1..$((${#L}-1))}; do
    #echo -n $(($L[$i]+1))--$L[$i+1]" "
    echo -n "<"$L[$i+1]" "
    n_close=$(cat $1 | awk 'length($3)==1 && length($4)==1' | gsort -k1,1 -k2,2n |python src/add_closest_in_same.py | awk '$NF<='$L[$i+1]'{ n+=1 } END {print n}')

    n_all=$(cat $1 | awk 'length($3)==1 && length($4)==1' | wc -l | awk '{print $1}')
    

    n_close_random=$(for j in {1..100}; do
    	cat $1 | awk 'length($3)==1 && length($4)==1' | gsort -k1,1 -k2,2n | awk '{print $1"_"$2,$5}' | python ~/Scripts/permute_pheno.py | gsed 's/_/\t/' | awk '{print $1,$2,"NA","NA",$3}' | python src/add_closest_in_same.py | awk 'BEGIN {n=0} $NF<='$L[$i+1]'{ n+=1 } END {print n}'
    done | awk 'BEGIN {sum=0;n=0} { sum+=$1;n+=1} END {print sum/n}')

    frac_obs=$(echo "scale=4;" $n_close / $n_all | bc -l)
    frac_random=$(echo "scale=4;" $n_close_random / $n_all | bc -l)
    frac_dependent=$(echo "scale=4; ("$n_close "-" $n_close_random")" "/" $n_all | bc -l)

    echo $n_close $n_all $frac_obs $n_close_random $frac_random $frac_dependent

done