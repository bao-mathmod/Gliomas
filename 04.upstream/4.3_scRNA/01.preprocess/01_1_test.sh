# for f in \
# /mnt/18T/chibao/gliomas/data/output_cell/scRNA_2/set1/PRJNA577146/*_features.tsv.gz
# do
#   echo ">>> $f"
#   # how many columns?
#   cols=$(zcat "$f" | awk -F'\t' 'NR==1{print NF; exit}')
#   echo "columns: $cols"
#   # show first 5 lines with visible tab markers
#   zcat "$f" | head -n 5 | sed 's/\t/[TAB]/g'
#   echo
# done

for f in \
/mnt/18T/chibao/gliomas/data/output_cell/snRNA/set2/PRJNA901214/*_features.tsv.gz
do
  cols=$(zcat "$f" | awk -F'\t' 'NR==1{print NF; exit}')
  echo ">>> $f  (cols=$cols)"
  if [ "$cols" -ge 2 ]; then
    # sample 20 rows from each of the first two columns
    echo "col1 sample:"; zcat "$f" | cut -f1 | head -n 20
    echo "col2 sample:"; zcat "$f" | cut -f2 | head -n 20
  else
    echo "Only 1 column present."
    zcat "$f" | head -n 20
  fi
  echo
done

