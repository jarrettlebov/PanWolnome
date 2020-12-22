# PanWolnome
This is a work-in-progress repository dedicated to exploring the Wolbachia pangenome. 

<!-- MarkdownTOC autolink="true" -->

- [Phylogeny Directories](#phylogeny-directories)
- [Test](#test)
- [Test2](#test2)
- [T3](#t3)
- [Isolate strain accession list from mugsy generated multiple alignment file \(.maf\)](#isolate-strain-accession-list-from-mugsy-generated-multiple-alignment-file-maf)
- [this set of commands makes the processed directory, removes any prior contents, and then reads each line of the accesssion.list file and creates separate .fna file for each item in accession.list, and then starts a fasta header within each .fna file that contains the file name.](#this-set-of-commands-makes-the-processed-directory-removes-any-prior-contents-and-then-reads-each-line-of-the-accesssionlist-file-and-creates-separate-fna-file-for-each-item-in-accessionlist-and-then-starts-a-fasta-header-within-each-fna-file-that-contains-the-file-name)
- [This command adds a parses the mugsy produced .maf file and adds aligned sequences corresponding to each strain to the .fna files produced in the previous step.](#this-command-adds-a-parses-the-mugsy-produced-maf-file-and-adds-aligned-sequences-corresponding-to-each-strain-to-the-fna-files-produced-in-the-previous-step)
- [concatenate aligned sequence fastas, change directories to the Mugsy_out directory, and invoke mothur](#concatenate-aligned-sequence-fastas-change-directories-to-the-mugsy_out-directory-and-invoke-mothur)
- [Run the following mothur commands](#run-the-following-mothur-commands)
- [change the name of the mothur output. concat_alignment.fna is the file needed for running iqtree](#change-the-name-of-the-mothur-output-concat_alignmentfna-is-the-file-needed-for-running-iqtree)
- [change fasta headers of concat_alignment.fna](#change-fasta-headers-of-concat_alignmentfna)
- [Create tree file from strain_concat_alignment.fna using IQTree](#create-tree-file-from-strain_concat_alignmentfna-using-iqtree)
- [Swap accessions for abreviated names in concat_alignment.fna.treefile using accession.hash](#swap-accessions-for-abreviated-names-in-concat_alignmentfnatreefile-using-accessionhash)
- [Create tree file from a subset of strains, including wBm, wBp, wCle, wMel, wOo, wOv](#create-tree-file-from-a-subset-of-strains-including-wbm-wbp-wcle-wmel-woo-wov)
	- [Subset records for wBm, wBp, wCle, wMel, wOo, wO](#subset-records-for-wbm-wbp-wcle-wmel-woo-wo)
- [Selection Analysis Directories](#selection-analysis-directories)
- [Run PanOCT](#run-panoct)
	- [Make genomes.list file \(option --genome_list_file\)](#make-genomeslist-file-option---genome_list_file)
	- [Make combined.fasta containing all nt fasta records for coding sequences in all genomes \(option: --combined_fasta\)](#make-combinedfasta-containing-all-nt-fasta-records-for-coding-sequences-in-all-genomes-option---combined_fasta)
		- [Concatenate all coding sequences](#concatenate-all-coding-sequences)
		- [Simplify headers from verbose_combined.fasta so that they only include locus tags and delete empty lines](#simplify-headers-from-verbose_combinedfasta-so-that-they-only-include-locus-tags-and-delete-empty-lines)
	- [Make combined.att containing attributes for all nt fasta records for all coding sequenes in all genomes \(option: --combined.att\). This is necessary because the gbk files can contain locus tags that have no corresponding sequence in the coding seq fasta.](#make-combinedatt-containing-attributes-for-all-nt-fasta-records-for-all-coding-sequenes-in-all-genomes-option---combinedatt-this-is-necessary-because-the-gbk-files-can-contain-locus-tags-that-have-no-corresponding-sequence-in-the-coding-seq-fasta)
		- [Make 1st column of combined.att \(grab only the genome accession without the ".\[0-9\]"\)](#make-1st-column-of-combinedatt-grab-only-the-genome-accession-without-the-0-9)
		- [Make 2nd column of combined.att \(grab only the locus tag\)](#make-2nd-column-of-combinedatt-grab-only-the-locus-tag)
		- [Make 3rd and 4th columns of combined.att \(coordinates for each feature reversing coordinates when genes are on "complementary" strand\)](#make-3rd-and-4th-columns-of-combinedatt-coordinates-for-each-feature-reversing-coordinates-when-genes-are-on-complementary-strand)
		- [Make 5th column of combined.att \(grab the protein annotation\)](#make-5th-column-of-combinedatt-grab-the-protein-annotation)
		- [Make 6th column of combined.att containing the name of the genome \(e.g. wBp\) of the appropriate strain on every line](#make-6th-column-of-combinedatt-containing-the-name-of-the-genome-eg-wbp-of-the-appropriate-strain-on-every-line)
			- [print the number of coading sequences in each strain](#print-the-number-of-coading-sequences-in-each-strain)
			- [make a file that has the name of each strain printed the same number of times as coading sequences it contains \(for each line in the pasted .list and .num files \[separated by a \t which can be specified with IFS=$'\t'\], the command reads field 1 and field 2, prints the value in field 1 the number of times specified in field 2, and appends output to col6.att\)](#make-a-file-that-has-the-name-of-each-strain-printed-the-same-number-of-times-as-coading-sequences-it-contains-for-each-line-in-the-pasted-list-and-num-files-separated-by-a-t-which-can-be-specified-with-ifs%24t-the-command-reads-field-1-and-field-2-prints-the-value-in-field-1-the-number-of-times-specified-in-field-2-and-appends-output-to-col6att)
		- [Paste columns together to creat combined.att](#paste-columns-together-to-creat-combinedatt)
		- [Split combined.att file by genome id in column 6](#split-combinedatt-file-by-genome-id-in-column-6)
	- [Run panoct nt pipeline](#run-panoct-nt-pipeline)
- [Aggregate list of core genes for each strain](#aggregate-list-of-core-genes-for-each-strain)
	- [Run R code to generate a strain-specific list of](#run-r-code-to-generate-a-strain-specific-list-of)
	- [Delete first row of each strain-specific list of core coding sequences \(all files contain x in their first row\)](#delete-first-row-of-each-strain-specific-list-of-core-coding-sequences-all-files-contain-x-in-their-first-row)
	- [Concatenate strain-specific lists of core coding sequences](#concatenate-strain-specific-lists-of-core-coding-sequences)
	- [Pull out sequences corresponding to contents of core_ortholist.txt](#pull-out-sequences-corresponding-to-contents-of-core_ortholisttxt)
	- [reorder all seqs in core_genseqs.fna according to the order in core_ortholist.txt](#reorder-all-seqs-in-core_genseqsfna-according-to-the-order-in-core_ortholisttxt)
		- [run reorder_fasta.py](#run-reorder_fastapy)
	- [Split ordered_core_geneseqs.fna by strain accoding to header information](#split-ordered_core_geneseqsfna-by-strain-accoding-to-header-information)
		- [Create genome look-up table](#create-genome-look-up-table)
		- [rename strain-specific core genome files](#rename-strain-specific-core-genome-files)
	- [Create strain-specific core genome files consisting of concatenated core genes](#create-strain-specific-core-genome-files-consisting-of-concatenated-core-genes)
		- [Add '#' to each record in each core gene fasta](#add--to-each-record-in-each-core-gene-fasta)
		- [Add a fasta header containing the file name to the first line of each .core.fna](#add-a-fasta-header-containing-the-file-name-to-the-first-line-of-each-corefna)
		- [Remove .core.fna extension from the file name header and remove all gene record headers](#remove-corefna-extension-from-the-file-name-header-and-remove-all-gene-record-headers)
		- [Concatenate all core genomes and format the concatenated file so that each line contains 60 nt](#concatenate-all-core-genomes-and-format-the-concatenated-file-so-that-each-line-contains-60-nt)
- [Run translatorx](#run-translatorx)
	- [Run R code](#run-r-code)
				- [PARALOG ANALYSIS](#paralog-analysis)
- [this command iterates over each column, and writes the contents of each column to a new file titled with the header of the column](#this-command-iterates-over-each-column-and-writes-the-contents-of-each-column-to-a-new-file-titled-with-the-header-of-the-column)
	- [delete first row in all .txt files in $Panoct_out/results/Nonparalog_core_clusters](#delete-first-row-in-all-txt-files-in-%24panoct_outresultsnonparalog_core_clusters)
	- [Create multi-fasta files containing the orthologous cluster nt sequences from each genome](#create-multi-fasta-files-containing-the-orthologous-cluster-nt-sequences-from-each-genome)
		- [This loop iterates over all .txt files containing the lists of ortholog locus tags per cluster and extracts the corresponding sequneces from the $Combined_fasta file](#this-loop-iterates-over-all-txt-files-containing-the-lists-of-ortholog-locus-tags-per-cluster-and-extracts-the-corresponding-sequneces-from-the-%24combined_fasta-file)
		- [Remove stop codons from the end of each fasta record in each multifasta file](#remove-stop-codons-from-the-end-of-each-fasta-record-in-each-multifasta-file)
	- [Execute translatorx for each](#execute-translatorx-for-each)
- [Analyze nematode clade](#analyze-nematode-clade)
	- [Run additional R code](#run-additional-r-code)
	- [Create multi-fasta files containing the orthologous cluster nt sequences from each genome](#create-multi-fasta-files-containing-the-orthologous-cluster-nt-sequences-from-each-genome-1)
		- [This loop iterates over all .txt files containing the lists of ortholog locus tags per cluster and extracts the corresponding sequneces from the $Combined_fasta file](#this-loop-iterates-over-all-txt-files-containing-the-lists-of-ortholog-locus-tags-per-cluster-and-extracts-the-corresponding-sequneces-from-the-%24combined_fasta-file-1)
		- [Remove stop codons from the end of each fasta record in each multifasta file](#remove-stop-codons-from-the-end-of-each-fasta-record-in-each-multifasta-file-1)
	- [Execute translatorx for each](#execute-translatorx-for-each-1)
- [Alphanumeritize fasta records](#alphanumeritize-fasta-records)
- [Convert all fasta alignments to phylip alignments](#convert-all-fasta-alignments-to-phylip-alignments)
- [Rename all sequence header labels so that labels match the treefile strains](#rename-all-sequence-header-labels-so-that-labels-match-the-treefile-strains)
- [Run PAML](#run-paml)
- [creat .ctl file in text editor, and upload to $TransX/Phylip_nt/Paml_run](#creat-ctl-file-in-text-editor-and-upload-to-%24transxphylip_ntpaml_run)

<!-- /MarkdownTOC -->
#Phylogeny Directories
Working_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/Mugsy
mkdir $Working_dir/Mugsy_out
Mugsy_out=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/Mugsy/Mugsy_out
Mothur_dir=/usr/local/packages/mothur-1.40.4
IQTree=/local/aberdeen2rw/julie/Jarrett/Programs/IQtree/iqtree-1.6.12-Linux/bin/iqtree

#Test

# Test2

#   T3

#Isolate strain accession list from mugsy generated multiple alignment file (.maf)	
awk '{print $2}' $Mugsy_out/mugsy.maf | grep -v "score" |sort | uniq | sed '/^$/d' | grep -v "ver" | sed 's/.*\.//g' > accession.list 

#this set of commands makes the processed directory, removes any prior contents, and then reads each line of the accesssion.list file and creates separate .fna file for each item in accession.list, and then starts a fasta header within each .fna file that contains the file name. 
mkdir $Mugsy_out/Processed
rm $Mugsy_out/Processed/*

while read line
do
	echo ">"$line"" > "$Mugsy_out"/Processed/"$(echo "$line")".fna
done < "$Mugsy_out"/accession.list

Accession_list="$Mugsy_out"/accession.list


# This command adds a parses the mugsy produced .maf file and adds aligned sequences corresponding to each strain to the .fna files produced in the previous step. 
while read line
do
  marker=$(echo "$line" | awk -F " " '{print $1}')
  if [ "$marker" = "a" ]
  then
    MULT=$(echo "$line" | awk -F " " '{print $4}')
  else
    if [ "$MULT" = "mult="$(wc -l "$Accession_list" | awk '{print $1}')"" ]
    then
      SPECIES=$(echo "$line" | awk -F " " '{print $2}' | sed "s/.*\.//g")
      echo "$line" | awk -F " " '{print $7}' >> "$Mugsy_out"/Processed/"$SPECIES".fna
    fi
  fi
done < "$Mugsy_out"/mugsy.maf
rm "$Mugsy_out"/Processed/*=*


while read line; do	echo "$line" | awk '{sub/$1/,$2)}1' accession.hash;	done < test


#concatenate aligned sequence fastas, change directories to the Mugsy_out directory, and invoke mothur
cat "$Mugsy_out"/Processed/* > "$Mugsy_out"/concat.fasta
cd "$Mugsy_out"
rm mothur*
"$Mothur_dir"/mothur
#Run the following mothur commands 
filter.seqs(fasta=concat.fasta, vertical=F, trump=-)
filter.seqs(fasta=concat.filter.fasta, vertical=F, trump=.)
quit()
#change the name of the mothur output. concat_alignment.fna is the file needed for running iqtree
mv "$Mugsy_out"/concat.filter.filter.fasta "$Mugsy_out"/concat_alignment.fna

#change fasta headers of concat_alignment.fna
seqkit replace -p "(.+)" -r '{kv}' -k "$Mugsy_out"/accession.hash "$Mugsy_out"/concat_alignment.fna > "$Mugsy_out"/strain_concat_alignment.fna 

"

# Create tree file from strain_concat_alignment.fna using IQTree
$IQTree -s "$Mugsy_out"/strain_concat_alignment.fna -nt 4 -bb 1000 -redo

#Swap accessions for abreviated names in concat_alignment.fna.treefile using accession.hash

#Create tree file from a subset of strains, including wBm, wBp, wCle, wMel, wOo, wOv
##Subset records for wBm, wBp, wCle, wMel, wOo, wO
awk 'BEGIN {RS=">"} /wMel|wBm|wBp|wCle|wOo|wOv/ {print ">"$0}' strain_concat_alignment.fna | sed '/^$/d' > fil_wMelOG_concat_alignment.fna
$IQTree -s "$Mugsy_out"/fil_wMelOG_concat_alignment.fna -nt 4 -bb 1000 -redo


#Selection Analysis Directories
Working_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol
Panoct_bin=/local/aberdeen2rw/julie/Jarrett/Programs/PanOCT/PanGenomePipeline-master/pangenome/bin
Gbks=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/GenBank/SimpleExten
Genes_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/Fasta/CodingSeqs/Nucleotides
Prots_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/Fasta/CodingSeqs/AminoAcids
Scripts=/local/aberdeen2rw/julie/Jarrett/Scripts
Panoct_out="$Working_dir"/PanOCTpl_out
Programs=/local/aberdeen2rw/julie/Jarrett/Programs
Paml_dir=/local/aberdeen2rw/julie/Jarrett/Programs/PAML/paml4.9j
IQTree=/local/aberdeen2rw/julie/Jarrett/Programs/IQtree/iqtree-1.6.12-Linux/bin/iqtree
Test=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/Transx_out/test/
TransX=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/No_stops_Fasta_seqs/Transx_out
NoStopsTest=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/No_stops_Fasta_seqs/Transx_out/Test2


#Run PanOCT
cd "$Working_dir"/
rm -r PanOCTpl_out
mkdir PanOCTpl_out
Panoct_out="$Working_dir"/PanOCTpl_out
##Make genomes.list file (option --genome_list_file)
ls -l $Genes_dir | awk '{print $NF}' | sed '1d' | sed 's/.fna//g' > $Panoct_out/genomes.list
##Make combined.fasta containing all nt fasta records for coding sequences in all genomes (option: --combined_fasta)
### Concatenate all coding sequences
cat "$Genes_dir"/* > "$Panoct_out"/verbose_combined.fasta
### Simplify headers from verbose_combined.fasta so that they only include locus tags and delete empty lines
sed 's/^.*\(locus_tag.*\)/>\1/g' "$Panoct_out"/verbose_combined.fasta | sed 's/]//g' | sed 's/locus_tag=//g' | awk '{print $1}' > "$Panoct_out"/combined.fasta
sed -i '/^$/d' "$Panoct_out"/combined.fasta
##Make combined.att containing attributes for all nt fasta records for all coding sequenes in all genomes (option: --combined.att). This is necessary because the gbk files can contain locus tags that have no corresponding sequence in the coding seq fasta.
### Make 1st column of combined.att (grab only the genome accession without the ".[0-9]")
awk -F "|" '{print $2}' "$Panoct_out"/verbose_combined.fasta | sed '/^$/d' | sed 's/\.[0-9].*//g' > "$Panoct_out"/col1.att
### Make 2nd column of combined.att (grab only the locus tag)
grep ">" "$Panoct_out"/verbose_combined.fasta | sed 's/^.*\(locus_tag.*\)/>\1/g' | sed 's/>locus_tag=//g' | awk -F " " '{print $1}' | sed 's/]//g' > "$Panoct_out"/col2.att
### Make 3rd and 4th columns of combined.att (coordinates for each feature reversing coordinates when genes are on "complementary" strand)
grep ">" "$Panoct_out"/verbose_combined.fasta | sed 's/.*\(location=\)//g' | awk '{print $1}' | tr -d "()[]join<>" | sed 's/\.\./\t/g' | awk -F "\t" '{print $1, $NF}' | sed 's/t/\t/g' | awk '{if ($1 == "cmpleme") print $3,$2; else print $1,$2}' | sed 's/ /\t/g' > "$Panoct_out"/col3-4.att 
### Make 5th column of combined.att (grab the protein annotation)
grep ">" "$Panoct_out"/verbose_combined.fasta | sed 's/.*\(protein=\)//g' | sed 's/\] \[/#/g' | awk -F "#" '{print $1}' > "$Panoct_out"/col5.att
### Make 6th column of combined.att containing the name of the genome (e.g. wBp) of the appropriate strain on every line
#### print the number of coading sequences in each strain
for f in "$Genes_dir"/*; do grep ">" $f | wc -l;  done > "$Panoct_out"/genes.num
#### make a file that has the name of each strain printed the same number of times as coading sequences it contains (for each line in the pasted .list and .num files [separated by a \t which can be specified with IFS=$'\t'], the command reads field 1 and field 2, prints the value in field 1 the number of times specified in field 2, and appends output to col6.att)
cd "$Panoct_out"
while read -r f1 f2;
do
	yes "$f1" | head -n "$f2"
done < <(paste genomes.list genes.num) > col6.att

### Paste columns together to creat combined.att
rm "$Panoct_out"/combined.att
cd "$Panoct_out"/
paste col1.att col2.att col3-4.att col5.att col6.att > "$Panoct_out"/combined.att
### Split combined.att file by genome id in column 6
cat "$Panoct_out"/genomes.list | while read line; do cat "$Panoct_out"/combined.att | grep $line > "$Panoct_out"/"$line.natt"; done
##Run panoct nt pipeline
cd "$Panoct_out"/
mkdir fasta_dir att_dir
mv combined.fasta ./fasta_dir
mv *.natt ./att_dir
perl "$Panoct_bin"/run_pangenome.pl -w "$Panoct_out" -g ./genomes.list  --no_grid --blast_local --panoct_local --use_nuc  --panoct_verbose

#Aggregate list of core genes for each strain 
mkdir "$Working_dir"/Core_genome_panoct

##Run R code to generate a strain-specific list of 
R 
setwd("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results")

anno=read.table("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/matchtable.txt", row.names = 1, header = F, sep = "\t")
bin=read.table("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/matchtable_0_1.txt", row.names = 1, header = F, sep = "\t")
anbin=as.data.frame(lapply(anno,gsub,pattern="----------",replacement=0,fixed=TRUE)) # This replaces blank cells in the annotation ortholog table with 0

headers=c("1-wAlbB","2-wAna2.1","3-wAu","4-wBm","5-wBp","6-wBt","7-wCauA","8-wCle","9-wFol","10-wHa","11-wInc_Cu","12-wIrr","13-wMau","14-wMeg","15-wMel","16-wNo","17-wOo","18-wOv","19-wPip","20-wRi","21-wTpre");headers;length(headers)

colnames(anno)=headers
colnames(bin)=headers
colnames(anbin)=headers

lab.core=anbin[rowSums(bin)==21, ] # This generates a table of orthologs that are present in all strains (each row has a value in ever cell)

length(lab.core$`1-wAlbB`)

setwd("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/Core_genome_panoct") # establish the output directory to write strain specific lists to

for (i in colnames(lab.core)) {
  raw_file <- paste0(i, ".txt")
  write.table(lab.core[,i], raw_file, quote = F, row.names = F, col.names = T, append = F)
} #this command iterates over each column, and writes the contents of each column to a new file titled with the header of the column

q()
y

##Delete first row of each strain-specific list of core coding sequences (all files contain x in their first row)
sed -i '1d' "$Working_dir"/Core_genome_panoct/*.txt
##Concatenate strain-specific lists of core coding sequences
cat "$Working_dir"/Core_genome_panoct/*.txt > "$Working_dir"/Core_genome_panoct/core_ortholist.txt
##Pull out sequences corresponding to contents of core_ortholist.txt
seqtk subseq "$Panoct_out"/fasta_dir/combined.fasta "$Working_dir"/Core_genome_panoct/core_ortholist.txt > "$Working_dir"/Core_genome_panoct/core_geneseqs.fna
##reorder all seqs in core_genseqs.fna according to the order in core_ortholist.txt
###run reorder_fasta.py
python "$Scripts_dir"/reorder_fasta.py "$Working_dir"/Core_genome_panoct/core_geneseqs.fna "$Working_dir"/Core_genome_panoct/core_ortholist.txt > ordered_core_geneseqs.fna
## Split ordered_core_geneseqs.fna by strain accoding to header information
cd "$Working_dir"/Core_genome_panoct/
awk '{if(substr($0,1,1) == ">"){split(substr($0,2,length($0)),a,/_/);filename=a[1]};print $0 > filename".core.fna"}' ordered_core_geneseqs.fna
'
### Create genome look-up table
awk -F "\t" -v OFS="\t" '{print $2, $6}' $Panoct_out/combined.att | sed 's/_.*\t/\t/g' | awk -v OFS="\t" '{print $2,$1}' | uniq > "$Working_dir"/Core_genome_panoct/genome.hash
### rename strain-specific core genome files
cd "$Working_dir"/Core_genome_panoct/
while read -r to from; do if [ -e "${from}.core.fna" ]; then mv "${from}.core.fna" "${to}.core.fna"; fi; done < genome.hash

##Create strain-specific core genome files consisting of concatenated core genes
### Add '#' to each record in each core gene fasta
sed -i 's/>/>#/g' *.core.fna
### Add a fasta header containing the file name to the first line of each .core.fna
cd "$Working_dir"/Core_genome_panoct/
for filename in $(ls *.core.fna); do sed "1s/^/>${filename} \n/" ${filename} > $filename.new; mv $filename.new $filename; done
### Remove .core.fna extension from the file name header and remove all gene record headers
cd "$Working_dir"/Core_genome_panoct/
for f in *.core.fna; do grep -v "#" $f | sed s'/.core.fna//g' > "$f.new"; mv $f.new $f; done
### Concatenate all core genomes and format the concatenated file so that each line contains 60 nt
cat *.core.fna > core_genome.fna
fasta_formatter -i core_genome.fna -o core_genome_w60.fna -w 60
mv core_genome_w60.fna core_genome.fna

#Run translatorx
cd $Panoct_out/results
mkdir $Panoct_out/results/Nonparalog_core_clusters
##Run R code
R
anno=read.table("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/matchtable.txt", row.names = 1, header = F, sep = "\t")
bin=read.table("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/matchtable_0_1.txt", row.names = 1, header = F, sep = "\t")
anbin=as.data.frame(lapply(anno,gsub,pattern="----------",replacement=0,fixed=TRUE)) # This replaces blank cells in the annotation ortholog table with 0

headers=c("1-wAlbB","2-wAna2.1","3-wAu","4-wBm","5-wBp","6-wBt","7-wCauA","8-wCle","9-wFol","10-wHa","11-wInc_Cu","12-wIrr","13-wMau","14-wMeg","15-wMel","16-wNo","17-wOo","18-wOv","19-wPip","20-wRi","21-wTpre");headers;length(headers)

colnames(anno)=headers
colnames(bin)=headers
colnames(anbin)=headers

lab.core=anbin[rowSums(bin)==21, ] # This generates a table of orthologs that are present in all strains (each row has a value in ever cell)

length(lab.core$`1-wAlbB`)
table(rowSums(bin)==21)

##### PARALOG ANALYSIS####
paramatch=read.table("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/matchtable_paralogs.txt", header = F, sep = "\t")

nonpara=paramatch[paramatch$V2<2,1]

orthonly=bin[nonpara,]
corthonly=orthonly[rowSums(orthonly)==21,]
length(corthonly$`1-wAlbB`)

corthanno=anbin[as.numeric(row.names(corthonly)),]

clustcol=t(corthanno)

setwd("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/Nonparalog_core_clusters")

for (i in colnames(clustcol)) {
  raw_file <- paste0(i, ".txt")
  write.table(clustcol[,i], raw_file, quote = F, row.names = F, col.names = T, append = F)
}
#this command iterates over each column, and writes the contents of each column to a new file titled with the header of the column

quit()

##delete first row in all .txt files in $Panoct_out/results/Nonparalog_core_clusters
sed -i '1d' $Panoct_out/results/Nonparalog_core_clusters/*.txt
##Create multi-fasta files containing the orthologous cluster nt sequences from each genome 
mkdir $Panoct_out/results/Nonparalog_core_clusters/Fasta_seqs
Combined_fasta=$Panoct_out/fasta_dir/combined.fasta
rm $Panoct_out/results/Nonparalog_core_clusters/Fasta_seqs/*
cd $Panoct_out/results/Nonparalog_core_clusters/
###This loop iterates over all .txt files containing the lists of ortholog locus tags per cluster and extracts the corresponding sequneces from the $Combined_fasta file
for f in ./*.txt;
	do seqtk subseq $Combined_fasta $f > $Panoct_out/results/Nonparalog_core_clusters/Fasta_seqs/$f".fa";
done

###Remove stop codons from the end of each fasta record in each multifasta file
cd $Panoct_out/results/Nonparalog_core_clusters/Fasta_seqs/
mkdir No_stops_Fasta_seqs
cp *.fa No_stops_Fasta_seqs
cd No_stops_Fasta_seqs
for f in ./*.fa; do sed 's/TAG$//g' $f | sed 's/TAA$//g' | sed 's/TGA$//g' > "$f.new"; mv $f.new $f; done 

##Execute translatorx for each 
$Programs/translatorx.pl -i 1.txt.fa -c 11

cd $Panoct_out/results/Nonparalog_core_clusters/Fasta_seqs
mkdir Transx_out
rm Transx_out/*
for f in ./*.fa;
	do $Programs/translatorx.pl -i $f -c 11 -o Transx_out/$f;
done

#Analyze nematode clade
cd $Panoct_out/results
mkdir $Panoct_out/results/Nonparalog_nemacore_clusters

##Run additional R code
R
nemabin=bin[,c(4,5,8,15,17,18)]#grab specific columns from the binary match table
nemaorth=nemabin[nonpara,]#Exclude paralogs
nema.core.nopara=nemaorth[rowSums(nemaorth)==6,]#include only orthologs common to all members in the group
nema.core.anbin=anbin[as.numeric(row.names(nema.core.nopara)),c(4,5,8,15,17,18)]#extract the group-core ortholog tags
nemaclustcol=t(nema.core.anbin)#transpose core ortholog table so that each column is a cluster

setwd("V:/readandwrite/julie/Jarrett/Ref_genomes/Wolbachia/PanWol/PanOCTpl_out/results/Nonparalog_nemacore_clusters")
for (i in colnames(nemaclustcol)) {
  raw_file <- paste0(i,".txt")
  write.table(nemaclustcol[,i], raw_file, quote = F, row.names = F, col.names = F, append = F)
}#write each column to a text file 
quit()

##Create multi-fasta files containing the orthologous cluster nt sequences from each genome 
mkdir $Panoct_out/results/Nonparalog_nemacore_clusters/Fasta_seqs
Combined_fasta=$Panoct_out/fasta_dir/combined.fasta
rm $Panoct_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/*
cd $Panoct_out/results/Nonparalog_nemacore_clusters/
###This loop iterates over all .txt files containing the lists of ortholog locus tags per cluster and extracts the corresponding sequneces from the $Combined_fasta file
for f in ./*.txt;
	do seqtk subseq $Combined_fasta $f > $Panoct_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/$f".fa";
done

###Remove stop codons from the end of each fasta record in each multifasta file
cd $Panoct_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/
mkdir No_stops_Fasta_seqs
cp *.fa No_stops_Fasta_seqs
cd No_stops_Fasta_seqs
for f in ./*.fa; do sed 's/TAG$//g' $f | sed 's/TAA$//g' | sed 's/TGA$//g' > "$f.new"; mv $f.new $f; done 

##Execute translatorx for each 
$Programs/translatorx.pl -i 1.txt.fa -c 11

cd $Panoct_out/results/Nonparalog_nemacore_clusters/Fasta_seqs/No_stops_Fasta_seqs
mkdir Transx_out
rm Transx_out/*
for f in ./*.fa;
	do $Programs/translatorx.pl -i $f -c 11 -o Transx_out/$f;
done

#Alphanumeritize fasta records
cd $TransX
for f in ./*.txt.fa.nt_ali.fasta; do awk '/^>/{key=$1} {print key, $0}' $f | sort -k1,1 -s | cut -d' ' -f2- > "$f.new"; mv $f.new $f; done 

'
#Convert all fasta alignments to phylip alignments
rm -r $TransX/Phylip_nt #remove Phylip_nt
mkdir $TransX/Phylip_nt #make Phylip_nt
for f in "$TransX"/*nt_ali.fasta;
do $Programs/Fasta2Phylip.pl $f $f".phy";
done #convert all nt fasta alinments to phylip format
mv $TransX/*.phy $TransX/Phylip_nt/ # move phylip alignments to the appropriate folder
cd $TransX/Phylip_nt/ #change directories to Phylip_nt/
rename txt.fa.nt_ali.fasta.phy nt_ali.phy.locus *.txt.fa.nt_ali.fasta.phy #change file extension of all files to nt_ali.phy

#Rename all sequence header labels so that labels match the treefile strains
rm $TransX/Phylip_nt/*.tips
for f in $TransX/Phylip_nt/*.phy.locus;
do sed 's/DEI.*\t/wAlbB  \n/g' $f | sed 's/EA6.*\t/wAna2_1  \n/g' | sed 's/WPA.*\t/wAu  \n/g' | sed 's/WBM.*\t/wBm  \n/g' | sed 's/WBP.*\t/wBp  \n/g' | sed 's/BBB.*\t/wBt  \n/g' | sed 's/wCau.*\t/wCauA  \n/g' | sed 's/WCLE.*\t/wCle  \n/g' | sed 's/ASM.*\t/wFol  \n/g' | sed 's/WHA.*\t/wHa  \n/g' | sed 's/WG6.*\t/wInc_Cu  \n/g' | sed 's/E04.*\t/wIrr  \n/g' | sed 's/EJB.*\t/wMau  \n/g' | sed 's/CAI.*\t/wMeg  \n/g' | sed 's/WD_.*\t/wMel  \n/g' | sed 's/WNO.*\t/wNo  \n/g' | sed 's/WOO.*\t/wOo  \n/g' | sed 's/WOV.*\t/wOv  \n/g' | sed 's/WP_.*\t/wPip  \n/g' | sed 's/WRI.*\t/wRi  \n/g' | sed 's/wTp.*\t/wTpre  \n/g' > $f".tips";
done
cd $TransX/Phylip_nt/
rename nt_ali.phy.locus.tips nt_ali.phy.tips *.nt_ali.phy.locus.tips 
sed -i 's/\t/\n/g' *.locus

#Run PAML
rm -r $TransX/Phylip_nt/Paml_run
mkdir $TransX/Phylip_nt/Paml_run
cd $TransX/Phylip_nt/Paml_run
cp $Mugsy_out/fil_wMelOG_concat_alignment.fna.treefile .
sed 's/wMel/4/g' *.treefile | sed 's/wBm/1/g' | sed 's/wBp/2/g' | sed 's/wOo/5/g' | sed 's/wOv/6/g' | sed 's/wCle/3/g' > num.tree
cat $TransX/Phylip_nt/*.locus > all.loci
#creat .ctl file in text editor, and upload to $TransX/Phylip_nt/Paml_run
$Paml_dir/bin/codeml mods12_all_loci.ctl #Report run time here:

