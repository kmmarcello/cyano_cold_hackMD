# project outline
**note:** change it into separating by language type: R, Bash, Python, etc.
:::info
**local**
1. GTDB found cyanobacteria with this criteria: ![Screenshot 2023-12-05 at 11.29.09 AM](https://hackmd.io/_uploads/ryfG9epHT.png)
2. downloaded the genomes and .csv file with the criteria to local
3. Opened the .sh folder provided into BBEdit (contains all lines of code to use NCBI datasets to download the genomes) `CTL + a` `CTL + c` 
:::

:::success
`#!/bin/bash`
4. `CTL + v` into FARM cyanobacteria/ncbi_datasets `ENTER`
5. this is where the info is coming [here](https://hackmd.io/SKhz08m-SBOzyE4ksgw1Bw?both#May-15-2023-Monday)
6. pull out actual genomes with for loop
```
#!/bin/bash

for i in *.zip; do
    mkdir "$i-dir"
    cd "$i-dir"
    unzip "../$i"
    for j in *; do
        mv "$j" "$i.${j##*.}"
    done
    cd ..
done
```
7. `mkdir zipped` 
8. `mkdir unzipped`
9. move files to **zipped** `mv *.zip-dir unzipped`
10. `mkdir cyano_fna`
11. `mv cyano_fna unzipped`
12. `find . -maxdepth 5 -name '*.fna' -exec mv {} ./cyano_fna/ \;` 614 genomes moved to **cyano_fna**
13. `mv unzipped/cyano_fna/ .` moved to ~/cyanobacteria/ncbi_datasets
14. move Snakefile, hmms, and create file_names.txt to directory
15. edit Snakefile to match directories
16. `snakemake -n` check
17. `snakemake`
18. `srun -p bmm -J anvio -t 5:00:00 --mem=32G --pty bash`
19. `conda activate anvio-7`
20. deleted these bc they were not good ssequences: GCA_003326215.1_ASM332621v1_genomic, GCA_003326195.1_ASM332619v1_genomic, GCA_002796835.1_ASM279683v1_genomic, GCA_003724315.1_ASM372431v1_genomic, GCA_022241785.1_ASM2224178v1_genomic
21. after deleting each one `touch /home/kmrcello/cyanobacteria/ncbi_datasets/data/*` the data and rerun `snakemake`: snakemake file will get cold tolerance genes. need to get Bacteria 71, photo genes, and TCS genes.
22. run photo_hmms and membrane_hmms

**photo_hmms**
`scp -r ./photo_hmms kmrcello@farm.cse.ucdavis.edu:/home/kmrcello/cyanobacteria/ncbi_datasets`
`srun -p bmm -J anvio -t 10:00:00 --mem=32G --pty bash`
`conda activate anvio-7`
```
#!/bin/bash

for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
    anvi-run-hmms -c $i.db -H /home/kmrcello/cyanobacteria/ncbi_datasets/photo_hmms
done
```

**membrane_hmms**
```
#!/bin/bash

for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
    anvi-run-hmms -c $i.db -H /home/kmrcello/cyanobacteria/ncbi_datasets/membrane_hmms
done
```

**bacteria_71**
```
#!/bin/bash

for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
    anvi-run-hmms -c $i.db -I Bacteria_71
done
```
23. Get sequences for each set of genes:
```
#!/bin/bash

# List of objects to process
# this list of objects is actually the cold tolerance genes
object_list=("COG2609.faa.final_tree.fa" "COG0508.faa.final_tree.fa" "COG1278.faa.final_tree.fa" "COG0513.faa.final_tree.fa" "COG3239.faa.final_tree.fa" "COG0593.faa.final_tree.fa" "COG0484.faa.final_tree.fa" "COG0443.faa.final_tree.fa" "COG0188.faa.final_tree.fa" "COG0776.faa.final_tree.fa" "COG0361.faa.final_tree.fa" "COG0532.faa.final_tree.fa" "COG0290.faa.final_tree.fa" "COG0195.faa.final_tree.fa" "COG0380.faa.final_tree.fa" "COG1185.faa.final_tree.fa" "COG0858.faa.final_tree.fa" "COG0468.faa.final_tree.fa" "COG0557.faa.final_tree.fa" "COG0544.faa.final_tree.fa" "COG1115.faa.final_tree.fa")
```
**example**

```
# Loop through each object
for object in "${object_list[@]}"; do
    echo "Processing object: $object"
    
    # Loop through each .db file
    for file in *.db; do
        echo "Processing file: $file"
        
        # Perform some action for each object and file combination
        # For example, here's a placeholder command to show how you might use 'object' and 'file':
        # echo "Processing object $object with file $file"
        # You can replace this echo command with your specific processing command
        # For example, run a command using 'object' and 'file':
        # some_command --object="$object" --file="$file" > "${object}_${file}_output.txt"
    done
done
```

**real deal**
```
# Loop through each object
for object in "${object_list[@]}"; do
    echo "Processing object: $object"
    
    # Loop through each .db file
    for file in *.db; do
        anvi-get-sequences-for-hmm-hits -c $file \
                                --hmm-source hmms \
                                --gene-names $object \
                                -o ${object}-fastas/cold/${object}-${file}.fa;
    done
done

# Loop through each object
for object in "${photo_genes[@]}"; do
    echo "Processing object: $object"
    
    # Loop through each .db file
    for file in *.db; do
        anvi-get-sequences-for-hmm-hits -c $file \
                                --hmm-source photo_hmms \
                                --gene-names $object \
                                -o ${object}-fastas/photo/${object}-${file}.fa;
    done
done

# Loop through each object
for object in "${membrane_genes[@]}"; do
    echo "Processing object: $object"
    
    # Loop through each .db file
    for file in *.db; do
        anvi-get-sequences-for-hmm-hits -c $file \
                                --hmm-source membrane_hmms \
                                --gene-names $object \
                                -o ${object}-fastas/membrane/${object}-${file}.fa;
    done
done
```

```

# made these bc there was an error but they should be made first
for object in "${object_list[@]}"; do
    mkdir gene_fastas/cold/${object}-fastas
done

for object in "${photo_genes[@]}"; do
    mkdir gene_fastas/photo/${object}-fastas
done

for object in "${membrane_genes[@]}"; do
    mkdir gene_fastas/membrane/${object}-fastas
done


```


```
#!/bin/bash

# List of objects to process

photo_genes=("COG1143.faa.final_tree.fa" "1G6I8.faa.final_tree.fa" "1G1ET.faa.final_tree.fa" "1FZXJ.faa.final_tree.fa" "2ZBN6.faa.final_tree.fa" "3313D.faa.final_tree.fa" "2Z87P.faa.final_tree.fa" "1G08A.faa.final_tree.fa" "2Z8JK.faa.final_tree.fa" "2Z7TN.faa.final_tree.fa" "2Z7VA.faa.final_tree.fa" "2Z7ZP.faa.final_tree.fa")

membrane_genes=("1TVTF.faa.final_tree.fa" "COG4585.faa.final_tree.fa" "COG3239.faa.final_tree.fa" "1G096.faa.final_tree.fa" "1FZVK.faa.final_tree.fa" "COG1398.faa.final_tree.fa" "1G100.faa.final_tree.fa" "1G2GY.faa.final_tree.fa" "COG5002.faa.final_tree.fa" "1FZWA.faa.final_tree.fa" "1G0YA.faa.final_tree.fa")

```
* note on how to look at these lists again:
```
# Access individual elements by index
echo "${cold_genes[0]}"  # Outputs the first element
echo "${cold_genes[1]}"  # Outputs the second element

# Print all elements using a loop
for gene in "${cold_genes[@]}"; do
    echo "$gene"
done

# Verify the entire array content
echo "${cold_genes[@]}"
```

24. change names (this is hard bc the names are different than they were the last time I did this)
```
for object in "${object_list[@]}"; do
	for i in *.fa; do 
	    [ -f "$i" ] || continue
	    mv "$i" "${i//.db/}"
	done
done
```
```
for f in *.fa; do sed -i "s/^>/>${f}_/" "$f"; done
```




25. concatenate the fastas from each genome into gene.fa
```
cat *.fa > PsbA-nuc.fa
```

26. scp them to local
```
scp kmrcello@farm.cse.ucdavis.edu:/home/kmrcello/cyanobacteria/new_wgs_data/not_duplicates/db/PsbA-nuc.fa/PsbA-nuc.fa .
```
:::












New snakefile
:::danger
This is a work in progress bc I'm learning. don't actually use this!!!!
```
FILENAMES ="file_names.txt"
#TEST_FILENAMES ="subset_file_names.txt"
#TEST_SAMPLE_LST= [x.strip().split(".fna")[0] for x in open(FILENAMES, 'r')]

EXTERNAL_GENOMES ="external_genomes.txt"    # no idea if this is ok

SAMPLE_LST= [x.strip().split(".fna")[0] for x in open(FILENAMES, 'r')]
#TEST_SAMPLE_LST= SAMPLE_LST[3:5]

HMM_DIR_NAME="hmms"
PHOTO_HMM_DIR_NAME="photo_hmms"
MEMBRANE_HMM_DIR_NAME="membrane_hmms"

rule all:
    input:
        #reformat_fasta
        #    expand("{sample_i}.fa", sample_i=SAMPLE_LST),
        #        tsv is key file for contigs
        #    expand("{sample_i}.tsv", sample_i=SAMPLE_LST)
        #    expand("{sample_i}.db", sample_i=SAMPLE_LST)
        #    expand("{sample_i}.stats.txt", sample_i=SAMPLE_LST)
        #    expand("{sample_i}.stats.photo.txt", sample_i=SAMPLE_LST)
        #    expand("{sample_i}.stats.membrane.txt", sample_i=SAMPLE_LST)
        #    expand("reports/seq_report/{sample_i}.hmm.sequence.txt", sample_i=SAMPLE_LST)
        #    expand("reports/seq_report/{sample_i}.photo.sequence.txt", sample_i=SAMPLE_LST)
        #    expand("reports/seq_report/{sample_i}.membrane.sequence.txt", sample_i=SAMPLE_LST)
        #    expand("")
        
    

rule reformat_fasta:
    input: 
        "mags_fna/{sample_i}.fna"
    output:
        seq="reformatted_mags/{sample_i}.fa",
        rpt="reports/reformat_reports/{sample_i}.tsv"
    shell: 
        """
        anvi-script-reformat-fasta {input} -o {output.seq}\
        -l 0 --simplify-names --report-file {output.rpt}
        """

rule create_contigs_db:
    input:
        "reformatted_mags/{sample_i}.fa"
    output:
        "db/{sample_i}.db"
    threads: 4
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output}\
        -T {threads} --project-name ctg_cyanos
        """


rule run_hmms:
    input:
        "db/{sample_i}.db" 
    output:
        "reports/hmm_hits/{sample_i}.stats.txt"
    shell:
        """
        anvi-run-hmms -c {input} --hmm-profile-dir {HMM_DIR_NAME}
        anvi-display-contigs-stats {input} --report-as-text -o {output}
        """
        
rule run_photo_hmms:
    input:
        "db/{sample_i}.db"
    output:
        "reports/hmm_hits/{sample_i}.stats.photo.txt"
    shell:
        """
        anvi-run-hmms -c {input} --hmm-profile-dir {PHOTO_HMM_DIR_NAME}
        anvi-display-contigs-stats {input} --report-as-text -o {output}
        """
        
rule run_membrane_hmms:
    input:
        "db/{sample_i}.db"
    output:
        "reports/hmm_hits/{sample_i}.stats.membrane.txt"
    shell:
        """
        anvi-run-hmms -c {input} --hmm-profile-dir {MEMBRANE_HMM_DIR_NAME}
        anvi-display-contigs-stats {input} --report-as-text -o {output}
        """

rule hmm_sequence_hits:
    input:
        db="db/{sample_i}.db",
        
        #we dont actually use these anymore, they are a placeholder for snakemake while .dbs are updated
        rpt="reports/hmm_hits/{sample_i}.stats.txt"
    output:
        "reports/seq_report/{sample_i}.hmm.sequence.txt"
    shell:
        """
        anvi-get-sequences-for-hmm-hits -c {input.db} --hmm-source {HMM_DIR_NAME} -o {output}
        """
        
rule photo_sequence_hits:
    input:
        db="db/{sample_i}.db",
        
        #we dont actually use these anymore, they are a placeholder for snakemake while .dbs are updated
        rpt="reports/hmm_hits/{sample_i}.stats.txt"
    output:
        "reports/seq_report/{sample_i}.photo.sequence.txt"
    shell:
        """
        anvi-get-sequences-for-hmm-hits -c {input.db} --hmm-source {PHOTO_HMM_DIR_NAME} -o {output}
        """
        
rule membrane_sequence_hits:
    input:
        db="db/{sample_i}.db",
        
        #we dont actually use these anymore, they are a placeholder for snakemake while .dbs are updated
        rpt="reports/hmm_hits/{sample_i}.stats.txt"
    output:
        "reports/seq_report/{sample_i}.membrane.sequence.txt"
    shell:
        """
        anvi-get-sequences-for-hmm-hits -c {input.db} --hmm-source {MEMBRANE_HMM_DIR_NAME} -o {output}
        """

rule hmm_matrix:
    input:
        "./{list_i}.txt"
    output:
    shell:
        """
        anvi-script-gen-hmm-hits-matrix-across-genomes -e external-genomes-filamentous-names-v2.txt \
                                               -o {output}.txt \
                                               --hmm-source {MEMBRANE_HMM_DIR_NAME}
        """
        

```
`# touch /home/kmrcello/cyanobacteria/cyano_MAGs/mags_fna/*`
:::



