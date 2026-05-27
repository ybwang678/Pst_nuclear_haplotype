
# rnaseq mapping
for sample in `cat $3 | head -n $1 | tail -n $2`
do
    name=$(echo $sample | sed 's/\r//')
    /dt1/share/software/00.pub_bin/STAR --runThreadN 10 --genomeDir /dt1/ybwang/data/raw_data/ref_genome/pst_134_pri --readFilesIn /dt2/share/data/004.pst_rnaseq_ybwang_231028/$name/${name}_1.fastq.gz /dt2/share/data/004.pst_rnaseq_ybwang_231028/$name/${name}_2.fastq.gz  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000 --outFileNamePrefix /dt2/share/data/020.pst_bam_ybwang_240429/pst134pri_rna/${name}/rust --outSAMmapqUnique 60  --readFilesCommand zcat
    echo "$name mapping finish"
    echo `date`
done

# WGS data clean and snp calling
for sample in `cat $3 | head -n $1 | tail -n $2`
do
        name=${sample##*\/}
        fastp -i $4/${name}.R1.fq.gz -o $4/${name}_clean.R1.fq.gz -I $4/${name}.R2.fq.gz -O $4/${name}_clean.R2.fq.gz
        bwa mem -M -t 60 $5 $4/${name}/${name}_1.fq.gz $4/${name}/${name}_2.fq.gz > $6/${name}.sam
        /dt1/share/software/samtools-1.18/samtools view -bS  $6/${name}.sam -o $6/$name.bam
        rm $6/$name.sam
        java -jar /dt1/ybwang/software/picard.jar  SortSam INPUT=$6/$name.bam OUTPUT=$6/$name.sorted.bam SORT_ORDER=coordinate
        java -jar /dt1/ybwang/software/picard.jar MarkDuplicates REMOVE_DUPLICATES= true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=$6/${name}.sorted.bam OUTPUT=$6/$name.repeatmark.bam METRICS_FILE=$6/$name.bam.metrics
        java -jar /dt1/ybwang/software/picard.jar AddOrReplaceReadGroups     -I $6/$name.repeatmark.bam     -O $6/$name.AddOrReplaceReadGroups.bam     --RGID $name     --RGLB $name     --RGPL illumina     --RGPU machine     --RGSM $name
        /dt1/share/software/samtools-1.18/samtools index $6/$name.AddOrReplaceReadGroups.bam
done
nohup /dt1/ybwang/software/freebayes-master/scripts/freebayes-parallel <(/dt1/ybwang/software/freebayes-master/scripts/fasta_generate_regions.py /dt1/ybwang/data/raw_data/ref_genome/pst_134_pri/GCF_021901695.1_Pst134E36_v1_pri_genomic.fna.fai 100000) 16 
    --use-best-n-alleles 4  --ploidy 2 -f /dt1/ybwang/ref_genome/pst_134_pri/GCF_021901695.1_Pst134E36_v1_pri_genomic.fna bam.list > pst_fb.vcf &


# mash
for sample in `cat $3 | head -n $1 | tail -n $2`
do
    name=${sample##*\/}
    name=$(echo $sample | sed 's/\r//')
    mash screen /dt1/ybwang/data/raw_data/ref_genome/pst104E_tam2025/pst104E.msh /dt2/share/data/005.pst_wgs_ybwang_231028/$name/${name}_1_clean.fastq.gz > /dt3/project/006.pst_revision_ybwang_20250927/004.info/mash_kmer/${name}_pst104E.tab
done

#vcf filt
bcftools view -i 'QUAL > 20 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & AC > 0' all_sample395_fb_2.vcf -o /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt.vcf #初步过滤freebayes call出的vcf
bcftools view -v snps /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt.vcf > /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp.vcf #选取snp
vcftools --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp.vcf --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --recode --recode-INFO-all --out /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp #再次过滤
vcftools --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp.vcf --thin 0.1 --recode --recode-INFO-all --out /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp #只选取10%位点进行群体分析
vcftools --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename.recode.vcf --extract-FORMAT-info GT --out /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename.recode #形成gt_data以便进行pca分析
vcftools --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename.recode.vcf --thin 0.1 --recode --recode-INFO-all --out /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename_thin0.1 #只选取10%位点进行群体分析
plink --bfile /dt1/ybwang/project/pst395_fb/pst395_fb_bed/all_sample395_filt_thin0.1-temporary --recode vcf --out /dt1/ybwang/project/pst395_fb/all_sample395_filt_thin0.1

#admixture
plink --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename.recode.vcf --make-bed --out /dt1/ybwang/project/pst395_fb/pst395_fb_bed/all_sample395_filt_thin0.1 --thin 0.1 --allow-extra-chr
plink --allow-extra-chr --noweb -file all_sample395_filt_chr1 --geno 0.05 --maf 0.05 --hwe 0.0001 --make-bed --out all_sample395_filt_chr1
nohup bash /dt1/ybwang/shell_script/adm.sh all_sample395_filt_thin0.1-temporary.bed &

#iqtree
python /dt1/ybwang/software/vcf2phylip-master/vcf2phylip.py --input /dt1/ybwang/project/pst395_fb/all_sample395_filt_thin0.1.vcf
nohup iqtree -s all_sample395_filt_thin0.1.min4.phy -m GTR+F+G4+ASC -fast -nt 50 -st DNA &

#fis（het）
vcftools --vcf /dt1/ybwang/project/pst395_fb/all_sample395_fb_filt_snp_rename.recode.vcf --het --out all_sample395_fb_filt_snp_rename

#pca
plink --allow-extra-chr --threads 20 -bfile /dt1/ybwang/project/bgt_fb/001.vcf/plink_vcf/bgt_all1086Iso_with_mis_forpopgenetics_LDfilt --pca 20 --out /dt1/ybwang/project/bgt_fb/002.info/bgt1086_pca

#fst、pi
vcftools --vcf  all_sample395_fb_chr1.vcf  --keep    pop_china.txt   --window-pi   50000  --window-pi-step  25000   --out out.txt #pi
vcftools --vcf all_sample395_fb_chr1.vcf --weir-fst-pop  pop_china.txt --weir-fst-pop pop_na.txt --fst-window-size 50000 --fst-window-step 25000 --out out_china_na.txt #fst
