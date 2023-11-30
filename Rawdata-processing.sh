## data clean
fastp -i ./$x.R1.fq.gz -I ./$x.R2.fq.gz -o ../cleandata/$x.r1.fq.gz -O ../cleandata/$x.r2.fq.gz --detect_adapter_for_pe -w 4 --compression 9 -h ../cleandata/$x.html -j ../cleandata/$x.json
## mapping & deduplication
bwa mem -M -t 4 /path/reference/cs/iwgsc1/bwa_index/cs ../cleandata/$x.r1.fq.gz ../cleandata/$x.r2.fq.gz > ../bwa/$x.sam &&
samtools view -bS -F 1804 -f 2 -q 30 ../bwa/$x.sam | samtools sort - | samtools rmdup -s - ../bwa/$x.raw.bam && java -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -Xmx8000M -jar ~/bin/picard.jar MarkDuplicates I=../bwa/$x.raw.bam O=../bwa/$x.rmdup.bam M=../bwa/$x.txt REMOVE_DUPLICATES=true
## peak calling
macs2 callpeak -t $x.bam -p 1e-3 -f BAM --keep-dup all -g 14600000000 -n ../narrowpeak/$x
macs2 callpeak -t $x.bam -f BAM --keep-dup all -g 14600000000 --broad --broad-cutoff 0.05 -n ../broadpeak/$x
