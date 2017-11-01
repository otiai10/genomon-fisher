#!/bin/bash
set -eux

# Validate file existence
ls /var/refs/${REFERENCE}
ls /var/data/in/${TUMOR_READ_PAIR}
ls /var/data/in/${CONTROL_READ_PAIR}

# Decompress Tumor files
mkdir /var/data/in/tumor
tar -zxvf /var/data/in/${TUMOR_READ_PAIR} -C /var/data/in/tumor
DIR_TUMOR="/var/data/in/tumor/`ls /var/data/in/tumor`"
TUMORS=(); i=0; for f in `ls ${DIR_TUMOR} | sort`; do
    TUMORS[$i]="${DIR_TUMOR}/$f"
    i=`expr $i + 1`
done

# Decompress Control files
mkdir /var/data/in/control
tar -zxvf /var/data/in/${CONTROL_READ_PAIR} -C /var/data/in/control
DIR_CONTROL="/var/data/in/control/`ls /var/data/in/control`"
CONTROLS=(); i=0; for f in `ls ${DIR_CONTROL} | sort`; do
    CONTROLS[$i]="${DIR_CONTROL}/$f"
    i=`expr $i + 1`
done

# Align provided TUMOR reads
/bin/bwa mem /var/refs/${REFERENCE} ${TUMORS[0]} ${TUMORS[1]} > /var/data/out/result.tumor.sam
# Bamify result TUMOR sam file
/bin/samtools view -S -b -h /var/data/out/result.tumor.sam > /var/data/out/result.tumor.bam
# Sort TUMOR bam
/bin/samtools sort /var/data/out/result.tumor.bam > /var/data/out/result.tumor.sorted.bam

# Align provided CONTROL reads
/bin/bwa mem /var/refs/${REFERENCE} ${CONTROLS[0]} ${CONTROLS[1]} > /var/data/out/result.control.sam
# Bamify result CONTROL sam file
/bin/samtools view -S -b -h /var/data/out/result.control.sam > /var/data/out/result.control.bam
# Sort CONTROL bam
/bin/samtools sort /var/data/out/result.control.bam > /var/data/out/result.control.sorted.bam

# Execute Fisher's test
/usr/bin/fisher comparison \
  -o /var/data/out/result.fisher.txt \
  --ref_fa /var/refs/$REFERENCE \
  -1 /var/data/result.tumor.sorted.bam \
  -2 /var/data/result.control.sorted.bam \
  --samtools_path /bin/samtools \
  --min_depth 8 \
  --base_quality 15 \
  --min_variant_read 4 \
  --max_allele_freq 0.1 \
  --fisher_value 0.1 \
  --samtools_params "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"

echo "[genomon-fisher] Congratulations! Everything has been completed successfully. Bye, see you again!"
