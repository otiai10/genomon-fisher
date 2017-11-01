#!/bin/bash

echo "[genomon-fisher][0] Checking input parameters (mostly provided as env variables)"
if [ ! -f "/var/refs/${REFERENCE}" ]; then
    echo "[genomon-fisher][ERROR] REFERENCE file not found: specified with '${REFERENCE}', searched at /var/refs"
    exit 1
fi
if [ ! -f "/var/data/in/${TUMOR_READ_PAIR}" ]; then
    echo "[genomon-fisher][ERROR] TUMOR_READ_PAIR file not found: specified with '${TUMOR_READ_PAIR}', searched at /var/data/in"
    exit 1
fi
if [ ! -f "/var/data/in/${CONTROL_READ_PAIR}" ]; then
    echo "[genomon-fisher][ERROR] CONTROL_READ_PAIR file not found: specified with '${CONTROL_READ_PAIR}', searched at /var/data/in"
    exit 1
fi

echo "[genomon-fisher][1] >>> Decompress Tumor file"
mkdir /var/data/in/tumor
tar -zxvf /var/data/in/${TUMOR_READ_PAIR} -C /var/data/in/tumor || exit $?
DIR_TUMOR="/var/data/in/tumor/`ls /var/data/in/tumor`"
TUMORS=(); i=0; for f in `ls ${DIR_TUMOR} | sort`; do
    TUMORS[$i]="${DIR_TUMOR}/$f"
    i=`expr $i + 1`
done

echo "[genomon-fisher][2] >>> Decompress Control file"
mkdir /var/data/in/control
tar -zxvf /var/data/in/${CONTROL_READ_PAIR} -C /var/data/in/control || exit $?
DIR_CONTROL="/var/data/in/control/`ls /var/data/in/control`"
CONTROLS=(); i=0; for f in `ls ${DIR_CONTROL} | sort`; do
    CONTROLS[$i]="${DIR_CONTROL}/$f"
    i=`expr $i + 1`
done

echo "[genomon-fisher][3] >>> Align provided TUMOR reads: '${TUMORS[0]}', '${TUMORS[1]}'"
/bin/bwa mem /var/refs/${REFERENCE} ${TUMORS[0]} ${TUMORS[1]}> /var/data/out/result.tumor.sam
echo "[genomon-fisher][4] >>> Bamify result TUMOR sam file: 'result.tumor.sam'"
/bin/samtools view -S -b -h /var/data/out/result.tumor.sam > /var/data/out/result.tumor.bam
echo "[genomon-fisher][5] >>> Sort TUMOR bam: 'result.tumor.bam'"
/bin/samtools sort /var/data/out/result.tumor.bam > /var/data/out/result.tumor.sorted.bam

echo "[genomon-fisher][6] >>> Align provided CONTROL reads: '${CONTROLS[0]}', '${CONTROLS[1]}'"
/bin/bwa mem /var/refs/${REFERENCE} ${CONTROLS[0]} ${CONTROLS[1]}> /var/data/out/result.control.sam
echo "[genomon-fisher][7] >>> Bamify result CONTROL sam file: 'result.control.sam'"
/bin/samtools view -S -b -h /var/data/out/result.control.sam > /var/data/out/result.control.bam
echo "[genomon-fisher][8] >>> Sort CONTROL bam: 'result.control.bam'"
/bin/samtools sort /var/data/out/result.control.bam > /var/data/out/result.control.sorted.bam


echo "[genomon-fisher][9] >>> Execute Fisher's test: 'result.tumor.sorted.bam' / 'result.control.sorted.bam'"
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
  --samtools_params "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP" \
  || exit $?

echo "[genomon-fisher][6] >>> Congratulations! Everything has been completed successfully. Bye, see you again!"
