#!/bin/bash


if [ "$#" -ne 2 ]; then
    echo "Использование: $0 <имя_файла без расширения> <референсный_геном>"
    exit 1
fi

fastqc "$1.fastq.gz"

./minimap2/minimap2 -a -x map-pb "$2" "$1".fastq > "$1".sam
if [ $? -ne 0 ]; then
    echo "Ошибка при выполнении Minimap2."
    exit 1
fi
echo "Minimap finished successfully."


samtools view -b -o "$1.bam" "$1.sam"
if [ $? -ne 0 ]; then
    echo "Ошибка при выполнении Samtools (конвертация)."
    exit 1
fi
echo "Samtools (конвертация) finished successfully."


samtools flagstat "$1.bam" > "$1.txt"
if [ $? -ne 0 ]; then
    echo "Ошибка при выполнении Samtools (статистика)."
    exit 1
fi
echo "Samtools (статистика) finished successfully."


mapped=$(grep -E "^([[:digit:]]+ \+ [[:digit:]]+ mapped \([[:digit:]]{,2}\.[[:digit:]]{,2}%)" "$1.txt" | awk '{print $5}' | tr -d "(%")

if (( $(echo "$mapped < 90" | bc) )); then
    echo "Mapping quality is below 90%. Aborting pipeline."
    exit 1
else
    echo "Mapping quality is sufficient."
fi


samtools sort -o "$1.sorted.bam" "$1.bam"
if [ $? -ne 0 ]; then
    echo "Ошибка при выполнении Samtools (сортировка)."
    exit 1
fi


freebayes -f "$2" "$1.sorted.bam" > "$1.vcf"
if [ $? -ne 0 ]; then
    echo "Ошибка при выполнении FreeBayes."
    exit 1
fi

echo "Pipeline finished successfully. Output VCF file: $1.vcf"

