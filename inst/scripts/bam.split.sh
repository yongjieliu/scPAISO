#!/bin/bash

if [ $# -eq 4 ]; then
  R1Bam=$1
  R2Bam=$2
  prefix=$3
  outpath=$4

  mkdir -p $outpath/chr.bam.R2/$prefix
  mkdir -p $outpath/chr.bam.R1/$prefix

  if [ ! -f ${R1Bam}.bai ];then
    samtools index $R1Bam
  fi
  if [ ! -f ${R2Bam}.bai ];then
    samtools index $R2Bam
  fi

  pids=()
  while read contig; do
      { nice samtools view -b $R2Bam $contig > $outpath/chr.bam.R2/$prefix/$contig.bam && samtools index $outpath/chr.bam.R2/$prefix/$contig.bam; } &
      pids+=($!)
      { nice samtools view -b $R1Bam $contig > $outpath/chr.bam.R1/$prefix/$contig.bam && samtools index $outpath/chr.bam.R1/$prefix/$contig.bam; } &
      pids+=($!)
  done < <(samtools view -H $R2Bam | grep -P '^@SQ' | cut -f 2 -d ':' | cut -f 1)
elif [ $# -eq 3 ]; then
  R2Bam=$1
  prefix=$2
  outpath=$3
  mkdir -p $outpath/chr.bam.R2/$prefix

  if [ ! -f ${R2Bam}.bai ];then
    samtools index $R2Bam
  fi

  pids=()
  while read contig; do
      { nice samtools view -b $R2Bam $contig > $outpath/chr.bam.R2/$prefix/$contig.bam && samtools index $outpath/chr.bam.R2/$prefix/$contig.bam; } &
      pids+=($!)
  done < <(samtools view -H $R2Bam | grep -P '^@SQ' | cut -f 2 -d ':' | cut -f 1)
fi

for pid in "${pids[@]}"; do
  wait $pid || { echo "samtools command failed"; exit 1; }
done
