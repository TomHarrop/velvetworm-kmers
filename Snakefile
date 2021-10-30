#!/usr/bin/env python3

from pathlib import Path


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return Path(x).resolve().as_posix()


###########
# GLOBALS #
###########

bbmap_container = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.3/container-bbmap.bbmap_38.90.sif'
merqury = 'https://github.com/deardenlab/container-merqury/releases/download/v1.3/container-merqury.v1.3.sif'

########
# MAIN #
########

#########
# RULES #
#########

rule target:
    input:
        'output/020_bbnorm/hist_out.txt',
        'output/030_merqury/illumina.meryl/merylIndex'


rule merqury_make_kmers:
    input:
        fq = 'output/010_read-prep/short_reads.fq'
    output:
        'output/030_merqury/illumina.meryl/merylIndex'
    params:
        wd = 'output/030_merqury',
        fq = lambda wildcards, input: resolve_path(input.fq)
    log:
        resolve_path('output/logs/merqury_make_kmers.log')
    threads:
        workflow.cores
    container:
        merqury
    shell:
        'cd {params.wd} || exit 1 ; '
        'bash -c \''
        'meryl '
        'k=21 '
        'cpus={threads} '
        'count '
        'output illumina.meryl '
        '{params.fq} '
        '\''
        '&> {log}'


rule bbnorm:
    input:
        fq = 'output/010_read-prep/short_reads.fq'
    output:
        fq_norm = 'output/020_bbnorm/short_reads.norm.fq.gz',
        fq_toss = 'output/020_bbnorm/short_reads.toss.fq.gz',
        hist = 'output/020_bbnorm/hist.txt',
        hist_out = 'output/020_bbnorm/hist_out.txt',
        peaks = 'output/020_bbnorm/peaks.txt'
    log:
        norm = 'output/logs/bbnorm.log'
    params:
        target = 60
    threads:
        workflow.cores
    container:
        bbmap_container
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min=5 '
        'peaks={output.peaks} '
        '2> {log.norm} '



rule combine_illumina:
    input:
        expand('output/010_read-prep/run{run}.trim.fastq',
               run=['1', '2'])
    output:
        'output/010_read-prep/short_reads.fq'
    container:
        bbmap_container
    shell:
        'cat {input} > {output}'

rule trim_illumina:
    input:
        'output/010_read-prep/run{run}.filter.fastq'
    output:
        p = temp('output/010_read-prep/run{run}.trim.fastq'),
        stats = 'output/010_read-prep/run{run}_trim.txt'
    log:
        'output/logs/trim_illumina.{run}.log'
    params:
        trim = '/adapters.fa'
    threads:
        (workflow.cores // 4) - 1
    container:
        bbmap_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.stats} '
        '>> {output.p} '
        '2> {log} '


rule filter_illumina:
    input:
        'output/010_read-prep/run{run}.repair.fastq'
    output:
        p = pipe('output/010_read-prep/run{run}.filter.fastq'),
        stats = 'output/010_read-prep/run{run}.filter.txt'
    log:
        'output/logs/filter_illumina.{run}_filter.log'
    params:
        filter = '/phix174_ill.ref.fa.gz'
    threads:
        (workflow.cores // 4) - 1
    container:
        bbmap_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.stats} '
        '>> {output.p} '
        '2> {log}'

rule repair_illumina:
    input:
        r1 = 'data/illumina_run{run}/CC481_R1.fq.gz',
        r2 = 'data/illumina_run{run}/CC481_R2.fq.gz',
    output:
        p = pipe('output/010_read-prep/run{run}.repair.fastq'),
    log:
        'output/logs/repair_illumina.{run}.log'
    container:
        bbmap_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '>> {output.p} '
        '2> {log}'
