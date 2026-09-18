#!/usr/bin/env python

rule index_reference:
    input:
        fasta = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN.fa'
    output:
        marker = touch('{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.done')
    conda: 
        '../envs/map_reads_env.yaml'
    threads:
        4
    resources:
        load = 50,
        mem = 10000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 1400 
    log:
        command = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.command',
        log = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.log',
        qerr = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.qerr',
        qout = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.qout'
    shell:
        '''
        command="bwa index {input.fasta}"
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''

def get_read_helper(sample, suffix):
    if "ISOG" in sample:
        readfile = seqstor_dict[sample + suffix]
    else:
        readfile = seqstor_dict[samples_new2old[sample] + suffix]
    if '*' in readfile:
        readsplit = readfile.split('*')
        if len(readsplit) != 2:
            raise Exception("If there is * in path, it should split in two...")
        readlist = list(Path(readsplit[0]).glob("*" + readsplit[1]))
        if len(readlist) != 1:
            raise Exception(f"There should one file matching {readfile}... but found {len(readlist)} ({readlist}), that shouldn't be.")
        readfile = readlist[0]

    return readfile


def give_me_time(wildcards, attempt):
    if wildcards.sample.startswith("TARA_"):
        time = 1400 if attempt == 1 else 7000
    else:
        time = 240 if attempt == 1 else 1400
    return time


rule map_reads:
    input:
        fasta = ancient('{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN.fa'),
        index = ancient('{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-bwa-index.done')
    output:
        marker = touch('{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.done'),
        # ------ Server version ------ 
        #bamm = temp('{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.m.filtered.sorted.bam'),
        #bams = temp('{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.s.filtered.sorted.bam'),
        #bamp = temp('{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.p.filtered.sorted.bam'),
        # ------ end version ------
        bam_all_filtered_sorted = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.all.filtered.sorted.bam'
    params:
        fasta = 'COMBINED_REF_GENOMES_FOR_POPGEN.fa',
        fqm = lambda wildcards: get_read_helper(wildcards.sample, '.m.fq.gz'),
        fq1 = lambda wildcards: get_read_helper(wildcards.sample, '.1.fq.gz'),
        fq2 = lambda wildcards: get_read_helper(wildcards.sample, '.2.fq.gz'),
        fqs = lambda wildcards: get_read_helper(wildcards.sample, '.s.fq.gz')
    conda:
        '../envs/map_reads_env.yaml'
    threads:
        24
    resources:
        load = 10,
        mem = 2500,#1000,#2500,
        scratch = 4000,
        time = 1400 
        #time = give_me_time 
    benchmark:
        '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.benchmark'
    log:
        command = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.command',
        log = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.log',
        qerr = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.qerr',
        qout = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.qout'
    shell:
        '''
        command="
        # ------ Euler version ------ 
        # sync to node's scratch
        cd $TMPDIR
        rsync -L {input.fasta}* . 
        if [ {wildcards.sample} == TARA_* ];
        then
            rsync -L {params.fqm} m.fq.gz
            rsync -L {params.fqs} s.fq.gz
            rsync -L {params.fq1} r1.fq.gz
            rsync -L {params.fq2} r2.fq.gz
        else
            bbrename.sh -Xmx4G threads={threads} pigz=t bgzip=f in1={params.fq1} in2={params.fq2} out1=r1.fq.gz out2=r2.fq.gz prefix={wildcards.sample}-p
            bbrename.sh -Xmx4G threads={threads} pigz=t bgzip=f in1={params.fqm} out1=m.fq.gz prefix={wildcards.sample}-m
            bbrename.sh -Xmx4G threads={threads} pigz=t bgzip=f in1={params.fqs} out1=s.fq.gz prefix={wildcards.sample}-s
        fi
        
        # Map merged reads
        bwa mem -a -t {threads} {params.fasta} m.fq.gz | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > m.bam 
        #
        # Map single reads
        bwa mem -a -t {threads} {params.fasta} s.fq.gz | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > s.bam 
        #    
        # Map paired reads
        bwa mem -a -t {threads} {params.fasta} r1.fq.gz r2.fq.gz | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > p.bam 
        #
        # Merge filtered and sorted bamfiles
        samtools merge {output.bam_all_filtered_sorted} m.bam s.bam p.bam 
        # ----- end version -----

        # ----- Server version -----
        # FIXME some reads will need renaming
        # Map merged reads
        #bwa mem -a -t {threads} {input.fasta} {params.fqm} | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > {{output.bamm}}
        #
        # Map single reads
        #bwa mem -a -t {threads} {input.fasta} {params.fqs} | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > {{output.bams}}
        #
        # Map paired reads
        #bwa mem -a -t {threads} {input.fasta} {params.fq1} {params.fq2} | samtools view -F 4 -h - | sushicounter filter -u -i 0.95 -c 0.8 -a 45 - - | samtools view -bh - | samtools sort -O bam -@ 4 -m 4G - > {{output.bamp}}
        #
        # Merge filtered and sorted bamfiles
        #samtools merge {output.bam_all_filtered_sorted} {{output.bamm}} {{output.bams}} {{output.bamp}}
        # ----- end verion -----

        # Index final bam file
        samtools index {output.bam_all_filtered_sorted}
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule anvio_profile:
    input:
        contigs_db = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db',
        contigs_list = '{path}/species_clusters/{species}/{reference}-contigs.txt',
        merge_bam = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.all.filtered.sorted.bam',
        marker = '{path}/refgenomes/bams/{sample}/{sample}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.done'
    output:
        marker = touch('{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}.anvio.done'),
        anvio_profile = directory('{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}-ANVIO_PROFILE')
    params:
        anvio_sample = lambda wildcards: wildcards.sample.replace("-", "_")
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        12#24#12
    resources:
        load = 10,
        scratch = 500,
        mem = 5000,#2500,
        #mem = lambda wildcards, attempt: attempt * 5000, 
        time = 240#1400 
        #time = lambda wildcards, attempt: 240 if attempt == 1 else 1400 
    log:
        command = '{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}.anvio.command',
        log = '{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}.anvio.log',
        qerr = '{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}.anvio.qerr',
        qout = '{path}/species_clusters/{species}/profiles/{sample}-vs-{reference}.anvio.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8

        anvi-profile --input-file {input.merge_bam} \
            --contigs-db {input.contigs_db} \
            --contigs-of-interest {input.contigs_list} \
            --output-dir {output.anvio_profile} \
            --sample-name {params.anvio_sample} \
            --min-coverage-for-variability 3 \
            --profile-SCVs \
            --num-threads {threads} 
            # actually not so good to use variability full because it doesn't do what I had hoped
            # --report-variability-full \
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule anvio_contigs_db:
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-CONTIGS.done'),
        contigs_db = '{path}/species_clusters/{species}/{reference}-CONTIGS.db'
    params:
        fasta = '{path}/species_clusters/{species}/{reference}.fa'
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        8
    resources:
        load = 50,
        mem = 2000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 240
    log:
        log = '{path}/species_clusters/{species}/{reference}-CONTIGS.log',
        command = '{path}/species_clusters/{species}/{reference}-CONTIGS.command',
        qerr = '{path}/species_clusters/{species}/{reference}-CONTIGS.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-CONTIGS.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8;

        anvi-gen-contigs-database --project-name {wildcards.reference} \
            --contigs-fasta {params.fasta} \
            --output-db-path {output.contigs_db} \
            -T {threads}
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule anvio_merge_profile:
    input:
        contigs_db = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db',
        profiles = lambda wildcards: input_anvio_merge[wildcards.species]
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.done')
    params:
        prefix = '{path}/species_clusters/{species}/{reference}',
        profiles = '{path}/species_clusters/{species}/profiles' 
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        6 
    resources:
        load = 100,
        mem = lambda wildcards, attempt: attempt * 10000, 
        scratch = 500,
        time = 240 
    log:
        command = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.command',
        log = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.log',
        qerr = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8

        # Merge profiles, note: Takes the profile as first positional argument
        anvi-merge {params.profiles}/*-vs-{wildcards.reference}-ANVIO_PROFILE/*PROFILE.db \
            --output-dir {params.prefix}-MERGED_PROFILE \
            --contigs-db {input.contigs_db} \
            --sample-name {wildcards.reference}

        # Just add a DEFAULT collection with a single bin containing all the contigs named GENOME
        anvi-script-add-default-collection \
            --pan-or-profile-db {params.prefix}-MERGED_PROFILE/PROFILE.db --bin-id {wildcards.reference} 
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule anvio_genes_consensus:
    input:
        contigs_db = ancient('{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db'),
        combined_gff = ancient('{path}/refgenomes/go_microbiomics-integrated-cpl50_ctn10-combined-prokka.gff'),
        combined_gene_calls = ancient('{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-external_gene_ids.tsv'),
        merged_profile_db = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE/PROFILE.db'
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.done')
    params:
        genes_names = '{path}/species_clusters/{species}/{reference}-GENES.names',
        genes_map = '{path}/species_clusters/{species}/{reference}-GENES-ANVIO.map',
        anvio_ids = '{path}/species_clusters/{species}/{reference}-GENES-ANVIO.ids',
        consensus = '{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.fasta'
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        1 
    resources:
        load = 200,
        mem = 50000, 
        scratch = 500,
        time = 240 
    log:
        command = '{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.command',
        log = '{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.log',
        qerr = '{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-GENES_CONSENSUS.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8
        
        grep -w {wildcards.reference} {input.combined_gff} | grep -w CDS | cut -f10 | sed 's/;.*//g' | cut -f2 -d'=' > {params.genes_names} 
        grep -wf {params.genes_names} {input.combined_gene_calls} > {params.genes_map} 
        cut -f1 {params.genes_map} > {params.anvio_ids} 

        anvi-gen-gene-consensus-sequences \
                -p {input.merged_profile_db} \
                -c {input.contigs_db} \
                -o {params.consensus} \
                --genes-of-interest {params.anvio_ids} \
                --quince-mode
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule anvio_variability:
    #input:
    #    contigs_db = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db',
    #    marker = '{path}/species_clusters/{species}/{reference}-MERGED_PROFILE.done'
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.done')
    params:
        contigs_db = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db',
        prefix = '{path}/species_clusters/{species}/{reference}'
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        32 
    resources:
        load = 50,
        mem = 70000,#lambda wildcards, attempt: attempt * 50000, 
        scratch = 500,
        time = 1400#lambda wildcards, attempt: attempt * 1400 
    log:
        command = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.command',
        log = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.log',
        qerr = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8

        anvi-gen-variability-profile \
            -c {params.contigs_db} \
            -p {params.prefix}-MERGED_PROFILE/PROFILE.db \
            -C "DEFAULT" \
            -b {wildcards.reference} \
            -o {params.prefix}-ANVIO_VARIABILITY-NT.tsv \
            --engine NT \
            --quince-mode \
            --include-contig-names

        anvi-gen-variability-profile \
            -c {params.contigs_db} \
            -p {params.prefix}-MERGED_PROFILE/PROFILE.db \
            -C "DEFAULT" \
            -b {wildcards.reference} \
            -o {params.prefix}-ANVIO_VARIABILITY-CDN.tsv \
            --engine CDN \
            --kiefl-mode \
            --include-contig-names
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule compute_pnps:
    input:
        codon_profile = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN.tsv',
        contigs_db = '{path}/refgenomes/COMBINED_REF_GENOMES_FOR_POPGEN-CONTIGS.db',
        marker = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.done'
    output:
        pnps = directory('{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS'),
        marker = touch('{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS.done')
    #conda:
    #    '../envs/anvio_env.yaml'
    threads:
        1 
    resources:
        load = 20,
        mem = 1000000,#lambda wildcards, attempt: attempt * 50000, 
        scratch = 500,
        time = 240 #lambda wildcards, attempt: attempt * 1400 
    log:
        command = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS.command',
        log = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS.log',
        qerr = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-CDN-PNPS.qout'
    shell:
        '''
        command="
        export LC_ALL=en_US.utf8

        anvi-get-pn-ps-ratio \
                -V {input.codon_profile} \
                -c {input.contigs_db} \
                -o {output.pnps} 
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule identify_populations:
    input:
        nucleotide_profile = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-NT.tsv',
        marker = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.done'
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-POPULATIONS.done')
    conda:
        '../envs/pylearn_env.yaml'
    threads:
        32 
    resources:
        load = 10,
        mem = 4000,#lambda wildcards, attempt: attempt * 50000, 
        scratch = 500,
        time = 7000#lambda wildcards, attempt: attempt * 1400 
    log:
        command = '{path}/species_clusters/{species}/{reference}-POPULATIONS.command',
        log = '{path}/species_clusters/{species}/{reference}-POPULATIONS.log',
        qerr = '{path}/species_clusters/{species}/{reference}-POPULATIONS.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-POPULATIONS.qout'
    shell:
        '''
        command="
        export OMP_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export MKL_NUM_THREADS={threads}
        export GOTO_NUM_THREADS={threads}
        export NUMBA_NUM_THREADS={threads}
         
        python -W ignore cluster_populations_from_nt.py -i {input.nucleotide_profile} -t {threads} 
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''


rule neutral_distances:
    input:
        nucleotide_profile = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY-NT.tsv',
        marker = '{path}/species_clusters/{species}/{reference}-ANVIO_VARIABILITY.done'
    output:
        marker = touch('{path}/species_clusters/{species}/{reference}-NEUTRAL.done')
    conda:
        '../envs/pylearn_env.yaml'
    threads:
        1 
    resources:
        load = 10,
        mem = 300000,#lambda wildcards, attempt: attempt * 50000, 
        scratch = 500,
        time = 240#lambda wildcards, attempt: attempt * 1400 
    log:
        command = '{path}/species_clusters/{species}/{reference}-NEUTRAL.command',
        log = '{path}/species_clusters/{species}/{reference}-NEUTRAL.log',
        qerr = '{path}/species_clusters/{species}/{reference}-NEUTRAL.qerr',
        qout = '{path}/species_clusters/{species}/{reference}-NEUTRAL.qout'
    shell:
        '''
        command="
        python neutral_from_nt.py -i {input.nucleotide_profile}
        "
        echo "$command" > {log.command}
        eval "$command" &> {log.log}
        '''

