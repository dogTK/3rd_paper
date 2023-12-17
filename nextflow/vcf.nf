params.inVcf = '../rawdata/PG5103_C01_var_eff.vcf.gz'
params.outVcf = 'filtered_tech.vcf'

process filterVcf {
    input:
    path inVcf

    output:
    path "${params.outVcf}"

    script:
    """
    java -jar ~/snpEff/SnpSift.jar filter "( GEN[PG5103_01_b].DP >= 25 ) & ( GEN[PG5103_01_b].AF <= 0.08 )" ${inVcf} > ${params.outVcf}
    """
}

workflow {
    inVcf = file(params.inVcf)
    filterVcf(inVcf)
}