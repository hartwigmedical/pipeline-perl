package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import htsjdk.variant.variantcontext.VariantContext.Type
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils

class VariantFilter extends QScript {
    // Create an alias 'qscript' to be able to access variables in the HardFilter.
    // 'qscript' is now the same as 'HardFilter.this'
    qscript =>

    // Required arguments.  All initialized to empty values.
    @Input(doc="The reference file", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="Raw vcf file", shortName="V")
    var rawVCF: File =_

    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maximum amount of memory", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Argument(doc="SNP filter types", shortName="snpType", required=true)
    var snpFilterTypes: List[String] = _

    @Argument(doc="INDEL filter types", shortName="indelType", required=true)
    var indelFilterTypes: List[String] = _

    @Argument(doc="A optional list of SNPfilter names.", shortName="snpFilterName", required=false)
    var snpFilterNames: List[String] = _

    @Argument(doc="An optional list of filter expressions.", shortName="snpFilterExpression", required=false)
    var snpFilterExp: List[String] = _

    @Argument(doc="An optional list of INDEL filter names.", shortName="indelFilterName", required=false)
    var indelFilterNames: List[String] = _

    @Argument(doc="An optional list of INDEL filter expressions.", shortName="indelFilterExpression", required=false)
    var indelFilterExp: List[String] = _

    @Argument(doc="The number of SNPs which make up a cluster.", shortName="cluster", required=false)
    var clusterSize: Int = 0

    @Argument(doc="The window size (in bases) in which to evaluate clustered SNPs.", shortName="window", required=false)
    var clusterWindowSize: Int = 0

    // This trait allows us set the variables below in one place and then reuse this trait on each CommandLineGATK function below.
    trait HF_Arguments extends CommandLineGATK {
      this.reference_sequence = referenceFile
      this.memoryLimit = maxMem
    }

    def script() {
			val selectSNP = new SelectVariants with HF_Arguments
			selectSNP.V = rawVCF
			for (snpType <- snpFilterTypes) {
        if(snpType == "SNP") {
            selectSNP.selectType :+= Type.SNP
        } else if(snpType == "MNP") {
            selectSNP.selectType :+= Type.MNP
        } else if(snpType == "NO_VARIATION") {
            selectSNP.selectType :+= Type.NO_VARIATION
        }
	    }
      selectSNP.out = qscript.out + ".raw_snps.vcf"

      val SNPfilter = new VariantFiltration with HF_Arguments
      SNPfilter.scatterCount = numScatters
      SNPfilter.V = selectSNP.out
      SNPfilter.out = qscript.out + ".filtered_snps.vcf"
      SNPfilter.filterExpression = snpFilterExp
      SNPfilter.filterName = snpFilterNames

      if( clusterSize != 0 && clusterWindowSize != 0 ){
          SNPfilter.clusterSize = clusterSize
          SNPfilter.clusterWindowSize = clusterWindowSize
      }

      val selectINDEL = new SelectVariants with HF_Arguments
      selectINDEL.V = rawVCF
      for (indelType <- indelFilterTypes) {
        if (indelType == "INDEL") {
          selectINDEL.selectType :+= Type.INDEL
        } else if(indelType == "MIXED") {
          selectINDEL.selectType :+= Type.MIXED
        }
      }
      selectINDEL.out = qscript.out + ".raw_indels.vcf"

      val INDELfilter = new VariantFiltration with HF_Arguments
      INDELfilter.scatterCount = numScatters
      INDELfilter.V = selectINDEL.out
      INDELfilter.out = qscript.out + ".filtered_indels.vcf"
      INDELfilter.filterExpression = indelFilterExp
      INDELfilter.filterName = indelFilterNames

      val CombineVars = new CombineVariants with HF_Arguments
      CombineVars.V :+= SNPfilter.out
      CombineVars.V :+= INDELfilter.out
      CombineVars.out = qscript.out + ".filtered_variants.vcf"
      CombineVars.assumeIdenticalSamples = true

      add(selectSNP, SNPfilter)
      add(selectINDEL, INDELfilter)
      add(CombineVars)
    }
}
