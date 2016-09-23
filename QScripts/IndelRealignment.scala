package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction

class IndelRealignment extends QScript {
    qscript =>

    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _ 

    @Input(doc="One or more bam files.", shortName="I",required=true)
    var bamFiles: List[File] = Nil

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of data threads", shortName="nt", required=true)
    var numDataThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Input(doc ="Input VCF file with known indels.", shortName="known", required=false)
    var knownIndelFiles: List[File] = Nil
    
    // This trait allows us set the variables below in one place,
    // and then reuse this trait on each CommandLineGATK function below.
    trait TCIR_Arguments extends CommandLineGATK {
			this.reference_sequence = qscript.referenceFile
			this.memoryLimit = maxMem
    }

    def script() {
			for (bamFile <- bamFiles) {
				val targetCreator = new RealignerTargetCreator with TCIR_Arguments
				val indelRealigner = new IndelRealigner with TCIR_Arguments

				targetCreator.input_file :+= bamFile
				if (knownIndelFiles != Nil) {
						targetCreator.known = knownIndelFiles
				}
				targetCreator.num_threads = numDataThreads
				targetCreator.scatterCount = numScatters
				targetCreator.out = swapExt(bamFile, "bam", "target.intervals.list")

				indelRealigner.targetIntervals = targetCreator.out
				indelRealigner.input_file +:= bamFile
				indelRealigner.scatterCount = numScatters
				indelRealigner.out = swapExt(bamFile, "bam", "realigned.bam")

				add(targetCreator, indelRealigner)
			}

    }
}