dna-sequence-alignment
======================

Part of a bioinformatics course homework this is a sequence alignment library for DNA or proteins.
Supports local and global alignments considering BLOSUM or PAM cost matrix.

NOTE: 
web interface [dna-sequence-alignment-web](https://github.com/arianpasquali/dna-sequence-alignment-web)
console interface [dna-sequence-alignment-cli](https://github.com/arianpasquali/dna-sequence-alignment-cli)


Build
-----

Once you clone this repo just run maven install and check test results:

    mvn install test	 
	 
Usage
-----

Declare it as a maven dependency:

    <dependency>
		<groupId>pt.fcup.bioinformatics</groupId>
		<artifactId>dna-sequence-alignment</artifactId>
		<version>0.1-SNAPSHOT</version>
	 </dependency>

Example of global alignment using BLOSUM:

	//Pyrococcus furiosus
	String sequenceA = "WKVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI"; 
	
	//Thermococcus sibiricus
    String sequenceB = "MKVKLDKDTCIGCGVCASICPDVFEMDDDGKAKVIMEETDLECAKEAAESCPTGSI"; 
    
	GlobalAlignment ga = new GlobalAlignment();
    AlignmentResult result = ga.align(new BlosumCostMatrix(),sequenceA, sequenceB);

    System.out.println(result);	
    
result:

    AlignmentResult{
 		method='global:blosum'
	 	score=14

		KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI
		KVKLDKDTCIGCGVCASICPDVFEMDDDGKA-KVIMEETDLECAKEA--A-ESCPTGSI
	}
    
	 	 
	 

 
