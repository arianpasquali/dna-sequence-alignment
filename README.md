dna-sequence-alignment
======================

Part of a bioinformatics course homework this is a sequence alignment library for DNA or proteins.
Supports local and global alignments considering BLOSUM or PAM cost matrix.

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

	String sequenceA = "GAATTCAGTTA";
    String sequenceB = "GGATCGA";
    
	GlobalAlignment ga = new GlobalAlignment();
    AlignmentResult result = ga.align(new BlosumCostMatrix(),sequenceA, sequenceB);
	
	System.out.println(result);    
    
result:

    AlignmentResult{
 		method='global:blosum'
	 	score=4

		AATTCAGTTA
		----GATCGA
	}
    
	 	 
	 

 
