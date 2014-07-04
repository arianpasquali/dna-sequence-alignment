package pt.fcup.bioinformatics;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import pt.fcup.bioinformatics.sequencealignment.AlignmentResult;
import pt.fcup.bioinformatics.sequencealignment.GlobalAlignment;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.BlosumCostMatrix;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.PamCostMatrix;

/**
 * Unit test for simple App.
 */
public class GlobalAlignmentTest
    extends TestCase
{

    private String sequenceA = "WKVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI"; //Pyrococcus furiosus
    private String sequenceB = "MKVKLDKDTCIGCGVCASICPDVFEMDDDGKAKVIMEETDLECAKEAAESCPTGSI"; //Thermococcus sibiricus

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public GlobalAlignmentTest(String testName)
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( GlobalAlignmentTest.class );
    }

    /**
     * Test simple alignment
     */
    public void testSimpleAlighment()
    {

        GlobalAlignment ga = new GlobalAlignment();

        AlignmentResult result = ga.align(sequenceA, sequenceB);

        System.out.println(result);

        assertEquals("global:simple",result.getMethodName());
        assertEquals(52, result.getScore());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVI-EDE-ELYNCAKEAMEACP-VSAI",result.getAlignedSequenceA());
        assertEquals("KVKLDKDTCIGCGVCASICPDVFEMDDDGKA--K--VIME-ETDL-ECAKEAAESCPTGS-I",result.getAlignedSequenceB());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVI-EDE-ELYNCAKEAMEACP-VSAI\n"+
                     "KVKLDKDTCIGCGVCASICPDVFEMDDDGKA--K--VIME-ETDL-ECAKEAAESCPTGS-I",result.getAlignment());
    }

    /**
     * Test alignment with BLOSUM cost matrix
     */
    public void testBlosumAlighment()
    {

        GlobalAlignment ga = new GlobalAlignment();

        AlignmentResult result = ga.align(new BlosumCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals("global:blosum",result.getMethodName());
        assertEquals(14,result.getScore());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI",result.getAlignedSequenceA());
        assertEquals("KVKLDKDTCIGCGVCASICPDVFEMDDDGKA-KVIMEETDLECAKEA--A-ESCPTGSI",result.getAlignedSequenceB());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI\n"+
                     "KVKLDKDTCIGCGVCASICPDVFEMDDDGKA-KVIMEETDLECAKEA--A-ESCPTGSI",result.getAlignment());
    }

    /**
     * Test alignment with PAM cost matrix
     */
    public void testPamAlighment()
    {

        GlobalAlignment ga = new GlobalAlignment();

        AlignmentResult result = ga.align(new PamCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals(15,result.getScore());
        assertEquals("global:pam",result.getMethodName());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI",result.getAlignedSequenceA());
        assertEquals("KVKLDKDTCIGCGVCASICPDVFEMDDDGKAKVIMEETDLECA-KEA--A-ESCPTGSI",result.getAlignedSequenceB());

        assertEquals("KVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI\n"+
                     "KVKLDKDTCIGCGVCASICPDVFEMDDDGKAKVIMEETDLECA-KEA--A-ESCPTGSI",result.getAlignment());
    }
}
