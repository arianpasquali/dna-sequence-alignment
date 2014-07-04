package pt.fcup.bioinformatics;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import pt.fcup.bioinformatics.sequencealignment.AlignmentResult;
import pt.fcup.bioinformatics.sequencealignment.LocalAlignment;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.BlosumCostMatrix;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.PamCostMatrix;

/**
 * Unit test for simple App.
 */
public class LocalAlignmentTest
    extends TestCase
{

    private String sequenceA = "WKVSVDQDTCIGDAICASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEAMEACPVSAI"; //Pyrococcus furiosus
    private String sequenceB = "MKVKLDKDTCIGCGVCASICPDVFEMDDDGKAKVIMEETDLECAKEAAESCPTGSI"; //Thermococcus sibiricus

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public LocalAlignmentTest(String testName)
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( LocalAlignmentTest.class );
    }

    /**
     * Test simple alignment
     */
    public void testSimpleAlighment()
    {

        LocalAlignment la = new LocalAlignment();

        AlignmentResult result = la.align(sequenceA, sequenceB);

        System.out.println(result);

        assertEquals("local:simple",result.getMethodName());
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

        LocalAlignment la = new LocalAlignment();

        AlignmentResult result = la.align(new BlosumCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals("local:blosum",result.getMethodName());
        assertEquals(14,result.getScore());

        assertEquals("ASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEA",result.getAlignedSequenceA());
        assertEquals("ASICPDVFEMDDDGKA-KVIMEETDLECAKEA--A",result.getAlignedSequenceB());

        assertEquals("ASLCPDVFEMNDEGKAQPKVEVIEDEELYNCAKEA\n"+
                     "ASICPDVFEMDDDGKA-KVIMEETDLECAKEA--A",result.getAlignment());
    }

    /**
     * Test alignment with PAM cost matrix
     */
    public void testPamAlighment()
    {

        LocalAlignment la = new LocalAlignment();

        AlignmentResult result = la.align(new PamCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals(15,result.getScore());
        assertEquals("local:pam",result.getMethodName());

        assertEquals("S--VDQDTCIGDAICAS",result.getAlignedSequenceA());
        assertEquals("AKVIMEETDLECAKEAA",result.getAlignedSequenceB());

        assertEquals("S--VDQDTCIGDAICAS\n"+
                     "AKVIMEETDLECAKEAA",result.getAlignment());
    }
}
