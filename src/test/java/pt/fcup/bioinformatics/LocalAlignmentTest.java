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

    private String sequenceA = "GAATTCAGTTA";
    private String sequenceB = "GGATCGA";

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
        assertEquals(6, result.getScore());

        assertEquals("ATTC-A",result.getAlignedSequenceA());
        assertEquals("A-TCGA",result.getAlignedSequenceB());

        assertEquals("ATTC-A\n"+
                     "A-TCGA",result.getAlignment());
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
        assertEquals(8,result.getScore());

        assertEquals("ATTCA",result.getAlignedSequenceA());
        assertEquals("ATCGA",result.getAlignedSequenceB());

        assertEquals("ATTCA\n"+
                     "ATCGA",result.getAlignment());
    }

    /**
     * Test alignment with PAM cost matrix
     */
    public void testPamAlighment()
    {

        LocalAlignment la = new LocalAlignment();

        AlignmentResult result = la.align(new PamCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals(14,result.getScore());
        assertEquals("local:pam",result.getMethodName());

        assertEquals("TCAGT",result.getAlignedSequenceA());
        assertEquals("ATCGA",result.getAlignedSequenceB());

        assertEquals("TCAGT\n"+
                     "ATCGA",result.getAlignment());
    }
}
