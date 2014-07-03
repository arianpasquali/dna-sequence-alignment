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

    private String sequenceA = "GAATTCAGTTA";
    private String sequenceB = "GGATCGA";

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
        assertEquals(5, result.getScore());

        assertEquals("AATTCAGTTA",result.getAlignedSequenceA());
        assertEquals("GA-TC-G--A",result.getAlignedSequenceB());

        assertEquals("AATTCAGTTA\n"+
                     "GA-TC-G--A",result.getAlignment());
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
        assertEquals(4,result.getScore());

        assertEquals("AATTCAGTTA",result.getAlignedSequenceA());
        assertEquals("----GATCGA",result.getAlignedSequenceB());

        assertEquals("AATTCAGTTA\n"+
                     "----GATCGA",result.getAlignment());
    }

    /**
     * Test alignment with PAM cost matrix
     */
    public void testPamAlighment()
    {

        GlobalAlignment ga = new GlobalAlignment();

        AlignmentResult result = ga.align(new PamCostMatrix(),sequenceA, sequenceB);

        System.out.println(result);

        assertEquals(12,result.getScore());
        assertEquals("global:pam",result.getMethodName());

        assertEquals("AATTCAGTTA",result.getAlignedSequenceA());
        assertEquals("--GA-TCGA-",result.getAlignedSequenceB());

        assertEquals("AATTCAGTTA\n"+
                     "--GA-TCGA-",result.getAlignment());
    }
}
