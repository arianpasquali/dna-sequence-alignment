package pt.fcup.bioinformatics.sequencealignment;

/**
 * Created with IntelliJ IDEA.
 * User: arian
 * Date: 03/07/14
 * Time: 20:53
 */
public class AlignmentResult {

    public AlignmentResult(String methodName) {
        this.methodName = methodName;
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public String getAlignedSequenceA() {
        return alignedSequenceA;
    }

    public void setAlignedSequenceA(String alignedSequenceA) {
        this.alignedSequenceA = alignedSequenceA;
    }

    public String getMethodName() {
        return methodName;
    }

    public void setMethodName(String methodName) {
        this.methodName = methodName;
    }

    private int score;
    private String alignedSequenceA;

    public String getAlignedSequenceB() {
        return alignedSequenceB;
    }

    public void setAlignedSequenceB(String alignedSequenceB) {
        this.alignedSequenceB = alignedSequenceB;
    }

    private String alignedSequenceB;
    private String methodName;

    @Override
    public String toString() {
        return "AlignmentResult{" + "\n" +
                " method='" + methodName + '\'' + "\n" +
                " score=" + score + "\n\n" +

                alignedSequenceA + "\n" +
                alignedSequenceB + "\n\n" +

                '}';
    }

    public String getAlignment() {
        return alignedSequenceA + "\n" +
                alignedSequenceB;
    }
}
