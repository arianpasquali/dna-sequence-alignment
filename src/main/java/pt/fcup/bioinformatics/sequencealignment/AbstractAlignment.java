package pt.fcup.bioinformatics.sequencealignment;

import pt.fcup.bioinformatics.sequencealignment.costmatrix.CostMatrix;

/**
 * Created with IntelliJ IDEA.
 * User: arian
 * Date: 03/07/14
 * Time: 15:21
 */
public abstract class AbstractAlignment {


    public boolean isDebugMode() {
        return debugMode;
    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    private boolean debugMode;

    protected int gap = -1;
    protected int match = 2;
    protected int dismatch = -1;

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    protected String type;


    public AbstractAlignment(String type,int gap, int match, int dismatch) {
        this.type = type;

        this.gap = gap;
        this.match = match;
        this.dismatch = dismatch;
    }

    public AbstractAlignment(String type) {
        this(type,-1,2,-1);
    }

    public int getDismatch() {
        return dismatch;
    }

    public void setDismatch(int dismatch) {
        this.dismatch = dismatch;
    }

    public int getGap() {
        return gap;
    }

    public void setGap(int gap) {
        this.gap = gap;
    }

    public int getMatch() {
        return match;
    }

    public void setMatch(int match) {
        this.match = match;
    }

    protected String fixSequencePrefix(String sequence){
        //Para nao estar a alterar os indices por causa da primeira linha e primeira coluna da matriz
        if(sequence.startsWith(" ")){
            sequence = " ".concat(sequence);
        }

        return sequence;
    }

    // score. 1 para igual, -1 para diferente
    protected int score(char a, char b) {
        if (a == b) {
            return match;
        }
        return dismatch;
    }

    protected int[][] prepareResultMatrix(String x, String y) {
        //inicializacao da matriz
        int[][] mat = new int[x.length()][y.length()];
        mat[0][0] = 0;

        //primeira linha
        for (int i = 1; i < x.length(); i++) {
            mat[i][0] = mat[i - 1][0] + gap;
        }

        //primeira coluna
        for (int j = 1; j < y.length(); j++) {
            mat[0][j] = mat[0][j - 1] + gap;
        }
        return mat;
    }

    public abstract AlignmentResult align(String sequenceA, String sequenceB);
    public abstract AlignmentResult align(CostMatrix costMatrix, String sequenceA, String sequenceB);

    protected abstract AlignmentResult simpleAlignment(int[][] mat, String x, String y);
    protected abstract AlignmentResult blosumBasedAlignment(CostMatrix bcm, int[][] mat, String x, String y);
    protected abstract AlignmentResult pamBasedAlignment(CostMatrix bcm, int[][] mat, String x, String y);
}
