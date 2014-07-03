package pt.fcup.bioinformatics.sequencealignment;

import pt.fcup.bioinformatics.sequencealignment.costmatrix.BlosumCostMatrix;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.CostMatrix;
import pt.fcup.bioinformatics.sequencealignment.costmatrix.PamCostMatrix;

// Needleman–Wunsch algorithm
public class GlobalAlignment extends AbstractAlignment {

    public GlobalAlignment(int gap, int match, int dismatch) {
        super("global",gap, match, dismatch);
    }

    public GlobalAlignment() {
        super("global");
    }

    @Override
    public AlignmentResult align(CostMatrix costMatrix, String sequenceA, String sequenceB){
        AlignmentResult result = null;

        int[][] mat = prepareResultMatrix(sequenceA, sequenceB);

        if(costMatrix == null){
            return simpleAlignment(mat, sequenceA, sequenceB);
        }
        else
            if(costMatrix instanceof BlosumCostMatrix){
                return  blosumBasedAlignment(costMatrix, mat, sequenceA, sequenceB);
            }else
                if(costMatrix instanceof PamCostMatrix){
                    return pamBasedAlignment(costMatrix , mat,sequenceA,sequenceB);
                }

        return result;
    }

    @Override
    public AlignmentResult align(String sequenceA, String sequenceB){
        sequenceA = fixSequencePrefix(sequenceA);
        sequenceB = fixSequencePrefix(sequenceB);

        return align(null, sequenceA,sequenceB);
    }

    @Override
    protected AlignmentResult simpleAlignment(int[][] mat, String x, String y) {

        AlignmentResult result = new AlignmentResult(getType() + ":simple");

        //ciclo principal
        for (int i = 1; i < x.length(); i++) {
            for (int j = 1; j < y.length(); j++) {
                int calc1 = mat[i - 1][j - 1] + score(x.charAt(i), y.charAt(j));
                int calc2 = mat[i - 1][j] + gap;
                int calc3 = mat[i][j - 1] + gap;
                mat[i][j] = Math.max(Math.max(calc1, calc2), calc3);
            }
        }

        //impressao da matriz
        if(this.isDebugMode()){
            for (int i = 0; i < x.length(); i++) {
                for (int j = 0; j < y.length(); j++) {
                    System.out.printf("%d ", mat[i][j]);
                }
                System.out.println();
            }
        }

        //score final
        int finalScore = mat[x.length() - 1][y.length() - 1];
        result.setScore(finalScore);

        String a = "";
        String b = "";
        int i = x.length() - 1;
        int j = y.length() - 1;
        while (i > 0 && j > 0) {
            int score = mat[i][j];
            int scoreDiag = mat[i - 1][j - 1];
            int scoreUp = mat[i][j - 1];
            int scoreLeft = mat[i - 1][j];

            if (score == (scoreDiag + score(x.charAt(i), y.charAt(j)))) {
                a = x.charAt(i) + a;
                b = y.charAt(j) + b;
                i = i - 1;
                j = j - 1;
            } else if (score == scoreLeft + gap) {
                a = x.charAt(i) + a;
                b = "-" + b;
                i = i - 1;
            } else {
                b = y.charAt(j) + b;
                a = "-" + a;
                j = j - 1;
            }
        }
        while (i > 0) {
            a = x.charAt(i) + a;
            b = "-" + b;
            i--;
        }
        while (j > 0) {
            a = "-" + a;
            b = y.charAt(j) + b;
            j--;
        }

        result.setAlignedSequenceA(a);
        result.setAlignedSequenceB(b);

        return result;
    }

    @Override
    protected AlignmentResult blosumBasedAlignment(CostMatrix bcm, int[][] mat, String x, String y) {

        AlignmentResult result = new AlignmentResult(getType() + ":" + bcm.getName());

        for (int i = 1; i < x.length(); i++) {
            for (int j = 1; j < y.length(); j++) {
                int calc1 = mat[i - 1][j - 1] + bcm.score(x.charAt(i), y.charAt(j));
                int calc2 = mat[i - 1][j] + gap;
                int calc3 = mat[i][j - 1] + gap;
                mat[i][j] = Math.max(Math.max(calc1, calc2), calc3);
            }
        }

        //impressao da matriz
        if(this.isDebugMode()){
            for (int i = 0; i < x.length(); i++) {
                for (int j = 0; j < y.length(); j++) {
                    System.out.printf("%d ", mat[i][j]);
                }
                System.out.println();
            }
        }

        //score final
        int finalScore = mat[x.length() - 1][y.length() - 1];

        result.setScore(finalScore);

        String a = "";
        String b = "";
        int i = x.length() - 1;
        int j = y.length() - 1;
        while (i > 0 && j > 0) {
            int score = mat[i][j];
            int scoreDiag = mat[i - 1][j - 1];
            int scoreUp = mat[i][j - 1];
            int scoreLeft = mat[i - 1][j];

            if (score == (scoreDiag + bcm.score(x.charAt(i), y.charAt(j)))) {
                a = x.charAt(i) + a;
                b = y.charAt(j) + b;
                i = i - 1;
                j = j - 1;
            } else if (score == scoreLeft + gap) {
                a = x.charAt(i) + a;
                b = "-" + b;
                i = i - 1;
            } else {
                a = "-" + a;
                b = y.charAt(j) + b;
                j = j - 1;
            }
        }
        while (i > 0) {
            a = x.charAt(i) + a;
            b = "-" + b;
            i--;
        }
        while (j > 0) {
            a = "-" + a;
            b = y.charAt(j) + b;
            j--;
        }

        result.setAlignedSequenceA(a);
        result.setAlignedSequenceB(b);

        return result;
    }

    @Override
    protected AlignmentResult pamBasedAlignment(CostMatrix pcm,int[][] mat , String x, String y) {

        AlignmentResult result = new AlignmentResult(getType() + ":" +pcm.getName());

        for (int i = 1; i < x.length(); i++) {
            for (int j = 1; j < y.length(); j++) {
                int calc1 = mat[i - 1][j - 1] + pcm.score(x.charAt(i), y.charAt(j));
                int calc2 = mat[i - 1][j] + gap;
                int calc3 = mat[i][j - 1] + gap;
                mat[i][j] = Math.max(Math.max(calc1, calc2), calc3);
            }
        }

        //impressao da matriz
        if(this.isDebugMode()){
            for (int i = 0; i < x.length(); i++) {
                for (int j = 0; j < y.length(); j++) {
                System.out.printf("%d ", mat[i][j]);
                }
                System.out.println();
            }
        }

        //score final
        int finalScore = mat[x.length() - 1][y.length() - 1];
        result.setScore(finalScore);

        String a = "";
        String b = "";
        int i = x.length() - 1;
        int j = y.length() - 1;
        while (i > 0 && j > 0) {
            int score = mat[i][j];
            int scoreDiag = mat[i - 1][j - 1];
            int scoreUp = mat[i][j - 1];
            int scoreLeft = mat[i - 1][j];

            if (score == (scoreDiag + pcm.score(x.charAt(i), y.charAt(j)))) {
                a = x.charAt(i) + a;
                b = y.charAt(j) + b;
                i = i - 1;
                j = j - 1;
            } else if (score == scoreLeft + gap) {
                a = x.charAt(i) + a;
                b = "-" + b;
                i = i - 1;
            } else {
                a = "-" + a;
                b = y.charAt(j) + b;
                j = j - 1;
            }
        }
        while (i > 0) {
            a = x.charAt(i) + a;
            b = "-" + b;
            i--;
        }
        while (j > 0) {
            a = "-" + a;
            b = y.charAt(j) + b;
            j--;
        }

        result.setAlignedSequenceA(a);
        result.setAlignedSequenceB(b);

        return result;
    }
}
