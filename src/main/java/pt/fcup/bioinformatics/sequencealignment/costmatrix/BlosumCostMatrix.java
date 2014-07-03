package pt.fcup.bioinformatics.sequencealignment.costmatrix;


public class BlosumCostMatrix extends CostMatrix {

    public BlosumCostMatrix() {
        super(27, "blosum");
    }

//    public static void main(String[] a) throws IOException {
//        new BlosumCostMatrix().loadMatrix();
//    }
}