package pt.fcup.bioinformatics.sequencealignment.costmatrix;

import java.io.*;
import java.util.Scanner;
import java.util.regex.MatchResult;

/**
 * Created with IntelliJ IDEA.
 * User: arian
 * Date: 03/07/14
 * Time: 15:33
 */
public abstract class CostMatrix {

    public int size;

    protected String name;
    protected int[][] matrix;

    protected CostMatrix(int size, String name) {
        this.size = size;
        this.name = name;

        this.initialise();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private void initialise() {
        if(this.name != null){
            try {
                InputStream is = this.getClass().getResourceAsStream(this.name + ".txt");
                this.matrix = readMatrixContentFromFile(is);
          } catch (IOException e) {
                e.printStackTrace();
      }
        }
        else{
            throw new IllegalStateException("matrix name must be informed (e.g. blosum, pam)");
        }
    }

    private String createPattern() {
       return "-?(\\d+)";
    }

    private int[][] readMatrixContentFromFile(InputStream path) throws IOException {
        //1,0,1,1,0,1,0,1,0,1
        int[][] result = new int[this.size][this.size];

        BufferedReader reader = new BufferedReader(new InputStreamReader(path));
        String line = null;
        Scanner scanner = null;
        line = reader.readLine();

        if(line == null) {
            return result;
        }
        String pattern = createPattern();
        int lineNumber = 0;
        MatchResult temp = null;
        while(line != null) {
            scanner = new Scanner(line);
            scanner.findInLine(pattern);
            temp = scanner.match();
            int count = temp.groupCount();
            for(int i=1;i<=count;i++) {
                result[lineNumber][i-1] = Integer.parseInt(temp.group(i));
            }
            lineNumber++;
            scanner.close();
            line = reader.readLine();
        }
        return result;
    }

    public int score(char a, char b) {
        try {
            return matrix[ a - 'A'][ b - 'A'];
        } catch (ArrayIndexOutOfBoundsException e) {
            String s = e.getMessage();
            return matrix[ 'J' - 'A'][ 'J' - 'A']; // by default, return unknown codon
        }
    }
}
