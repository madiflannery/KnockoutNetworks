import org.apache.commons.io.filefilter.DirectoryFileFilter;

import java.io.*;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Created by Madison on 25/07/2015.
 */
public class DiffExpProcess {
    private static final String DIFF_EXP_FILE = "C:\\Users\\Madison\\Desktop\\Data Mapping 2\\zDiff_MDA468_Control_vs_EGF_counts.txt";
    private static final String OFF_GENES = "C:\\Users\\Madison\\Desktop\\Data Mapping 2\\out_off_genes.txt";
    private static final int ensembl_col = 0,
                        ctrl1_col = 1,
                        ctrl2_col = 2,
                        ctrl3_col = 3;
    private static final int MAXREADS = 5;

    public static void main(String[] args) {
        HashSet<String> ids = readDiffExpFile();
        //writeToFile(ids);

    }

    private static void writeToFile(HashSet<String> ids) {
        Iterator keys = ids.iterator();
        PrintWriter out = null;
        try {
            out = new PrintWriter(new BufferedWriter(new FileWriter(OFF_GENES, true)));
        } catch (IOException e) {
            e.printStackTrace();
        }
        while(keys.hasNext()) {
            out.println(keys.next());
        }
        out.close();
    }

    public static HashSet<String> readDiffExpFile() {
        HashSet<String> ensembl_ids = new HashSet<>();

        //read all lines from file
        File file = new File(DIFF_EXP_FILE);
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            PrintWriter out = null;
            try {
                out = new PrintWriter(new BufferedWriter(new FileWriter(OFF_GENES)));
            } catch (IOException e) {
                e.printStackTrace();
            }

            String line;
            line = br.readLine();
            int num = 0;
            //split on '=' to separate variable and value
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("\t");
                if ((Integer.parseInt(splitted[ctrl1_col].trim()) == 0) || (Integer.parseInt(splitted[ctrl2_col].trim()) == 0) || (Integer.parseInt(splitted[ctrl3_col].trim()) == 0)) {
                    //at least one val is 0
                    if((Integer.parseInt(splitted[ctrl1_col].trim()) < MAXREADS) && (Integer.parseInt(splitted[ctrl2_col].trim()) < MAXREADS) && (Integer.parseInt(splitted[ctrl3_col].trim()) < MAXREADS)) {
                        //and all less than 5
                        System.out.println(splitted[ensembl_col]);
                        num++;
                        out.println(splitted[ensembl_col].trim()); // + "\t" + splitted[ctrl1_col].trim() +"\t"+ splitted[ctrl2_col].trim() +"\t"+ splitted[ctrl3_col].trim());
                    }
                }
            }
            out.close();
            br.close();
            System.out.println(num);
        } catch (IOException e) {
            //exit on error
            e.printStackTrace();
            System.exit(1);
        }
        return ensembl_ids;
    }


}
