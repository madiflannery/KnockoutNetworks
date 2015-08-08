import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Created by Madison on 28/07/2015.
 */
public class ProcessNetworkSpeciesFile {
    private static final String NETWORK_FILE = "species.txt";
    private static final int RDF_COL = 1; //0 indexed
    private static final int INPUT_COL = 5;
    //public static ArrayList<String[]> network = dataMapping.readFile(NETWORK_FILE, false);

    public static HashSet<String> findInputRDFs(ArrayList<String[]> network) {
        int count = 0;
        HashSet<String> result = new HashSet<>();
        for(String[] netLine : network) {
            if(netLine[INPUT_COL].contains("INPUT")) {
                if(netLine[RDF_COL].contains("Complex") || netLine[RDF_COL].contains("Protein")) {
                    result.add(netLine[RDF_COL]);
                    count++;
                }
            }
        }
        System.out.println("Number of INPUT's: "+ count);
        return result;
    }

    public static ArrayList<ArrayList<String>> switchFormat(ArrayList<String[]> network) {
        ArrayList<ArrayList<String>> result = new ArrayList<>();
        for(String[] s : network) {
            result.add(new ArrayList<String>(Arrays.asList(s)));
        }
        return result;
    }


}
